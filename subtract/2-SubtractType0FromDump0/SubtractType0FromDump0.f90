!Subtract positions,a-b=c
!
!Input : a.input 
!Input : b.input 
!Output : c.output
!
!lwh
!2023.4.21
!----------------------------------------------------------------------

program main
implicit none

!-----------state variables
!-----varis:GLOBAL
integer :: i,j,k,read1_flag,read2_flag,flag
!-----varis:PART-1,read in dump.0
integer :: read1_num,read1_typenum
double precision :: xlo1,xhi1,ylo1,yhi1,zlo1,zhi1
integer,allocatable,dimension(:) :: read1_id,read1_type
double precision,allocatable,dimension(:,:) :: read1_pos
!-----varis:PART-2,read in awaydislocations.output
integer :: read2_num,read2_typenum
double precision :: xlo2,xhi2,ylo2,yhi2,zlo2,zhi2
integer,allocatable,dimension(:) :: read2_id,read2_type,read2_clu,read2_clu_size
double precision,allocatable,dimension(:,:) :: read2_pos
!-----varis:PART-3,a-b
integer :: c_num
double precision :: x1,y1,z1,x2,y2,z2
integer,allocatable,dimension(:) :: c_id,c_type
double precision,allocatable,dimension(:,:) :: c_pos


!----------GLOBAL INIT
read1_flag = 2 !1 for data,2 for dump
read2_flag = 1 !1 for data,2 for dump
!----------PART-1(READ-IN-1):read in dump.0 , lammps-DUMP-format
if(read1_flag==1)then
open(1,file="a.input",status="old")
read(1,*)
read(1,*)read1_num !number-of-atoms
read(1,*)read1_typenum
read(1,*)xlo1,xhi1
read(1,*)ylo1,yhi1
read(1,*)zlo1,zhi1
read(1,*)
read(1,*)
read(1,*)
allocate(read1_pos(1:read1_num,1:3))
allocate(read1_id(1:read1_num))
allocate(read1_type(1:read1_num))
do i = 1 , read1_num
    read(1,*)read1_id(i),read1_type(i),read1_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read1_pos:",read1_num
endif
if(read1_flag==2)then
open(1,file="a.input",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)read1_num !number-of-atoms
read(1,*)
read1_typenum = 2
read(1,*)xlo1,xhi1
read(1,*)ylo1,yhi1
read(1,*)zlo1,zhi1
read(1,*)
allocate(read1_pos(1:read1_num,1:3))
allocate(read1_id(1:read1_num))
allocate(read1_type(1:read1_num))
do i = 1 , read1_num
    read(1,*)read1_id(i),read1_type(i),read1_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read1_pos:",read1_num
endif
!----------PART-2(READ-IN-2):read in type0.data , lammps-data-format
if(read2_flag==1)then
open(1,file="b.input",status="old")
read(1,*)
read(1,*)read2_num !number-of-atoms
read(1,*)read2_typenum
read(1,*)xlo2,xhi2
read(1,*)ylo2,yhi2
read(1,*)zlo2,zhi2
allocate(read2_pos(1:read2_num,1:3))
allocate(read2_id(1:read2_num))
allocate(read2_type(1:read2_num))
allocate(read2_clu(1:read2_num))
read(1,*)
read(1,*)
read(1,*)
do i = 1 , read2_num
    read(1,*)read2_id(i),read2_type(i),read2_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read2_pos:",read2_num
endif
if(read2_flag==2)then
open(1,file="b.input",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)read2_num !number-of-atoms
read(1,*)
read2_typenum = 2
read(1,*)xlo2,xhi2
read(1,*)ylo2,yhi2
read(1,*)zlo2,zhi2
allocate(read2_pos(1:read2_num,1:3))
allocate(read2_id(1:read2_num))
allocate(read2_type(1:read2_num))
allocate(read2_clu(1:read2_num))
allocate(read2_clu_size(1:read2_num))
read(1,*)
do i = 1 , read2_num
    read(1,*)read2_id(i),read2_type(i),read2_pos(i,1:3),read2_clu(i),read2_clu_size(i)
enddo
close(1)
write(*,*)"-->The num of atoms in read2_pos:",read2_num
endif
!----------PART-3():a-b
allocate(c_id(1:read1_num))
allocate(c_type(1:read1_num))
allocate(c_pos(1:read1_num,1:3))
c_num = 0
do i = 1 , read1_num
    x1 = read1_pos(i,1)
    y1 = read1_pos(i,2)
    z1 = read1_pos(i,3)
    flag = 0
    do j = 1 , read2_num
        x2 = read2_pos(j,1)
        y2 = read2_pos(j,2)
        z2 = read2_pos(j,3)
        if(x1==x2.and.y1==y2.and.z1==z2)then !delete
            flag = 1
            exit
        endif
    enddo
    if(flag==0)then !save
        c_num = c_num + 1
        c_id(c_num) = read1_id(i)
        c_type(c_num) = read1_type(c_num)
        c_pos(c_num,1:3) = read1_pos(i,1:3)
    endif
enddo
write(*,*)"-->The num of atoms subtracted:",read1_num-c_num
write(*,*)"-->The num of atoms after subtracted:",c_num
!-----OUTPUT-1:DATA-format
open(13,file='c.data.output',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)c_num , "atoms"
write(13,*)"2 atom types"
write(13,*) xlo1,xhi1," xlo xhi"
write(13,*) ylo1,yhi1," ylo yhi"
write(13,*) zlo1,zhi1," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , c_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') c_id(i), c_type(i) , c_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DUMP-format
open(13,file='c.dump.output',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) c_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo1,xhi1
write(13,*) ylo1,yhi1
write(13,*) zlo1,zhi1
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , c_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') c_id(i), c_type(i) , c_pos(i,1:3)
enddo
close(13)
end program main        
