!This program is a part of csp operation
!
!Input : dump.0
!Input : type0.data
!
!Output : good.data/dump
!Output : bad.data/dump
!
!lwh
!2023.5.12
!----------------------------------------------------------------------

program main
use omp_lib
implicit none

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum,nthreads
double precision :: g_alat,g_7,g_8,g_9

!-----------varis:globle
integer :: i,j,k,readin1_flag,readin2_flag
double precision :: alat,alat_nearest,alat_2nd_nearest,alo,ahi,ahi2
!-----------varis:PART-1,read in dump.0
integer :: read1_num,read1_typenum
double precision :: xlo1,xhi1,ylo1,yhi1,zlo1,zhi1
integer,allocatable,dimension(:) :: read1_id,read1_type
double precision,allocatable,dimension(:,:) :: read1_pos
!-----------varis:PART-2,read in type0.data
integer :: read2_num,read2_typenum
double precision :: xlo2,xhi2,ylo2,yhi2,zlo2,zhi2
integer,allocatable,dimension(:) :: read2_id,read2_type
double precision,allocatable,dimension(:,:) :: read2_pos
!-----------varis:Part-3,select
double precision :: x1,y1,z1,x2,y2,z2,dist,flag
!-----------varis:PART-4,output matrix
integer :: away_num,near_num
integer,allocatable,dimension(:) :: away_id,away_type,near_id,near_type
double precision,allocatable,dimension(:,:) :: away_pos,near_pos




!----------PART-0:GLOBAL INIT
open(1,file= '../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_alat    ! 1 
read(1,*)g_typenum ! 2 
read(1,*)          ! 3
read(1,*)          ! 4  
read(1,*)          ! 5
read(1,*)          ! 6   
read(1,*)          ! 7
read(1,*)          ! 8
read(1,*)g_7       ! 9 nearest
read(1,*)g_8       ! 10 2-ed nearest
read(1,*)g_9       ! 11 cutoff when selecting atoms near type0
read(1,*)          ! 12
read(1,*)          ! 13
read(1,*)nthreads  ! 14-th ! 12-th para

close(1)



alat = g_alat
!alat_nearest = sqrt(3.0)*alat*0.5
alat_nearest = alat*g_7
alat_2nd_nearest = alat*g_8
alo = alat*0.90
ahi = alat*1.10
ahi2 = alat*g_9 ! cutoff of selecting atoms near type0

readin1_flag = 2 !1 for data,2 for dump
readin2_flag = 1
!----------PART-1(READ-IN-1):read in dump.0 , lammps-DUMP-format
if(readin1_flag==1)then
open(1,file="dump.0",status="old")
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
if(readin1_flag==2)then
open(1,file="dump.0",status="old")
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
!----------PART-1(READ-IN-1):read in dump.0 , lammps-DUMP-format
if(readin2_flag==1)then
open(1,file="wrapped.type0.data",status="old")
read(1,*)
read(1,*)read2_num !number-of-atoms
read(1,*)read2_typenum
read(1,*)xlo2,xhi2
read(1,*)ylo2,yhi2
read(1,*)zlo2,zhi2
read(1,*)
read(1,*)
read(1,*)
allocate(read2_pos(1:read2_num,1:3))
allocate(read2_id(1:read2_num))
allocate(read2_type(1:read2_num))
do i = 1 , read2_num
    read(1,*)read2_id(i),read2_type(i),read2_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read1_pos:",read2_num
endif
if(readin2_flag==2)then
open(1,file="type0.data",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)read2_num !number-of-atoms
read(1,*)
read2_typenum = 2
read(1,*)xlo2,xhi2
read(1,*)ylo2,yhi2
read(1,*)zlo2,zhi2
read(1,*)
allocate(read2_pos(1:read2_num,1:3))
allocate(read2_id(1:read2_num))
allocate(read2_type(1:read2_num))
do i = 1 , read2_num
    read(1,*)read2_id(i),read2_type(i),read2_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read1_pos:",read2_num
endif
!----------PART-2():select
away_num = 0
allocate(away_pos(1:read1_num,1:3))
allocate(away_id(1:read1_num))
allocate(away_type(1:read1_num))
near_num = 0
allocate(near_pos(1:read1_num,1:3))
allocate(near_id(1:read1_num))
allocate(near_type(1:read1_num))




write(*,*)"!!!! Start of opm!"
call omp_set_num_threads(nthreads)
!$OMP PARALLEL PRIVATE(i,j,x1,y1,z1,x2,y2,z2,flag,dist)&
!$OMP& SHARED(away_num,away_pos,away_id,away_type,near_num,near_type,near_id,near_pos)
!$OMP DO


do i = 1 , read1_num !loop dump.0
        x1 = read1_pos(i,1)
        y1 = read1_pos(i,2)
        z1 = read1_pos(i,3)
        flag = 0 !0 for 
        do j = 1 , read2_num !type0.data
                x2 = read2_pos(j,1)
                y2 = read2_pos(j,2)
                z2 = read2_pos(j,3)
                dist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                if(dist<ahi2)then
                        flag = 1 !near to type0
                        exit
                endif
        enddo
        if(flag==0)then !away from type0
                !$OMP CRITICAL
                away_num = away_num + 1
                away_pos(away_num,1:3) = read1_pos(i,1:3)
                away_id(away_num) = read1_id(i)
                away_type(away_num) = read1_type(i)
                !$OMP END CRITICAL
        elseif(flag==1)then
                !$OMP CRITICAL
                near_num = near_num + 1
                near_pos(near_num,1:3) = read1_pos(i,1:3)
                near_id(near_num) = read1_id(i)
                near_type(near_num) = read1_type(i)
                !$OMP END CRITICAL
        endif
enddo


!$OMP END DO
!$OMP END PARALlEL
write(*,*)"!!!! End of opm!"




write(*,*)"-->The num of atoms away from type0:",away_num
write(*,*)"-->The num of atoms near to type0:",near_num


!-----OUTPUT-1:DATA-format
open(13,file='away.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)away_num , "atoms"
write(13,*)read2_typenum," atom types"
write(13,*) xlo2,xhi2," xlo xhi"
write(13,*) ylo2,yhi2," ylo yhi"
write(13,*) zlo2,zhi2," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , away_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') away_id(i), away_type(i) , away_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DUMP-format
open(13,file='away.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) away_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo2,xhi2
write(13,*) ylo2,yhi2
write(13,*) zlo2,zhi2
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , away_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') away_id(i), away_type(i) , away_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DATA-format
open(13,file='near.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)near_num , "atoms"
write(13,*)read2_typenum," atom types"
write(13,*) xlo2,xhi2," xlo xhi"
write(13,*) ylo2,yhi2," ylo yhi"
write(13,*) zlo2,zhi2," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , near_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') near_id(i), near_type(i) , near_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DUMP-format
open(13,file='near.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) near_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo2,xhi2
write(13,*) ylo2,yhi2
write(13,*) zlo2,zhi2
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , near_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') near_id(i), near_type(i) , near_pos(i,1:3)
enddo
close(13)



end program main
