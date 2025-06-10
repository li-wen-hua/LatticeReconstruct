!This program is a part of csp operation
!
!Input : wrapped.data.output
!Input : near.data
!
!Output : good.data/dump
!Output : bad.data/dump
!
!lwh
!2023.5.15
!----------------------------------------------------------------------

program main
use omp_lib
implicit none

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum,nthreads
double precision :: g_alat,g_7,g_8,g_10,g_11


!-----------varis:globle
integer :: i,j,k,readin1_flag,readin2_flag
double precision :: alat,alat_nearest,alat_2nd_nearest,alo,ahi
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
integer :: flag
double precision :: x0,y0,z0,x1,y1,z1,x2,y2,z2,dist,para1,para2
double precision :: dx,dy,dz,flag_x,flag_y,flag_z,flag_mirror
double precision,allocatable,dimension(:) :: read2_para1,read2_para2
!-----------varis:PART-4,output matrix
integer :: good_num,bad_num
integer,allocatable,dimension(:) :: good_id,good_type,bad_id,bad_type
double precision,allocatable,dimension(:) :: good_para,bad_para
double precision,allocatable,dimension(:,:) :: good_pos,bad_pos
!-----------varis:PART-5,na
integer,allocatable,dimension(:) :: na_num
double precision,allocatable,dimension(:,:,:) :: na_list




!----------PART-0:GLOBAL INIT
open(1,file= '../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_alat
read(1,*)g_typenum
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_7 ! nearest
read(1,*)g_8 ! 2-ed nearest
read(1,*) ! cutoff when selecting atoms near type0
read(1,*)g_10 ! nearest cutoff of selecting common atoms
read(1,*)g_11 ! farest cutoff of selecting common atoms
read(1,*)nthreads

close(1)





alat = g_alat
alat_nearest = alat * g_7
alat_2nd_nearest = alat * g_8
alo = alat * g_10
ahi = alat * g_11

readin1_flag = 1 !1 for data,2 for dump
readin2_flag = 1
!----------PART-1(READ-IN-1):read in dump.0 , lammps-DUMP-format
if(readin1_flag==1)then
open(1,file="wrapped.data.output",status="old")
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
open(1,file="near.data",status="old")
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
write(*,*)"-->The num of atoms in read2_pos:",read2_num
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
write(*,*)"-->The num of atoms in read2_pos:",read2_num
endif
!----------PART-2():ca-list
allocate(na_num(1:read2_num)) !how many ca in atom i which is na_num(i) 
allocate(na_list(1:read2_num,1:24,1:3)) !the coordinates of ca
na_num = 0
na_list = 0



write(*,*)"!!!! Start of omp!"
call omp_set_num_threads(nthreads)
!$OMP PARALLEL PRIVATE(i,j,x2,y2,z2,x1,y1,z1,dist)&
!$OMP& SHARED(na_num,na_list)
!$OMP DO



do i = 1 , read2_num !loop near atoms(type0 almost), slower
        ! write(*,*)'Neighbor atom:',i,read2_num
        x2 = read2_pos(i,1)
        y2 = read2_pos(i,2)
        z2 = read2_pos(i,3)
        do j = 1 , read1_num !dump0
                x1 = read1_pos(j,1)
                y1 = read1_pos(j,2)
                z1 = read1_pos(j,3)
                dist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                if(dist>alo.and.dist<ahi)then !nearest or 2nd_nearest atoms/na
                        !$OMP CRITICAL
                        na_num(i) = na_num(i) + 1
                        na_list(i,na_num(i),1:3) = read1_pos(j,1:3)
                        !$OMP END CRITICAL
                endif
        enddo
enddo
!write(*,*)"-->The num of na atoms(TEST):",na_num(1:4)




!$OMP END DO
!$OMP END PARALLEL
write(*,*)"!!!! End of omp!"





!----------PART-3():selection
good_num = 0
allocate(good_pos(1:read2_num,1:3))
allocate(good_id(1:read2_num))
allocate(good_type(1:read2_num))
bad_num = 0
allocate(bad_pos(1:read2_num,1:3))
allocate(bad_id(1:read2_num))
allocate(bad_type(1:read2_num))

allocate(read2_para1(1:read2_num))
allocate(read2_para2(1:read2_num))

do i = 1 , read2_num ! i-th atom, NEAR atom
        x0 = read2_pos(i,1)
        y0 = read2_pos(i,2)
        z0 = read2_pos(i,3)
        para1 = 0 !count how many no-matching atom(atoms with no mirror atom)
        para2 = 0 !count matched atoms
        do j = 1 , na_num(i) !j-th na
                x1 = na_list(i,j,1)
                y1 = na_list(i,j,2)
                z1 = na_list(i,j,3)
                flag = 0 !0 for do not find mirror atom
                if(x1>10000) cycle
                do k = 1 , na_num(i) !k-th na
                        x2 = na_list(i,k,1)
                        y2 = na_list(i,k,2)
                        z2 = na_list(i,k,3)
                        if(x2>10000) cycle
                        flag_x = abs(2*x0-x1-x2)
                        flag_y = abs(2*y0-y1-y2)
                        flag_z = abs(2*z0-z1-z2)
                        flag_mirror = flag_x+flag_y+flag_z
                        !write(*,*)'MONITOR:',flag_mirror
                        if(flag_mirror<=0.7)then !a pair of mirror atom
                                na_list(i,j,1:3)=(/10001.0,0.0,0.0/)
                                na_list(i,k,1:3)=(/10001.0,0.0,0.0/)
                                flag = 1 ! find mirror atom
                                exit
                        endif
                enddo
                if(flag==0)then !not find mirror atom
                        para1 = para1 + 1
                elseif(flag==1)then
                        para2 = para2 + 1
                endif
        enddo
        read2_para1(i) = para1
        read2_para2(i) = para2
        !write(*,*)'!!!',para1,para2
        if(para1==0.and.para2==7)then !all have mirror atom and perf crystal
                good_num = good_num + 1
                good_pos(good_num,1:3) = read2_pos(i,1:3)
                good_id(good_num) = read2_id(i)
                good_type(good_num) = read2_type(i)
        else
                bad_num = bad_num + 1
                bad_pos(bad_num,1:3) = read2_pos(i,1:3)
                bad_id(bad_num) = read2_id(i)
                bad_type(bad_num) = read2_type(i)
        endif
enddo

write(*,*)"-->The num of good/bad atoms(TEST):",good_num,bad_num
!-----OUTPUT-1:DATA-format
open(13,file='good.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)good_num , "atoms"
write(13,*)read2_typenum," atom types"
write(13,*) xlo2,xhi2," xlo xhi"
write(13,*) ylo2,yhi2," ylo yhi"
write(13,*) zlo2,zhi2," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , good_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') good_id(i), good_type(i) ,good_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DUMP-format
open(13,file='good.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) good_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo2,xhi2
write(13,*) ylo2,yhi2
write(13,*) zlo2,zhi2
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , good_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x),2x,f19.8)') good_id(i), good_type(i) ,good_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DATA-format
open(13,file='bad.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)bad_num , "atoms"
write(13,*)read2_typenum," atom types"
write(13,*) xlo2,xhi2," xlo xhi"
write(13,*) ylo2,yhi2," ylo yhi"
write(13,*) zlo2,zhi2," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , bad_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') bad_id(i), bad_type(i) , bad_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DUMP-format
open(13,file='bad.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) bad_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo2,xhi2
write(13,*) ylo2,yhi2
write(13,*) zlo2,zhi2
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , bad_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x),2x,f19.8)') bad_id(i), bad_type(i),bad_pos(i,1:3)
enddo
close(13)
!!-----OUTPUT-1:DUMP-format
!open(13,file='test.para.dump',status='replace')
!write(13,'(A)') "ITEM: TIMESTEP"
!write(13,'(A)') "0"
!write(13,'(A)') "ITEM: NUMBER OF ATOMS"
!write(13,*) read2_num
!write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
!write(13,*) xlo2,xhi2
!write(13,*) ylo2,yhi2
!write(13,*) zlo2,zhi2
!write(13,'(A)') "ITEM: ATOMS id type x y z para1 para2"
!do i = 1 , read2_num
!  write(13,'(I8,2x,I2,2x,3(f19.8,1x),2x,f19.8)') read2_id(i), read2_type(i),read2_pos(i,1:3),read2_para1(i),read2_para2(i)
!enddo
!close(13)


end program main
