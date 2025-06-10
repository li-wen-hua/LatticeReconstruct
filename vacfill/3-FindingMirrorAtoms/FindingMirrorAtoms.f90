! Lattice constant;
! BCC;
!
!This program is used to fill the vacancies(step-1)
!
!Input1 : wrapped.dump.output-1 :: DUMP-format,from lammps. 
!
!Input2 : awaydislocations.data.output :: DATA-format, from ovito
!
!Output : 
!
!lwh
!2023.4.20
!----------------------------------------------------------------------

program main
use omp_lib
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!-----varis:MODULE
!!!!integer :: 
!!!!double precision :: 
!!!!integer,allocatable,dimension(:) :: 
!!!!double precision,allocatable,dimension(:,:) :: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum,nthreads
double precision :: g_alat,g_2,g_3,g_5

!-----varis:GLOBAL
integer :: i,j,k,l,m,n
double precision :: alat,alo,ahi,distmin
!-----varis:PART-1,read in dump.0
integer :: read1_num
double precision :: xlo1,xhi1,ylo1,yhi1,zlo1,zhi1
integer,allocatable,dimension(:) :: read1_id,read1_type
double precision,allocatable,dimension(:,:) :: read1_pos
!-----varis:PART-2,read in awaydislocations.output
integer :: read2_num,typenum2
double precision :: xlo2,xhi2,ylo2,yhi2,zlo2,zhi2
integer,allocatable,dimension(:) :: read2_id,read2_type
double precision,allocatable,dimension(:,:) :: read2_pos
!-----varis:PART-3,finding mirror atoms
integer :: mr_num,flag_rep
double precision :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,dist1,dist2
integer,allocatable,dimension(:) :: mr_id,mr_type
double precision,allocatable,dimension(:,:) :: mr_pos


!----------INI
open(1,file= '../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_alat
read(1,*)g_typenum
read(1,*) ! near/away from dis-lines
read(1,*)g_2 ! associated atoms, nearest
read(1,*)g_3 ! associated atoms, farest
read(1,*) ! shell thickness
read(1,*)g_5 ! if it is a vacancy 
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)nthreads

close(1)




alat = g_alat
alo = alat * g_2 !nearest 0.707/a range which determines a pair of atoms is associated like CNA
ahi = alat * g_3
distmin = alat * g_5 !approprate smaller. if it is a vacancy
!----------PART-1(READ-IN-1):read in dump.0 , lammps-DUMP-format
open(1,file="input_mirror.dump",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)read1_num !number-of-atoms
read(1,*)
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
write(*,*)"-->The num of atoms in read1_pos(wrapped.dump.output):",read1_num
!----------PART-2(READ-IN-2):read in type0.data , lammps-data-format
open(1,file="input_atomsself.data",status="old")
read(1,*)
read(1,*)read2_num !number-of-atoms
read(1,*)typenum2
read(1,*)xlo2,xhi2
read(1,*)ylo2,yhi2
read(1,*)zlo2,zhi2
allocate(read2_pos(1:read2_num,1:3))
allocate(read2_id(1:read2_num))
allocate(read2_type(1:read2_num))
read(1,*)
read(1,*)
read(1,*)
do i = 1 , read2_num
    read(1,*)read2_id(i),read2_type(i),read2_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in read2_pos(type0.data):",read2_num
!----------PART-3():finding mirror atoms
mr_num = 0 !"mr" = "mirror"
allocate(mr_id(1:read2_num*8))
allocate(mr_type(1:read2_num*8))
allocate(mr_pos(1:read2_num*8,1:3))


write(*,*)"!!!! Start of omp!"
call omp_set_num_threads(nthreads)
!$OMP PARALLEL PRIVATE(i,j,k,x1,y1,z1,x2,y2,z2,dist1,x0,y0,z0,flag_rep,x3,y3,z3,dist2) SHARED(mr_num,mr_id,mr_type,mr_pos)
!$OMP DO


do i = 1 , read2_num !loop atoms-away-from-dislines(type==0)
    x2 = read2_pos(i,1)
    y2 = read2_pos(i,2)
    z2 = read2_pos(i,3)
    do j = 1 , read1_num !loop wrapped-box
        x1 = read1_pos(j,1)
        y1 = read1_pos(j,2)
        z1 = read1_pos(j,3)
        dist1 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if(dist1>alo.and.dist1<=ahi)then
            x0 = 2*x2-x1 !symmetry atom
            y0 = 2*y2-y1 !symmetry atom
            z0 = 2*z2-z1 !symmetry atom
            flag_rep = 0 !if this atom is a vacancy or repeated
            do k = 1 , read1_num !loop wrapped-box
                x3 = read1_pos(k,1)
                y3 = read1_pos(k,2)
                z3 = read1_pos(k,3)
                dist2 = sqrt((x0-x3)**2+(y0-y3)**2+(z0-z3)**2)
                if(dist2<distmin)then !this atom is repeated
                    flag_rep = 1
                    exit
                endif
            enddo
            if(flag_rep==0)then !this position is a vacancy
                    !$OMP CRITICAL
                    mr_num = mr_num + 1
                    mr_id(mr_num) = mr_num
                    mr_type(mr_num) = 2
                    mr_pos(mr_num,1) = x0
                    mr_pos(mr_num,2) = y0
                    mr_pos(mr_num,3) = z0
                    !$OMP END CRITICAL
            endif
        endif
    enddo
enddo



!$OMP END DO
!$OMP END PARALlEL
write(*,*)"!!!! End of omp!"



write(*,*)"-->The num of mirror atoms in mr-matrix:",mr_num
!----------PART-4():output mirror atoms
!-----OUTPUT-1:DATA-format
open(13,file='output_mirroratoms.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)mr_num , "atoms"
write(13,*)"2 atom types"
write(13,*) xlo1,xhi1," xlo xhi"
write(13,*) ylo1,yhi1," ylo yhi"
write(13,*) zlo1,zhi1," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , mr_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') mr_id(i), mr_type(i) , mr_pos(i,1:3)
enddo
close(13)
!-----OUTPUT-1:DATA-format
open(13,file='output_mirroratoms.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) mr_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo1,xhi1
write(13,*) ylo1,yhi1
write(13,*) zlo1,zhi1
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , mr_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') mr_id(i), mr_type(i) , mr_pos(i,1:3)
enddo
close(13)

end program main
