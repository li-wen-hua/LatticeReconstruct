! Lattice constant;
!
!This program is used to create a "shell" of a simulation box , 
!so it can satisfy the periodic condition.
!
!Input : atoms.input :: Data fileformat,id-type-x-y-z. Contains 
!atoms we are about to classify them into different clusters.
!
!Output : wrapped.output :: Data file format.id-type-x-y-z.
!
!lwh
!2023.4.17
!----------------------------------------------------------------------

program main
implicit none

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum
double precision :: g_alat,g_2,g_3,g_4

!-----------varis:globle
integer :: i,j,k,l,m,n,readin_flag
double precision :: alat,a_min,a_max,ahi

!-----------varis:readin 
integer :: n0,typenum0
double precision :: xlo0,xhi0,ylo0,yhi0,zlo0,zhi0
double precision,allocatable,dimension(:,:) :: readin_pos
integer,allocatable,dimension(:) :: readin_id,readin_type

!-----------varis:mirror/mr-matrix 
integer :: mr_num
integer,allocatable,dimension(:) :: mr_id,mr_type
double precision,allocatable,dimension(:,:) :: mr_pos

!-----------varis:find mirror atoms
integer :: flag,flag1,flag2,flag3,flag4,flag5,flag6
double precision :: dx,dy,dz,x1,y1,z1,d_xlo,d_xhi,d_ylo,d_yhi,d_zlo,d_zhi

!-----------varis:warning
integer :: warning_num
!-----------ini set
open(1,file= '../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_alat
read(1,*)g_typenum
read(1,*) ! near/away from dis-lines
read(1,*)g_2 ! associated atoms, nearest
read(1,*)g_3 ! associated atoms, farest
read(1,*)g_4 ! shell thickness
read(1,*) ! if it is a vacancy 

close(1)



readin_flag = 2 !1-->data-format,2-->dump-format
alat = g_alat
a_min = alat * g_2 !nearest 0.707/a range which determines a pair of atoms is associated like CNA
a_max = alat * g_3
ahi = a_max * g_4 !the mirror operation criterion ##### i.e. THE SHELL THICKNESS #####

!-----------END ini set
!----------DATA FILE:readin inp -- lammps data file format
if(readin_flag==1)then
open(1,file="input.data",status="old")
read(1,*)
read(1,*)n0 !number-of-atoms
read(1,*)typenum0 !number-of-atom-types
read(1,*)xlo0,xhi0
read(1,*)ylo0,yhi0
read(1,*)zlo0,zhi0
read(1,*)
read(1,*)
read(1,*)
allocate(readin_pos(1:n0,1:3))
allocate(readin_id(1:n0))
allocate(readin_type(1:n0))
do i = 1 , n0
    read(1,*)readin_id(i),readin_type(i),readin_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in readin_pos(inp):",n0
endif
!----------END of readin inp -- lammps data file format
!----------DUMP FILE:readin inp -- lammps data file format
if(readin_flag==2)then
open(1,file="input.dump",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)n0 !number-of-atoms
read(1,*)
typenum0 = 1 !default 1
read(1,*)xlo0,xhi0
read(1,*)ylo0,yhi0
read(1,*)zlo0,zhi0
read(1,*)
allocate(readin_pos(1:n0,1:3))
allocate(readin_id(1:n0))
allocate(readin_type(1:n0))
do i = 1 , n0
    read(1,*)readin_id(i),readin_type(i),readin_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in readin_pos(inp):",n0
endif
!----------END of readin inp -- lammps data file format

!----------create a "shell" to accomadate the periodic condition
warning_num = 0
mr_num = 0
allocate(mr_id(1:3*n0))
allocate(mr_type(1:3*n0))
allocate(mr_pos(1:3*n0,1:3))

dx = abs(xlo0-xhi0)
dy = abs(ylo0-yhi0)
dz = abs(zlo0-zhi0)
do i = 1 , n0 !i-th atom in readin-matrix
        flag = 0 !if the atom position is changed
        flag1 = 0
        flag2 = 0
        flag3 = 0
        flag4 = 0
        flag5 = 0
        flag6 = 0
        x1 = readin_pos(i,1)
        y1 = readin_pos(i,2)
        z1 = readin_pos(i,3)
        d_xlo = x1-xlo0 !distance_with_xlo
        d_xhi = xhi0-x1
        d_ylo = y1-ylo0
        d_yhi = yhi0-y1
        d_zlo = z1-zlo0
        d_zhi = zhi0-z1
        if(d_xlo<0.or.d_xhi<0.or.d_ylo<0.or.d_yhi<0.or.d_zlo<0.or.d_zhi<0)then
                warning_num = warning_num + 1
        endif
        if (d_xlo<=ahi) then !means the atom near this edge
                x1 = x1 + dx !change its position to its mirror position
                flag1 = flag1 + 1
                flag = flag + 1
        endif
        if (d_xhi<=ahi) then
                x1 = x1 - dx
                flag2 = flag2 + 1
                flag = flag + 1
        endif
        if (d_ylo<=ahi) then
                y1 = y1 + dy
                flag3 = flag3 + 1
                flag = flag + 1
        endif
        if (d_yhi<=ahi) then
                y1 = y1 - dy
                flag4 = flag4 + 1
                flag = flag + 1
        endif
        if (d_zlo<=ahi) then
                z1 = z1 + dz
                flag5 = flag5 + 1
                flag = flag + 1
        endif
        if (d_zhi<=ahi) then
                z1 = z1 - dz
                flag6 = flag6 + 1
                flag = flag + 1
        endif
        !Then we find mirror atoms and add them into mr-matrix
        if (flag1==1) then !means this atom's position changed/faceX1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif
        if (flag2==1) then !faceX2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif
        if (flag3==1) then !faceY1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif
        if (flag4==1) then !faceY2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif
        if (flag5==1) then !faceZ1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif
        if (flag6==1) then !faceZ2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif
        if (flag1==1.and.flag3==1) then !edgeXY1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif           
        if (flag1==1.and.flag4==1) then !edgeXY2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif           
        if (flag2==1.and.flag3==1) then !edgeXY3
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif           
        if (flag2==1.and.flag4==1) then !edgeXY4
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = readin_pos(i,3)
        endif           
        if (flag3==1.and.flag5==1) then !edgeYZ1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag3==1.and.flag6==1) then !edgeYZ2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag4==1.and.flag5==1) then !edgeYZ3
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag4==1.and.flag6==1) then !edgeYZ4
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = readin_pos(i,1)
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag5==1) then !edgeXZ1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag6==1) then !edgeXZ2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag5==1) then !edgeXZ3
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag6==1) then !edgeXZ4
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = readin_pos(i,2)
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag3==1.and.flag5==1) then !corner1
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag3==1.and.flag6==1) then !corner2
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag4==1.and.flag5==1) then !corner3
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag1==1.and.flag4==1.and.flag6==1) then !corner4
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag3==1.and.flag5==1) then !corner5
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag3==1.and.flag6==1) then !corner6
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag4==1.and.flag5==1) then !corner7
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
        if (flag2==1.and.flag4==1.and.flag6==1) then !corner8
                mr_num = mr_num + 1 !"mr"=="mirror"
                mr_id(mr_num) = readin_id(i)
                mr_type(mr_num) = readin_type(i)
                mr_pos(mr_num,1) = x1
                mr_pos(mr_num,2) = y1
                mr_pos(mr_num,3) = z1
        endif           
enddo
if(warning_num>0)then
        write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*)"!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*)"!!!!!!!!!!Atom Outside of The Box!!!!!!!!!!!!!!"
        write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*)"-->WARNING ATOMS: ",warning_num
endif
write(*,*)"-->The number of atoms in mirror-matrix(the shell) : ",mr_num

!----------END of create a "shell" to accomadate the periodic condition

!!---------------TEST:output mirror-matix atoms/The "shell"
!open(13,file='shell.test',status='replace')
!write(13,*)"Fill the vacancies , lammps data file , by fortran"
!write(13,*)mr_num , "atoms"
!write(13,*)typenum0,"atom types"
!write(13,*) xlo0,xhi0," xlo xhi"
!write(13,*) ylo0,yhi0," ylo yhi"
!write(13,*) zlo0,zhi0," zlo zhi"
!write(13,*)
!write(13,'(A)') "Atoms"
!write(13,*)
!do i = 1 , mr_num
!  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') mr_id(i), mr_type(i) , mr_pos(i,1:3)
!enddo
!close(13)
!write(*,*)"-->The number of atoms in shell.test : ",mr_num
!
!---------------DATA FORMAT:Output the box that wrapped a shell
open(13,file='output.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)mr_num+n0 , "atoms"
write(13,*)typenum0,"atom types"
write(13,*) xlo0,xhi0," xlo xhi"
write(13,*) ylo0,yhi0," ylo yhi"
write(13,*) zlo0,zhi0," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , n0
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') readin_id(i), readin_type(i) , readin_pos(i,1:3)
enddo
do i = n0+1 , n0+mr_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') mr_id(i-n0), mr_type(i-n0) , mr_pos(i-n0,1:3)
!NOTE : The id of atoms in mr-matrix must be the same with atoms in
!readin-matrix, i.e. its original atom.Becuase we will use the id in the next step.
enddo
close(13)
write(*,*)"-->The number of atoms in output : ",mr_num+n0
!---------------DUMP FORMAT:Output the box that wrapped a shell
open(13,file='output.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) mr_num+n0
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo0,xhi0
write(13,*) ylo0,yhi0
write(13,*) zlo0,zhi0
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , n0
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') readin_id(i), readin_type(i) , readin_pos(i,1:3)
enddo
do i = n0+1 , n0+mr_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') mr_id(i-n0), mr_type(i-n0) , mr_pos(i-n0,1:3)
!NOTE : The id of atoms in mr-matrix must be the same with atoms in
!readin-matrix, i.e. its original atom.Becuase we will use the id in the next step.
enddo
close(13)

write(*,*)"WRAP A SHELL IS DONE!"
end program main
