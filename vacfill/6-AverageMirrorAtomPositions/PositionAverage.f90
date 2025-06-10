! NONE.
!
!This program is used to average the repeated atoms(step-2)
!
!Input1  : cluster.dump.output :: cluster analysis, dump file
!format,id-type-x-y-z-clusters).Sort by size
!
!Note  : the cutoff radius shouls set "1.0" the cluster analysis, so
!that we can calculate the average positions. 
!
!Input2  : dump.0 :: ini box,lammps data format file,from lammps
!write_data,if remove the mass part
!
!Output : filledbox.output(final filled box matrix)
!
!lwh
!2023.3.31
!----------------------------------------------------------------------

program main
implicit none

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum

!-----------varis:globle
integer :: i,j,k,l,m,n

!-----------varis:readin repclu.dat 
integer :: rep_num
double precision :: xlo,xhi,ylo,yhi,zlo,zhi
integer,allocatable,dimension(:) :: rep_id,rep_type,rep_clu
double precision,allocatable,dimension(:,:) :: rep_pos

!----------varis:counting/seperate int_pos into clusters
integer :: clu_num
integer,allocatable,dimension(:) :: clu_tc

!---------varis:seperating/......
integer :: clu_size_max
integer,allocatable,dimension(:,:) :: clu_type
double precision,allocatable,dimension(:,:,:) :: clu_pos

!--------varis:average rep_pos-->add_pos
integer :: add_num
double precision :: possum_x,possum_y,possum_z
integer,allocatable,dimension(:) :: add_id,add_type
double precision,allocatable,dimension(:,:) :: add_pos

!--------varis:read-in cg.dat
integer :: n0,typenum0
double precision :: xlo0,xhi0,ylo0,yhi0,zlo0,zhi0
double precision,allocatable,dimension(:,:) :: readin_pos
integer,allocatable,dimension(:) :: readin_id,readin_type

!--------varis:periodic
integer :: p_count
double precision :: p_cut,lx,ly,lz,px,py,pz,x1,y1,z1,x2,y2,&
z2,dx,dy,dz
!-------varis:average
integer :: fillatom_type

!-------varis:output modi_clu.dat
integer :: count1,lammps_type



open(1,file= '../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_typenum

close(1)




!-----------read in the LAMMPS dump file,containing repeated clusters
open(12,file='input_cluster.dump',status='old')
do i = 1,3
  read(12,*)
enddo
read(12,*)rep_num
allocate(rep_id(1:rep_num))
allocate(rep_type(1:rep_num))
allocate(rep_pos(1:rep_num,1:3))
allocate(rep_clu(1:rep_num))
read(12,*)
read(12,*)xlo,xhi
read(12,*)ylo,yhi
read(12,*)zlo,zhi
read(12,*)
do i = 1,rep_num
  read(12,*)rep_id(i),rep_type(i),rep_pos(i,1:3),rep_clu(i)
enddo
close(12)
write(*,*)"-->Total number of atoms in repclu.dat:",rep_num

!----------seperate rep_pos into clusters/counting
clu_num = maxval(rep_clu)
write(*,*)"-->MAX number of clusters:",clu_num

allocate(clu_tc(1:clu_num)) !cluster_type&count
clu_tc = 0
do i = 1 , rep_num
        clu_tc(rep_clu(i)) = clu_tc(rep_clu(i)) + 1
enddo
!write(*,*)"-->The size of those clusters (NOT sorted by size):"
!write(*,*)clu_tc(1:clu_num)

!----------seperate int_pos into clusters/seperating
clu_size_max = maxval(clu_tc)
write(*,*)"-->The num of atoms of the largest cluster:",clu_size_max

allocate(clu_type(1:clu_num,1:clu_size_max))
allocate(clu_pos(1:clu_num,1:clu_size_max,1:3))

clu_tc = 0
do i = 1 , rep_num
                clu_tc(rep_clu(i)) = clu_tc(rep_clu(i)) + 1
                clu_type(rep_clu(i),clu_tc(rep_clu(i))) = rep_type(i)
                clu_pos(rep_clu(i),clu_tc(rep_clu(i)),1) = rep_pos(i,1)
                clu_pos(rep_clu(i),clu_tc(rep_clu(i)),2) = rep_pos(i,2)
                clu_pos(rep_clu(i),clu_tc(rep_clu(i)),3) = rep_pos(i,3)
enddo

!----------readin ini box:dump.0
open(1,file="input_added.dump",status="old")
read(1,*)
read(1,*)
read(1,*)
read(1,*)n0 !number-of-atoms
typenum0 = 2 !number-of-atom-types
read(1,*)
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
write(*,*)"-->The num of atoms in readin_pos(dump.0):",n0

!----------if the cluster is seperated by periodic boundary
p_cut = 0.8
lx = xhi-xlo
ly = yhi-ylo
lz = zhi-zlo
px = lx*p_cut !indicator of periodic
py = ly*p_cut
pz = lz*p_cut
p_count = 0 !count the number of clusters which seperated by boundary
do i = 1 , clu_num !i-th cluster
    if(clu_tc(i)==0)cycle !take size-0 cluster out
    do j = 1 , clu_tc(i)-1 !j-th atom of i-th cluster
        x1 = clu_pos(i,j,1)
        y1 = clu_pos(i,j,2)
        z1 = clu_pos(i,j,3)
        do k = j+1 , clu_tc(i) !loop radin_pos
            x2 = clu_pos(i,k,1)
            y2 = clu_pos(i,k,2)
            z2 = clu_pos(i,k,3)
            dx = abs(x1-x2)
            dy = abs(y1-y2)
            dz = abs(z1-z2)
            if (dx>px)then
                p_count = p_count + 1
                if(x1<0)then
                    x1 = x1 + lx
                    clu_pos(i,j,1) = x1
                endif
                if(x2<0)then
                    x2 = x2 + lx
                    clu_pos(i,k,1) = x2
                endif
            endif
            if (dy>py)then
                p_count = p_count + 1
                if(y1<0)then
                    y1 = y1 + ly
                    clu_pos(i,j,2) = y1
                endif
                if(y2<0)then
                    y2 = y2 + ly
                    clu_pos(i,k,2) = y2
                endif
            endif
            if (dz>pz)then
                p_count = p_count + 1
                if(z1<0)then
                    z1 = z1 + lz
                    clu_pos(i,j,3) = z1
                endif
                if(z2<0)then
                    z2 = z2 + lz
                    clu_pos(i,k,3) = z2
                endif
            endif
        enddo
    enddo
enddo
write(*,*)"-->The number of cluster atoms which seperated by periodic boundaries:",p_count

!----------average positions/add_pos matrix
fillatom_type = g_typenum   !#### SHOW filled atoms DIFFERENT COLOR IN OVITO ####
lammps_type = 1
add_num = 0
allocate(add_id(1:n0))
allocate(add_type(1:n0))
allocate(add_pos(1:n0,1:3))
do i = 1 , clu_num !i-th cluster
    if(clu_tc(i)==0)cycle !take size-0 cluster out
    possum_x = 0
    possum_y = 0
    possum_z = 0
    do j = 1 , clu_tc(i) !j-th atom of i-th cluster
        possum_x = possum_x + clu_pos(i,j,1)
        possum_y = possum_y + clu_pos(i,j,2)
        possum_z = possum_z + clu_pos(i,j,3)
    enddo
    add_num = add_num + 1
    add_id(add_num) = n0 + add_num
    add_type(add_num) = fillatom_type !set add atom type, 1 for simulation, 2 for test 
    add_pos(add_num,1) = possum_x/clu_tc(i)
    add_pos(add_num,2) = possum_y/clu_tc(i)
    add_pos(add_num,3) = possum_z/clu_tc(i)
enddo
write(*,*)"-->The num of atoms in add_pos:",add_num
!---------------confine atoms inner box
do i = 1 , add_num
        if(add_pos(i,1)<xlo0)add_pos(i,1)=add_pos(i,1)+abs(xlo0-xhi0)
        if(add_pos(i,1)>xhi0)add_pos(i,1)=add_pos(i,1)-abs(xlo0-xhi0)
        if(add_pos(i,2)<ylo0)add_pos(i,2)=add_pos(i,2)+abs(ylo0-yhi0)
        if(add_pos(i,2)>yhi0)add_pos(i,2)=add_pos(i,2)-abs(ylo0-yhi0)
        if(add_pos(i,3)<zlo0)add_pos(i,3)=add_pos(i,3)+abs(zlo0-zhi0)
        if(add_pos(i,3)>zhi0)add_pos(i,3)=add_pos(i,3)-abs(zlo0-zhi0)
enddo
write(*,*)"-->@@@@The num of atoms in add_pos(confined):",add_num

!!---------------output modified repclu.dat,i.e.clu_pos-->clu_pos
!count1 = 0 !count how many atom in this matrix
!do i = 1 , clu_num
!        if(clu_tc(i)==0)cycle !take size-0 cluster out
!        count1 = count1 + clu_tc(i)
!enddo
!
!open(13,file='modi_clu.dat',status='replace')
!write(13,*)"Fill the vacancies , lammps data file , by fortran"
!write(13,*) count1 , "atoms"
!write(13,*)"2 atom types"
!write(13,*) xlo,xhi," xlo xhi"
!write(13,*) ylo,yhi," ylo yhi"
!write(13,*) zlo,zhi," zlo zhi"
!write(13,*)
!write(13,'(A)') "Atoms"
!write(13,*)
!count1 = 0
!do i = 1 , clu_num
!        if(clu_tc(i)==0)cycle !take size-0 cluster out
!        do j = 1 , clu_tc(i)
!                count1 = count1 + 1
!                write(13,'(I8,2x,I2,2x,3(f19.8,1x))') count1, clu_type(i,j) ,clu_pos(i,j,1:3)
!        enddo
!enddo
!close(13)
!
!!---------------output add_pos
!open(13,file='fillatoms.test',status='replace')
!write(13,*)"Fill the vacancies , lammps data file , by fortran"
!write(13,*)add_num , "atoms"
!write(13,*)"2 atom types"
!write(13,*) xlo,xhi," xlo xhi"
!write(13,*) ylo,yhi," ylo yhi"
!write(13,*) zlo,zhi," zlo zhi"
!write(13,*)
!write(13,'(A)') "Atoms"
!write(13,*)
!do i = 1 , add_num
!  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') add_id(i), add_type(i) , add_pos(i,1:3)
!enddo
!close(13)
!write(*,*)"-->The num of atoms in fillatoms.dat file:",add_num
!

!---------------output filled box, finally
open(13,file='filledbox.data',status='replace')
write(13,*)"Fill the vacancies , lammps data file , by fortran"
write(13,*)n0+add_num , "atoms"
write(13,*)g_typenum,"atom types"
write(13,*) xlo,xhi," xlo xhi"
write(13,*) ylo,yhi," ylo yhi"
write(13,*) zlo,zhi," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1 , n0
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') i, readin_type(i) , readin_pos(i,1:3)
enddo
do i = 1 , add_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') n0+i, add_type(i) , add_pos(i,1:3)
enddo
close(13)
write(*,*)"-->The num of atoms in filledbox.output file:",add_num+n0
!----------------OUTPUT-1:DUMP-format
open(13,file='filledbox.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) n0+add_num
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo,xhi
write(13,*) ylo,yhi
write(13,*) zlo,zhi
write(13,'(A)') "ITEM: ATOMS id type x y z"
do i = 1 , n0
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') i, readin_type(i) , readin_pos(i,1:3)
enddo
do i = 1 , add_num
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') n0+i, add_type(i) , add_pos(i,1:3)
enddo
close(13)
write(*,*)"-->The num of atoms in filledbox.dump file:",add_num+n0
!!---------------output filled box, finally-->LAMMPS file
!open(13,file='filledbox.lmp',status='replace')
!write(13,*)"Fill the vacancies , lammps data file , by fortran"
!write(13,*)n0+add_num , "atoms"
!write(13,*)g_typenum,"atom types"
!write(13,*) xlo,xhi," xlo xhi"
!write(13,*) ylo,yhi," ylo yhi"
!write(13,*) zlo,zhi," zlo zhi"
!write(13,*)
!write(13,'(A)') "Atoms"
!write(13,*)
!do i = 1 , n0
!  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') i, lammps_type , readin_pos(i,1:3)
!enddo
!do i = 1 , add_num
!  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') n0+i, lammps_type , add_pos(i,1:3)
!enddo
!close(13)
!write(*,*)"-->The num of atoms in filledbox.output file:",add_num+n0

write(*,*)"POSITION AVERAGE is DONE!"
end program main
