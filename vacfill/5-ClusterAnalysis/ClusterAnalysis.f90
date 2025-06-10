! Lattice Constant;
!
!This program is a cna operation to classify clusters
!
!Note : Taking the periodic condition into our consideration.
!
!Input1 : mirroratoms.data.output/atoms.data.input :: Data fileformat,id-type-x-y-z.
! Containing atoms we are about to classify them into different clusters.
!
!Input2 : wrapped.data.output(*for-convenient*)(wrap atoms.input). Data format file.
!
!Output : cluster.dump.output :: Dump file format.id-type-x-y-z(for next step).
!
!lwh
!2023.4.18
!----------------------------------------------------------------------

program main
implicit none

!-----------state variables
!ADDED in 2024.4.27
integer :: g_typenum
double precision :: g_alat,g_2,g_3,g_6

!-----------varis:globle
integer :: i,j,k,l,m,n,distmin
double precision :: alat,a_min,a_max

!-----------varis:readin 
integer :: n0,typenum0
double precision :: xlo0,xhi0,ylo0,yhi0,zlo0,zhi0
double precision,allocatable,dimension(:,:) :: readin_pos
integer,allocatable,dimension(:) :: readin_id,readin_type

!-----------varis:mirror/mr-matrix 
integer :: wrap_num
integer,allocatable,dimension(:) :: wrap_id,wrap_type
double precision,allocatable,dimension(:,:) :: wrap_pos


!-----------varis:cn-finding
double precision :: x1,y1,z1,x2,y2,z2,dist
integer,allocatable,dimension(:) :: cn
integer,allocatable,dimension(:,:) :: cnid

!----------varis:classify clusters
integer :: cluster_num,id1,id2,id3,flag
integer,allocatable,dimension(:) :: cluster_atom_num
integer,allocatable,dimension(:,:) :: cluster_atom_id,cluster_atom_order
integer :: asso_num
integer,allocatable,dimension(:) :: asso_cluster
!----------varis:clear empty cluster
integer :: non_empty_cluster
!----------varis:output
integer :: od


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
read(1,*) ! shell thickness
read(1,*) ! if it is a vacancy 
read(1,*)g_6 ! if it is a vacancy 

close(1)




alat = g_alat
a_min = alat * g_2 !nearest 0.707/a range which determines a pair of atoms is associated like CNA
a_max = alat * g_3

distmin = g_6 !cluster distance cutoff
!-----------END ini set
!----------readin inp -- lammps data file format
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
write(*,*)"-->The num of atoms in input:",n0
!----------END of readin inp -- lammps data file format

!-----------read in wrapped.dat,which wrapped a shell on inp to
!-----------satisfy periodic boundary
open(1,file="input_wrapped.data",status="old")
read(1,*)
read(1,*)wrap_num !number-of-atoms
read(1,*) !wrap_type_num !number-of-atom-types
read(1,*) !xlo0,xhi0
read(1,*) !ylo0,yhi0
read(1,*) !zlo0,zhi0
read(1,*)
read(1,*)
read(1,*)
allocate(wrap_pos(1:wrap_num,1:3))
allocate(wrap_id(1:wrap_num))
allocate(wrap_type(1:wrap_num))
do i = 1 , wrap_num
    read(1,*)wrap_id(i),wrap_type(i),wrap_pos(i,1:3)
enddo
close(1)
write(*,*)"-->The num of atoms in input_wrapped:",wrap_num
!-----------END of reading wrapped.dat

!-----------find the "id" of common neighbor atom list of each atom

allocate(cn(1:n0)) !"cn" = the number of "common neighbor" of i-th atom
allocate(cnid(1:n0,1:wrap_num)) !"cnid" = the "id" of "common neighbor"
cn(1:n0) = 0
cnid(1:n0,1:wrap_num) = 0
do i = 1 , n0 !i-th atom in readin_pos
        x1 = readin_pos(i,1)
        y1 = readin_pos(i,2)
        z1 = readin_pos(i,3)
        do j = 1 , wrap_num !j-th atom in wrap_pos
                if(i==j)cycle
                x2 = wrap_pos(j,1)
                y2 = wrap_pos(j,2)
                z2 = wrap_pos(j,3)
                dist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                !if(dist>=a_min.and.dist<a_max)then !used for cna
                if(dist<distmin)then !used for cluster
                        cn(i) = cn(i) + 1 !how many cn-atoms around i-th atom
                        cnid(i,cn(i)) = wrap_id(j) !ids of those cn-atoms
                endif
        enddo
enddo
write(*,*)"-->Cluster Analysis is done!"
!-----------END find the common neighbor atom list of each atom

!-----------classify atoms into clusters
asso_num = 0
allocate(asso_cluster(1:n0))
asso_cluster(1:n0) = 0

allocate(cluster_atom_num(1:n0))
allocate(cluster_atom_id(1:n0,1:wrap_num))
allocate(cluster_atom_order(1:n0,1:wrap_num))
cluster_num = 0
cluster_atom_num(1:n0) = 0
cluster_atom_id(1:n0,1:wrap_num) = 0
cluster_atom_order(1,n0) = 0

!cluster_num = 1
!cluster_atom_num(cluster_num) = 1 !initialization
!cluster_atom_id(cluster_num,cluster_atom_num(cluster_num)) = cnid(1,1) 
!cluster_atom_order(1,1) = 1

do i = 1 , n0 !loop readin-matrix or cnid-matrix
        id1 = readin_id(i) !atom1-id
        asso_num = 0
        asso_cluster(1:n0) = 0
        do j = 1 , cluster_num !loop cluster_atom_id matrix
                flag = 0 !if the two atoms are in the same cluster
                do k = 1 , cluster_atom_num(j)
                        id2 = cluster_atom_id(j,k) !atom2-id / k-th atom in j-th cluster
                        !now we see if the two atoms are in the same cluster,i.e. atom1
                        !is a cn atom of atom2
                        do l = 1 , cn(i) !loop i-th atom's cna-list ids
                                id3 = cnid(i,l) !atom1's cn atom ids
                                if(id3==id2) then
                                        flag = 1
                                        exit
                                endif
                        enddo
                        if(flag==1)then
                                exit
                        endif
                enddo
                if (flag==1)then !same cluster
                        asso_num = asso_num + 1 !the number of cluster owning atom1
                        asso_cluster(asso_num) = j !id of cluster
                endif
        enddo
        if(asso_num==1)then !only one cluster is associated with atom1
                cluster_num = cluster_num !don't change
                cluster_atom_num(asso_cluster(asso_num)) = cluster_atom_num(asso_cluster(asso_num)) + 1
                cluster_atom_id(asso_cluster(asso_num),cluster_atom_num(asso_cluster(asso_num))) = id1
                cluster_atom_order(asso_cluster(asso_num),cluster_atom_num(asso_cluster(asso_num))) = i
        endif
        if(asso_num==0)then !atom1 belongs to no cluster,creadte a new one
                cluster_num = cluster_num + 1 !add 1
                cluster_atom_num(cluster_num) = 1
                cluster_atom_id(cluster_num,1) = id1
                cluster_atom_order(cluster_num,1) = i
        endif
        if(asso_num>1)then !combine clusters
                cluster_num = cluster_num !don't change
                cluster_atom_num(asso_cluster(1)) = cluster_atom_num(asso_cluster(1)) + 1
                cluster_atom_id(asso_cluster(1),cluster_atom_num(asso_cluster(1))) = id1
                cluster_atom_order(asso_cluster(1),cluster_atom_num(asso_cluster(1))) = i
                do n = 2 , asso_num
                        do m = 1 , cluster_atom_num(asso_cluster(n))
                                cluster_atom_num(asso_cluster(1)) =&
cluster_atom_num(asso_cluster(1)) + 1
                                cluster_atom_id(asso_cluster(1),cluster_atom_num(asso_cluster(1)))&
= cluster_atom_id(asso_cluster(n),m)
                                cluster_atom_order(asso_cluster(1),cluster_atom_num(asso_cluster(1)))&
= cluster_atom_order(asso_cluster(n),m)
                        enddo
                        cluster_atom_num(asso_cluster(n)) = 0
                enddo
        endif
enddo
write(*,*)"-->Initially , the number of cluster is :",cluster_num
!-----------END classify atoms into clusters

!-----------clear empty clusters
non_empty_cluster = 0
do i = 1 , cluster_num
        if(cluster_atom_num(i)/=0) non_empty_cluster=non_empty_cluster + 1
enddo
write(*,*)"-->***Finally , the number of non-empty-cluster is :",non_empty_cluster

!-----------END clear empty clusters
!-----------print max cluster size
write(*,*)"-->***Max cluster:",maxval(cluster_atom_num(1:cluster_num))
!write(*,*)"-->Clusters:"
!write(*,*)"-->",cluster_atom_num(1:cluster_num)
write(*,*)"Cluster Analysis is done!:"
!-----------END printing

!---------------TEST:output mirror-matix atoms-DUMP file format
open(13,file='output_cluster.dump',status='replace')
write(13,'(A)') "ITEM: TIMESTEP"
write(13,'(A)') "0"
write(13,'(A)') "ITEM: NUMBER OF ATOMS"
write(13,*) n0
write(13,'(A)') "ITEM: BOX BOUNDS pp pp pp"
write(13,*) xlo0,xhi0
write(13,*) ylo0,yhi0
write(13,*) zlo0,zhi0
write(13,'(A)') "ITEM: ATOMS id type x y z Cluster Order num-of-NAs"
do i = 1 , cluster_num
        do j = 1 , cluster_atom_num(i)
                od = cluster_atom_order(i,j)
                write(13,*) readin_id(od), readin_type(od) , readin_pos(od,1:3),i,od,cn(od)
        enddo
enddo
close(13)

end program main
