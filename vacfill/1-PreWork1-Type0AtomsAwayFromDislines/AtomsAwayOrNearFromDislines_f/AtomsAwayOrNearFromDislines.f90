! None.
!
!This program aims to delete atoms around the dislocations. Inputs are the ovito
!CA format file contains vertexes of dis-lines, specified format see example
!file, and the file contains interstitials in the ovito dump file format. Then 
!the output files are the file contains atoms away from dis-lines.
!
!FILES:
!Input: type0.input, interstitials , ovito lammps dumpfile . Note that this file
!need DXA analysis and as well as the WS analysis first , especially the WS
!analysis . Or maybe error !
!Input: disline.input, part of the CA file , dis-line vertexes
!Output: awaydislocations.output, atom away from dis-lines , we desire
!Output: neardislocations.output, atom near dis-lines , we want to delete
!Output: output.dat, dis-line vertexes
!Output: outputpro.dat, vertexes with inserted positions in the gaps
!
!Note : in order to form a approximate cyinder, or pipe, we add many positions
!into the gaps of dis-lines. 
!
!lwh 
!2022.12.9  
!modified : 2022.12.12
!changed cutoff and file name : 2023.3.1

program main
implicit none
!------------------------------state the variables-------------------
!ADDED in 2024.4.27
integer :: g_typenum
double precision :: g_alat,g_1

!added variables
integer :: lines,l0
character(len=20) :: char1
!--------------------------------part 1 varias-----------------------
integer,allocatable,dimension(:) :: LineNum , StruType ,NumLinepos !the order-number & 
!of dislocation lines , structure type , vertex quantity
double precision,allocatable,dimension(:,:) :: Pos , LineBergers    !atom positions , &
!the bergers of dislocation lines 
integer :: Num , N  !the quantity of lines , the quantity of positions
integer :: i , j , k , k2  !omit
character(len=20) :: char0    !omit
integer,allocatable,dimension(:) :: itype   !atom type

!-------------------------------part 2 varias------------------------
integer :: IntNum , FN , near , FN2  !interstitials , quantity of Final ints
double precision ::  xlo , xhi , ylo , yhi , zlo , zhi !box region
double precision,allocatable,dimension(:,:) :: IntPos , Fpos,Fpos2  !positons of &
!interstitials , final positions
integer,allocatable,dimension(:) :: Id , itype2 , Occ , itypepro
double precision :: cutoff , dx , dy , dz , dr2 , dr !distance cutoff , distance

!-----------------------------part 3 varias-----------------------------
integer :: extra,type00   !extra atoms inserted in dis-lines
double precision :: dist   !distance , omit
double precision :: a0(3), b0(3), c0(3)  !basic vecters
double precision,allocatable,dimension(:,:) :: ProPos   !newly defined matrix for&
!inserted atoms and dis-line atoms
integer,allocatable,dimension(:) :: itype3 
integer :: l , pronum , extraonece   !omit,position quantity in ProPos matrix,&
!extra atom quantity in one gap
integer :: extraN  !atoms of original and inserted
!---------------------------------initialization-----------------------
open(1,file= '../../repair.params' ,status="old") !READ IN GLOBAL PARAMETERS
read(1,*)
read(1,*)
read(1,*)
read(1,*)g_alat
read(1,*)g_typenum
read(1,*)g_1

close(1)

cutoff = g_1*g_alat   !selection threshold

!---------------------------------read in words.dat-------------------
open(1,file= 'words.dat' ,status="old") !access how many lines in this file
read(1,*)lines
close(1)
write(*,*)"-->Total lines:",lines
!which line is the start line
open(1,file= 'disline.input' ,status="old") !dislocations start at line-i
do i = 1 , lines
        read(1,*)char1
        if(char1=='DISLOCATIONS')then
                l0 = i !line starts
                write(*,*)char1,'start at line',l0
        endif
enddo
close(1)
write(*,*)"-->Start line:",l0

!------------------------------read data from input file,except Pos------------
open(10,file='disline.input',status='old')
do i = 1 , l0-1 !omit first 'l0-1' lines
        read(10,*)
enddo
read(10,*)char0,Num
allocate(LineNum(1:Num))
allocate(StruType(1:Num))
allocate(NumLinepos(1:Num))
allocate(LineBergers(1:Num,1:3))

do i = 1,Num
  read(10,*)LineNum(i)
  read(10,*)LineBergers(i,1:3)
  read(10,*)StruType(i)
  read(10,*)NumLinepos(i)
  do j = 1,NumLinepos(i)
    read(10,*)
  enddo
enddo

close(10)

!-----------------------allocate and load Pos-------------------------- 
N = sum(NumLinepos(1:Num))
allocate(Pos(N,1:3))
allocate(itype(1:N))
k = 0

open(10,file='disline.input',status='old')
do i = 1 , l0-1 !omit first 'l0-1' lines
        read(10,*)
enddo
read(10,*)

do i = 1,Num
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  do j = 1,NumLinepos(i)
    read(10,*)Pos(k+j,1:3)
    itype(k+j) = 1
  enddo
  k = k + NumLinepos(i)
enddo

close(10)
!--------------------------cal insert atoms-----------------------
k = 0
dist = 0
extra = 0
do i = 1,Num  !one line by one line
  do j = 1,NumLinepos(i)
    if(j==NumLinepos(i))then
      exit
    else
      dist = sqrt( (Pos(k+j,1)-Pos(k+j+1,1))**2 + (Pos(k+j,2)-Pos(k+j+1,2))**2 + (Pos(k+j,3)-Pos(k+j+1,3))**2 )
      if(dist>cutoff/2)then
        extra = extra + floor(2*dist/cutoff)
      endif
    endif
  enddo
  k = k + NumLinepos(i)
enddo
!-------------------------insert atoms into dis-lines-------------
a0 = (/1.0d0, 0.0d0, 0.0d0/)
b0 = (/0.0d0, 1.0d0, 0.0d0/)
c0 = (/0.0d0, 0.0d0, 1.0d0/)
extraN = extra + N  !atom number of inserted plus original
allocate(ProPos(1:extraN,1:3))
allocate(itype3(1:extraN))
extraonece = 0
dist = 0
k = 0
pronum = 0   !represents how many atoms in the ProPos matrix at present
do i = 1,Num  !one line by one line
  do j = 1,NumLinepos(i)
    if(j==NumLinepos(i))then   !if this atom the last atom of the dis-line
      ProPos(pronum+1,1:3) = Pos(k+j,1:3)
      itype3(pronum+1) = 1
      pronum = pronum + 1
      exit
    else
      ProPos(pronum+1,1:3) = Pos(k+j,1:3)   !the first atom of this gap
      itype3(pronum+1) = 1
      pronum = pronum + 1
      dist = sqrt( (Pos(k+j,1)-Pos(k+j+1,1))**2 + (Pos(k+j,2)-Pos(k+j+1,2))**2 + (Pos(k+j,3)-Pos(k+j+1,3))**2 )
      if(dist>cutoff/2)then
        extraonece = floor(2*dist/cutoff)  !number of inserted atoms in this gap
        do l = 1,extraonece
          ProPos(pronum+l,1:3) = Pos(k+j,1:3) + (Pos(k+j+1,1)-Pos(k+j,1))/(extraonece+1)*l*a0 + &
(Pos(k+j+1,2)-Pos(k+j,2))/(extraonece+1)*l*b0 + (Pos(k+j+1,3)-Pos(k+j,3))/(extraonece+1)*l*c0 
          itype3(pronum+l) = 1
        enddo
        pronum = pronum + extraonece
      endif
    endif
  enddo
  k = k + NumLinepos(i)
enddo
      
!---------------------------read the interstitial positions , LAMMPS dump file format------------
open(12,file='input.data',status='old')
read(12,*)
read(12,*)IntNum
allocate(Id(1:IntNum))
allocate(itype2(1:IntNum))
allocate(IntPos(1:IntNum,1:3))
allocate(Occ(1:IntNum))
read(12,*)type00
read(12,*)xlo,xhi
read(12,*)ylo,yhi
read(12,*)zlo,zhi
read(12,*)
read(12,*)
read(12,*)
do i = 1,IntNum
  read(12,*)Id(i),itype2(i),IntPos(i,1:3)
enddo

close(12)
!------------------------modify Pos matrix into box region----------
do i = 1,N
  if (Pos(i,1)<xlo)then
    Pos(i,1) = Pos(i,1) + (xhi-xlo)
  endif
  if (Pos(i,1)>xhi)then
    Pos(i,1) = Pos(i,1) - (xhi-xlo)
  endif
  if (Pos(i,2)<ylo)then
    Pos(i,2) = Pos(i,2) + (yhi-ylo)
  endif
  if (Pos(i,2)>yhi)then
    Pos(i,2) = Pos(i,2) - (yhi-ylo)
  endif
  if (Pos(i,3)<zlo)then
    Pos(i,3) = Pos(i,3) + (zhi-zlo)
  endif
  if (Pos(i,3)>zhi)then
    Pos(i,3) = Pos(i,3) - (zhi-zlo)
  endif
enddo
!------------------------modify ProPos matrix into box region----------
do i = 1,extraN
  if (ProPos(i,1)<xlo)then
    ProPos(i,1) = ProPos(i,1) + (xhi-xlo)
  endif
  if (ProPos(i,1)>xhi)then
    ProPos(i,1) = ProPos(i,1) - (xhi-xlo)
  endif
  if (ProPos(i,2)<ylo)then
    ProPos(i,2) = ProPos(i,2) + (yhi-ylo)
  endif
  if (ProPos(i,2)>yhi)then
    ProPos(i,2) = ProPos(i,2) - (yhi-ylo)
  endif
  if (ProPos(i,3)<zlo)then
    ProPos(i,3) = ProPos(i,3) + (zhi-zlo)
  endif
  if (ProPos(i,3)>zhi)then
    ProPos(i,3) = ProPos(i,3) - (zhi-zlo)
  endif
enddo
!----------------------Test Tor----------------------------

!write(*,*)"test:Ndis-line,Ndis-line atom,loop aom,inserted atom,all atom"
!write(*,*)Num , N , k , extra , extraN

!!------------------------------output position data, Pos matrix----------------
!open(11,file='output.dat',status='replace')
!
!write(11,*)"dislocation atom positions"
!write(11,*)N , "atoms"
!write(11,*)"1 atom types"
!write(11,*) "0.0 0.0 xlo xhi"
!write(11,*) "0.0 0.0 ylo yhi"
!write(11,*) "0.0 0.0 zlo zhi"
!write(11,*)
!write(11,'(A)') "Atoms"
!write(11,*)
!do i = 1,N
!  write(11,'(I8,2x,I2,2x,3(f19.8,1x))') i, itype(i) , Pos(i,1:3)
!enddo
!
!close(11)
!
!!------------------------------output position data, ProPos matrix----------------
!allocate(itypepro(1:extraN))
!itypepro = 1
!open(14,file='outputpro.dat',status='replace')
!
!write(14,*)"dislocation atom positions"
!write(14,*)extraN , "atoms"
!write(14,*)"1 atom types"
!write(14,*) "0.0 0.0 xlo xhi"
!write(14,*) "0.0 0.0 ylo yhi"
!write(14,*) "0.0 0.0 zlo zhi"
!write(14,*)
!write(14,'(A)') "Atoms"
!write(14,*)
!do i = 1,extraN
!  write(14,'(I8,2x,I2,2x,3(f19.8,1x))') i, itypepro(i) , ProPos(i,1:3)
!enddo
!
!close(11)

!------------------------count atoms away from dis-lines-----------------

FN = 0  !atoms away from dis-lines
FN2 = 0
dr = 0
!cutoff = 5.0
do i = 1,IntNum !loop interstitial atoms
  near = 0
  do j = 1,extraN   !loop dis-line atoms
    dx = IntPos(i,1)-ProPos(j,1)
    dy = IntPos(i,2)-ProPos(j,2)
    dz = IntPos(i,3)-ProPos(j,3)
    dr2 = dx**2+dy**2+dz**2
    dr = sqrt(dr2)
    if (dr<cutoff) then
      near = near + 1
      exit 
    else
    continue
    endif
  enddo
  if (near == 0) then
    FN = FN + 1
  else
    FN2 = FN2 + 1
  endif
enddo

!---------------------select atoms away from dis-lines---------
k = 0  !count atoms away from dis-lines
k2 = 0  !count atoms near dis-lines
allocate(Fpos(1:FN,1:3))
allocate(Fpos2(1:FN2,1:3))
do i = 1,IntNum
  near = 0
  do j = 1,extraN
    dx = IntPos(i,1)-ProPos(j,1)
    dy = IntPos(i,2)-ProPos(j,2)
    dz = IntPos(i,3)-ProPos(j,3)
    dr2 = dx**2+dy**2+dz**2
    dr = sqrt(dr2)
    if (dr<cutoff) then
      near = near + 1
      exit
    else
      continue
    endif
  enddo
  if (near == 0) then !means away from dis-lines
    k = k + 1
    Fpos(k,1:3) = IntPos(i,1:3)
  else  !means near the dis-lines
    k2 = k2 + 1
    Fpos2(k2,1:3) = IntPos(i,1:3)
  endif
enddo

!------------------output final/selected atom positions,away----------------------
open(13,file='output_away.data',status='replace')

write(13,*)"selected interstitial atom positions, away"
write(13,*)FN , "atoms"
write(13,*)g_typenum," atom types"
write(13,*) xlo,xhi," xlo xhi"
write(13,*) ylo,yhi," ylo yhi"
write(13,*) zlo,zhi," zlo zhi"
write(13,*)
write(13,'(A)') "Atoms"
write(13,*)
do i = 1,FN
  write(13,'(I8,2x,I2,2x,3(f19.8,1x))') i, itype2(i) , Fpos(i,1:3)
enddo

close(13)

!------------------output final/selected atom positions,near----------------------
open(15,file='output_near.data',status='replace')

write(15,*)"selected interstitial atom positions, near"
write(15,*)FN2 , "atoms"
write(15,*)g_typenum," atom types"
write(15,*) "0.0 0.0 xlo xhi"
write(15,*) "0.0 0.0 ylo yhi"
write(15,*) "0.0 0.0 zlo zhi"
write(15,*)
write(15,'(A)') "Atoms"
write(15,*)
do i = 1,FN2
  write(15,'(I8,2x,I2,2x,3(f19.8,1x))') i, itype2(i) , Fpos2(i,1:3)
enddo

close(15)

!----------------------test for----------------------------

!write(*,*)'test:vertexes,N of ints,atoms away from dis-line,atoms near dis-line,Nall'
!write(*,*)N , IntNum , FN , FN2 , FN+FN2

!----------------------moniter datas-----------------------

write(*,*)"-----------Data moniter:-----------"
write(*,*)'@@@Cutoff:                          ',cutoff
write(*,*)'Number of dislocation lines:     ',Num
write(*,*)'Number of vertexes:              ',N
write(*,*)'Number of inserted vertexes:     ',extra
write(*,*)'Total vertexes:                  ',extraN
write(*,*)'Number of interstitials:         ',IntNum
write(*,*)'@@@Atoms away from dis-lines',FN
write(*,*)'@@@Atoms near dis-lines     ',FN2
write(*,*)'Total ints                       ',FN+FN2

!----------------------end-----------------------------------
end
