program main

implicit none
!Event number
integer::ie,nevent
!Line number
integer                                 ::il
!String event number
character(4)::xstring
!Indices
integer::i,j,k,m,ii
!First arrival travel times (nz,nx,ny)
real,allocatable,dimension(:,:,:)::tcalp,tcals
!Grid size
integer                          ::nz,nx,ny
!Mesh spacing h : dz=dx=dy
real::h
!Number of stations
integer::nstatot,nstatp,nstats
!Station names
character(3),allocatable,dimension(:)::namestatot,namestatp,namestats
!Station indices
integer,allocatable,dimension(:)::numberstatp,numberstats
!Observed arrival times
real,allocatable,dimension(:,:)::tobsp,tobss
!Result tables
real,allocatable,dimension(:,:,:)::d,res
real,allocatable,dimension(:,:,:)::resmc
!Weights sum
real::q,q2
!Quality parameters
real::a,b
!Results
integer::loczn,locxn,locyn
real::resn,t0n
integer                                 ::locmczn,locmcxn,locmcyn
real                                    ::t0mcn,resmcn
!Read
integer::IOstatus
!Event data
character(4)                            ::temp
character(15)                           ::temp1,temp2,temp4,temp5
character(21)::temp3
!Event-end lines
integer,dimension(9999)                 ::num
!Result
real,allocatable,dimension(:,:)                  ::output
!Store res matrices
character(len=8) :: fmt ! format descriptor
character(5) xie

!Write( xstring, '(i0)') l
! set size of 3D model
! Grille de 25 x 21 x 9km
nx=251
ny=211
nz=91

!Mesh step
h=0.100

!Number of stations
nstatot=20

! Allocate array 
! -> declare la taille des tableaux
allocate(namestatot(nstatot))
allocate(namestatp(nstatot))
allocate(namestats(nstatot))
allocate(numberstatp(nstatot))
allocate(numberstats(nstatot))
allocate(tobsp(nstatot,2))
allocate(tobss(nstatot,2))
allocate(tcalp(nz,nx,ny))
allocate(tcals(nz,nx,ny))
allocate(d(nz,nx,ny))
allocate(res(nz,nx,ny))
allocate(resmc(nz,nx,ny))

num=0
open(12,file=('R98in'),status='unknown')  ! lit les records
IOstatus=0  ! gérer les messages d'erreur
ie=1  ! indice event
il=1  ! indice ligne
do while(IOstatus .eq. 0)
  temp='0'
  read(12,*,IOSTAT=IOstatus)temp  ! lit la ligne
  if (temp .eq. '10') then  ! si 10 -> new event
    num(ie+1)=il  ! pour savoir où commence chaque event
    ie=ie+1
  endif
  il=il+1
enddo

nevent=ie-1  ! car on a augmenté pour le dernier
il=il-2
close(12)  ! ferme le fichier de records
allocate(output(nevent,13))  ! nb lignes, 13 colonnes


open(10,file=('stautm98.t'),status='unknown')  ! lit le nom des stations
do j=1,nstatot	
	read(10,'(A3)')namestatot(j)  ! lit 3 caractères dans le fichier 10 -> puis l'écrit dans le tableau des stations
enddo
close(10)

ie=1
open(42,file='R98in',status='unknown')  ! relit les records
do while(ie .le. nevent)  ! Grosse boucle sur les events
  nstatp=0
  nstats=0
  namestatp='000'
  namestats='000'
  numberstatp=0
  numberstats=0
  tobsp=0.
  tobss=0.
  tcalp=0.
  tcals=0.
  do ii=1,num(ie+1)-num(ie)-1  ! combien de lignes par event
    read(42,'(A5,A4,A21,A7,A3)')temp1,temp2,temp3,temp4,temp5  ! décrit la ligne
    !cas ligne P
!    	write(*,*)ii,'p',temp5
    if(temp2(1:1) .eq. 'P') then
      nstatp=nstatp+1  ! pour connaitre le nombre de stations P
!	write(*,*)ii,'p',temp1(2:4),temp3   
      ! Récupérer le nom de la station (' xxx' ou 'xxx1')
      if (temp1(1:1) .eq. ' ') then
      namestatp(nstatp)=temp1(2:4)
      else
      namestatp(nstatp)=temp1(1:3)
      endif	
        tobsp(nstatp,2)=ichar(temp2(3:3))-48  ! convertit le facteur de qualité en entier 
        ! On gère les jours (en rajoutant +1j ou +2j en secondes car on est le 7 pour les 1ers)
        if(temp3(6:6) .eq. '8')then
          tobsp(nstatp,1)=86400
       endif
       if(temp3(6:6) .eq. '9')then
         tobsp(nstatp,1)=172800
       endif
	write(*,*) 
      ! convertit tout en secondes : AAMMJJHHMMSS.ss
       tobsp(nstatp,1)=tobsp(nstatp,1)+(10*(ichar(temp3(7:7))-48)+ &
(ichar(temp3(8:8))-48))*3600+(10*(ichar(temp3(9:9))-48)+ &
(ichar(temp3(10:10))-48))*60+10*(ichar(temp3(11:11))-48)+ &
(ichar(temp3(12:12))-48)+0.1*(ichar(temp3(14:14))-48)+ &
0.01*(ichar(temp3(15:15))-48)
    else
      ! stations S
      nstats=nstats+1
 !    write(*,*)ii,'s',temp1,temp2,temp3,temp4,temp5
      if (temp1(1:1) .eq. ' ') then
      namestats(nstats)=temp1(2:4)
      else
      namestats(nstats)=temp1(1:3)
      endif
!write(*,*)namestats(nstats)
      tobss(nstats,2)=ichar(temp5(3:3))-48
      if(temp3(6:6) .eq. '8')then
        tobss(nstats,1)=86400
      endif
      if(temp3(6:6) .eq. '9')then
        tobss(nstats,1)=172800
      endif
      tobss(nstats,1)=tobss(nstats,1)+(10*(ichar(temp3(7:7))-48)+ &
(ichar(temp3(8:8))-48))*3600+(10*(ichar(temp3(9:9))-48)+ &
(ichar(temp3(10:10))-48))*60+100*(ichar(temp4(1:1))-48)+ &
10*(ichar(temp4(2:2))-48)+(ichar(temp4(3:3))-48) &
+0.1*(ichar(temp4(5:5))-48)+0.01*(ichar(temp4(6:6))-48)
    endif
  enddo
  read(42,*)
  ! Verifie qu'on a bien toutes les stations
  if (nstats+nstatp .ne. num(ie+1)-num(ie)-1) then
    write(*,*)'erreur1'
  endif

!Indice station à la place du nom car dans `grids/` on utilise les stations
  do i=1,nstatp  ! vient de R98in
    do j=1,nstatot  ! vient de stautm98
      if(namestatp(i) .eq. namestatot(j)) then  ! qd on a le match
        numberstatp(i)=j
 !    write(*,*)numberstatp(i)
      endif
    enddo
  enddo

  ! pareil pour les S
  do i=1,nstats
    do j=1,nstatot
      if(namestats(i) .eq. namestatot(j)) then
        numberstats(i)=j
!     write(*,*)numberstats(i)
      endif
    enddo
  enddo

  ! On convertit les facteurs de qualité (3 = mauvais -> 0.125 vs 0. bien -> 1 pour bien prendre en compte le résidu dans le solveur)
!q total weight
  do i=1,nstatp
    if (tobsp(i,2)==0) then
      tobsp(i,2)=1
    else if (tobsp(i,2)==1) then
      tobsp(i,2)=0.5
    else if (tobsp(i,2)==2) then
      tobsp(i,2)=0.25
    else if (tobsp(i,2)==3) then
      tobsp(i,2)=0.125
    else if (tobsp(i,2)==4) then
      tobsp(i,2)=0.  ! observations avec 0 en facteur de qualité ne compte même pas dans le solveur
    else
      write(*,*)'erreur2 p',ie,nstatp
    endif
  enddo
  do i=1,nstats
    if (tobss(i,2)==0) then
      tobss(i,2)=1
    else if (tobss(i,2)==1) then
      tobss(i,2)=0.5
    else if (tobss(i,2)==2) then
      tobss(i,2)=0.25
    else if (tobss(i,2)==3) then
      tobss(i,2)=0.125
    else if (tobss(i,2)==4) then
      tobss(i,2)=0.
    else
      write(*,*)'erreur2 s',ie,nstats
    endif
  enddo
  ! somme des facteurs de qualité -> pour normaliser en norme L1 (moins sensible aux outliers)
  q=(sum(tobsp(:,2))+sum(tobss(:,2)))
  ! somme des facteurs de qualité au carré -> pour normaliser en norme L2
  q2=0.
  do i=1,nstatp
    q2=q2+tobsp(i,2)**2
  enddo
  do i=1,nstats
    q2=q2+tobss(i,2)**2
  enddo

  ! compute sum of tobs - tpropa -> c'est le t0 estimé pour chaque station
  ! Pour chaque station et chaque event, on calcule t0 = tobs(event) - tcal(station)(point) en chaque point de grille
  d(:,:,:)=0.
  ! on parcourt les stations P
  do i=1,nstatp
    !open grid tcalp(i)
    j=numberstatp(i)  ! trouver le numéro pour ouvrir la grille des tcal
    Write( xstring, '(i0)') j	
    open(30,file='./grids/gridp'//trim(xstring)//'.bin',form='unformatted', &
status='unknown', access='direct', recl=4*ny*nx*nz)
    ! On récupère les temps de parcours calculés pour cette station
    read(30,rec=1)tcalp
    close(30)
    !tobsp(:,2) quality factor
    d=d+tobsp(i,2)*(tobsp(i,1)-tcalp(:,:,:))
    ! dans d, en chaque point, on somme pondérée les t0 = tobs - tcal trouvé dans chaque station
  end do
  ! On fait pareil pour les stations S
  do i=1,nstats
    !open grid tcals(i)
    j=numberstats(i)
    Write( xstring, '(i0)') j	
	open(30,file='./grids/grids'//trim(xstring)//'.bin',form='unformatted', &
status='unknown', access='direct', recl=4*ny*nx*nz)
    read(30,rec=1)tcals
    close(30)
    !quality factor
    d=d+tobss(i,2)*(tobss(i,1)-tcals(:,:,:))
  end do
  ! On normalise la somme pondérée (divise par la somme des facteurs de qualité)
  d=d/q  ! -> chaque point de d contient le t0 de cet event estimé par moyenne pondérée sur les stations

  ! On calcule les résidus cette fois: sum(tobs(i) - (t_0 + t_cal(i))) (au carré pour L2) avec t_0 estimé ci-dessus
  res  =0.  ! norme L1
  resmc=0.  ! norme L2
  do m=1,nstatp
    j=numberstatp(m)
    Write( xstring, '(i0)') j	
    open(30,file='./grids/gridp'//trim(xstring)//'.bin',form='unformatted', &
status='unknown',access='direct', recl=4*ny*nx*nz)
    read(30,rec=1)tcalp
    close(30)
    res=res+tobsp(m,2)*abs(tobsp(m,1)-(tcalp(:,:,:)+d))
    resmc=resmc+(tobsp(m,2)**2)*((tobsp(m,1)-(tcalp(:,:,:)+d))**2)
  end do
  do m=1,nstats
    j=numberstats(m)
    Write( xstring, '(i0)') j	
    open(30,file='./grids/grids'//trim(xstring)//'.bin',form='unformatted', &
status='unknown',access='direct', recl=4*ny*nx*nz)
    read(30,rec=1)tcals
    close(30)
    res=res+tobss(m,2)*abs(tobss(m,1)-(tcals(:,:,:)+d))
    resmc=resmc+(tobss(m,2)**2)*((tobss(m,1)-(tcals(:,:,:)+d))**2)
  end do
  ! Normalise
  res=res/q
  resmc=resmc/q2

! Pour enregistrer la matrice des résidus
  if (mod(ie, 50) .eq. 0) then ! pour un event sur 100
    fmt = '(I4.4)' ! an integer of width 3 with zeros at the left
    write (xie,fmt) ie ! converting integer to string using a 'internal file'
    open(110,file='./results/res/resL1_'//trim(xie)//'.bin',form='unformatted',status='unknown',access='direct', recl=4*ny*nx*nz)
    write(110,rec=1)res
    close(110)
    open(120,file='./results/res/resL2_'//trim(xie)//'.bin',form='unformatted',status='unknown',access='direct', recl=4*ny*nx*nz)
    write(120,rec=1)resmc
    close(120)
  end if

! On cherche les résidus minimaux
!res2 unicity
  ! on part de grands résidus
  resn=999.
  resmcn=999.
  do k=1,ny
  do j=1,nx
  do i=1,nz
    ! si nouveau résidu plus petit, on stocke la valeur et là où atteint
    if(res(i,j,k) .lt. resn) then
      resn = res(i,j,k)
      locxn = j
      locyn = k
      loczn = i
      t0n= d(i,j,k)
    end if
    if(resmc(i,j,k) .lt. resmcn) then
      resmcn = resmc(i,j,k)
      locmcxn = j
      locmcyn = k
      locmczn = i
      t0mcn= d(i,j,k)
    end if

  end do
  end do
  end do

  output(ie,1)=ie
  output(ie,2)=t0n
  ! Calculs pour se replacer en coordonnées UTM
  output(ie,3)=(loczn-1)*h-2.7
  output(ie,4)=(locxn-1)*h+352
  output(ie,5)=(locyn-1)*h-2361
  output(ie,6)=resn
  output(ie,7)=q
  output(ie,13)=t0mcn
  output(ie,8)=(locmczn-1)*h-2.7
  output(ie,9)=(locmcxn-1)*h+352
  output(ie,10)=(locmcyn-1)*h-2361
  output(ie,11)=resmcn
  output(ie,12)=q2
  CALL system('clear')
  write(*,*)100*ie/nevent,'%'  ! barre de progression

  ie=ie+1
enddo
open(18,file=('./results/output'),status='unknown')
do i=1,nevent
  write(18,*) output(i,1),output(i,2),output(i,3),output(i,4),output(i,5), &
output(i,6),output(i,7),output(i,13),output(i,8),output(i,9),output(i,10), &
output(i,11),output(i,12)
enddo
  close(18)
call exit
end program main

