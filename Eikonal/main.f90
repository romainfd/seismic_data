program main

implicit none

real,allocatable,dimension(:,:,:)::velp,vels,tcalp,tcals
integer				 ::nx,ny,nz,ierr,i,j,k,nsweep,nstat
real				 ::dzin,dxin,dyin
real,allocatable,dimension(:)  	::zr,xr,yr
!Station names
character(3),allocatable,dimension(:)	::namestat
real,allocatable,dimension(:,:)  	::coordstat
character (2)				::xstring

nz=91
nx=251
ny=213
dzin=100.
dxin=100.
dyin=100.
nsweep=1
nstat=20

allocate(tcalp(nz,nx,ny),stat=ierr)
allocate(tcals(nz,nx,ny),stat=ierr)
allocate(velp(nz-1,nx-1,ny-1),stat=ierr)
allocate(vels(nz-1,nx-1,ny-1),stat=ierr)
allocate(namestat(nstat))
allocate(coordstat(nstat,3))
allocate(zr(nstat))
allocate(xr(nstat))
allocate(yr(nstat))

!Velocity model
velp(1:27,:,:)=3000.
velp(28:37,:,:)=3400.
velp(38:47,:,:)=3800.
velp(48:57,:,:)=4250.
velp(58:77,:,:)=5000.
velp(78:90,:,:)=6000.
vels=velp/1.75


open(10,file=('stautm98.t'),status='unknown')
do j=1,nstat	
read(10,*) namestat(j),coordstat(j,2),coordstat(j,3),coordstat(j,1)
zr(j)=(coordstat(j,1)+2.7)*1000.	
xr(j)=(coordstat(j,2)-352.)*1000.
yr(j)=(coordstat(j,3)+2361.)*1000.
write(*,*) zr(j),xr(j),yr(j)
enddo
close(10)


!	open(130,file='velpsyn.bin',form='unformatted',status='unknown',access='direct', &
!	     recl=4*(nx-1)*(ny-1)*(nz-1))
!	write(130,rec=1)vels
!	close(130)



do i=1,nstat
	Write( xstring, '(i0)') i
write(*,*) xstring
	call FTeik3d_2(velp,tcalp,nz,nx,ny,zr(i),xr(i),yr(i),dzin,dxin,dyin,nsweep,5.)
write(*,*) "Eik OK"
	open(30,file='./grids/gridp'//trim(xstring)//'.bin',form='unformatted',&
status='unknown',access='direct', recl=4*ny*nx*nz)
	write(30,rec=1)tcalp
	close(30)
	call FTeik3d_2(vels,tcals,nz,nx,ny,zr(i),xr(i),yr(i),dzin,dxin,dyin,nsweep,5.)
	open(20,file='./grids/grids'//trim(xstring)//'.bin',form='unformatted',&
status='unknown',access='direct', recl=4*ny*nx*nz)
	write(20,rec=1)tcals
	close(20)
end do

call exit
end program main
