program compute_coulomb
implicit none
integer i, j, nqm, nmm, ntot
double precision dist, E_elec, E_field2, E_pol
double precision, dimension(:,:), allocatable :: r
double precision, dimension(:), allocatable :: q, a


!Read input files r.in and q.in

open (unit=55, file="r.in", status='old', action='read')
read(55,*) nqm, nmm, ntot
allocate(r(3,ntot),q(ntot),a(nqm))
do i=1,ntot
  read(55,*) (r(j,i),j=1,3)
end do
close (unit=55)

!r=r/0.529177d0

open (unit=66, file="q.in", status='old', action='read')
do i=1,ntot
  read(66,*) q(i)
end do
close (unit=66)

open (unit=77, file="a.in", status='old', action='read')
do i=1,nqm
  read(77,*) a(i)
end do
close (unit=77)

E_elec=0.d0

do i=1,nqm
  do j=nqm+1,ntot
     E_elec=E_elec + q(i)*q(j)/dist(r(1,i),r(2,i),r(3,i),r(1,j),r(2,j),r(3,j))
  end do
end do

E_pol=0.d0
do i=1,nqm
  E_pol=E_pol+a(i)*E_field2(r,q,nqm,ntot,i)
end do
E_pol=-0.5d0*E_pol
!E_pol=-E_pol

write(*,*) "----------------------------------------------"
write(*,*) "EQMMM_elec+EQMMM_nuc (Hartree)"
write(*,*) "----------------------------------------------"
write(*,*) "without polarization: ", E_elec
write(*,*) "polarization: ", E_pol
write(*,*) "summ: ", E_pol+E_elec
!write(*,*) "LIO", -0.033594d0
!write(*,*) "estimated/LIO", (E_pol+E_elec)/(-0.033594d0)
write(*,*) "---------------------------------------------"

end program compute_coulomb


function dist(x1,y1,z1,x2,y2,z2)
implicit none
double precision dist,x1,y1,z1,x2,y2,z2
dist = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
dist = sqrt(dist)
return
end function dist

function E_field2(r,q,nqm,ntot,i)
implicit none
double precision, dimension(3,ntot) :: r
double precision, dimension(ntot) :: q
double precision, dimension(3) :: rij, E_field
double precision distij, distij2, E_field2
integer nqm, ntot, i, j, k

E_field=0.d0
do j=nqm+1,ntot
  rij(1) = r(1,j)-r(1,i)
  rij(2) = r(2,j)-r(2,i)
  rij(3) = r(3,j)-r(3,i)
  distij2 = rij(1)**2 + rij(2)**2 + rij(3)**2
  distij = DSQRT(distij2)
  rij = rij/distij
  do k = 1,3
    E_field(k)=E_field(k)+(q(j)/distij2)*rij(k)
  end do
end do

E_field2=E_field(1)**2+E_field(2)**2+E_field(3)**2

end function
