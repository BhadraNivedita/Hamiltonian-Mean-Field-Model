	! integrates the driven hmf model by rk4 method

	implicit none

	integer, parameter :: mm=2000,nn=2*mm
	integer, parameter :: nh=1,tran=9000,iter=10000
	double precision   :: xx(nn),yy(nn),dy(nn)
	double precision   :: tau,dtau,r,phi,grnd
	double precision   :: KK,a,w
	integer            :: i,j,k,l,m,n,tran1,itlen
	double precision   :: KE,PE,H,H0,U0,pRange,Uavg,Tavg,Ravg

	common KK

	call sgrnd(101)

	open (unit=1,file='qss_psection.dat')

	tau=0.0d0; dtau=0.1d0
	tran1=tran+1; itlen=iter-tran

	a=1.0d0

	do m = 1, 1	! drive frequency

	w=0.0d0

	do n = 1, 1	! initial energy

	U0=0.88d0; H0=U0*mm
	Uavg=0.0d0; Tavg=0.0d0; Ravg=0.0d0

	do k = 1, nh	! history

	KK=a

	do i = 1, mm
	yy(i) = (grnd()-0.5)*3.4	!r=0.6
	end do

	PE = 0.0d0
	do i = 1, mm
	do j = 1, mm
	PE = PE + dcos(yy(i)-yy(j))
	end do
	end do
	PE = (KK/nn)*(mm**2-PE)

	pRange = 2.0d0*dsqrt((12.0d0/nn)*(H0-PE))

	do i = mm+1, nn
        yy(i) = pRange*(grnd()-0.5d0)
        end do

	KE = 0.0d0
	do i = mm+1, nn
	KE = KE + yy(i)**2
	end do
	KE = KE/2.0d0

	do i = mm+1, nn
	yy(i) = yy(i) * dsqrt((H0-PE)/KE)
	end do

	call order_parameter(yy,r,phi)

	write(*,*) r

	do l = 1, tran	! transients

	KK = abs(a*dcos(w*l*dtau))

	call derivs(tau,yy,dy)
        call rk4(yy,dy,nn,tau,dtau,xx,derivs)

	yy = xx

	end do	! transients

	do l = tran1, iter	! iterations

	KK = abs(a*dcos(w*l*dtau))

	call derivs(tau,yy,dy)
        call rk4(yy,dy,nn,tau,dtau,xx,derivs)

	yy = xx

	call order_parameter(xx,r,phi)

	KE = 0.0d0
	do i = mm+1, nn
	KE = KE + yy(i)**2
	end do
	KE = KE/2.0d0

	PE = 0.0d0
	do i = 1, mm
	do j = 1, mm
	PE = PE + dcos(yy(i)-yy(j))
	end do
	end do
	PE = (KK/nn)*(mm**2-PE)

	H = PE + KE

	Uavg = Uavg + H/mm
	Tavg = Tavg + 2*KE/mm
	Ravg = Ravg + r

	write(1,'(10000(xf16.10))') r,(xx(mm+i*50),xx(i*50)-phi,i=1,40)

	end do	! iterations

        end do	! history

	Uavg = Uavg/(nh*itlen)
	Tavg = Tavg/(nh*itlen)
	Ravg = Ravg/(nh*itlen)

	end do	! initial energy

	write(*,*) Ravg

	end do	! drive frequency

	end program

	! calculates the time derivative 

        subroutine derivs(tau,yy,dy)

        integer, parameter :: mm=2000, nn=2*mm
	double precision   :: yy(nn), dy(nn)
	double precision   :: tau, dtau
	integer            :: i, j, k, l, m, n
        double precision   :: KK, r, phi

	common KK

	call order_parameter(yy,r,phi)

        do i = 1, mm
        dy(i) = yy(mm+i)
        end do
	do i = mm+1, nn
        dy(i) = KK*r*dsin(phi-yy(i-mm))
        end do

        return
        end subroutine

	! finds the magnitude and phase of the order parameter

	subroutine order_parameter(theta,r,phi)

        integer, parameter :: mm=2000, nn=2*mm
        double precision   :: theta(nn)
        double precision   :: r, phi, real_sum, imag_sum

        real_sum = 0.0d0
        imag_sum = 0.0d0

        do j=1,mm
        real_sum = real_sum + dcos(theta(j))
        imag_sum = imag_sum + dsin(theta(j))
        enddo

        real_sum = real_sum/mm
        imag_sum = imag_sum/mm
        r = dsqrt((real_sum)**2 + (imag_sum)**2)
        phi = dacos(real_sum/r)

        return
        end

	SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)

	INTEGER n,NMAX
	REAL*8 h,x,dydx(n),y(n),yout(n)
	EXTERNAL derivs
	PARAMETER (NMAX=4000)
	INTEGER i
	REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)

	hh=h*0.5
	h6=h/6.
	xh=x+hh

	do i=1,n
        yt(i)=y(i)+hh*dydx(i)
	end do

	call derivs(xh,yt,dyt)
	do i=1,n
        yt(i)=y(i)+hh*dyt(i)
	end do

	call derivs(xh,yt,dym)
	do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
	end do

	call derivs(x+h,yt,dyt)
	do i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
	end do

	return
	END

	include 'mt.f'
