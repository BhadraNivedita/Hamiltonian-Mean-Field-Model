	! integrates the driven hmf model by rk4 method

	implicit none

	integer, parameter :: mm=1000,nn=2*mm
	integer, parameter :: nh=20,tran=90000,iter=100000
	double precision   :: xx(nn),yy(nn),dy(nn)
	double precision   :: tau,dtau,r,phi,grnd
	double precision   :: KK,a,w
	integer            :: i,j,k,l,m,n,tran1,itlen
	double precision   :: KE,PE,H,H0,U0,pRange,Uavg,Tavg,Ravg

	common KK

	call sgrnd(201)

	open (unit=1,file='amp0.5a.dat')

	tau=0.0d0; dtau=0.01d0
	tran1=tran+1; itlen=iter-tran

	a=0.50d0

	do m = 1, 1	! drive frequency

	w=.5d0

	do n = 1, 10	! initial energy

	U0=n*0.01d0; H0=U0*mm
	Uavg=0.0d0; Tavg=0.0d0; Ravg=0.0d0

	do k = 1, nh	! history

	KK=a

	do i = 1, mm
	yy(i) = 0.01*(grnd()-0.5d0)
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

	!write(1,'(41(xf12.6))') H/mm, 2*KE/mm, r

	Uavg = Uavg + H/mm
	Tavg = Tavg + 2*KE/mm
	Ravg = Ravg + r

	
	end do	! iterations

        end do	! history

	Uavg = Uavg/(nh*itlen)
	Tavg = Tavg/(nh*itlen)
	Ravg = Ravg/(nh*itlen)

	write(1,'(100(xf12.6))') U0, Uavg, Tavg, Ravg

	end do	! initial energy

	write(1,'(100(xf12.6))')

	end do	! drive frequency

	end program

	! calculates the time derivative 

        subroutine derivs(tau,yy,dy)

        integer, parameter :: mm=1000, nn=2*mm
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

        integer, parameter :: mm=1000, nn=2*mm
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
	PARAMETER (NMAX=2000)
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
