	! integrates hmf model by rk4 method and finds the variation of temperature and order parameter (magnetization) with energy per particle

	implicit none

	integer, parameter :: mm=100, nn=2*mm
	integer, parameter :: nh=50, tran=10000, iter=20000
	double precision   :: xx(nn), yy(nn), dy(nn)
	double precision   :: tau, dtau, amp, omega, r, phi
 	double precision   :: pi, pi2, grnd, omega_span(5)
	integer            :: i,j,k,l,m,n,tran1,itlen
	double precision   :: KE, PE, H, U, a, Tavg, Ravg, Uavg

	common amp, omega

	pi   = 4.0d0*datan(1.0d0); pi2 = 2.0d0*pi!defining paramter values
	tau  = 0.0d0; dtau = 0.01d0
	tran1=tran+1; itlen=iter-tran

	amp=1.0d0
	omega_span(1)=5.0d0; omega_span(2)=10.0d0; omega_span(3)=50.0d0;!diff combination of frequency w
	omega_span(4)=1.0d0; omega_span(5)=100.0d0;

	call sgrnd(148)

	open (unit=1,file='Khmfsmall.dat',status='unknown')!opening a file to write

	do l = 1, 5	! omega

	omega = omega_span(l)

	do i = 1, 50	! energy per particle

	U = i*0.1d0; H = U*mm
	Tavg=0.0d0; Ravg=0.0d0; Uavg=0.0d0

	do j = 1, nh	! history

	do n = 1, mm
	yy(n) = 0.1*(grnd()-0.5d0)
	end do

	PE = 0.0d0
	do m = 2, mm
	do n = 1, mm
	PE = PE + dcos(yy(m)-yy(n))
	end do
	PE=PE+abs(amp*cos(omega*tau))*dcos(yy(1)-yy(n))
	end do
	PE = (1/nn)*(mm**2-PE)
	
	

	a = 2.0d0*sqrt((12.0d0/nn)*(H-PE))

	do n = mm+1, nn
        yy(n) = a*(grnd()-0.5d0)
        end do

	do k = 1, tran	! transients

	call derivs(tau,yy,dy)
        call rk4(yy,dy,nn,tau,dtau,xx,derivs)

	yy = xx

	end do	! transients

	do k = tran1, iter	! iterations

	call derivs(tau,yy,dy)
        call rk4(yy,dy,nn,tau,dtau,xx,derivs)

	yy = xx

	call order_parameter(xx,r,phi)!calling order paramater

	KE = 0.0d0
	do n = mm+1, nn
	KE = KE + yy(n)**2
	end do
	KE = KE/2.0d0
	
        PE = 0.0d0
	do m = 2, mm
	do n = 1, mm
	PE = PE + dcos(yy(m)-yy(n))
	end do
	PE=PE+abs(amp*cos(omega*tau))*dcos(yy(1)-yy(n))
	end do
	PE = (1/nn)*(mm**2-PE)
	
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

	write(1,'(41(xf10.6))') omega, Uavg, Tavg, Ravg!writing average values to datafile

	end do	! energy per particle

	write(1,'(41(xf10.6))')

	end do	! amplitude

	end program

!----------------------------------------------------------------------------- calculates the time derivative ------------------------------------------

        subroutine derivs(tau,yy,dy)

        integer, parameter :: mm=100, nn=2*mm
	double precision   :: yy(nn), dy(nn)
	double precision   :: pi, tau, dtau
	integer            :: i, j, k, l, m, n
        double precision   :: kk, amp, omega, r, phi

	common amp, omega

	call order_parameter(yy,r,phi)

	kk = abs(amp*dcos(omega*tau))
	

        do i = 1, mm
        dy(i) = yy(mm+i) 
        end do



	do i = mm+2, nn
        dy(i) = r*dsin(phi-yy(i-mm))
	dy(mm+1) = kk*r*dsin(phi-yy(i-mm))
        end do




        return
        end subroutine

!---------------------------------------------------------------- finds the magnitude and phase of the order parameter--------------------------------------

	subroutine order_parameter(theta,r,phi)

        integer, parameter :: mm=100, nn=2*mm
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

!---------------------------------------------------------------------------Subroutine for RK$ algorithm---------------------------------------------------------

	SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)

	INTEGER n,NMAX
	REAL*8 h,x,dydx(n),y(n),yout(n)
	EXTERNAL derivs
	PARAMETER (NMAX=200)
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
