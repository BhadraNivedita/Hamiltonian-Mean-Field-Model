! Timeseries for different values of w

	implicit none

	integer, parameter :: mm=100,nn=2*mm
	integer, parameter :: tran=48000,iter=50000
	double precision   :: xx(nn),yy(nn),dy(nn)
	double precision   :: tau,dtau,r,phi,rx,ry,grnd
	double precision   :: KK,a,w
	integer            :: i,j,k,l,m,n,p,tran1,itlen
	double precision   :: KE,PE,H,H0,U0,pRange,ravg

	common KK

	call sgrnd(143)

	open (unit=1,file='Mvsw.dat')

	tau=0.0d0; dtau=0.01d0
	tran1=tran+1 ; itlen=iter-tran1

	a=1.0d0; 
	!w=1.0d0

	do p = 1, 1 ! w loop

	w = 0.5*p

	U0=0.4d0; H0=U0*mm;ravg=0

	KK=a

	do i = 1, mm
	yy(i) = 0.1*(grnd()-0.5d0)
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

	call order_parameter(xx,r,phi,rx,ry)

	write(1,'(41(xf12.6))') l*dtau, r

	ravg=ravg+r

	end do	! iterations

	ravg=ravg/itlen

	!write(1,'(41(xf12.6))') w, ravg	

	end do   !w loop 


	end program

	! calculates the time derivative 

        subroutine derivs(tau,yy,dy)

        integer, parameter :: mm=100, nn=2*mm
	double precision   :: yy(nn), dy(nn)
	double precision   :: tau, dtau
	integer            :: i, j, k, l, m, n
        double precision   :: KK, r, phi, rx, ry

	common KK

	call order_parameter(yy,r,phi,rx,ry)

        do i = 1, mm
        dy(i) = yy(mm+i)
        end do
	do i = mm+1, nn
        dy(i) = KK*r*dsin(phi-yy(i-mm))
        end do

        return
        end subroutine

	! finds the magnitude and phase of the order parameter

	subroutine order_parameter(theta,r,phi,rx,ry)

        integer, parameter :: mm=100, nn=2*mm
        double precision   :: theta(nn)
        double precision   :: r, phi, rx, ry

        rx = 0.0d0
        ry = 0.0d0

        do j=1,mm
        rx = rx + dcos(theta(j))
        ry = ry + dsin(theta(j))
        enddo

        rx = rx/mm
        ry = ry/mm
        r = dsqrt((rx)**2 + (ry)**2)
        phi = dacos(rx/r)

        return
        end

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
