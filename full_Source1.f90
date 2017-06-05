     program main
!             3-D COMPRESSIBLE MHD CODE
!    ************************************************
!
! 1. Scheme  ---Modify 2 step Lax-Wendroff. 
! 2. History --- completed in June, 1991.
! 3. Author  --- Zhiwei Ma (Geophysical Institute, UAF, AK 99775)
! 4. Project ---- Asymmetric Local magnetic reconnection.
!```````````````````````````````````````````````````````````
! 5. revised by cjxiao on 2004-04-20 (in subroutine init):
!    to chang X-line not perpendicular to MR plane. 
!
!
!                  Introduction
! The present version of this code simulate one half of the
! simulation box in x-y-z direction. In the six boundary,two boundaries
! at x=Lx and x=-Lx use inflow ones,three at y=Ly, y=-Ly , and z=Lz
! free (outgoing) ones,and one at z=-Lz symmetric (or antisymmetric ) 
! ones.
!
!     **********************************************
!-----------------
!     Basic variable definitions:
!     x(mx,my,mz,1) to x(mx,my,mz,8) represents, respectively,
!     rho, rho*vx, rho*vy, rho*vz, bx, by, bz and energy.
!-----------------
!
      include 'ma3ds1.for'
      include 'ma3ds2.for'
!
      dtime=10.         ! concernd about continuing running
      dnstep=1.         ! concernd about continuing running
      nstop=49          ! maximum time-step
	  t0=0.0    
	  nst2=5            ! 给nst 或者nst1 赋值的
	  cont=2    
	          
	  call input               ! time=0, nstep=0, nst1=1
!     call initia

! cwm add:start
      call harris_initia
!     call fluc_at_bndry
      call fluc_at_neutral_line
! cwm add:end
!
!     if(.not.lrstrt)then
!     call recrd
!     call facur
!     end if

  	  if(lrstrt) then 
	  call readin(nst2,cont,dtime)
	  end if

	  if(time.eq.0) then
      call recrd
      call recrd1
      call facur
      end if

! cwm add: start
      open(unit=112,file="rhoVx.dat",status="unknown",form="formatted")
      open(unit=117,file="Bz.dat",status="unknown",form="formatted")
      open(unit=120,file="steptotime.dat",status="unknown",form="formatted")
      write(112,9) (((x(jx,jy,jz,2),jx=1,mx),jy=1,my),jz=1,mz)      
      write(117,9) (((x(jx,jy,jz,7),jx=1,mx),jy=1,my),jz=1,mz)
    9 format(11767(1x,e22.15))      ! 11767=mx*my*mz
! cwm add: end

  200 continue
      call setdt

      call stepon
	
      nstep=nstep+1
      time=time+dt
      		 
!	  write(*,*)'nst=',nst
!     open(unit=16,file='stepnm',status='unknown',form='formatted')
!     write(16,99) nstep,time
!  99 format(i5,f9.5)
!     close(16)

!zxg to continue
! 	  if(abs(mod(time,2.d0)).le.dt)then
!	  open(unit=17,file='continue',status="unknown",form="unformatted")
!	  write(17)ncase,nstep,time,nst
!	  write(17)x
!	  close(17)
!	  else
!	  endif
!zxg continue end

	  if(abs(time-nst*dtime).le.dt.and.&
	  (nstep-nstp(nst)).ge.dnstep) then
      nst=nst+nint
      nstp(nst)=nstep
      call recrd
      call recrd1
!     call facur
      end if

!!yg----------------            !!!!!!!!  底下的没有考虑到新程序中
!     if ((time.gt.t0).and.((time-t0).le.dt)) then
!	  call incident_plasma(x,xi,time,t0,0)
!!    call initialSF(x,xi,time,t0)   
!     end if !
!     if ((time.gt.t0).and.((time-t0).gt.dt)) then
!	  call incident_plasma(x,xi,time,t0,1)
!     end if
!!yg---------------

! cwm add: start      
      if( mod(nstep,12) .eq. 0 ) then
      write(112,9) (((x(jx,jy,jz,2),jx=1,mx),jy=1,my),jz=1,mz)      
      write(117,9) (((x(jx,jy,jz,7),jx=1,mx),jy=1,my),jz=1,mz)
      endif
      write(120,*) nstep,' ','time=',time,'','dt=',dt
! cwm add: end

      if(nstep.gt.nstop) goto 300
      if(nst.lt.nend) goto 200

  300 continue

! cwm add: start
      close(112)
      close(117)
! cwm add: end

      stop
      end



      subroutine input
! --------------------------
!   This routine inputs parameters and define basic variables
!   LRSTRT: =.f., starting from t=0; =.t., continueing from
!           steps as given by NST.
!   NEND:   the final steps intended for the current run ,
!           including the steps from the previous run.
!   NSP:    time step interval for diagnostic plots.
!   DELTAV: thickness of the initial velocity shear.
!   DELTAC: thickness of the initial current sheet (for
!           normalization coventions , see routine INITIA.
!   RATIO:  ratio of dy to dx, or the ratio of the width to
!           length of the box, if the number of grid points is such
!           that my=mx.
!   NSMTHX:  the number of rows starting from incoming
!           boundary for the smoothing region in x-direction.
!   NSMTHY:  the number of rows starting from incoming
!           boundary for the smoothing region in y-direction.
!   ETA:    exact inverse of magnetic Renolds number.
!   GAMMA:  adiabatic constant.
!   BETAS:   ratio of kinetic to magnetic pressure at magnetosheath side
!   MACHS:   Alfven mach number of incoming flow  at magnetosheath side.
!   MACHM:   Alfven mach number of incoming flow at magnetopause side.
!   ANGS:    bz=b*cos(ang) and by=b*sin(ang), at magnetosheath side.
!   ANGM:    bz=b*cos(ang) and by=b*sin(ang), at magnetopause side.
!   TS0:   initia magnetosheath temperature
!   Tm0:   initia magnetopause temperature.
!   BS0:   initia magnetosheath magnetic field strength
!   BM0:   initia magnetopause magnetic field strength.
!   NCASE:  case number of the run.
!   NST:    beginning data file number for the current run.
!   NINT:   data file number increment.
! --------------------------
!
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      include 'ma3ds3.for'
!
! Coordinates system
!                A z     /
!                l     / 
!                l   /
!                l /
!  <-------------*--------------
!  x            /l
!             /  l
!           /    l
!         /      l
!        y
!
      nxp1=mx
      nyp1=my
      nzp1=mz
      nx=mx-1
      ny=my-1
      nz=mz-1
      call gridpnt
      open(unit=11,file='grid.dat',status='unknown',form='formatted')
      write(11,99)(xx(jx),jx=1,mx),(zz(jz),jz=1,mz)
   99 format(5(1x,e11.4))
!
      time=0.
      nstep=0
      nst1=1
!
      return
      end
!
!
!
      subroutine initia
!
!----------------
! Defines coordinates system and specifies initial configuration.
! Normlization convention:
!   1. Density --- normalised to asymtotic value, i.e., rho=1
!   2. Magnetic field --- normalised to asymtotic value, i.e.,
!                         b0=1.
!   3. Velocity --- normalised to asymtotic Alfven speed, VA=1, a
!                   natural result of 1. and 2.
!   4. Length --- normalised to a=10*dx, i.e., dx=0.1
!   5. Time --- normalised to a/VA.
!---------------
!
      include 'ma3ds1.for'
      include 'ma3ds2.for'
!
! assign asymmetric quantities:
      pi     = 3.1415926
      rhom   = 1.0d0
      rhos   = 1.0d0
      pp     = 0.
      gaminv = 1./gamma
      phirad = 2.*pi*phi/360.
      delbz  = 0.5*sqrt( bm0**2+bs0**2-2.*bm0*bs0*cos(phirad) )
      by0    = 0.5*bs0*sin(phirad)/delbz
      bz0    = 0.25*(bm0**2 - bs0**2)/delbz
      pmsp   = betam*(0.5*bm0**2)             ! pressure at magnetopause
      pmsh   = pmsp + 0.5*(bm0**2 - bs0**2)    ! magnetopause and magnetosheath pressure balance
      betas  = pmsh/(0.5*bs0**2)
      p0     = pmsp + 0.5*(bm0**2 - bs0**2)
      bzm    = bz0 + delbz
      bzs    = bz0 - delbz
      ymax   = -yy(1)
!
      signy=1
      if(phi.lt.180) signy=-1

!cjx  the angle between the X-line(Y-direction) and MR plane
      thita_Xline=0.0
      thita_Xline_rad = 2.*pi*thita_Xline/360.
      tan_Xline=tan(thita_Xline_rad)
!cjx  -------------------------------------------------	   
      do 3 jz=1,nzp1
      do 3 jy=1,nyp1
      do 3 jx=1,nxp1
     x(jx,jy,jz,1) = 0.5*(rhom+rhos) +0.5*(rhom-rhos)*&          !rhom=rhos
                        tanh((xx(jx)+yy(jy)*tan_Xline)/aw)
      x(jx,jy,jz,5)  = -bm0*3*xx(jx)*zz(jz)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)       ! 注意这里的 5/2 =2 !!!!!!! 确定要这样子做？

	  if  ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  x(jx,jy,jz,5)=x(jx-1,jy,jz,5) 
	  end if

      x(jx,jy,jz,6)  = -bm0*3*yy(jy)*zz(jz)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)
	  
	  if   ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  x(jx,jy,jz,6)=x(jx-1,jy,jz,6) 
	  end if

!	x(jx,jy,jz,6)  =0.1
      x(jx,jy,jz,7)  = -bm0*(2*zz(jz)**2-xx(jx)**2-yy(jy)**2)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)
	  
	  if   ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  x(jx,jy,jz,7)=x(jx-1,jy,jz,7) 
	  end if

      xi(jx,jy,jz,1) = 0.5*(rhom+rhos) &
        +0.5*(rhom-rhos)*tanh(((xx(jx)+yy(jy)*tan_Xline)+dx/2.)/aw)
      xi(jx,jy,jz,5)  =-bm0*3*xx(jx)*zz(jz)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)
	  
	  if   ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  xi(jx,jy,jz,5)=xi(jx-1,jy,jz,5) 
	  end if

      xi(jx,jy,jz,6)  = -bm0*3*yy(jy)*zz(jz)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)

	  	  if   ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0) then
	  xi(jx,jy,jz,6)=xi(jx-1,jy,jz,6) 
	  end if
!	xi(jx,jy,jz,6)  =0.1
      xi(jx,jy,jz,7)  = -bm0*(2*zz(jz)**2-xx(jx)**2-yy(jy)**2)/(xx(jx)**2+yy(jy)**2+zz(jz)**2)**(5/2)

	  	  if   ((xx(jx)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  xi(jx,jy,jz,7)=xi(jx-1,jy,jz,7) 
	  end if
    3 continue
!
      do 4 jz=1,nzp1
      do 4 jy=1,nyp1
      do 4 jx=1,nxp1
      x(jx,jy,jz,2)=  0.d0
      x(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*(1.-tanh((xx(jx)+yy(jy)*tan_Xline)/aw))*x(jx,jy,jz,1)/bs0
      x(jx,jy,jz,4)=0.5*v0*x(jx,jy,jz,7)*(1.-tanh((xx(jx)+yy(jy)*tan_Xline)/aw))*x(jx,jy,jz,1)/bs0
      xi(jx,jy,jz,2)= 0.d0
      xi(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*(1.-tanh(((xx(jx)+yy(jy)*tan_Xline)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
      xi(jx,jy,jz,4)  = 0.5*v0*x(jx,jy,jz,7)*(1.-tanh(((xx(jx)+yy(jy)*tan_Xline)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
    4 continue
!cjx----------------------------------------------------------------
!      do 3 jz=1,nzp1
!      do 3 jy=1,nyp1
!      do 3 jx=1,nxp1
!      x(jx,jy,jz,1) = 0.5*(rhom+rhos) +0.5*(rhom-rhos)*tanh(xx(jx)/aw)
!      x(jx,jy,jz,5)  = 0.0
!      x(jx,jy,jz,6)  = signy*sqrt(by0*by0+epsilon*delbz*delbz
!     1              /cosh(xx(jx)/aw)**2)
!      x(jx,jy,jz,7)  = bz0 - delbz*tanh(xx(jx)/aw)
!      xi(jx,jy,jz,1) = 0.5*(rhom+rhos) 
!     1   +0.5*(rhom-rhos)*tanh((xx(jx)+dx/2.)/aw)
!      xi(jx,jy,jz,5)  = 0.0
!      xi(jx,jy,jz,6)  = signy*sqrt(by0*by0+epsilon*delbz*delbz
!     1              /cosh((xx(jx)+dx/2.)/aw)**2)
!      xi(jx,jy,jz,7)  = bz0 - delbz*tanh((xx(jx)+dx/2.)/aw)
!    3 continue

!      do 4 jz=1,nzp1
!      do 4 jy=1,nyp1
!      do 4 jx=1,nxp1
!      x(jx,jy,jz,2)=0.
!      x(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*(1.-tanh(xx(jx)/aw))
!     1                    *x(jx,jy,jz,1)/bs0
!      x(jx,jy,jz,4)=0.5*v0*x(jx,jy,jz,7)*(1.-tanh(xx(jx)/aw))
!     1                    *x(jx,jy,jz,1)/bs0
!      xi(jx,jy,jz,2)=0.
!      xi(jx,jy,jz,3)=0.5*v0*x(jx,jy,jz,6)*
!     1        (1.-tanh((xx(jx)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
!      xi(jx,jy,jz,4)  = 0.5*v0*x(jx,jy,jz,7)*
!     1        (1.-tanh((xx(jx)+dx/2.)/aw))*x(jx,jy,jz,1)/bs0
!    4 continue
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      presc=0.5*bm0**2+pmsp
      do 5 jz=1,nzp1
      do 5 jy=1,nyp1
      do 5 jx=1,nxp1
      hb2=.5*(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2)
      hv2=.5*(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2+x(jx,jy,jz,4)**2)/x(jx,jy,jz,1)
      pr(jx,jy,jz)=presc-hb2
      x(jx,jy,jz,8)=hv2+hb2+pr(jx,jy,jz)/(gamma-1)
      hb2=.5*(xi(jx,jy,jz,5)**2+xi(jx,jy,jz,6)**2+xi(jx,jy,jz,7)**2)
      hv2=.5*(xi(jx,jy,jz,2)**2+xi(jx,jy,jz,3)**2+xi(jx,jy,jz,4)**2)/xi(jx,jy,jz,1)
      pri=presc-hb2
      xi(jx,jy,jz,8)=hv2+hb2+pri/(gamma-1)
    6 continue
    5 continue
!
      call current(xi,1)
      do 71 jz=1,nzp1
      do 71 jy=1,nyp1
      do 71 jx=1,nxp1
      xi(jx,jy,jz,3)=xi(jx,jy,jz,3)+vyi0*di*w0(jx,jy,jz,2)    ! vyi0=0
   71 continue
      call current(x,1)
      do 72 jz=1,nzp1
      do 72 jy=1,nyp1
      do 72 jx=1,nxp1
      x(jx,jy,jz,3)=x(jx,jy,jz,3)+vyi0*di*w0(jx,jy,jz,2)
   72 continue
      do 8 m=1,8
      do 8 jx=1,nxp1
      fx(jx,m)=x(jx,1,1,m)
      fxi(jx,m)=xi(jx,1,1,m)
    8 continue
      do 9 m=1,8
      do 9 jz=1,nzp1
      do 9 jy=1,nyp1
      do 9 jx=1,nxp1
      xm(jx,jy,jz,m)=0.
      xi(jx,jy,jz,m)=0.
    9 continue

      return
      end
!
! cwm add: start
      subroutine harris_initia
      include 'ma3ds1.for'
      include 'ma3ds2.for'
! assign asymmetric quantities:      
      do 3 jz=1,nzp1
      do 3 jy=1,nyp1
      do 3 jx=1,nxp1
      x(jx,jy,jz,1) = (1./cosh(xx(jx)/normlambda))**2+rhoinfinity
      x(jx,jy,jz,7) = tanh(xx(jx)/normlambda)
	  x(jx,jy,jz,6) = 0.
      x(jx,jy,jz,5) = 0.

      xi(jx,jy,jz,1) = (1./cosh((xx(jx)+dx/2.)/normlambda))**2+rhoinfinity
      xi(jx,jy,jz,7) = tanh((xx(jx)+dx/2.)/normlambda)
      xi(jx,jy,jz,6) = 0.
      xi(jx,jy,jz,5) = 0.
    3 continue
!
      do 4 jz=1,nzp1
      do 4 jy=1,nyp1
      do 4 jx=1,nxp1
      x(jx,jy,jz,2) = 0.
      x(jx,jy,jz,3) = x(jx,jy,jz,1)*balcoeff*2./normlambda
      x(jx,jy,jz,4) = 0.
      xi(jx,jy,jz,2) = 0.
      xi(jx,jy,jz,3) = xi(jx,jy,jz,1)*balcoeff*2./normlambda
      xi(jx,jy,jz,4) = 0.
    4 continue

      do 5 jz=1,nzp1
      do 5 jy=1,nyp1
      do 5 jx=1,nxp1
      hb2=.5*(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2)
      hv2=.5*(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2+x(jx,jy,jz,4)**2)/x(jx,jy,jz,1)
      pr(jx,jy,jz)=balcoeff*x(jx,jy,jz,1)
      x(jx,jy,jz,8)=hv2+hb2+pr(jx,jy,jz)/(gamma-1)
      hb2=.5*(xi(jx,jy,jz,5)**2+xi(jx,jy,jz,6)**2+xi(jx,jy,jz,7)**2)
      hv2=.5*(xi(jx,jy,jz,2)**2+xi(jx,jy,jz,3)**2+xi(jx,jy,jz,4)**2)/xi(jx,jy,jz,1)
      pri=balcoeff*xi(jx,jy,jz,1)
      xi(jx,jy,jz,8)=hv2+hb2+pri/(gamma-1)
    6 continue
    5 continue

      do 8 m=1,8
      do 8 jx=1,nxp1
      fx(jx,m)=x(jx,1,1,m)
      fxi(jx,m)=xi(jx,1,1,m)
    8 continue
      do 9 m=1,8
      do 9 jz=1,nzp1
      do 9 jy=1,nyp1
      do 9 jx=1,nxp1
      xm(jx,jy,jz,m)=0.
      xi(jx,jy,jz,m)=0.
    9 continue

      return
      end

      subroutine fluc_at_bndry
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      do 9 jz=1,nzp1
      do 9 jy=1,nyp1
      x(1,jy,jz,5) = x(1,jy,jz,5)-kz*fluc*sin(kz*zz(jz))
      x(2,jy,jz,5) = x(2,jy,jz,5)-kz*fluc*sin(kz*zz(jz))
      x(nx,jy,jz,5) = x(nx,jy,jz,5)-kz*fluc*sin(kz*zz(jz))
      x(nxp1,jy,jz,5) = x(nxp1,jy,jz,5)-kz*fluc*sin(kz*zz(jz))
    9 continue
      end

      subroutine fluc_at_neutral_line
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      do 9 jz=1,nzp1
      do 9 jy=1,nyp1
      do 9 jx=1,nxp1
      if (abs(xx(jx))<normlambda) then
      x(jx,jy,jz,5) = x(jx,jy,jz,5)-kz*fluc*cos(kx*xx(jx))*sin(kz*zz(jz))
      x(jx,jy,jz,7) = x(jx,jy,jz,7)+kx*fluc*sin(kx*xx(jx))*cos(kz*zz(jz))
      endif
    9 continue
      end
! cwm add: end


      subroutine stepon
!
!     This routine time-advances X's by 2 step Lax-Wendroff scheme
!     note: X is always the up-to-date value while Xm being the
!           intermediate value.
!
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      dimension xhelp(216,12)
!
! 1. first step
!
! 1.1 Calculate fluxes

! 1.2 Advance the first step
      hdt=0.5*dt                       ! hdt for half_dt
      call current(x,1)
      call foreta(time,1)
      call pressure(x,1)
      do 1000 m=1,8
      call flux(x,fs,gs,hs,mx,my,mz,m,1)
      do 100 jz=1,nz
      do 100 jy=1,ny
      do 100 jx=1,nx
      fdx=-hdt*(fs(jx+1,jy+1,jz+1)-fs(jx,jy+1,jz+1) &
         +fs(jx+1,jy+1,jz)-fs(jx,jy+1,jz)&
         +fs(jx+1,jy,jz+1)-fs(jx,jy,jz+1)&
         +fs(jx+1,jy,jz)-fs(jx,jy,jz))/(4.*dx)
      gdy=-hdt*(gs(jx+1,jy+1,jz+1)-gs(jx+1,jy,jz+1)&
         +gs(jx+1,jy+1,jz)-gs(jx+1,jy,jz)&
         +gs(jx,jy+1,jz+1)-gs(jx,jy,jz+1)&
         +gs(jx,jy+1,jz)-gs(jx,jy,jz))/(4.*dy)
      hdz=-hdt*(hs(jx+1,jy+1,jz+1)-hs(jx+1,jy+1,jz) &
         +hs(jx+1,jy,jz+1)-hs(jx+1,jy,jz)&
         +hs(jx,jy+1,jz+1)-hs(jx,jy+1,jz)&
         +hs(jx,jy,jz+1)-hs(jx,jy,jz))/(4.*dz)
      ppx=(xi(jx,jy,jz,m)+xi(jx+1,jy,jz,m)+xi(jx,jy+1,jz,m)&
         +xi(jx,jy,jz+1,m)+xi(jx+1,jy+1,jz,m)+xi(jx,jy+1,jz+1,m)&
         +xi(jx+1,jy,jz+1,m)+xi(jx+1,jy+1,jz+1,m))/8.0
      xm(jx,jy,jz,m)=fxi(jx,m)+ppx+fdx+gdy+hdz
  100 continue
 1000 continue
!
! 2. Second step
!
! 2.1 Calculate fluxes
! 2.2 Advance the second time step
      call current(xm,2)
      call foreta(time,2)
      call pressure(xm,2)
      do 2000 m=1,8
      call flux(xm,fs,gs,hs,nx,ny,nz,m,2)
      do 200 jz=2,nz
      do 200 jy=2,ny
      do 200 jx=2,nx
      fdx=-dt*(fs(jx,jy,jz)-fs(jx-1,jy,jz)&
         +fs(jx,jy,jz-1)-fs(jx-1,jy,jz-1)&
         +fs(jx,jy-1,jz)-fs(jx-1,jy-1,jz)&
         +fs(jx,jy-1,jz-1)-fs(jx-1,jy-1,jz-1))/(4.*dx)
      gdy=-dt*(gs(jx,jy,jz)-gs(jx,jy-1,jz)&
         +gs(jx,jy,jz-1)-gs(jx,jy-1,jz-1)&
         +gs(jx-1,jy,jz)-gs(jx-1,jy-1,jz)&
         +gs(jx-1,jy,jz-1)-gs(jx-1,jy-1,jz-1))/(4.*dy)
      hdz=-dt*(hs(jx,jy,jz)-hs(jx,jy,jz-1)&
         +hs(jx,jy-1,jz)-hs(jx,jy-1,jz-1)&
         +hs(jx-1,jy,jz)-hs(jx-1,jy,jz-1)&
         +hs(jx-1,jy-1,jz)-hs(jx-1,jy-1,jz-1))/(4.*dz)
      x(jx,jy,jz,m)=x(jx,jy,jz,m)+fdx+gdy+hdz
  200 continue
 2000 continue
!
      call bndry(x,1)
!
!      if(mod(nstep,10).eq.0) call cleanb
      caf0=1.-0.25*tanh(time/100.)            ! concered about smthf(xi,caf0)
      do 250 m=1,8
      do 250 jz=1,mz
      do 250 jy=1,my
      do 250 jx=1,mx
      xi(jx,jy,jz,m)=x(jx,jy,jz,m)-fx(jx,m)     !  这里会改 xi 的值   ！注意！注意！注意！
  250 continue
!      call smthf(xi,caf0)                                  !   这个是另一种平滑方式没有改写到c++程序中
      if(mod(nstep,3).eq.0) then
      call smthxyz(xi,0,1)    ! Smooth x y z ? and call bndry(x,2)
      call avrg1(xi,0.996d0)  ! m=12348 6-point average call bndry(x,2)
      call avrg2(xi,0.996d0)  ! m=5,6,7  call bndry(x,2)
      endif
      do 260 m=1,8
      do 260 jz=1,mz
      do 260 jy=1,my
      do 260 jx=1,mx
      x(jx,jy,jz,m)=xi(jx,jy,jz,m)+fx(jx,m)
  260 continue
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 从这开始一直到结束的部分没有被改写。
      if(mod(nstep,10).eq.0) then
      call current(x,1)
      do 300 jz=1,mz
      do 300 jy=1,my
      do 300 jx=1,mx
      xm(jx,jy,jz,1)=sqrt(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2)
      xm(jx,jy,jz,2)=(w0(jx,jy,jz,1)*x(jx,jy,jz,5)+w0(jx,jy,jz,2)&
        *x(jx,jy,jz,6)+w0(jx,jy,jz,3)*x(jx,jy,jz,7))/xm(jx,jy,jz,1)
  300 continue
      open(unit=14,file='facmax.dat',status='unknown',form='formatted')
      write(14,12)time,nstep,ncase
      k=0
  400 k=k+1
      jz=(mz*k/9.)+1
      cmax=-10.
      cmin=10.
      do 310 jy=1,my
      do 310 jx=1,mx-2
! Find xm(x,y,z,2)==J \dot B / |B|, that is field aligned current, its maximum value and corresponding coordinate point
      if(xm(jx,jy,jz,2).gt.cmax) then
      cmax=xm(jx,jy,jz,2)
      xcmax=xx(jx)
      ycmax=yy(jy)
      zcmax=zz(jz)
      endif
! Find its minimum value and cooresponding location
      if(xm(jx,jy,jz,2).lt.cmin) then
      cmin=xm(jx,jy,jz,2)
      xcmin=xx(jx)
      ycmin=yy(jy)
      zcmin=zz(jz)
      endif
!
  310 continue
      write(14,13)xcmin,ycmin,zcmin,xcmax,ycmax,zcmax,cmin,cmax
      if(k.lt.8) goto 400
   12 format(1x,f8.2,2(2x,i5))
   13 format(6(1x,f6.1),2(1x,e11.3))
!      if(mod(nstep,30).eq.0) then
!      open(unit=15,file='sat.dat',status='unknown',form='formatted')
!      write(15,12)(time,nstep,ncase)
!      write(15,17)((xm(jx,jy,6,2),jx=11,nx,10),jy=1,my,4)
!   17 format(8(1x,e9.3))
!      endif
      endif
!      if(mod(nstep,20).eq.0) call energy  !w0 should be current
!
      return
      end
!                       !!!!!!!!!!!!!!!!!!!!!!!!!!!!! facur 这个子程序也没有被改写
      subroutine facur
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      character*8 output
      character*3 cn     ! count number from readin
      dimension fact(mx,my,mz),facp(mx,my,mz),faci(mx,my,mz)
!
      call current(x,1)
      call pressure(x,1)
      do 1 jz=1,mz
      do 1 jy=1,my
      do 1 jx=1,mx
      xm(jx,jy,jz,1)=sqrt(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2)
      xm(jx,jy,jz,2)=(w0(jx,jy,jz,1)*x(jx,jy,jz,5)+w0(jx,jy,jz,2)*x(jx,jy,jz,6)+w0(jx,jy,jz,3)*x(jx,jy,jz,7))/xm(jx,jy,jz,1)**2
    1 continue
!
      do 2 jz=2,nz
      do 2 jy=2,ny
      do 2 jx=2,nx
      pp1=w0(jx,jy,jz,1)*(pr(jx+1,jy,jz)-pr(jx-1,jy,jz))/(2.*dx)
      pp2=w0(jx,jy,jz,2)*(pr(jx,jy+1,jz)-pr(jx,jy-1,jz))/(2.*dy)
      pp3=w0(jx,jy,jz,3)*(pr(jx,jy,jz+1)-pr(jx,jy,jz-1))/(2.*dz)
      pp4=(pr(jx+1,jy,jz)-pr(jx-1,jy,jz))/(2.*dx)&
       *(x(jx,jy,jz,7)*(xm(jx,jy+1,jz,1)-xm(jx,jy-1,jz,1))/(2.*dy)&
       -x(jx,jy,jz,6)*(xm(jx,jy,jz+1,1)-xm(jx,jy,jz-1,1))/(2.*dz))
      pp5=(pr(jx,jy+1,jz)-pr(jx,jy-1,jz))/(2.*dy)&
       *(x(jx,jy,jz,5)*(xm(jx,jy,jz+1,1)-xm(jx,jy,jz-1,1))/(2.*dz)&
       -x(jx,jy,jz,7)*(xm(jx+1,jy,jz,1)-xm(jx-1,jy,jz,1))/(2.*dx))
      pp6=(pr(jx,jy,jz+1)-pr(jx,jy,jz-1))/(2.*dz)&
       *(x(jx,jy,jz,6)*(xm(jx+1,jy,jz,1)-xm(jx-1,jy,jz,1))/(2.*dx)&
       -x(jx,jy,jz,5)*(xm(jx,jy+1,jz,1)-xm(jx,jy-1,jz,1))/(2.*dy))
      facp(jx,jy,jz)=-(pp1+pp2+pp3)/xm(jx,jy,jz,1)**2 &
                    +2.*(pp4+pp5+pp6)/xm(jx,jy,jz,1)**3
    2 continue
!
      do 3 jz=2,nz
      do 3 jy=2,ny
      do 3 jx=2,nx
      fact(jx,jy,jz)=x(jx,jy,jz,5)*(xm(jx+1,jy,jz,2)-xm(jx-1,jy,jz,2))/(2.*dx)&
        +x(jx,jy,jz,6)*(xm(jx,jy+1,jz,2)-xm(jx,jy-1,jz,2))/(2.*dy)&
        +x(jx,jy,jz,7)*(xm(jx,jy,jz+1,2)-xm(jx,jy,jz-1,2))/(2.*dz)
      faci(jx,jy,jz)=fact(jx,jy,jz)-facp(jx,jy,jz)
    3 continue
      call bndry1(fact,1)
      call bndry1(facp,1)
      call bndry1(faci,1)
      call vorticity(x)       ! Calculate vorticity form velocity variables
!
      output='fac'//cn(nst)
      open(unit=8,file=output,status="unknown",form="formatted")
      write(8,9)(((fact(jx,jy,jz),facp(jx,jy,jz),faci(jx,jy,jz),&
         pr(jx,jy,jz),(w0(jx,jy,jz,m),m=1,3),jx=1,mx),jy=1,my)&
          ,jz=1,mz)
    9 format(7(1x,e11.4))
      close(8)
      return
      end
!
!      subroutine energy
!      include 'ma3ds1.for'
!      include 'ma3ds2.for'
!      dimension wyz(my,mz),fyz(my,mz),gyz(my,mz),hyz(my,mz)
!      dimension wz(mz),fz(mz),gz(mz),hz(mz),cr1(mz),cr2(mz)
!      dimension nk1(mz),nk2(mz)
!!
!!  define statement functions
!!  d2fc= d2 f / dx2   with central difference
!!      d2fc(fm,f0,fp,xm1,x0,xp1)=
!!     1 2.*((fp-f0)/(xp1-x0)-(f0-fm)/(x0-xm1))/(xp1-xm1)  fp: f plus, fm: f minus
!!  d1fc= d f / dx  with  central difference
!!      d1fc(fm,f0,fp,xm1,x0,xp1)=((xm1-x0)/(xp1-x0)*(fp-f0)-(xp1-x0)/(xm1-x0)*(fm-f0))/(xm1-xp1)
!!
!      do 1 jz=2,mz-1
!      do 1 jy=2,my-1
!      do 1 jx=2,mx-1
!      fs(jx,jy,jz)=-d1fc(pr(jx-1,jy,jz),pr(jx,jy,jz),pr(jx+1,jy,jz)&
!                       ,xx(jx-1),xx(jx),xx(jx+1))&
!        +w0(jx,jy,jz,2)*x(jx,jy,jz,7)-w0(jx,jy,jz,3)*x(jx,jy,jz,6)
!      gs(jx,jy,jz)=-d1fc(pr(jx,jy-1,jz),pr(jx,jy,jz),pr(jx,jy+1,jz)&
!                              ,yy(jy-1),yy(jy),yy(jy+1))&
!        +w0(jx,jy,jz,3)*x(jx,jy,jz,5)-w0(jx,jy,jz,1)*x(jx,jy,jz,7)
!      hs(jx,jy,jz)=-d1fc(pr(jx,jy,jz-1),pr(jx,jy,jz),pr(jx,jy,jz+1)&
!                              ,zz(jz-1),zz(jz),zz(jz+1))&
!        +w0(jx,jy,jz,1)*x(jx,jy,jz,6)-w0(jx,jy,jz,2)*x(jx,jy,jz,5)
!    1 continue
!      call bndry1(fs,0)
!      call bndry1(gs,0)
!      call bndry1(hs,0)
!      do 2 jz=1,mz
!      do 2 jy=1,my
!      do 2 jx=1,mx
!      xm(jx,jy,jz,1)=fs(jx,jy,jz)**2+gs(jx,jy,jz)**2+hs(jx,jy,jz)**2
!      fs(jx,jy,jz)=0.5*(x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2&
!                 +x(jx,jy,jz,7)**2)
!      gs(jx,jy,jz)=0.5*(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2&
!                  +x(jx,jy,jz,4)**2)/x(jx,jy,jz,1)
!      hs(jx,jy,jz)=pr(jx,jy,jz)/(gamma-1.)
!      xm(jx,jy,jz,2)=(w0(jx,jy,jz,1)*x(jx,jy,jz,5)+&
!          w0(jx,jy,jz,2)*x(jx,jy,jz,6)+w0(jx,jy,jz,3)*&
!           x(jx,jy,jz,7))/sqrt(2.*fs(jx,jy,jz))
!! w0 should be current, thus xm(2)=J \dot B/|B|
!    2 continue
!! Accumulation begins, nk1 is for positive value counting numer; 
!! cr1 is accumulated result. nk2 and cr2 is for negative
!      do 3 jz=1,mz,mz/10
!      cr1(jz)=0.
!      cr2(jz)=0.
!      nk1(jz)=1
!      nk2(jz)=1
!      do 3 jy=1,my
!      do 3 jx=1,mx
!      crj=xm(jx,jy,jz,2)
!      if(crj.ge.0.) then
!      nk1(jz)=nk1(jz)+1
!      cr1(jz)=cr1(jz)+crj
!      else
!      nk2(jz)=nk2(jz)+1
!      cr2(jz)=cr2(jz)+crj
!      endif
!    3 continue
!! 4 for integration in x-direction
!! 5 for integration in y-direction
!! 6 for integration in z-direction
!! integrating xm(1) to wt meaning forces
!!             fs    to ft for Magnetic energy
!!             gs    to gt for kinetic energy
!!             hs    to ht for thermal energy
!      do 4 jz=1,mz
!      do 4 jy=1,my
!      call integ(xm(1,jy,jz,1),aa,xx,mx)
!      call integ(fs(1,jy,jz),bb,xx,mx)
!      call integ(gs(1,jy,jz),cc,xx,mx)
!      call integ(hs(1,jy,jz),dd,xx,mx)
!      wyz(jy,jz)=aa
!      fyz(jy,jz)=bb
!      gyz(jy,jz)=cc
!      hyz(jy,jz)=dd
!    4 continue
!      do 5 jz=1,mz
!      call integ(wyz(1,jz),aa,yy,my)
!      call integ(fyz(1,jz),bb,yy,my)
!      call integ(gyz(1,jz),cc,yy,my)
!      call integ(hyz(1,jz),dd,yy,my)
!      wz(jz)=aa
!      fz(jz)=bb
!      gz(jz)=cc
!      hz(jz)=dd
!    5 continue
!      call integ(wz,wt,zz,mz)
!      call integ(fz,ft,zz,mz)
!      call integ(gz,gt,zz,mz)
!      call integ(hz,ht,zz,mz)
!      open(unit=11,file='energy.dat',status='unknown',form='formatted')
!!      write(11,9)('Force$','ME$','KE$','TE$','Time$')
!      write(11,6)wt,ft,gt,ht,time
!    6 format(5(1x,e13.3))
!      open(unit=12,file='facur.dat',status='unknown',form='formatted')
!      write(12,7)(time)
!      write(12,11)(nk2(jz),nk2(jz),jz=1,mz,mz/10)
!      write(12,8)(cr1(jz),cr2(jz),jz=1,mz,mz/10)
!    7 format(f9.3)
!    8 format(8(1x,e9.3))
!   11 format(8(1x,i8))
!!    9 format(5(1x,a13))
!      return
!      end
!
      subroutine integ(fin,fout,x,mx)
      double precision fin(mx),x(mx),fout
      fout=0
      do 1 jx=2,mx
      fout=fout+(fin(jx-1)+fin(jx))*(x(jx)-x(jx-1))/2.
    1 continue
      return
      end
!
      subroutine setdt
      include 'ma3ds1.for'
      include 'ma3ds2.for'      
!
      dimension temp(mx)
      call foreta(time,1)
      call pressure(x,1)
      dtmin=1000.      ! A new variables
      dxyz=.5*dx*dy*dz/sqrt((dx**2*dy**2+dy**2*dz**2+dx**2*dz**2))
      do 1 jz=1,mz
      do 1 jy=1,my
      do 2 jx=1,mx
!      x(jx,jy,jz,1)=cvmgm(0.0001,x(jx,jy,jz,1),x(jx,jy,jz,1))
      temp(jx)=dxyz/( sqrt(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2&
              +x(jx,jy,jz,4)**2)/x(jx,jy,jz,1)&
              +sqrt((x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2&
             +x(jx,jy,jz,7)**2&
              +gamma*pr(jx,jy,jz))/x(jx,jy,jz,1)))
    2 continue
      test=1000.
      do 3 jx=1,nxp1
      if(temp(jx).lt.test) test=temp(jx)
    3 continue
      dtmin=dmin1(test,dtmin) ! d-min-1 pick minimal value
    1 continue
      dt=0.5*dtmin
      return
      end
!
!
      subroutine bndry(x,nlt)
!
!------------------------------
! Set the boundaries for X's
!------------------------------

! nlt seems unused!!!
!
      include 'ma3ds1.for'
!
      dimension x(mx,my,mz,8)
	lsz=1.0
!
! magnetosheath and -pause b.!.
! inflow boundary
      do 1 m=1,8
      do 1 jz=2,nz
      do 1 jy=2,ny
      x(1,jy,jz,m)=x(2,jy,jz,m)
      x(mx,jy,jz,m)=x(nx,jy,jz,m)
!cyg------------------------------
    
	
!cyg----------------------------	
    1 continue
!
      if(halfx) then
      do 2 jz=2,nz
      do 2 jy=2,ny
      x(mx,jy,jz,1)=x(nx-1,my-jy+1,jz,1)
      x(mx,jy,jz,2)=-x(nx-1,my-jy+1,jz,2)
      x(mx,jy,jz,3)=-x(nx-1,my-jy+1,jz,3)
      x(mx,jy,jz,4)=x(nx-1,my-jy+1,jz,4)
      x(mx,jy,jz,5)=x(nx-1,my-jy+1,jz,5)
      x(mx,jy,jz,6)=x(nx-1,my-jy+1,jz,6)
      x(mx,jy,jz,7)=-x(nx-1,my-jy+1,jz,7)
      x(mx,jy,jz,8)=x(nx-1,my-jy+1,jz,8)
    2 continue
      endif
!
! out flowing B.C.
!
!
      if(periody) then
      do 3 m=1,8
      do 3 jz=2,nz
      do 3 jx=1,mx
      x(jx,1,jz,m)=x(jx,ny,jz,m)
      x(jx,my,jz,m)=x(jx,2,jz,m)
    3 continue
      else
      do 4 m=1,8
      do 4 jz=2,nz
      do 4 jx=1,mx
      x(jx,1,jz,m)=x(jx,2,jz,m)
      x(jx,my,jz,m)=x(jx,ny,jz,m)
    4 continue
      endif
!
      do 5 m=1,8
      do 5 jy=1,my
      do 5 jx=1,mx
      x(jx,jy,1,m)=x(jx,jy,2,m)
      x(jx,jy,mz,m)=x(jx,jy,nz,m)
    5 continue
      if(halfz) then
      do 6 jy=1,my
      do 6 jx=1,mx
      x(jx,jy,mz,1)=x(jx,jy,nz-1,1)
      x(jx,jy,mz,2)=x(jx,jy,nz-1,2)
      x(jx,jy,mz,3)=x(jx,jy,nz-1,3)
      x(jx,jy,mz,4)=-x(jx,jy,nz-1,4)
      x(jx,jy,mz,5)=-x(jx,jy,nz-1,5)
      x(jx,jy,mz,6)=-x(jx,jy,nz-1,6)
      x(jx,jy,mz,7)=x(jx,jy,nz-1,7)
      x(jx,jy,mz,8)=x(jx,jy,nz-1,8)
    6 continue
      endif
!
      return
      end
!
!
      subroutine bndry1(x,nlt)
!
!------------------------------
! Set the boundaries for X's, 
! nlt is for equivalent extrapolation
! nlt==0 for ositive dot-symmetry
! nlt.neq.0 for negative dot-symmytry
!------------------------------
!
      include 'ma3ds1.for'
!
      dimension x(mx,my,mz)
!
      if(nlt.eq.0) then
! magnetosheath and -pause b.!.
      do 1 jz=2,nz
      do 1 jy=2,ny
      x(1,jy,jz)=x(2,jy,jz)
      x(nxp1,jy,jz)=x(nx,jy,jz)
    1 continue
!
      if(halfx) then
      do 2 jz=2,nz
      do 2 jy=2,ny
      x(nxp1,jy,jz)=x(nx-1,my-jy+1,jz)
    2 continue
      endif
!
      if(periody) then
      do 3 jz=2,nz
      do 3 jx=1,nxp1
      x(jx,1,jz)=x(jx,ny,jz)
      x(jx,nyp1,jz)=x(jx,2,jz)
    3 continue
      else
      do 4 jz=2,nz
      do 4 jx=1,nxp1
      x(jx,1,jz)=x(jx,2,jz)
      x(jx,nyp1,jz)=x(jx,ny,jz)
    4 continue
      endif
!
      do 5 jy=1,nyp1
      do 5 jx=1,nxp1
      x(jx,jy,1)=x(jx,jy,2)
      if(halfz) then
      x(jx,jy,nzp1)=x(jx,jy,nz-1)
      else
      x(jx,jy,mz)=x(jx,jy,nz)
      endif
    5 continue
!
      else
! magnetosheath and -pause b.!.
      do 6 jz=2,nz
      do 6 jy=2,ny
      x(1,jy,jz)=x(2,jy,jz)
      x(nxp1,jy,jz)=x(nx,jy,jz)
    6 continue
! 
      if(halfx) then
      do 66 jz=2,nz
      do 66 jy=2,ny
      x(nxp1,jy,jz)=-x(nx-1,my-jy+1,jz)
   66 continue
      endif
!
      if(periody) then
      do 7 jz=2,nz
      do 7 jx=1,nxp1
      x(jx,1,jz)=x(jx,ny,jz)
      x(jx,nyp1,jz)=x(jx,2,jz)
    7 continue
      else
      do 8 jz=2,nz
      do 8 jx=1,nxp1
      x(jx,1,jz)=x(jx,2,jz)
      x(jx,nyp1,jz)=x(jx,ny,jz)
    8 continue
      endif
!
      do 9 jy=1,nyp1
      do 9 jx=1,nxp1
      x(jx,jy,1)=x(jx,jy,2)
      if(halfz) then
      x(jx,jy,mz)=-x(jx,jy,nz-1)
      else
      x(jx,jy,mz)=x(jx,jy,nz)
      endif
    9 continue
      endif
!
      return
      end
!
!
      subroutine flux(x,fs,gs,hs,nnx,nny,nnz,m,mm)
!
!-------------------------------
!  Calculate fluxes
!  Notations: X1    X2     X3     X4     X5  X6  x7  
!             rho   rhovx  rhovy  rhovz  bx  by  bz  e
!-------------------------------
!
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8),fs(mx,my,mz)
      dimension hs(mx,my,mz),gs(mx,my,mz)
!
      if(m.eq.1) then
      do 1 jz=1,nnz
      do 1 jy=1,nny
      do 1 jx=1,nnx
!
! [1] Continuity eq.
!
      fs(jx,jy,jz)=x(jx,jy,jz,2)
      gs(jx,jy,jz)=x(jx,jy,jz,3)
      hs(jx,jy,jz)=x(jx,jy,jz,4)
    1 continue
!
      else
      if(m.eq.2) then
! [2] Momentum eq.
      do 2 jz=1,nnz
      do 2 jy=1,nny
      do 2 jx=1,nnx
!
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      fs(jx,jy,jz)=x(jx,jy,jz,2)**2/x(jx,jy,jz,1)+pr(jx,jy,jz)&
                    +.5*b2-x(jx,jy,jz,5)**2
      gs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,3)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,6)*x(jx,jy,jz,5)
      hs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,4)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,7)*x(jx,jy,jz,5)
    2 continue
      else
      if(m.eq.3) then
      do 3 jz=1,nnz
      do 3 jy=1,nny
      do 3 jx=1,nnx
!
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      fs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,3)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,6)*x(jx,jy,jz,5)
      gs(jx,jy,jz)=x(jx,jy,jz,3)**2/x(jx,jy,jz,1)+pr(jx,jy,jz)&
                    +.5*b2-x(jx,jy,jz,6)**2
      hs(jx,jy,jz)=x(jx,jy,jz,4)*x(jx,jy,jz,3)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,6)*x(jx,jy,jz,7)
!
    3 continue
      else
      if(m.eq.4) then
!
      do 4 jz=1,nnz
      do 4 jy=1,nny
      do 4 jx=1,nnx
!
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      fs(jx,jy,jz)=x(jx,jy,jz,2)*x(jx,jy,jz,4)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,5)*x(jx,jy,jz,7)
      gs(jx,jy,jz)=x(jx,jy,jz,3)*x(jx,jy,jz,4)/x(jx,jy,jz,1)&
                    -x(jx,jy,jz,6)*x(jx,jy,jz,7)
      hs(jx,jy,jz)=x(jx,jy,jz,4)**2/x(jx,jy,jz,1)+pr(jx,jy,jz)&
                    +.5*b2-x(jx,jy,jz,7)**2
    4 continue
      else
      if(m.eq.5) then
      call current(x,mm)
      call foreta(time,mm)
      do 5 jz=1,nnz
      do 5 jy=1,nny
      do 5 jx=1,nnx
!
! [3] Magnetic induction eq.
!      vcrb-z: V-J cross B's z component
      vcrbz=((x(jx,jy,jz,2)-di*w0(jx,jy,jz,1))*x(jx,jy,jz,6)&
           -(x(jx,jy,jz,3)-di*w0(jx,jy,jz,2))*x(jx,jy,jz,5))&
             /x(jx,jy,jz,1)
      vcrby=((x(jx,jy,jz,4)-di*w0(jx,jy,jz,3))*x(jx,jy,jz,5)&
           -(x(jx,jy,jz,2)-di*w0(jx,jy,jz,1))*x(jx,jy,jz,7))&
             /x(jx,jy,jz,1)
!
      fs(jx,jy,jz)=0.
      gs(jx,jy,jz)=-vcrbz+etaf(jx,jy,jz)*w0(jx,jy,jz,3)
      hs(jx,jy,jz)=vcrby-etaf(jx,jy,jz)*w0(jx,jy,jz,2)
    5 continue
      else
      if(m.eq.6) then
      call current(x,mm)
      call foreta(time,mm)
      do 6 jz=1,nnz
      do 6 jy=1,nny
      do 6 jx=1,nnx
!
! [3] Magnetic induction eq.
!
!      vcrb-z: V-J cross B's z component
!
      vcrbz=((x(jx,jy,jz,2)-di*w0(jx,jy,jz,1))*x(jx,jy,jz,6)&
           -(x(jx,jy,jz,3)-di*w0(jx,jy,jz,2))*x(jx,jy,jz,5))&
             /x(jx,jy,jz,1)
      vcrbx=((x(jx,jy,jz,3)-di*w0(jx,jy,jz,2))*x(jx,jy,jz,7)&
           -(x(jx,jy,jz,4)-di*w0(jx,jy,jz,3))*x(jx,jy,jz,6))&
             /x(jx,jy,jz,1)
!
      fs(jx,jy,jz)=vcrbz-etaf(jx,jy,jz)*w0(jx,jy,jz,3)
      gs(jx,jy,jz)=0.
      hs(jx,jy,jz)=-vcrbx+etaf(jx,jy,jz)*w0(jx,jy,jz,1)
    6 continue
      else
      if(m.eq.7) then
      call current(x,mm)
      call foreta(time,mm)
      do 7 jz=1,nnz
      do 7 jy=1,nny
      do 7 jx=1,nnx
!
! [3] Magnetic induction eq.
!
      vcrby=((x(jx,jy,jz,4)-di*w0(jx,jy,jz,3))*x(jx,jy,jz,5)&
           -(x(jx,jy,jz,2)-di*w0(jx,jy,jz,1))*x(jx,jy,jz,7))&
             /x(jx,jy,jz,1)
      vcrbx=((x(jx,jy,jz,3)-di*w0(jx,jy,jz,2))*x(jx,jy,jz,7)&
           -(x(jx,jy,jz,4)-di*w0(jx,jy,jz,3))*x(jx,jy,jz,6))&
             /x(jx,jy,jz,1)
!
      fs(jx,jy,jz)=-vcrby+etaf(jx,jy,jz)*w0(jx,jy,jz,2)
      gs(jx,jy,jz)=vcrbx-etaf(jx,jy,jz)*w0(jx,jy,jz,1)
      hs(jx,jy,jz)=0.
    7 continue
      else
      if(m.eq.8) then
      do 8 jz=1,nnz
      do 8 jy=1,nny
      do 8 jx=1,nnx
!
! [4] Energy eq.
!
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      bdotv=(x(jx,jy,jz,5)*x(jx,jy,jz,2)+x(jx,jy,jz,6)*x(jx,jy,jz,3)&
            +x(jx,jy,jz,7)*x(jx,jy,jz,4))/x(jx,jy,jz,1)
      eng=x(jx,jy,jz,8)+pr(jx,jy,jz)+.5*b2
!
      fs(jx,jy,jz)=eng*x(jx,jy,jz,2)/x(jx,jy,jz,1)&
                 -bdotv*x(jx,jy,jz,5)&
                 +etaf(jx,jy,jz)*(w0(jx,jy,jz,2)*x(jx,jy,jz,7)&
                -w0(jx,jy,jz,3)*x(jx,jy,jz,6))
      gs(jx,jy,jz)=eng*x(jx,jy,jz,3)/x(jx,jy,jz,1)&
                 -bdotv*x(jx,jy,jz,6)&
                 +etaf(jx,jy,jz)*(w0(jx,jy,jz,3)*x(jx,jy,jz,5)&
                 -w0(jx,jy,jz,1)*x(jx,jy,jz,7))
      hs(jx,jy,jz)=eng*x(jx,jy,jz,4)/x(jx,jy,jz,1)&
                 -bdotv*x(jx,jy,jz,7)&
                 +etaf(jx,jy,jz)*(w0(jx,jy,jz,1)*x(jx,jy,jz,6)&
                 -w0(jx,jy,jz,2)*x(jx,jy,jz,5))
!
    8 continue
      else
      endif
      endif
      endif
      endif
      endif
      endif
      endif
      endif
!
      return
      end
!
!
      subroutine current(x,mm)
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8)
!
!
      if(mm.eq.1) then
      do 10 jz=2,nz
      do 10 jy=2,ny
      do 10 jx=2,nx
      w0(jx,jy,jz,1)=.5*(x(jx,jy+1,jz,7)-x(jx,jy-1,jz,7))/dy&
                 -.5*(x(jx,jy,jz+1,6)-x(jx,jy,jz-1,6))/dz
      w0(jx,jy,jz,2)=-.5*(x(jx+1,jy,jz,7)-x(jx-1,jy,jz,7))/dx&
                  +.5*(x(jx,jy,jz+1,5)-x(jx,jy,jz-1,5))/dz
      w0(jx,jy,jz,3)=.5*(x(jx+1,jy,jz,6)-x(jx-1,jy,jz,6))/dx&
                -.5*(x(jx,jy+1,jz,5)-x(jx,jy-1,jz,5))/dy
   10 continue
! boundary at jx=1,nxp1
!
      do 20 m=1,3
      do 20 jz=2,nz
      do 20 jy=2,ny
      w0(1,jy,jz,m)=2.*w0(2,jy,jz,m)-w0(3,jy,jz,m)
      w0(mx,jy,jz,m)=2.*w0(nx,jy,jz,m)-w0(nx-1,jy,jz,m)
   20 continue
!
      if(halfx) then
      do 25 jz=2,nz
      do 25 jy=2,ny
      w0(mx,jy,jz,1)=w0(nx-1,my-jy+1,jz,1)
      w0(mx,jy,jz,2)=w0(nx-1,my-jy+1,jz,2)
      w0(mx,jy,jz,3)=-w0(nx-1,my-jy+1,jz,3)
   25 continue
      endif
!
! boundary at jy=1,nyp1
      if(periody) then
      do 30 m=1,3
      do 30 jz=2,nz
      do 30 jx=1,mx
      w0(jx,1,jz,m)=w0(jx,ny,jz,m)
      w0(jx,nyp1,jz,m)=w0(jx,2,jz,m)
   30 continue
      else
      do 35 m=1,3
      do 35 jz=2,nz
      do 35 jx=1,mx
      w0(jx,1,jz,m)=w0(jx,2,jz,m)
      w0(jx,nyp1,jz,m)=w0(jx,ny,jz,m)
   35 continue
      endif
! b.!. at jz=1
      do 40 m=1,3
      do 40 jy=1,nyp1
      do 40 jx=1,nxp1
      w0(jx,jy,1,m)=w0(jx,jy,2,m)
      w0(jx,jy,mz,m)=w0(jx,jy,nz,m)
   40 continue
! b.!. at jz=nzp1
      if(halfz) then
      do 50 jy=1,nyp1
      do 50 jx=1,nxp1
      w0(jx,jy,nzp1,1)=w0(jx,jy,nz-1,1)
      w0(jx,jy,nzp1,2)=w0(jx,jy,nz-1,2)
      w0(jx,jy,nzp1,3)=-w0(jx,jy,nz-1,3)
   50 continue
      endif
      else
      do 60 jz=2,nz-1
      do 60 jy=2,ny-1
      do 60 jx=2,nx-1
      w0(jx,jy,jz,1)=.5*(x(jx,jy+1,jz,7)-x(jx,jy-1,jz,7))/dy&
                 -.5*(x(jx,jy,jz+1,6)-x(jx,jy,jz-1,6))/dz
      w0(jx,jy,jz,2)=-.5*(x(jx+1,jy,jz,7)-x(jx-1,jy,jz,7))/dx&
                 +.5*(x(jx,jy,jz+1,5)-x(jx,jy,jz-1,5))/dz
      w0(jx,jy,jz,3)=.5*(x(jx+1,jy,jz,6)-x(jx-1,jy,jz,6))/dx&
                -.5*(x(jx,jy+1,jz,5)-x(jx,jy-1,jz,5))/dy
   60 continue
! boundary at jx=1,nx
!
      do 70 m=1,3
      do 70 jz=2,nz-1
      do 70 jy=2,ny-1
      w0(1,jy,jz,m)=w0(2,jy,jz,m)
      w0(nx,jy,jz,m)=w0(nx-1,jy,jz,m)
   70 continue
!
      if(halfx) then
      do 75 jz=2,nz
      do 75 jy=2,ny
      w0(nx,jy,jz,1)=w0(nx-1,ny-jy+1,jz,1)
      w0(nx,jy,jz,2)=w0(nx-1,ny-jy+1,jz,2)
      w0(nx,jy,jz,3)=-w0(nx-1,ny-jy+1,jz,3)
   75 continue
      endif
!
! boundary at jy=1,ny
      if(periody) then
      do 80 m=1,3
      do 80 jz=2,nz-1
      do 80 jx=1,nx
      w0(jx,1,jz,m)=w0(jx,ny-1,jz,m)
      w0(jx,ny,jz,m)=w0(jx,2,jz,m)
   80 continue
      else
      do 85 m=1,3
      do 85 jz=2,nz-1
      do 85 jx=1,nx
      w0(jx,1,jz,m)=w0(jx,2,jz,m)
      w0(jx,ny,jz,m)=w0(jx,ny-1,jz,m)
   85 continue
      endif
      do 90 m=1,3
      do 90 jy=1,ny
      do 90 jx=1,nx
      w0(jx,jy,1,m)=w0(jx,jy,2,m)
      w0(jx,jy,nz,m)=w0(jx,jy,nz-1,m)
   90 continue
! b.!. at jz=nz
      if(halfz) then
      do 100 jy=1,ny
      do 100 jx=1,nx
      w0(jx,jy,nz,1)=w0(jx,jy,nz-1,1)
      w0(jx,jy,nz,2)=w0(jx,jy,nz-1,2)
      w0(jx,jy,nz,3)=-w0(jx,jy,nz-1,3)
  100 continue
      endif
      endif
!
      return
      end
!                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! vorticity 没有被改写
      subroutine vorticity(x)                    !!! 被 facur 子程序 调用
! Calculate vorticity form Velocity variables
! w0=\nabla \times \vec{V}
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8),xm(mx,my,mz,8)
!
!
      do 5 jz=1,mz
      do 5 jy=1,my
      do 5 jx=1,mx
      xm(jx,jy,jz,5)=x(jx,jy,jz,2)/x(jx,jy,jz,1) ! Stands for Vx
      xm(jx,jy,jz,6)=x(jx,jy,jz,3)/x(jx,jy,jz,1) ! Vy
      xm(jx,jy,jz,7)=x(jx,jy,jz,4)/x(jx,jy,jz,1) ! Vz
    5 continue
!
      do 10 jz=2,nz
      do 10 jy=2,ny
      do 10 jx=2,nx
      w0(jx,jy,jz,1)=.5*(xm(jx,jy+1,jz,7)-xm(jx,jy-1,jz,7))/dy&
                 -.5*(xm(jx,jy,jz+1,6)-xm(jx,jy,jz-1,6))/dz
      w0(jx,jy,jz,2)=-.5*(xm(jx+1,jy,jz,7)-xm(jx-1,jy,jz,7))/dx&
                  +.5*(xm(jx,jy,jz+1,5)-xm(jx,jy,jz-1,5))/dz
      w0(jx,jy,jz,3)=.5*(xm(jx+1,jy,jz,6)-xm(jx-1,jy,jz,6))/dx&
               -.5*(xm(jx,jy+1,jz,5)-xm(jx,jy-1,jz,5))/dy
   10 continue
! boundary at jx=1,nxp1
!
      do 20 m=1,3
      do 20 jz=2,nz
      do 20 jy=2,ny
      w0(1,jy,jz,m)=2.*w0(2,jy,jz,m)-w0(3,jy,jz,m)
      w0(mx,jy,jz,m)=2.*w0(nx,jy,jz,m)-w0(nx-1,jy,jz,m)
   20 continue
!
      if(halfx) then
      do 25 jz=2,nz
      do 25 jy=2,ny
      w0(mx,jy,jz,1)=w0(nx-1,my-jy+1,jz,1)
      w0(mx,jy,jz,2)=w0(nx-1,my-jy+1,jz,2)
      w0(mx,jy,jz,3)=-w0(nx-1,my-jy+1,jz,3)
   25 continue
      endif
!
! boundary at jy=1,nyp1
      if(periody) then
      do 30 m=1,3
      do 30 jz=2,nz
      do 30 jx=1,mx
      w0(jx,1,jz,m)=w0(jx,ny,jz,m)
      w0(jx,nyp1,jz,m)=w0(jx,2,jz,m)
   30 continue
      else
      do 35 m=1,3
      do 35 jz=2,nz
      do 35 jx=1,mx
      w0(jx,1,jz,m)=w0(jx,2,jz,m)
      w0(jx,nyp1,jz,m)=w0(jx,ny,jz,m)
   35 continue
      endif
! b.!. at jz=1
      do 40 m=1,3
      do 40 jy=1,nyp1
      do 40 jx=1,nxp1
      w0(jx,jy,1,m)=w0(jx,jy,2,m)
      w0(jx,jy,mz,m)=w0(jx,jy,nz,m)
   40 continue
! b.!. at jz=nzp1
      if(halfz) then
      do 50 jy=1,nyp1
      do 50 jx=1,nxp1
      w0(jx,jy,nzp1,1)=-w0(jx,jy,nz-1,1)
      w0(jx,jy,nzp1,2)=-w0(jx,jy,nz-1,2)
      w0(jx,jy,nzp1,3)=w0(jx,jy,nz-1,3)
   50 continue
      endif
!
      return
      end
!                             !!!!!!!!!!!!!!! readin 没有被改写
!      subroutine readin
!      include 'ma3ds1.for'
!      include 'ma3ds2.for'
!      character*8 contin
!	character*3 cn

!      contin='m3d'//cn(nst)
 !     open(unit=8,file=contin,status="unknown",form="unformatted")
!	open(unit=8,file='continue',status="unknown",form="unformatted")
!      read(8)ncase,nstep,time,nst
!      read(8)x
 !     close(8)
!	nst=1
!      return
!      end
!
!cyg-------------------------------
      subroutine readin(nst3,cont2,dtime1)              !!!!!!!!!  readin  没有被改写
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      character*8 contin
	  character*3 cn
	  nst=nst3

	if (cont2.eq.1) then
	open(unit=8,file='continue',status="unknown",form="unformatted")
	read(8)ncase,nstep,time,nst
	read(8)x
      close(8)

	else 

      contin='m3d'//cn(nst)
      open(unit=8,file=contin,status="unknown",form="unformatted")
      read(8)ncase,nstep,time
      read(8)x
	nst=ceiling(time)/dtime1+1
      close(8)
	end if

      return
      end

!cyg----------------------
      subroutine recrd
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      character*8 output
      character*3 cn
      output='m3d'//cn(nst)
      open(unit=8,file=output,status="unknown",form="unformatted")
      write(8)ncase,nstep,time
      write(8)x
      close(8)
      return
      end
!
      subroutine recrd1
      include 'ma3ds1.for'
      include 'ma3ds2.for'
      character*8 output
      character*3 cn                     ! It's a function!!!!
      output='m3ds'//cn(nst)
      call current(x,1)
      open(unit=8,file=output,status="unknown",form="formatted")
      write(8,9)((((x(jx,jy,jz,m),m=1,8),(w0(jx,jy,jz,m),m=1,3),jx=1,mx),jy=1,my),jz=1,mz)
    9 format(11(1x,e11.4))
      close(8)
      return
      end
!
!
      character*3 function cn(n)                     ! 把数字变成字母
!
!-----assume that n is no greater than 999
!
!
!-----separate the digits
!
      n1=n/100
      n2=(n-100*n1)/10
      n3=n-100*n1-10*n2
!
!-----stick together cn using char function
!
      n1=n1+48
      n2=n2+48
      n3=n3+48
      cn(1:1)=char(n1)
      cn(2:2)=char(n2)
      cn(3:3)=char(n3)
!
      return
      end
!
!
      subroutine gridpnt
      include 'ma3ds1.for'
      dimension  xxh(mx),yyh(my),zzh(mz)
      nnxh=(mx+1)/2
      nnx=2*(nnxh-1)
      nnxq=nnxh/2
      nnyh=(my+1)/2
      nny=2*(nnyh-1)
      nnyq=nnyh/2
      nnzh=(mz+1)/2
      nnz=2*(nnzh-1)
      nnzq=nnzh/2
      dx=(xmax-xmin)/float(nnx)
      do 1 jx=1,mx
      xx(jx)=xmin+(jx-1)*dx
    1 continue
!
      dy=(ymax-ymin)/float(nny)
      do 2 jy=1,my
      yy(jy)=ymin+(jy-1)*dy
    2 continue
!
      dz=(zmax-zmin)/float(nnz)
      do 3 jz=1,mz
      zz(jz)=zmin+(jz-1)*dz
    3 continue
! 
      return
      end
!
      subroutine smthxyz(x,lsmth,kk)     ! lsmth seemed unused!!!
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8)
      do 10 k=1,kk
      do 11 m=1,8
      do 12 jz=2,nz
      do 12 jy=2,ny
      do 12 jx=2,nx
      w0(jx,jy,jz,1)=((x(jx+1,jy,jz,m)+x(jx-1,jy,jz,m))&
                    +(x(jx,jy+1,jz,m)+x(jx,jy-1,jz,m))&
              +(x(jx,jy,jz+1,m)+x(jx,jy,jz-1,m))-6.d0*x(jx,jy,jz,m))
   12 continue
      do 13 jz=2,nz
      do 13 jy=2,ny
      do 13 jx=2,nsmthx
      theta=3.1415926d0*(jx-2)/(nsmthx-3)
      x(jx,jy,jz,m)=x(jx,jy,jz,m)+(1.d0/96.0d0)*(.5d0*(1.d0+cos(theta)))*w0(jx,jy,jz,1)
   13 continue
      if(.not.halfx) then
      do 23 jz=2,nz
      do 23 jy=2,ny
      do 23 jx=mx-nsmthx+1,nx
      theta=3.1415926d0*(mx-jx-1)/(nsmthx-3)
      x(jx,jy,jz,m)=x(jx,jy,jz,m)+(1.d0/96.0d0)*(.5d0*(1.d0+cos(theta)))*w0(jx,jy,jz,1)
   23 continue
      endif
      if(.not.periody) then
      do 33 jz=2,nz
      do 33 jy=2,nsmthy
      do 33 jx=2,nx
      theta=3.1415926d0*(jy-2)/(nsmthy-3)
      x(jx,jy,jz,m)=x(jx,jy,jz,m)+(1.d0/96.0d0)*(.5d0*(1.d0+cos(theta)))*w0(jx,jy,jz,1)
   33 continue
      do 35 jz=2,nz 
      do 35 jy=nyp1-nsmthy+1,ny
      do 35 jx=2,nx 
      theta=3.1415926d0*(nyp1-jy-1)/(nsmthy-3)
      x(jx,jy,jz,m)=x(jx,jy,jz,m)+(1.d0/96.0d0)*(.5d0*(1.d0+cos(theta)))*w0(jx,jy,jz,1)
   35 continue
      endif
      do 43 jz=2,nz
      do 43 jy=2,ny
      do 43 jx=2,nx
!      theta=2.*3.1415926*zz(jz)/zmin
      x(jx,jy,jz,m)=x(jx,jy,jz,m) +(1.d0/48.0d0)*w0(jx,jy,jz,1)
!     1     +(1./48.0)*(2.+cos(theta))/3.*w0(jx,jy,jz,1)
   43 continue
!      if(.not.halfz) then
!      do 53 jz=nzp1-nsmthz+1,nz
!      do 53 jy=2,ny
!      do 53 jx=2,nx
!      theta=3.1415926*(nzp1-jz-1)/(nsmthz-3)
!      x(jx,jy,jz,m)=x(jx,jy,jz,m)
!     1     +(1./96.0)*(.5*(1.+cos(theta)))*w0(jx,jy,jz,1)
!   53 continue
!      endif
   11 continue
      call bndry(x,2)
!
   10 continue
!
      return
      end
!
      subroutine avrg1(x,caf1)
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8)
      do 11 m=1,8
      if(m.eq.5.or.m.eq.6.or.m.eq.7) goto 11
      do 12 jz=2,nz
      do 12 jy=2,ny
      do 12 jx=2,nx
      w0(jx,jy,jz,1)=caf1*x(jx,jy,jz,m)+(1-caf1)*&
              ((x(jx+1,jy,jz,m)+x(jx-1,jy,jz,m))&
              +(x(jx,jy+1,jz,m)+x(jx,jy-1,jz,m))&
              +(x(jx,jy,jz+1,m)+x(jx,jy,jz-1,m)))/6.0d0
   12 continue
      do 13 jz=2,nz
      do 13 jy=2,ny
      do 13 jx=2,nx
      x(jx,jy,jz,m)=w0(jx,jy,jz,1)
   13 continue
   11 continue
      call bndry(x,2)
!
      return
      end
!
      subroutine avrg2(x,caf1)
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8)
      do 11 m=5,7
      do 12 jz=2,nz
      do 12 jy=2,ny
      do 12 jx=2,nx
      w0(jx,jy,jz,1)=caf1*x(jx,jy,jz,m)+(1-caf1)*&
              ((x(jx+1,jy,jz,m)+x(jx-1,jy,jz,m))&
              +(x(jx,jy+1,jz,m)+x(jx,jy-1,jz,m))&
              +(x(jx,jy,jz+1,m)+x(jx,jy,jz-1,m)))/6.0d0
   12 continue
      do 13 jz=2,nz
      do 13 jy=2,ny
      do 13 jx=2,nx
      x(jx,jy,jz,m)=w0(jx,jy,jz,1)
   13 continue
   11 continue
      call bndry(x,2)
!
      return
      end
!
      subroutine smthf(x,caf1)                  !   这个还没改写！！！！！！
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8), wh(mx,my,mz,3)
!
      do 11 m=1,8
      do 12 jz=2,nz
      do 12 jx=2,nx
      do 12 jy=2,ny
      wh(jx,jy,jz,1)=(x(jx+1,jy,jz,m)-x(jx,jy,jz,m))*&
                    (x(jx,jy,jz,m)-x(jx-1,jy,jz,m))
      wh(jx,jy,jz,2)=(x(jx,jy+1,jz,m)-x(jx,jy,jz,m))*&
                    (x(jx,jy,jz,m)-x(jx,jy-1,jz,m))
      wh(jx,jy,jz,3)=(x(jx,jy,jz+1,m)-x(jx,jy,jz,m))*&
                    (x(jx,jy,jz,m)-x(jx,jy,jz-1,m))
   12 continue
      call bndry1(wh(1,1,1,1),0)
      call bndry1(wh(1,1,1,2),0)
      call bndry1(wh(1,1,1,3),0)
      do 13 jz=2,nz
      do 13 jx=2,nx
      do 13 jy=2,ny
         if((wh(jx,jy,jz,1).lt.0.and.&
          (wh(jx+1,jy,jz,1).lt.0.or.wh(jx-1,jy,jz,1).lt.0))&
      .or.(wh(jx,jy,jz,2).lt.0.and.&
          (wh(jx,jy+1,jz,2).lt.0.or.wh(jx,jy-1,jz,2).lt.0))&
      .or.(wh(jx,jy,jz,3).lt.0.and.&
          (wh(jx,jy,jz+1,3).lt.0.or.wh(jx,jy,jz-1,3).lt.0))) then
      w0(jx,jy,jz,1)=caf1*x(jx,jy,jz,m)+(1-caf1)*&
              ((x(jx+1,jy,jz,m)+x(jx-1,jy,jz,m))&
              +(x(jx,jy+1,jz,m)+x(jx,jy-1,jz,m))&
             +(x(jx,jy,jz+1,m)+x(jx,jy,jz-1,m)))/6.0
      else
      w0(jx,jy,jz,1)=x(jx,jy,jz,m)
      endif
   13 continue
      do 14 jz=2,mz-1
      do 14 jx=2,nx
      do 14 jy=2,ny
      x(jx,jy,jz,m)=w0(jx,jy,jz,1)
   14 continue
   11 continue
      call bndry(x,2)
      return
      end
!
      subroutine pressure(x,mm)
      include 'ma3ds1.for'
      dimension x(mx,my,mz,8)
      if(mm.eq.1) then
      do 1 jz=1,mz
      do 1 jy=1,my
      do 1 jx=1,mx
      rv2=(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2+x(jx,jy,jz,4)**2)&
          /x(jx,jy,jz,1)
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      pr(jx,jy,jz)=(gamma-1.)*(x(jx,jy,jz,8)-.5*rv2-.5*b2)
    1 continue
      else
      do 2 jz=1,nz
      do 2 jy=1,ny
      do 2 jx=1,nx
      rv2=(x(jx,jy,jz,2)**2+x(jx,jy,jz,3)**2+x(jx,jy,jz,4)**2)&
          /x(jx,jy,jz,1)
      b2=x(jx,jy,jz,5)**2+x(jx,jy,jz,6)**2+x(jx,jy,jz,7)**2
      pr(jx,jy,jz)=(gamma-1.)*(x(jx,jy,jz,8)-.5*rv2-.5*b2)
    2 continue
      endif
      call positive(pr,1.d-5)
      return
      end
!
!
      subroutine positive(fn,c)
      include 'ma3ds1.for'
!	include 'ma3ds2.for'
      dimension fn(mx,my,mz)
	character*15 out
      do 1 jz=1,mz
      do 1 jy=1,my
      do 1 jx=1,mx

	if (fn(jx,jy,jz).lt.0.0)then
	out='finaltime.txt'
	open(unit=8,file=out,status="unknown",form="formatted")
	write(8,20)time
  20	format(9(1x,e11.4))
  200	format('finaltime=',i5)
!	write(*,*)'finaltime=',time
!	stop
	endif

      if(fn(jx,jy,jz).lt.c) then
      fn(jx,jy,jz)=c
      end if
    1 continue
      return
      end  
          
!
      subroutine foreta(t,mm)
! Set Etaf
      include 'ma3ds1.for'
      dimension etam(my)
!
      etac=0.5d0
      
!cjx  attentiion      
      etal=0.000d0
      etab=0.005d0
      alpha0=2.0d0
!cjx  attentiion      
      xlen=15.0d0
      xtrig=0.0d0
!cjx  attentiion      
      ztrig=0.0d0
!cjx  attentiion      
      awx=aw
      awz=2.*awx
!
      if(mm.eq.1) then
      do 1 jy=1,my
      if(abs(yy(jy)).le.xlen) then
      etam(jy)=1.
      else
      etam(jy)=(1.-tanh(abs(yy(jy))-xlen)**2)
      endif
    1 continue
!
      do 2 jz=1,mz
      do 2 jy=1,my
      do 2 jx=1,mx
! nonlinear resistivity
!
! localized resistivity perturbation
      if(abs(xx(jx)-xtrig).le.0.5) then
      etax=1.0d0
      else
      etax=1.-tanh((xx(jx)-xtrig)/awx)**2
      endif
      if(abs(abs(zz(jz))-ztrig).le.1.) then
      etaz=1.0d0
      else
      etaz=1.-tanh((abs(zz(jz))-ztrig)/awz)**2
      endif
      etaf(jx,jy,jz)=etab+etam(jy)*etal*etax*etaz
    2 continue
      else                         ! if mm .neq. 1
      do 3 jy=1,ny
      if(abs(yy(jy)+dy/2.).le.xlen) then
      etam(jy)=1.d0
      else
      etam(jy)=(1.-tanh(abs(yy(jy)+dy/2.)-xlen)**2)
      endif
    3 continue
!
      do 4 jz=1,nz
      do 4 jy=1,ny
      do 4 jx=1,nx
! localized resistivity perturbation
      if(abs(xx(jx)+dx/2.-xtrig).le.0.5) then
      etax=1.0d0
      else
      etax=1.-tanh((xx(jx)+dx/2.-xtrig)/awx)**2
      endif
      if(abs(abs(zz(jz)+dz/2.)-ztrig).le.1.) then
      etaz=1.0d0
      else
      etaz=1.-tanh((zz(jz)+dz/2.-ztrig)/awz)**2             !??????? 确定不需要 tanh((abs(zz(jz)+dz/2.)-ztrig)/aw2) 因为ztrig=0， 加不加abs都一样
      endif
      etaf(jx,jy,jz)=etab+etam(jy)*etal*etax*etaz
    4 continue
      endif
!
      return
      end

!cyg---------------------------------------                    !!!!!!!!!!!!!! 没有改写这个程序
	subroutine incident_plasma(x,xi,t,t1,Io)  !Io代表是否是初始，初始则为0,
	                                  !不是则为1
	include 'ma3ds1.for'
!	include 'ma3ds2.for'
	dimension x(mx,my,mz,8)
	dimension xi(mx,my,mz,8)
!	jxh=nxp1/2+1
!	real t
!	real time
!      include 'ma3ds2.for'
!	include 'ma3ds3.for'

	vs0=0.6
	lsz=1.  !剪切的半宽度
	lsx=1.
	lsy=1.   !两端随y衰减因子
	tao=100   !随时间变化因子
	t0=t1
    	pi=3.1415926
    xx0=0.
    xx1=xx0
    
	yy0=0.
    yy1=-yy0
    
    zz0=2.
    zz1=-zz0
    
    theta_in=45
    theta_in=theta_in*pi/180
    
    phi_in=0
    phi_in=phi_in*pi/180
    
    vx0_in=vs0*sin(theta_in)*cos(phi_in)
    vy0_in=vs0*sin(theta_in)*sin(phi_in)
    vz0_in=vs0*cos(theta_in)
    
    vx1_in=vx0_in
    vy1_in=vy0_in
    vz1_in=-vz0_in
    
    
	
	if (Io.eq.0) then
	a=1     !a初始值前面系数，初始值时则导数部分为零
	b=0     
	end if

	if (Io.eq.1) then
	a=1     !a初始值前面系数，初始值后时则初始值部分为零
	b=0
	end if

	if (t.gt.40) then
	do jz=1,mz
           do jy=1,my
	do jx=1,mx
	
	
!	write(*,*)'*****1'
	   x(jx,jy,jz,2)=0.0*x(jx,jy,jz,2) &
		  +x(jx,jy,jz,1) &
             *(vx0_in*exp(-(((xx(jx)-xx0)/lsx)**2+((yy(jy)-yy0)/lsy)**2+((zz(jz)-zz0)/lsz)**2)) &
	+ vx1_in*exp(-(((xx(jx)-xx1)/lsx)**2+((yy(jy)-yy1)/lsy)**2+((zz(jz)-zz1)/lsz)**2))) &
	
	*(a*1.0+b*1.0/tao*1.0/cosh(t/tao-t0/tao)**2)
          
          
          	  x(jx,jy,jz,3)=0.0*x(jx,jy,jz,3) &
		  +x(jx,jy,jz,1) &
             *(vy0_in*exp(-(((xx(jx)-xx0)/lsx)**2+((yy(jy)-yy0)/lsy)**2+((zz(jz)-zz0)/lsz)**2)) &
	+ vy1_in*exp(-(((xx(jx)-xx1)/lsx)**2+((yy(jy)-yy1)/lsy)**2+((zz(jz)-zz1)/lsz)**2))) &
	*(a*1.0+b*1.0/tao*1.0/cosh(t/tao-t0/tao)**2)
          
          	  x(jx,jy,jz,4)=0.0*x(jx,jy,jz,4) &
		  +x(jx,jy,jz,1) &
             *(vz0_in*exp(-(((xx(jx)-xx0)/lsx)**2+((yy(jy)-yy0)/lsy)**2+((zz(jz)-zz0)/lsz)**2)) &
	+ vz1_in*exp(-(((xx(jx)-xx1)/lsx)**2+((yy(jy)-yy1)/lsy)**2+((zz(jz)-zz1)/lsz)**2))) &
	*(a*1.0+b*1.0/tao*1.0/cosh(t/tao-t0/tao)**2)

	end do 

	ix_temp=mx/2+1
	if (x(ix_temp,jy,jz,5).lt.(-bm0*3*xx(ix_temp)*zz(jz)/(xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)**(5/2))) then

	 x(ix_temp,jy,jz,5)  =x(ix_temp,jy,jz,5)+tanh((t-40)/tao)* &
(-bm0*3*xx(ix_temp)*zz(jz)/(xx(1)**2+yy(jy)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,5))

	  if  ((xx(1)**2+yy(jy)**2+zz(jz)**2)==0)  then
	  x(ix_temp,jy,jz,5)  =x(ix_temp,jy,jz,5)+tanh((t-40)/tao)* &
(-bm0*3*xx(ix_temp)*zz(jz)/(xx(ix_temp)**2+yy(jy-1)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,5))
	  end if


	end if


	if (x(ix_temp,jy,jz,6).lt.-bm0*3*yy(jy)*zz(jz)/(xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)**(5/2)) then
	x(ix_temp,jy,jz,6)=x(ix_temp,jy,jz,6)+tanh((t-40)/tao)* &
(-bm0*tanh((t-60)/40)*3*yy(jy)*zz(jz)/(xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,6))

	if  ((xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)==0)  then

	x(ix_temp,jy,jz,6)=x(ix_temp,jy,jz,6)+tanh((t-40)/tao)* &
(-bm0*tanh((t-60)/40)*3*yy(jy-1)*zz(jz)/(xx(ix_temp)**2+yy(jy-1)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,6))

	end if

	end if


	if (x(ix_temp,jy,jz,7).lt.(-bm0*(2*zz(jz)**2-xx(ix_temp)**2-yy(jy)**2)/(xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)**(5/2))) then
	x(ix_temp,jy,jz,7)  =x(ix_temp,jy,jz,7)+tanh((t-40)/tao)* &
((2*zz(jz)**2-xx(ix_temp)**2-yy(jy)**2)/(xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,7))

	if  ((xx(ix_temp)**2+yy(jy)**2+zz(jz)**2)==0)  then

	x(ix_temp,jy,jz,7)=x(ix_temp,jy,jz,7)+tanh((t-40)/tao)* &
((2*zz(jz)**2-xx(1)**2-yy(jy-1)**2)/(xx(ix_temp)**2+yy(jy-1)**2+zz(jz)**2)**(5/2)-x(ix_temp,jy,jz,7))

	end if
	



	end if	

             	
	end do
	end do
	end if
	return
	end 
!cyg-------------------------
