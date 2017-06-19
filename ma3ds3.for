      lrstrt=.false.
      uniformx=.true.         ! not used temporarily
      uniformy=.true.         ! not used temporarily
      uniformz=.false.        ! not used temporarily
      periody=.true.  !.false.
!      midxmidx=.false.       ! not used temporarily
      midy=.true.             ! not used temporarily
      midz=.true.             ! not used temporarily
      halfx=.false.
      halfy=.false.           ! not used temporarily
      halfz=.false.       
      cbndry=.true.           ! not used temporarily
      ncase=1
      nend=20              ! maximum number of files
      nst=1
      nst1=1
      nint=1 
      dxmin=0.05d0         ! not used temporarily ! d0 is for double precision
      dymin=0.5d0          ! not used temporarily
      dzmin=0.5d0          ! not used temporarily
      xmin=-5.d0   ! -1.
      xmax=5.d0    ! 11.
      ymin=1.d0    ! -5.
      ymax=7.d0    ! 5.
      zmin=-5.d0 
      zmax=5.d0
      nsmthx=2*mx/3
      nsmthy=3*my/4
      nsmthz=2*mz/3
      if(halfz) nsmthz=mz-1
!      gatma=1.66667d0        ! not used temporarily ! should be gamma
      gamma=1.66667d0        ! This is added by cwm
      cj=1.d0                ! not used temporarily
      dis=1.d0               ! not used temporarily
      betam= 0.01d0
    
      di= 1.d0   !0.0d0          ! di is used in initializing rhoVy, and computing VcrossB in flux
                           ! togather with vyi0 
      phi=180.d0
 
 
      epsilon=0.d0
      bm0=1.d0
      bs0=1.d0
      v0=0.d0    
      vyi0=0.d0      
      apx=1.d0                ! not used temporarily
      apy=1.d0                ! not used temporarily
      aw=1.d0                ! width of rho
      alamda=0.01d0
      alpha=0.1d0             ! not used temporarily
      caf=0.997d0

! cwm add :start
! harris-current and fluctuation
      rhoinfinity=0.2d0
      balcoeff=1.d0
      normlambda=0.5d0
      fluc=0.2d0
      kx=4*3.1415926d0/xmax
      kz=2*3.1415926d0/zmax
! cwm add: end