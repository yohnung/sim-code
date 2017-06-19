      parameter (mx=41,my=7,mz=41)
      implicit double precision (a-h,o-z)
      logical lrstrt,uniformx,uniformy,uniformz,periody, &
      midx,midy,midz,cbndry,halfx,halfy,halfz
      dimension w0(mx,my,mz,3),flx(my),fx(mx,8),fxi(mx,8)
      dimension pr(mx,my,mz),etaf(mx,my,mz)
      dimension xx(mx), yy(my),zz(mz)
      dimension nstp(250)
! cwm add: start
      double precision rhoinfinity, balcoeff, normlambda, fluc, kx, kz
      common /fluctuation/ rhoinfinity, balcoeff, normlambda, &
      fluc, kx, kz
! cwm add: end
      common /index/ nx, nxp1, ny, nyp1, nz2, nz2p1,nz,nzp1, nstep
      common /para/  xmax,xmin,ymax,ymin,zmax,zmin,dxmin,dymin,dzmin
      common /cst/   time, dx, dy, dz, dt, ds, di, vyi0
      common /phy1/  gamma, cj, dis, alpha0
      common /phy2/  betam, phi, epsilon, v0, bm0, bs0, pmsp
      common /num/   nsmthx, nsmthy, nsmthz
      common /run/   nend,nsp,ncase,nst,nst1,nint,npt,nstp, nk
      common /pert/  time0, time1, alamda, alpha, caf, apx, apy 
      common /var1/  lrstrt,uniformx,uniformy,uniformz,periody, &
      midx,midy,midz,cbndry,halfx,halfy,halfz
      common /var11/ xx, yy, zz, aw, xp, yp, zp
      common /var12/ flx,pr,etaf,fx,fxi
      common /wrk/   w0,fw0,fw1
      common /bf/    spsi,scur,spr,Bxi,Bzi
