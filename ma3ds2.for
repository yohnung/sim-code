      double precision x(mx,my,mz,8),xi(mx,my,mz,8)
      double precision xm(mx,my,mz,8)
      double precision fs(mx,my,mz), gs(mx,my,mz) ,hs(mx,my,mz)

      common /var21/ x,xi     ! x for variables, and xi for variables
                              ! evaluated at x+dx/2.
      common /var22/ xm
      common /var24/ fs
      common /var25/ gs
      common /var26/ hs
