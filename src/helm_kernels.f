      subroutine slp(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk

      zk = zpars(1)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      u = zs*h0
      

      return
      end
c
c
c
c
c
      subroutine dlp(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk,h1,gx,gy

      zk = zpars(1)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      gx = -zs*zk*h1*(srcinfo(1)-targinfo(1))/rr
      gy = -zs*zk*h1*(srcinfo(2)-targinfo(2))/rr
      
      u = gx*srcinfo(3) + gy*srcinfo(4)

      return
      end
c
c
c
c
c
      subroutine comb(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk,h1,gx,gy

      zk = zpars(1)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      gx = -zs*zk*h1*(srcinfo(1)-targinfo(1))/rr
      gy = -zs*zk*h1*(srcinfo(2)-targinfo(2))/rr
      
      u = zpars(3)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(2)*zs*h0


      return
      end
