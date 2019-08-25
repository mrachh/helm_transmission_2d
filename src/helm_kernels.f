      subroutine slp(srcinfo,targinfo,dpars,zpars,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk,h1

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
c
      subroutine sprime(srcinfo,targinfo,dpars,zpars,ipars,u)
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
      
      u = -(gx*targinfo(3) + gy*targinfo(4))

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
c
c
c        transmission matrices
c
c
      subroutine transmission_dir(srcinfo,targinfo,dpars,zpars,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1} + beta S_{k2} + gamma D_{k1} + delta D_{k2}
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk,h1,gx,gy
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      gx = -zs*zk*h1*(srcinfo(1)-targinfo(1))/rr
      gy = -zs*zk*h1*(srcinfo(2)-targinfo(2))/rr
      
      u = zpars(5)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(3)*zs*h0

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      gx = -zs*zk2*h1*(srcinfo(1)-targinfo(1))/rr
      gy = -zs*zk2*h1*(srcinfo(2)-targinfo(2))/rr
      
      u = u+zpars(6)*(gx*srcinfo(3) + gy*srcinfo(4)) + zpars(4)*zs*h0


      return
      end
c
c
c
c
c
      subroutine transmission_neu(srcinfo,targinfo,dpars,zpars,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1}' + beta S_{k2}' + gamma D_{k1}' + delta D_{k2}'
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars(*)
      real *8 dpars(*),srcinfo(*),targinfo(*)
      complex *16 zpars(*),u,h0,ima,zs,z,zk,h1,gx,gy,h2,zk2
      complex *16 d2gdx2,d2gdy2,d2gdxdy
      complex *16 gd0,gs0,gd1,gs1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)

      xd = targinfo(1) - srcinfo(1)
      yd = targinfo(2) - srcinfo(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk

      d2gdx2 = (-h1 + h2*xd*xd)*zk
      d2gdxdy = h2*xd*yd*zk
      d2gdy2 = (-h1+h2*yd*yd)*zk

      gd0 = -zs*(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(3) + srcinfo(4)*targinfo(4)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      gx = zs*zk*h1*(srcinfo(1)-targinfo(1))/rr
      gy = zs*zk*h1*(srcinfo(2)-targinfo(2))/rr
      
      gs0 = gx*targinfo(3) + gy*targinfo(4)


      gd0 = 0

      u = zpars(3)*gs0 + zpars(5)*gd0



      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk2

      d2gdx2 = (-h1 + h2*xd*xd)*zk2
      d2gdxdy = h2*xd*yd*zk2
      d2gdy2 = (-h1+h2*yd*yd)*zk2

      gd1 = -zs*(d2gdx2*srcinfo(3)*targinfo(3) +
     1    d2gdxdy*(srcinfo(3)*targinfo(3) + srcinfo(4)*targinfo(4)) + 
     2    d2gdy2*srcinfo(4)*targinfo(4))

      gx = zs*zk2*h1*(srcinfo(1)-targinfo(1))/rr
      gy = zs*zk2*h1*(srcinfo(2)-targinfo(2))/rr
      
      gs1 = gx*targinfo(3) + gy*targinfo(4)

      gd1 = 0

      u = u+zpars(4)*gs1 + zpars(6)*gd1



      return
      end
c
c
c
c
c
