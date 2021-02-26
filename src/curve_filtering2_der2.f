

      subroutine simple_curve_resampler_mem_der2(n,xy,nb,eps,nmax,
     1   nlarge,nout,lsave,lused,ier)
      implicit real *8 (a-h,o-z)
      integer n,nb
      real *8 xy(2,n)
      integer, intent(out) :: nlarge,nout

      real *8, allocatable :: x0(:),y0(:),work(:),ts(:)
      real *8, allocatable :: xver(:),yver(:)

      allocate(x0(n+1),y0(n+1))
      allocate(ts(n+1))
      allocate(xver(n+1),yver(n+1))

      do i=1,n
        x0(i) = xy(1,i)
        y0(i) = xy(2,i)
      enddo

      x0(n+1) = x0(1)
      y0(n+1) = y0(1)

      lused = 0
      lsave = 0
      ier = 0

      epscheck = max(eps,1.0d-12)
 
      do ii=1,nmax
        nlarge = 64*2**(ii)*n
        lenw = 10*nlarge*n + 10000
        allocate(work(lenw))
        curvelen = 0
        derr = 0
        nbl = 0
        iert = 0
        work = 0
        call rsblcurve(iert,x0,y0,n,nb,nlarge,
     1       curvelen,nbl,derr,ts,work,lenw,lsave,lused)
        nout = nlarge


        if(iert.ne.0) goto 1111
        erra = 0
        ra = 0
        do i=1,n
          xver(i) = 0
          yver(i) = 0
          call eval_curve_der2(ier,ts(i),work,xver(i),yver(i),dxt,dyt,
     1     d2xdt2,d2ydt2)
          erra = erra + (xver(i)-xy(1,i))**2
          ra = ra + xy(1,i)**2
          erra = erra + (yver(i)-xy(2,i))**2
          ra = ra + xy(2,i)**2
        enddo
        erra = sqrt(erra/ra)

        if(erra.lt.epscheck) goto 1000 

 1111   continue        
        deallocate(work)
      enddo

      if(iert.ne.0) ier = 4
      if(iert.eq.0) ier = 2
 1000 continue      


      return
      end
c
c
c
c
c
      subroutine simple_curve_resampler_guru_der2(n,xy,nb,
     1   nlarge,lsave,lused,nout,srcinfoout,hout,curvelen,wsave,ts,ier)
      implicit real *8 (a-h,o-z)
      integer n,nout,lsave,lused
      real *8 xy(2,n),srcinfoout(6,nout),curvelen
      real *8 ts(n+1)
      real *8 wsave(lsave)
      real *8, allocatable :: x0(:),y0(:),work(:)

      allocate(x0(n+1),y0(n+1))

      do i=1,n
        x0(i) = xy(1,i)
        y0(i) = xy(2,i)
      enddo

      x0(n+1) = x0(1)
      y0(n+1) = y0(1)
      lenw = lused + 1000
      allocate(work(lenw))
      nbl = 0
      ier = 0

      
      call rsblcurve(ier,x0,y0,n,nb,nlarge,curvelen,nbl,
     1   derr,ts,work,lenw,lsave0,lused)

      
      wsave(1:lsave) = work(1:lsave)

      hout = curvelen/(nout+0.0d0)
      do i=1,nout
        t = (i-1)*hout
        call eval_curve_der2(ier,t,wsave,srcinfoout(1,i),
     1    srcinfoout(2,i),srcinfoout(3,i),srcinfoout(4,i),
     2    srcinfoout(5,i),srcinfoout(6,i))
      enddo

      return
      end
c
c
c
c
c
      subroutine eval_curve_multi_der2(n,ts,lsave,wsave,binfo)
      implicit real *8 (a-h,o-z)
      integer n,lsave
      real *8 ts(n),wsave(lsave),binfo(6,n)

      do i=1,n
        call eval_curve_der2(ier,ts(i),wsave,binfo(1,i),binfo(2,i),
     1   binfo(3,i),binfo(4,i),binfo(5,i),binfo(6,i))
      enddo

      return
      end
c
c
c
c
c
        subroutine eval_curve_der2(ier,t,w,x,y,dxdt,dydt,d2xdt2,
     1     d2ydt2)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1),w2(1)
        external funcurv_fast_der2
c
        nlarge=w(1)
        iw=w(2)
        iw2=w(3)
        its=w(4)
        iwr=w(5)
        iwv=w(6)
        nints=w(7)
        h=w(8)
        eps=w(9)
        m=w(10)
        ixgs=w(11)
        iwhts=w(12)
        curvelen = w(iw2+3-1)
c
        call anaptbl_der2(ier,t,nlarge,h,w(its),funcurv_fast_der2,
     1       w(iw),w(iw2),
     1       eps,m,w(ixgs),w(iwhts),w(iwr),w(iwv),nints,curvelen,tout,
     2       x,y,dxdt,dydt,d2xdt2,d2ydt2)
c
        return
c
        end
c
c
c
c
c
        subroutine funcurv_fast_der2(t,w,w2,x,y,dxdt,dydt,d2xdt2,
     1     d2ydt2)
c
        implicit real *8 (a-h,o-z)
        real *8 w(1),w2(1),tang(10),curv,z(10),curv0,der22(2)
        real *8 d2xdt2,d2ydt2
c
        n0=w2(1)
        nlarge=w2(2)
        curvelen=w2(3)
        tn1=w2(4)
        iw=5
c
c       undo the shift introduced in pertgauss that
c       ensures that t=0 corresponds to the x0(1),y0(1)
c       where x0,y0 is the user input data to rsblcurve
c
        t2=t+tn1*curvelen
c
c       evaluate the curve
c
        curv = 0
        call rsrespnt(t2,z,nlarge,tang,der22,w)
c
        x=z(1)
        y=z(2)
        dxdt=tang(1)
        dydt=tang(2)

        dtn = dxdt**2 + dydt**2
        dtn = sqrt(dtn)
        dxdt = dxdt/dtn
        dydt = dydt/dtn
        d2xdt2 = der22(1)
        d2ydt2 = der22(2)

c
c       evaluate the guassian perturbations
c
        call rsrespnt(t2,z,nlarge,tang,der22,w2(iw))


        x=x+z(1)
        y=y+z(2)
        dxdt=dxdt+tang(1)
        dydt=dydt+tang(2)

        d2xdt2 = d2xdt2 + der22(1)
        d2ydt2 = d2ydt2 + der22(2)

c
        return
c
        end
c
c
c
c 
c 
         subroutine anaptbl_der2(ier,x,n,h,ts,funcurve,par1,par2,eps,
     1       mm,xgs,whts,wright,wvals,nints,rl,tout,
     2       xout,yout,dxdtout,dydtout,d2xdt2out,d2ydt2out)
         implicit real *8 (a-h,o-z)
         real *8 ts(1),xgs(1),whts(1),wright(1),wvals(1)
         real *8 tz(10),xz(10),cz(10)
         external funcurve
c 
c       This subroutine finds the location (both on the
c       parametrizing interval [0,rl] and in R^2) of the
c       point specified by the user on the curve. The
c       curve has to have been preprocessed by a preceding
c       call to the subroutine anarsbl (see), and the point
c       is specified by its distance (in terms of arc-length)
c       from the beginning of the curve. PLEASE NOTE THAT THIS
C       SUBROUTINE DEPENDS ON THE SUBROUTINE ANARSBL FOR THE
C       PARAMETERS H,TS; THIS SUBROUTINE HAS NO FUNCTION AS
C       A STAND-ALONE DEVICE.
c 
c                        Input parameters:
c
c  x - the location (in terms of arc-length) on the curve of the
c       point to be found
c  n - the number of nodes in the equispaced discretization of the
c       curve produced by a preceding call to anarsbl
c  h - the sampling interval along the arc of the resampling returned
c       in array ts by the subroutine anarsbl
c  ts - the equispaced discretization of the curve produced by the
c       subroutine anarsbl
c  funcurve - the subroutine providing a parametrization
c        of the curve. the calling sequence of funcurve
c        must be
c 
c                  funcurve(t,par1,par2,x,y,dxdt,dydt,d2xdt2,d2ydt2),
c 
c        with:
c 
c       t - the parameter by which the user's curve is parametrizaed
c       x - the x-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       y - the y-coordinate of the point on the curve corresponding
c           to the values t of the parameter
c       dxdt - the derivative with respect to t of the  x-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c       dydt - the derivative with respect to t of the  y-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c       d2xdt2 - the 2nd derivative with respect to t of the  x-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c       d2ydt2 - the 2nd derivative with respect to t of the  y-coordinate
c           of the point on the curve corresponding to the values t
c           of the parameter
c  eps - the precision to which the point is to be found. Generally, it
c       is a good idea to make this parameter the same as in the
c       preceding call to anarsbl
c 
c                         Output parameters:
c 
c  tout - the value of the curve parameter corresponding to the location
c       x along the arc-length
c  xout, yout - the coordinates in R^2 of the point x on the curve
c  dxdt, dydt - the coordinates of the tangent (of length 1) to the
c       curve at the point (xout,yout)
c

c 
c       . . . find the two nodes (in  arc-length) between which the
c       user-supplied point is located
c 
        ier=0
        d=x/h
        m=d
c 
c       if the point is outside the user-supplied curve - bomb
c 
        if(m .gt. n) ier=512
        if(m .gt. n) return

        if(m.eq.n) tk = rl
        if(m.eq.n) goto 3000

        t1 = ts(m+1)
        x1 = m*h
c
c       use cubic interpolation to find the approximate location
c       of the curve parameter (as opposed to the arc-length)
c       corresponding to the user-supplied arc-length
c
        if (m.ge.1) tz(1)=ts(m)
        if (m.eq.0) tz(1)=ts(n)-ts(n+1)
        tz(2)=ts(m+1)
        tz(3)=ts(m+2)
        if(m.le.n-2) tz(4) = ts(m+3)
        if(m.eq.n-1) tz(4) = ts(n+1)+ts(2)
        xz(1)=(m-1)*h
        xz(2)=m*h
        xz(3)=(m+1)*h
        xz(4)=(m+2)*h
c
        cz(1)=1
        cz(2)=1
        cz(3)=1
        cz(4)=1

c
        do 100 i=1,4
        do 200 j=1,4
        if (i.eq.j) goto 200
        cz(i)=cz(i)*(x-xz(j))/(xz(i)-xz(j))
 200    continue
 100    continue
c
        t0=0
        do 300 i=1,4
        t0=t0+cz(i)*tz(i)
 300    continue
c 
c       use Newton to find the location of the curve parameter
c       (as opposed to the arc-length) corresponding to the
c       user-supplied arc-length
c 
        tk=t0
        tkold=t0
c 
cccc        epsgauss=eps*1000
        epsgauss=eps
        epsnewt=dsqrt(eps)/1000
        maxit=20
        fk=h
c
        ij0=1
        call findint(wright,nints,t1,ij0,ija)
c 
        do 1400 i=1,maxit
c 
c       evaluate the function whose root we are seeking
c 
        step=1
        do 1200 j=1,20
c
        call findint(wright,nints,tk,ija,ijb)
c
c
        if (ija.eq.ijb) then
        call anarsblg2_der2(t1,tk,funcurve,
     1       par1,par2,xgs,whts,mm,rltot)
        goto 2300
        end if
c
        call anarsblg2_der2(t1,wright(ija),funcurve,
     1       par1,par2,xgs,whts,mm,rltot)
c
        do 2200 kk=ija+1,ijb-1
        rltot=rltot+wvals(kk)
 2200   continue
c
        call anarsblg2_der2(wright(ijb-1),tk,funcurve,
     1       par1,par2,xgs,whts,mm,rltot2)
        rltot=rltot+rltot2
c
 2300   continue
c 
        fkold=fk
        fk=rltot+x1-x
c 
c       if the target function failed to decrease
c       - activate step-length control
c 
        if(abs(fk) .le. abs(fkold)) goto 1100
c 
        ier=1
        step=step/2
        tk=tkold-fkold/dldt *step
        goto 1200
c 
 1100 continue
c 
c       determine the dl/dt at that point and make a
c       newton step
c 
        call funcurve(tk,par1,par2,xk,yk,dxdt,dydt,d2xdt2,d2ydt2)
        dldt=dsqrt(dxdt**2+dydt**2)
c 
        tkold=tk
        tk=tk-fk/dldt *step
c 
c        if the required precision has been achieved
c        - terminate newton
c 
        if(abs(fk) .lt. epsnewt) goto 2000
c 
        goto 1400
 1200 continue
c 
        ier=64
        return
c 
 1400 continue
 2000 continue

 3000 continue
c
        tout=tk
        curvout = 0
        call funcurve(tout,par1,par2,xout,yout,dxdtout,dydtout,
     1    d2xdt2out,d2ydt2out)
c 
c 
        return
c
        end
c 
c 
c 
c 
c 
c 
        subroutine anarsblg2_der2(a,b,fun,par1,par2,t,w,m,rint)
        implicit real *8 (a-h,o-z)
        save
        dimension t(1),w(1),par1(1),par2(1)
c 
c       integrate the function fun on the interval [a,b]
c 
        rint=0
        u=(b-a)/2
        v=(b+a)/2
c
        do 1200 i=1,m
c
        tt=u*t(i)+v
c
        call fun(tt,par1,par2,x,y,dxdt,dydt,d2xdt2,d2ydt2)
c
        rint=rint+w(i)*dsqrt(dxdt**2+dydt**2)
c
 1200   continue
c
        rint=rint*u
c
        return
c
        end
