c
c     This code solves the exterior Dirichlet problem governed by the 
c     Helmholtz equation in 2D using Alpert quadratures.
c    
c     We form a boundary integral equation using either the SLP or the 
c     Combined field equation (u = DLP - i eta SLP).
c
      program axihelm_neu
      implicit real *8 (a-h,o-z)

      real *8, allocatable :: srcinfo(:,:)
      complex *16, allocatable :: amat(:,:),soln(:),rhs(:)
      complex *16 zk,alpha,beta,ima
      complex *16 zpars(10),ssum
      complex *16 uex,u
      real *8 dpars
      integer ipars
      real *8 targ(2),src(2)

      integer, allocatable :: ipvt(:)
      character *1 transa
      external slp,dlp,comb


      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0

      src(1) = 0.01d0
      src(2) = -0.07d0

      targ(1) = 12.1d0
      targ(2) = 5.2d0

      nmax = 8000
      allocate(srcinfo(5,nmax))
c
c       geometry params for larry cup
c

      zk = 1.0d1
      zk = 0.1d0
      zpars(1) = zk
      norder = 16
      ppw = 60

c
c       geometry is an ellipse, major axis=a, minor axis=b
c  


c
c      estimate length of curve 
c
      ninit = 200
      a = 1
      b = 0.3

      rl_approx = 0

      h = 2*pi/ninit
      call prin2('h=*',h,1)
      do i=1,ninit
        thet = (i-1)*h
        rl_approx = rl_approx+sqrt((a*sin(thet))**2+(b*cos(thet))**2)*h
      enddo

      n = max(ceiling(rl_approx*abs(zk)/(2*pi)*ppw),100)

      n = 200
      h = 2*pi/n

      call prinf('npts=*',n,1)
c
c       initialize the geometry
c
      do i=1,n
        thet = (i-1)*h
        srcinfo(1,i) = a*cos(thet)
        srcinfo(2,i) = b*sin(thet)
        srcinfo(5,i) = sqrt((a*sin(thet))**2 +(b*cos(thet))**2)
        srcinfo(3,i) = b*cos(thet)/srcinfo(5,i)
        srcinfo(4,i) = a*sin(thet)/srcinfo(5,i)
      enddo

      call slp(src,targ,dpars,zpars,ipars,uex)

      allocate(amat(n,n),rhs(n),soln(n))
c
c     form system matrix 
c
c
      zpars(2) = -ima*zk
      zpars(2) = 1
      zpars(3) = 1
      zpars(3) = 0
      call prinf('n=*',n,1)
      call prin2('h=*',h,1)
      call prin2('zpars=*',zpars,6)
      call prinf('norder=*',norder,1)
      call formmatbac(amat,norder,n,srcinfo,h,slp,dpars,zpars,ipars)
  
      print *, zpars(3)
      do i=1,n
        amat(i,i) = amat(i,i) + 0.5d0*zpars(3)
      enddo

      do i=1,n
        write(32,*) srcinfo(1,i),srcinfo(2,i)
      enddo

      do i=1,n
        call slp(src,srcinfo(1,i),dpars,zpars,ipars,rhs(i))
        soln(i) = rhs(i)
      enddo
      call prin2('rhs=*',rhs,24)

      allocate(ipvt(n))
      info = 0
      call zgetrf(n,n,amat,n,ipvt,info)
      
      transa = 'n'
      info = 0
cc      call prin2('soln=*',soln,2*n)
cc      call prin2('amat=*',amat,2*n*n)

      nn = 1
      call zgetrs(transa,n,nn,amat,n,ipvt,soln,n,info)
      call prin2('soln=*',soln,12)
      call prinf('info=*',info,1)

c
c     evaluate solution at target point     
c
   
      ssum = 0.0d0
      call prin2('targ=*',targ,2)
      call prin2('zpars=*',zpars,6)
      call prin2('h=*',h,1)
      do i = 1,n
        w = srcinfo(5,i)*h
        call slp(srcinfo(1,i),targ,dpars,zpars,ipars,u)
        ssum = ssum + soln(i)*u*w
      enddo

      write(6,*) ' ssum is ',ssum
      write(13,*) ' ssum is ',ssum
      write(6,*) ' uex is ',uex
      write(6,*) ' ssum/uex is ',(ssum/uex)
      write(13,*) ' ssum/uex is ',(ssum/uex)

c
      stop
      end
c
