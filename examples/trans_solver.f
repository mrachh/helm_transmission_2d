c
c     This code solves the exterior Dirichlet problem governed by the 
c     Helmholtz equation in 2D using Alpert quadratures.
c    
c
c     We assume the representation takes the form:
c
c     u0 = (-1/b0) S_k0 sigma + (1/b0) D_k0 mu  (exterior)
c
c     u = (-1/b) S_k sigma + (1/b) D_k mu     (interior)
c
c     so that
c
c     [au] = (1/2)[a0/b0 + a/b]mu + [(a0/b0) D^*_k0 - (a/b)D^*_k] mu 
c             - [ a0/b0 S_k0 -  a/b S_k] sigma = f1
c     [b dudn] = [D'_k0 - D'_k] mu  + sigma - [S'_k0 - S'_k] = f2
c
c     or (with q = (1/2)[a0/b0+a/b] )
c
c     [au]/q = mu + [(a0/b0) D^*_k0 - (a/b)D^*_k]/q mu 
c             - [ a0/b0 S_k0 -  a/b S_k]/q sigma = f1/q
c
      program axihelm_neu
      implicit real *8 (a-h,o-z)

      real *8, allocatable :: srcinfo(:,:)
      complex *16, allocatable :: amattmp(:,:)
      complex *16, allocatable :: amat(:,:),soln(:),rhs(:)
      complex *16 zk(2),alpha,beta,ima
      complex *16 zpars(10),ssum,ssum1,ssum2
      complex *16 uex,u,uin,dudnin,uout,dudnout,u1,u2
      real *8 dpars
      integer ipars
      real *8 targ(2,2),src(2,2)
      complex *16 a(2),b(2),q,zz

      integer, allocatable :: ipvt(:)
      character *1 transa
      external slp,dlp,comb,transmission_dir,transmission_neu


      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0

c
c       setup sources and targets
c       
c       the first source/target is in the interior
c       and the second one is in the exterior
c

      src(1,1) = 0.01d0
      src(2,1) = -0.07d0

      src(2,1) = 10.5d0
      src(2,2) = -3.1d0

      targ(1,1) = -0.02d0
      targ(2,1) = 0.03d0
      
      targ(1,2) = 12.1d0
      targ(2,2) = 5.2d0
c
c
c       initialize wave speeds, and jump parameters
c
      zk(1) = 1.0d0
      zk(2) = 2.0d0

      a(1) = 2.0d0
      a(2) = 1.0d0
      b(1) = 4.0d0
      b(2) = 3.0d0
 
      nmax = 8000
      allocate(srcinfo(5,nmax))

c
c       geometry params for larry cup
c
      zkm = max(abs(zk(1)),abs(zk(2)))
      norder = 16
      ppw = 20

c
c       geometry is an ellipse, major axis=a, minor axis=b
c  

c
c      estimate length of curve 
c
      ninit = 200
      a0 = 1.1d0
      b0 = 1.3d0

      rl_approx = 0

      h = 2*pi/ninit
      do i=1,ninit
        thet = (i-1)*h
        rl_approx = rl_approx+sqrt((a0*sin(thet))**2+
     1     (b0*cos(thet))**2)*h
      enddo

      n = max(ceiling(rl_approx*zkm/(2*pi)*ppw),100)
      n = 200

      h = 2*pi/n

c
c       initialize the geometry
c
      do i=1,n
        thet = (i-1)*h
        srcinfo(1,i) = a0*cos(thet)
        srcinfo(2,i) = b0*sin(thet)
        srcinfo(5,i) = sqrt((a0*sin(thet))**2 +(b0*cos(thet))**2)
        srcinfo(3,i) = b0*cos(thet)/srcinfo(5,i)
        srcinfo(4,i) = a0*sin(thet)/srcinfo(5,i)
      enddo
      nsys = 2*n



      allocate(amat(nsys,nsys),rhs(nsys),soln(nsys))
      allocate(amattmp(n,n))
c
c
c        the unknowns are organized as
c    \sigma_{1},\mu_{1},\sigma_{2},\mu_{2}...
c

c
c     form a/b sk - a0/b0 sk0 system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zk(1)
      zpars(2) = zk(2)
      zpars(3) = a(1)/b(1)
      zpars(4) = -a(2)/b(2)
      zpars(5) = 0
      zpars(6) = 0
      q = 0.5d0*(a(1)/b(1) + a(2)/b(2))

      call prin2('q=*',q,2)
      call prin2('zpars=*',zpars,12)
      call prinf('norder=*',norder,1)
      call prin2('h=*',h,1)
      call formmatbac(amattmp,norder,n,srcinfo,h,transmission_dir,
     1    dpars,zpars,ipars)

      

  
      do i=1,n
        do j=1,n
          zz = amattmp(i,j)/q
          amat(2*i-1,2*j-1) = zz 
        enddo
      enddo

c
c     form -a/b dk + a0/b0 dk0 system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zk(1)
      zpars(2) = zk(2)
      zpars(3) = 0
      zpars(4) = 0
      zpars(5) = -a(1)/b(1)
      zpars(6) = a(2)/b(2)
      call formmatbac(amattmp,norder,n,srcinfo,h,transmission_dir,
     1    dpars,zpars,ipars)
 
      do i=1,n
        do j=1,n
          amat(2*i-1,2*j) = amattmp(i,j)/q
          write(53,*) i,j,real(amattmp(i,j)/q),imag(amattmp(i,j)/q)
        enddo
        amat(2*i-1,2*i) = amat(2*i-1,2*i)+1 
      enddo

c
c     form sk' - sk0' system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zk(1)
      zpars(2) = zk(2)
      zpars(3) = 1
      zpars(4) = -1
      zpars(5) = 0
      zpars(6) = 0

      do i=1,n 
        do j=1,n
          amattmp(i,j) = 0
        enddo
      enddo

      call prinf('norder=*',norder,1)
      call prinf('n=*',n,1)
      call formmatbac(amattmp,norder,n,srcinfo,h,transmission_neu,
     1    dpars,zpars,ipars)
  
      do i=1,n
        do j=1,n
          write(56,*) i,j,real(amattmp(i,j)),imag(amattmp(i,j))
          amat(2*i,2*j-1) = amattmp(i,j) 
        enddo
      enddo

      stop

c
c     form -dk + dk0 system matrix (k is interior k0 is 
c      exterior
c
      zpars(1) = zk(1)
      zpars(2) = zk(2)
      zpars(3) = 0
      zpars(4) = 0
      zpars(5) = -1
      zpars(6) = 1
      call formmatbac(amattmp,norder,n,srcinfo,h,transmission_neu,
     1    dpars,zpars,ipars)
 
      do i=1,n
        do j=1,n
          amat(2*i,2*j) = amattmp(i,j) 
        enddo
        amat(2*i,2*i-1) = amat(2*i,2*i-1)+1 
      enddo

c
c       end of generating matrix
c

c
c      now generate right hand side
c

      
      do i=1,n
        call slp(src(1,2),srcinfo(1,i),dpars,zk(1),ipars,uin)
        call slp(src(1,1),srcinfo(2,i),dpars,zk(2),ipars,uout)
        call sprime(src(1,2),srcinfo(1,i),dpars,zk(1),ipars,dudnin)
        call sprime(src(1,1),srcinfo(2,i),dpars,zk(2),ipars,dudnout)

        rhs(2*i-1) = a(2)*uout - a(1)*uin
        rhs(2*i) = b(2)*dudnout - b(1)*dudnin
      enddo

      do i=1,nsys
        soln(i) = rhs(i)
      enddo

      allocate(ipvt(nsys))
      info = 0
      call zgetrf(nsys,nsys,amat,nsys,ipvt,info)
      
      transa = 'n'
      info = 0

      nn = 1
      call zgetrs(transa,nsys,nn,amat,nsys,ipvt,soln,nsys,info)

c
c     evaluate solution at target point a target point in 
c     interior and exterior
c
   
      ssum1 = 0.0d0
      ssum2 = 0.0d0
      do i = 1,n
        w = srcinfo(5,i)*h
        call slp(srcinfo(1,i),targ(1,1),dpars,zk(1),ipars,u1)
        call dlp(srcinfo(1,i),targ(1,1),dpars,zk(1),ipars,u2)

        ssum1 = ssum1 + (soln(2*i)*u2-soln(2*i-1)*u1)*w/b(1)

        call slp(srcinfo(1,i),targ(1,2),dpars,zk(2),ipars,u1)
        call dlp(srcinfo(1,i),targ(1,2),dpars,zk(2),ipars,u2)

        ssum2 = ssum2 + (soln(2*i)*u2-soln(2*i-1)*u1)*w/b(2)
      enddo

      call slp(src(1,2),targ(1,1),dpars,zk(1),ipars,u1)
      call slp(src(1,1),targ(2,1),dpars,zk(2),ipars,u2)

      write(6,*) ssum1,ssum2
      write(6,*) u1,u2


c
      stop
      end
c
