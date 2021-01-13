c
c     This code computes the solution to a transmission problem
c     and an impedance problem on the unit disk and compares it
c     to an analytic series solution.
c
c     The two problems are 
c     u0 - u1 = -uinc
c     du/dn0 - alpha du/dn1 = -dudnin
c
c     du/dn +  Z*u = 0
c    
c
      program test_circ
      implicit real *8 (a-h,o-z)

      real *8, allocatable :: srcinfo(:,:)
      complex *16, allocatable :: amattmp(:,:)
      complex *16, allocatable :: zplanew(:),zplanew_ser(:)
      complex *16, allocatable :: amat(:,:),soln(:),rhs(:)
      complex *16, allocatable :: amatt(:,:),solnt(:),rhst(:)
      complex *16, allocatable :: amati(:,:),solni(:),rhsi(:) 
      complex *16 zk(2),alpha,beta,ima
      complex *16 zpars(10),ssum,ssum1,ssum2
      complex *16 uex,u,uin,dudnin,uout,dudnout,u1,u2,z
      complex *16, allocatable :: fjvals(:),fjder(:)
      complex *16, allocatable :: fjvals2(:),fjder2(:)
      complex *16, allocatable :: fhvals(:),fhder(:)
      complex *16, allocatable :: zbeta(:),zgamma(:)
      real *8 dpars
      integer ipars
      real *8 targ(2,2),src(2,2)
      complex *16 a(2),b(2),q,zz,zimp

      integer, allocatable :: ipvt(:)
      character *1 transa
      external slp,dlp,comb,transmission_dir,transmission_neu



      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      ima = dcmplx(0.0d0,1.0d0)
      zs = ima/4.0d0

      allocate(fjvals(0:1000),fjder(0:1000))
      allocate(fjvals2(0:1000),fjder2(0:1000))
      allocate(fhvals(0:1000),fhder(0:1000))

c
c       setup sources and targets
c       
c       the first source/target is in the interior
c       and the second one is in the exterior
c

      src(1,1) = 0.1d0
      src(2,1) = 0.2d0

      src(1,2) = 12.1d0
      src(2,2) = 5.2d0

      targ(1,1) = -0.02d0
      targ(2,1) = 0.03d0
      
      targ(1,2) = 3.1d0
      targ(2,2) = 0.0d0
c
c
c       initialize wave speeds, and jump parameters
c
      zk(1) = 1.0d0
      zk(2) = 2.0d0

      a(1) = 1.0d0
      a(2) = 1.0d0
      b(1) = 1.2d0
      b(2) = 1.0d0
 
      nmax = 8000
      allocate(srcinfo(6,nmax))

      zkm = max(abs(zk(1)),abs(zk(2)))
      norder = 16
      ppw = 30

      zimp = ima*zk(1) + 1

c
c       geometry is an ellipse, major axis=a, minor axis=b
c  

c
c      estimate length of curve 
c
      ninit = 200
      a0 = 1.0d0
      b0 = 1.0d0

      rl_approx = 0

      h = 2*pi/ninit
      do i=1,ninit
        thet = (i-1)*h
        rl_approx = rl_approx+sqrt((a0*sin(thet))**2+
     1     (b0*cos(thet))**2)*h
      enddo

      n = max(ceiling(rl_approx*zkm/(2*pi)*ppw),100)

      h = 2*pi/n
      allocate(zplanew(n),zplanew_ser(n))

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
        d2xdt2 = -a0*cos(thet)
        d2ydt2 = -b0*sin(thet)
        dxdt = -a0*sin(thet)
        dydt = b0*cos(thet)
        srcinfo(6,i) = (dxdt*d2ydt2 - dydt*d2xdt2)/srcinfo(5,i)**3
        zplanew(i) = exp(ima*zk(2)*srcinfo(1,i))
      enddo
      nsys = 2*n

      kmax = 20
      rscale = 1.0d0
      ifder = 1 
      call jbessel2d(kmax,zk(2),rscale,fjvals,ifder,fjder)
      erra = 0
      ra = 0
      do i=1,n
        zplanew_ser(i) = 0
        thet = atan2(srcinfo(2,i),srcinfo(1,i))
        do j=-kmax,kmax
          zplanew_ser(i) = zplanew_ser(i) + ima**abs(j)*fjvals(abs(j))*
     1      exp(ima*j*thet)
        enddo
        erra = erra + abs(zplanew_ser(i)-zplanew(i))**2*h
        ra = ra + abs(zplanew(i))**2*h
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in plane wave=*',erra,1)




      allocate(amat(nsys,nsys),rhs(nsys),soln(nsys),rhsi(n),solni(n))

      do i=1,nsys
        do j=1,nsys
          amat(i,j) = 0
        enddo
      enddo
c
c
c        the unknowns are organized as
c    \sigma_{1},\mu_{1},\sigma_{2},\mu_{2}...
c

      call trans_mat(n,norder,h,srcinfo,zk,a,b,amat)

      allocate(amattmp(n,n),amati(n,n))
      do i=1,n
        do j=1,n
          amati(i,j) = 0
          amattmp(i,j) = 0
        enddo
      enddo

      call sprime_ext_mat(n,norder,h,srcinfo,zk(2),amati)
      call slp_mat(n,norder,h,srcinfo,zk(2),amattmp)

      do i=1,n
        do j=1,n
          amati(i,j) = amati(i,j) + zimp*amattmp(i,j)
        enddo
      enddo

      q = 0.5d0*(a(1)/b(1) + a(2)/b(2))
c
c      now generate right hand side
c
      do i=1,n
        rhs(2*i-1) = -zplanew(i)/q
        rhs(2*i) = -zplanew(i)*ima*zk(2)*srcinfo(3,i)
        rhsi(i) = -(zimp*zplanew(i) + zplanew(i)*ima*zk(2)*srcinfo(3,i))
      enddo

      do i=1,nsys
        soln(i) = rhs(i)
      enddo

      do i=1,n
        solni(i) = rhsi(i)
      enddo



      allocate(ipvt(nsys))
      info = 0
      call zgetrf(nsys,nsys,amat,nsys,ipvt,info)
      
      transa = 'n'
      info = 0

      nn = 1
      call zgetrs(transa,nsys,nn,amat,nsys,ipvt,soln,nsys,info)




      info = 0
      call zgetrf(n,n,amati,n,ipvt,info)
      info = 0
      nn = 1
      call zgetrs(transa,n,nn,amati,n,ipvt,solni,n,info)

c
c     evaluate solution at target point a target point an 
c     exterior point
c
   
      ssum1 = 0.0d0
      ssum2 = 0.0d0
      do i = 1,n
        w = srcinfo(5,i)*h

        call slp(srcinfo(1,i),targ(1,2),dpars,zk(2),ipars,u1)
        call dlp(srcinfo(1,i),targ(1,2),dpars,zk(2),ipars,u2)

        ssum2 = ssum2 + (soln(2*i)*u2-soln(2*i-1)*u1)*w/b(2)
        ssum1 = ssum1 + solni(i)*w*u1
      enddo

      u1 = 0

      allocate(zbeta(-kmax:kmax))
      allocate(zgamma(-kmax:kmax))
      ifder = 1
      call jbessel2d(kmax,zk(2),rscale,fjvals2,ifder,fjder2)
      call jbessel2d(kmax,zk(1),rscale,fjvals,ifder,fjder)

      call h2dall(kmax,zk(2),rscale,fhvals,ifder,fhder) 

      call prin2('b(1)=*',b(1),2)

      do j=-kmax,kmax
        jj = abs(j)
        zbeta(j) = ima**j*(-b(1)*zk(1)*fjder(jj)*fjvals2(jj) + 
     1     zk(2)*fjvals(jj)*fjder2(jj))/(-zk(2)*fhder(jj)*fjvals(jj)+
     2     b(1)*zk(1)*fhvals(jj)*fjder(jj))
        zgamma(j) = -ima**j*(zk(2)*fjder2(jj) + zimp*fjvals2(jj))/
     1      (zk(2)*fhder(jj) + zimp*fhvals(jj))
      enddo

      call prin2('zbeta=*',zbeta(-7),36)
      
      r = sqrt(targ(1,2)**2 + targ(2,2)**2)
      thet = atan2(targ(2,2),targ(1,2))
      call prin2('thet=*',thet,1)
      z = zk(2)*r
      call prin2('r=*',r,1)
      call h2dall(kmax,z,rscale,fhvals,ifder,fhder)
      u2 = 0
      do j=-kmax,kmax
        if(j.lt.0) then
          u2 = u2 + (-1.0d0)**j*zbeta(j)*fhvals(-j)*exp(ima*j*thet)
          u1 = u1 + (-1.0d0)**j*zgamma(j)*fhvals(-j)*exp(ima*j*thet)
        else
          u2 = u2 + zbeta(j)*fhvals(j)*exp(ima*j*thet)
          u1 = u1 + zgamma(j)*fhvals(j)*exp(ima*j*thet)
        endif
      enddo


      write(6,*) ssum2
      write(6,*) u2
      write(6,*) u2/ssum2

      write(6,*) " "
      write(6,*) " "

      write(6,*) ssum1
      write(6,*) u1
      write(6,*) u1/ssum1


c
      stop
      end
c
