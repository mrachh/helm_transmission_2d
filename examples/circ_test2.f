      implicit real *8 (a-h,o-z)
      complex *16 zk(2),alpha,beta,ima
      complex *16 zpars(10),ssum,ssum1,ssum2
      complex *16 uex,u,uin,dudnin,uout,dudnout,u1,u2,z
      complex *16, allocatable :: fjvals(:),fjder(:)
      complex *16, allocatable :: fjvals2(:),fjder2(:)
      complex *16, allocatable :: fhvals(:),fhder(:)
      complex *16, allocatable :: zbeta(:),zgamma(:)

      complex *16, allocatable :: zvals1(:),zvals2(:)
      
      real *8 dpars
      integer ipars
      real *8 targ(2,2),src(2,2)
      complex *16 a(2),b(2),q,zz,zimp

      integer, allocatable :: ipvt(:)
      character *1 transa
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)
      done = 1
      pi = atan(done)*4

      allocate(fjvals(0:1000),fjder(0:1000))
      allocate(fjvals2(0:1000),fjder2(0:1000))
      allocate(fhvals(0:1000),fhder(0:1000))

      zk(1) = (1.2d0+0.1d0*ima)*15
      zk(2) = (1.0d0)*15

      a(1) = 1.0d0
      a(2) = 1.0d0
      b(1) = 1.0d0/1.2d0
      b(2) = 1.0d0

      zimp = b(1)*(ima*zk(1) + 1.0d0)
      kmax = 100 

      allocate(zbeta(-kmax:kmax),zgamma(-kmax:kmax))
      rscale = 1.0d0
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

      call prin2('zbeta=*',zbeta(0),14)
      call prin2('zgamma=*',zgamma(0),14)

      nang = abs(zk(2))*40
      r = 400.1d0
      allocate(zvals1(nang),zvals2(nang))
      z = zk(2)*r
      call prin2('r=*',r,1)
      call h2dall(kmax,z,rscale,fhvals,ifder,fhder)
      call prinf('nang=*',nang,1)
      do i=1,nang
        thet = 2.0d0*pi*(i-1)/(nang+0.0d0)
        u1 = 0
        u2 = 0
        do j=-kmax,kmax
          if(j.lt.0) then
cc            u1 = u1 + (-1.0d0)**j*zbeta(j)*fhvals(-j)*exp(ima*j*thet)
cc            u2 = u2 + (-1.0d0)**j*zgamma(j)*fhvals(-j)*exp(ima*j*thet)
              u1 = u1 + (-ima)**j*exp(ima*j*thet)*zbeta(j)
              u2 = u2 + (-ima)**j*exp(ima*j*thet)*zgamma(j)
          else
cc            u1 = u1 + zbeta(j)*fhvals(j)*exp(ima*j*thet)
cc            u2 = u2 + zgamma(j)*fhvals(j)*exp(ima*j*thet)
              u1 = u1 + (-ima)**j*exp(ima*j*thet)*zbeta(j)
              u2 = u2 + (-ima)**j*exp(ima*j*thet)*zgamma(j)
          endif
        enddo
        zvals1(i) = u1
        zvals2(i) = u2
        write(77,*) real(zvals1(i)),imag(zvals1(i)),real(zvals2(i)),
     1    imag(zvals2(i))
      enddo

      call prin2('zvals1=*',zvals1,24)
      call prin2('zvals2=*',zvals2,24)



      stop
      end
