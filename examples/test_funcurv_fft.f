      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x(:),y(:),dxdt(:),dydt(:)
      real *8, allocatable :: xex(:),yex(:),dxex(:),dyex(:)
      real *8, allocatable :: x0(:),y0(:),dx0dt(:),dy0dt(:)
      complex *16, allocatable :: z(:),zf(:),dzdt(:),dzdtf(:)
      complex *16, allocatable :: zh(:),zhf(:),zdhdt(:),zdhdtf(:)
      complex *16, allocatable :: par1(:)
      real *8, allocatable :: wsave1(:),wsave2(:)
      real *8, allocatable :: work(:),tts(:)
      complex *16 ima
      real *8 par2(100)
      data ima/(0.0d0,1.0d0)/
      external funcurv_fft

      call prini(6,13)
      done = 1
      pi = atan(done)*4


      n1 = 121
      allocate(x0(n1),y0(n1),dx0dt(n1),dy0dt(n1),z(n1),zf(n1))
      allocate(dzdtf(n1),dzdt(n1))

      n2 = 101

      allocate(zh(n2),zdhdt(n2),zhf(n2),zdhdtf(n2))

      zhf = 0
      zdhdtf = 0

      nosch = 10
      r3 = 0.03d0

      nosc = 5
      r1 = 1.0d0
      r2 = 0.17d0

      h = 2*pi/(n1+0.0d0)

      rl = 0

      do i=1,n1
        t = (i-1)*h
        r = r1 + r2*cos(nosc*t) 
        drdt = -r2*nosc*sin(nosc*t)
        x0(i) = r*cos(t)
        y0(i) = r*sin(t)
        dx0dt(i) = -r*sin(t) + drdt*cos(t)
        dy0dt(i) = r*cos(t) + drdt*sin(t)

        rl = rl + sqrt(dx0dt(i)**2 + dy0dt(i)**2)*h
      enddo
      call prin2('x=*',x0,24)
      call prin2('initial length of curve=*',rl,1)

      hnew = rl/(n1+0.0d0)

      rc = 2*pi/rl
      rl2 = 0
      do i=1,n1
        t = (i-1)*hnew
        targ = t*rc
        r = r1 + r2*cos(nosc*targ)
        drdt = -r2*nosc*rc*sin(nosc*targ)
        x0(i) = r*cos(targ)
        y0(i) = r*sin(targ)
        z(i) = x0(i) + ima*y0(i)
        dx0dt(i) = -r*sin(targ)*rc + drdt*cos(targ)
        dy0dt(i) = r*cos(targ)*rc + drdt*sin(targ)

        dzdt(i) =  dx0dt(i) + ima*dy0dt(i)
        rl2 = rl2 + sqrt(dx0dt(i)**2 + dy0dt(i)**2)*hnew
      enddo


      h2 = 2*pi/(n2+0.0d0)
      print *, "rc=",rc
      do i=1,n2
        t = (i-1)*h2
        r = r3*sin(nosch*t)
        drdt = nosch*r3*cos(nosch*t)
        zh(i) = r
        zdhdt(i) = drdt*rc 
      enddo

      call prin2('x=*',x0,24)
      call prin2('new length of curve=*',rl2,1)
      erra = abs(rl-rl2)
      call prin2('error in length of curve=*',erra,1)

      allocate(wsave1(10*n1+100))
      allocate(wsave2(10*n2+100))
      call zffti(n1,wsave1)
      call zffti(n2,wsave2)

      do i=1,n1
        zf(i) = z(i)
        dzdtf(i) = dzdt(i)
      enddo

      do i=1,n2
        zhf(i) = zh(i)
        zdhdtf(i) = zdhdt(i)
      enddo
      

      call zfftf(n1,zf,wsave1)
      call zfftf(n1,dzdtf,wsave1)
      
      do i=1,n1
        zf(i) = zf(i)/n1
        dzdtf(i) = dzdtf(i)/n1
      enddo

      call zfftf(n2,zhf,wsave2)
      call zfftf(n2,zdhdtf,wsave2)

      do i=1,n2
        zhf(i) = zhf(i)/n2
        zdhdtf(i) = zdhdtf(i)/n2
      enddo

      call prin2('done with ffts*',i,0)
      allocate(par1(2*n1+2*n2+100))

      par1(1:n1) = zf
      par1(n1+1:2*n1) = dzdtf
      
      par1(2*n1+1:2*n1+n2) = zhf
      par1(2*n1+n2+1:2*n1+2*n2) = zdhdtf

      par2(1) = n1
      par2(2) = n2
      par2(3) = n1
      par2(4) = rl

      ntest = 100
      allocate(x(ntest),y(ntest),dxdt(ntest),dydt(ntest))
      allocate(xex(ntest),yex(ntest),dxex(ntest),dyex(ntest))

      erra = 0
      ra = 0
      call prin2('rl=*',rl,1)
      call prin2('rc=*',rc,1)
      do i=1,ntest
        t = hkrand(0)*rl
        call funcurv_fft(t,par1,par2,x(i),y(i),dxdt(i),dydt(i))
        targ = t*rc
        r = r1 + r2*cos(nosc*targ) 
        drdt = -r2*nosc*rc*sin(nosc*targ) 
        d2rdt2 = -r2*nosc*nosc*rc*rc*cos(nosc*targ)
        xex(i) = r*cos(targ)  
        yex(i) = r*sin(targ)
        dxex(i) = -r*sin(targ)*rc + drdt*cos(targ)
        dyex(i) = r*cos(targ)*rc + drdt*sin(targ)

        d2xdt2 = -r*cos(targ)*rc**2 - 2*drdt*rc*sin(targ) +
     1     d2rdt2*cos(targ)
        d2ydt2 = -r*sin(targ)*rc**2 + 2*drdt*rc*cos(targ) + 
     1     d2rdt2*sin(targ)
        print *, d2xdt2,d2ydt2

        dst = sqrt(dxex(i)**2 + dyex(i)**2)

        xex(i) = xex(i) - r3*sin(nosch*targ)*dyex(i)/dst
        yex(i) = yex(i) + r3*sin(nosch*targ)*dxex(i)/dst

        dy0 = dyex(i)
        dx0 = dxex(i)

        rdot = dx0*d2xdt2 + dy0*d2ydt2
        print *, "rdot=",rdot
        print *, "dst=",dst

        dxex(i) = dxex(i) - r3*sin(nosch*targ)*d2ydt2/dst - 
     1     r3*rc*nosch*cos(nosch*targ)*dy0/dst + 
     2     r3*sin(nosch*targ)*dy0*rdot/dst**3
        dyex(i) = dyex(i) + r3*sin(nosch*targ)*d2xdt2/dst + 
     1     r3*rc*nosch*cos(nosch*targ)*dx0/dst -
     2     r3*sin(nosch*targ)*dx0*rdot/dst**3
        
        erra = erra + (xex(i)-x(i))**2
        erra = erra + (yex(i)-y(i))**2
        erra = erra + (dxex(i)-dxdt(i))**2
        erra = erra + (dyex(i)-dydt(i))**2

        if(i.le.2) print *, x(i)/xex(i)
        if(i.le.2) print *, dyex(i),dydt(i)

        ra = ra + xex(i)**2
      enddo
      erra = sqrt(erra/ra)

      call prin2('x=*',x,24)
      call prin2('xex=*',xex,24)

      call prin2('dydt=*',dydt,24)
      call prin2('dyex=*',dyex,24)

      call prin2('dxdt=*',dxdt,24)
      call prin2('dxex=*',dxex,24)

      call prin2('erra=*',erra,1)

      lw = 100000
      allocate(work(lw))

      nnew = 100
      ier = 0
      eps = 1.0d-13
      lsave = 0 
      allocate(tts(nnew+1))
      print *, "here"
      call prin2('par1=*',par1,24)
      call prinf('n1=*',n1,1)
      call prinf('n2=*',n2,1)
      call prin2('rl=*',rl,1)
      work = 0
      call curve_resampler_guru(ier,n1,n2,par1,rl,nnew,eps,tts,
     1  h,rltot,work,lw,lsave)

      call prinf('lsave=*',lsave,1)
      call prin2('rltot=*',rltot,1)
      call prin2('tts=*',tts,24)
      call prin2('h=*',h,1)

      do i=1,3
        call funcurv_fft(tts(i),par1,par2,x3,y3,dxdt3,dydt3)
        ttest = (i-1)*h
        call anafpt(ier,ttest,nnew,h,tts,funcurv_fft,par1,par2,
     1    eps,tout,xout,yout,dxdtout,dydtout,work)
        dst = sqrt(dxdt3**2 + dydt3**2)
        dxdt3 = dxdt3/dst
        dydt3 = dydt3/dst
        print *, "i=",i
        print *, x3, xout
        print *, y3, yout
        print *, dxdt3, dxdtout
        print *, dydt3, dydtout
        print *, " "
        print *, " "
      enddo
      


      
      


      stop
      end
