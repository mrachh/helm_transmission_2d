      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x0(:),y0(:),wsave(:),xy(:,:)
      real *8, allocatable :: srcinfo(:,:)
      real *8, allocatable :: ts(:),xver(:),yver(:)
      real *8, allocatable :: ts0(:),wts0(:)

      call prini(6,13)

      n = 316
      allocate(x0(n+1),y0(n+1),ts(n+1),xver(n+1),yver(n+1))
      allocate(xy(2,n))

      open(unit=87,file='bunny10out')
      do i=1,n
        read(87,*) xtmp,ytmp
        x0(i) = xtmp
        y0(i) = ytmp
        xy(1,i) = x0(i)
        xy(2,i) = y0(i)
      enddo
      close(87)

cc      open(unit=87,file='bunny10out')
cc      do i=1,n
cc        write(87,*) x0(i),y0(i)
cc      enddo

      


      nb = 400
      eps = 1.0d-12

      ierm = 0
      nmax = 4
      call simple_curve_resampler_mem(n,xy,nb,eps,nmax,nlarge,
     1   nout,lsave,lused,ierm)
      if(ierm.eq.4) then
        call 
     1     prinf('nb too small, resulting curve self intersecting*',i,0)
        call prinf('increase nb*',i,0)
      else if(ierm.eq.2) then
        call prinf('desired inteprolation accuracy not reached*',i,0)
        call prinf('try increasing nmax*',i,0)
      endif
     

      call prinf('nb=*',nb,1)
      call prinf('nlarge=*',nlarge,1)
      call prinf('nout=*',nout,1)
      print *, "lsave=",lsave
      print *, "lused=",lused



      allocate(wsave(lsave))
      ier = 0
      allocate(srcinfo(6,nout))
      call simple_curve_resampler_guru(n,xy,nb,nlarge,lsave,lused,nout,
     1   srcinfo,h,curvelen,wsave,ts,ier)
      call prinf('ier=*',ier,1)
      call prinf('nout=*',nout,1)
      call prinf('nlarge=*',nlarge,1)
      call prinf('nb=*',nb,1)
      call prin2('ts=*',ts,24)
      call prin2('curvelen=*',curvelen,1)
      call prin2('curvature=*',srcinfo(6,1:nout),24)
      
      erra = 0
      ra = 0
      do i=1,n
        xver(i) = 0
        yver(i) = 0
        call eval_curve(ier,ts(i),wsave,xver(i),yver(i),dxt,dyt,curv)
        erra = erra + (xver(i)-xy(1,i))**2
        ra = ra + xy(1,i)**2
        erra = erra + (yver(i)-xy(2,i))**2
        ra = ra + xy(2,i)**2
      enddo
      erra = sqrt(erra/ra)
      print *, "erra=",erra
      print *, "ra=",ra

      call prin2('error in curve interpolation=*',erra,1)

      open(unit=87,file='out.dat')

      do i=1,nout
        write(87,*) srcinfo(1,i),srcinfo(2,i)
      enddo
c
c  test curvature here
c
      
      rstart = 0.5d0
      rend = curvelen*(0.33d0*(1-rstart) + rstart)

      k = 2*n
      call prinf('k=*',k,1)
      itype = 1
      allocate(ts0(k),wts0(k))
      call legeexps(itype,k,ts0,umat,vmat,wts0)

      curvint = 0
      r1 = 0
      do i=1,k
        teval = (ts0(i)+1)/2*rend
        call eval_curve(ier,teval,wsave,xtmp,ytmp,dxt,dyt,curv)
        curvint = curvint + wts0(i)*rend/2*curv
        r1 = r1 + wts0(i)*rend/2
      enddo

      call prin2('rend=*',rend,1)
      call prin2('r1=*',r1,1)

      call prin2('curvint=*',curvint,1)
      call eval_curve(ier,rend,wsave,xtmp,ytmp,dxt,dyt,curv0)

      print *, dxt,dyt,sqrt(dxt**2 + dyt**2)

      r0 = 0.0d0
      call eval_curve(ier,r0,wsave,xtmp,ytmp,dxt0,dyt0,curv0)

      print *, dxt0,dyt0,dxt,dyt
      r1 = dxt0*dxt + dyt0*dyt
      r2 = cos(curvint)

      call prin2('r1=*',r1,1)
      call prin2('r2=*',r2,1)

      erra = abs(r2-r1)
      call prin2('error=*',erra,1)
      


      

      stop
      end
