      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x0(:),y0(:),wsave(:),xy(:,:)
      real *8, allocatable :: srcinfo(:,:)
      real *8, allocatable :: ts(:),xver(:),yver(:)
      real *8, allocatable :: ts0(:),wts0(:)
      real *8, allocatable :: work(:)

      call prini(6,13)

      n = 316
      allocate(x0(n+1),y0(n+1),ts(n+1),xver(n+1),yver(n+1))
      allocate(xy(2,n))

      open(unit=87,file='bunny10out')
c
c   input vertices
c
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

      

c
c  nb = band limit
c
      nb = 200
c
c  interpolation tolerance
c
      eps = 1.0d-9

      ierm = 0
      nmax = 1
      call simple_curve_resampler_mem_der2(n,xy,nb,eps,nmax,nlarge,
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

      print *, "Here"



      allocate(wsave(lsave))
      ier = 0
      allocate(srcinfo(6,nout))
      call simple_curve_resampler_guru_der2(n,xy,nb,nlarge,lsave,
     1   lused,nout,
     1   srcinfo,h,curvelen,wsave,ts,ier)
      call prinf('ier=*',ier,1)
      call prinf('nout=*',nout,1)
      call prinf('nlarge=*',nlarge,1)
      call prinf('nb=*',nb,1)
      call prin2('ts=*',ts,24)
      call prin2('curvelen=*',curvelen,1)
      call prin2('srcinfo=*',srcinfo,24)


      
      erra = 0
      ra = 0
      do i=1,n
        xver(i) = 0
        yver(i) = 0
        call eval_curve_der2(ier,ts(i),wsave,xver(i),yver(i),dxt,dyt,
     1      d2xdt2,d2ydt2)
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


      stop
      end
