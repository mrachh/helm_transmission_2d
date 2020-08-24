      implicit real *8 (a-h,o-z)
      real *8, allocatable :: x0(:),y0(:),wsave(:),xy(:,:)
      real *8, allocatable :: srcinfo(:,:)
      real *8, allocatable :: ts(:),xver(:),yver(:)

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

      ierm = 0
      call simple_curve_resampler_mem(n,xy,nb,nlarge,nout,lsave,
     1   lused,ierm)
      if(ierm.ne.0) then
        call 
     1     prinf('nb too small, resulting curve self intersecting*',i,0)
        call prinf('increase nb*',i,0)
      endif
     

      call prinf('nb=*',nb,1)
      call prinf('nlarge=*',nlarge,1)
      call prinf('nout=*',nout,1)
      call prinf('lsave=*',lsave,1)
      call prinf('lused=*',lused,1)


      allocate(wsave(lsave))
      ier = 0
      allocate(srcinfo(5,nout))
      call simple_curve_resampler_guru(n,xy,nb,nlarge,lsave,lused,nout,
     1   srcinfo,h,curvelen,wsave,ts,ier)
      call prinf('ier=*',ier,1)
      call prinf('nout=*',nout,1)
      call prinf('nlarge=*',nlarge,1)
      call prinf('nb=*',nb,1)
      call prin2('ts=*',ts,24)
      call prin2('curvelen=*',curvelen,1)
      
      erra = 0
      ra = 0
      do i=1,n
        xver(i) = 0
        yver(i) = 0
        call eval_curve(ier,ts(i),wsave,xver(i),yver(i),dxt,dyt)
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
