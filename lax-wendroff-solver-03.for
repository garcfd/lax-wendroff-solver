
      program main
      implicit none

      integer, parameter :: ni=400,nj=400,nk=1
      integer live(0:ni+1,0:nj+1)

      real xn(ni+1),yn(nj+1),zn(nk+1)
      real r(ni,nj),u(ni,nj),v(ni,nj)
      real p(ni,nj),t(ni,nj)
      real Q(0:ni+1,0:nj+1,4)
      real E(0:ni+1,0:nj+1,4)
      real F(0:ni+1,0:nj+1,4)
      real rhs(4),src(4)
      real resid(4),resid1(4)

      integer i,j,k,n,it,sa,nits,nsav
      integer ccln,ii,jj,nout,iloc,jloc
      integer im1,jm1,ip1,jp1

      real dx,dy,dz,dt
      real uu,vv,rr,pp,tt,ee,eps
      real rinf,uinf,vinf,pinf,tinf
      real rgas,cp,cv,kth,ga,ma,cs
      real ptot,ttot,vm,cfl
 
      real unrm,vnrm,unew,vnew,tweak
      real unew1,vnew1,unew2,vnew2
      real usum,vsum,rsum,tsum,psum
      real nrm1,nrm2,mag1,mag2,dotp
      real diff1,diff2,uold,vold
      real urf,xx,yy

      character(30) outfile
      character(3) string
      CHARACTER(LEN=1)  :: lf
      CHARACTER(LEN=10) :: str1,str2,str3,str4
      lf = char(10)



C-----specify constants
      dx = 0.1      ! spacing
      dy = 0.1      ! spacing
      dz = 0.1      ! spacing
      eps = 8.0     ! wiggle
      urf = 1.0     ! smooth
      cfl = 0.1     ! Ncfl


C-----air properties
      rgas = 287.1
      cp   = 1004.5
      cv   = 717.5
      ga   = 1.4
      kth  = 0.026 

      rinf = 1.225 
      uinf = 100.0
      vinf = 0.0
      pinf = 101325.0
      tinf = 288.0
      ma = 0.9

c      ptot = 107587.8847 ! 100.0
c      ttot =    292.9776 ! 100.0

c      ptot = 128034.843 ! 200.0
c      ttot =    307.910 ! 200.0

      ptot = pinf*(1.0 + 0.2*ma*ma)**3.5
      ttot = tinf*(1.0 + 0.2*ma*ma)

      nsav = 200
      nits = 1000
      nout = 100


C-----specify live/dead regions (if needed)
      live = 1    ! everything live
      iloc = 100  ! blockage centre
      jloc = nj/2 ! blockage centre
       
      do i = 1, ni
      do j = 1, nj

        xx = real(i-iloc)
        yy = real(j-jloc)
        rr = sqrt(xx*xx + yy*yy)

        if (rr.le.15) then
          live(i,j) = 0 ! dead region
        endif

      enddo
      enddo


C-----create uniform mesh nodes
      xn(1) = 0.0 ! start point
      yn(1) = 0.0 ! start point
      zn(1) = 0.0 ! start point
      zn(2) = dz  ! for 2d case

      do i = 2, ni+1
        xn(i) = xn(i-1) + dx
      enddo
      do j = 2, nj+1
        yn(j) = yn(j-1) + dy
      enddo


C-----initialise arrays
      do i = 1, ni
      do j = 1, nj

        r(i,j) = rinf 
        u(i,j) = uinf 
        v(i,j) = vinf
        p(i,j) = pinf 
        t(i,j) = tinf 

      enddo
      enddo

C-----residial data
      OPEN(UNIT=10, FILE="conv.dat")




C-----start iterations
      do sa = 1,nsav

C==============
!$acc data copy(r,u,v,p,t,Q,E,F,live)
C==============

      do it = 1,nits

!$acc  parallel loop collapse(2)
!$acc& default(none)
!$acc& private(src,rhs)
!$acc& reduction(+:resid)


C-----step 1 - compute vectors at cell centres
      do i = 0,ni+1
      do j = 0,nj+1

C-------bottom wall
        if ((j.eq.0).and.(i.ge.1).and.(i.le.ni)) then
          rr = r(i,1)
          uu = u(i,1)
          vv = 0.0 ! v(i,1) ! wall
          pp = p(i,1)
          tt = t(i,1)
        endif

C-------top wall
        if ((j.eq.nj+1).and.(i.ge.1).and.(i.le.ni)) then
          rr = r(i,nj)
          uu = u(i,nj)
          vv = 0.0 ! v(i,nj) ! wall
          pp = p(i,nj)
          tt = t(i,nj)
        endif

C-------left wall
        if ((i.eq.0).and.(j.ge.1).and.(j.le.nj)) then

          uu = u(1,j)
          vv = v(1,j)
          vm = sqrt(uu*uu + vv*vv)       ! Vmag
          ma = vm/sqrt(ga*rgas*t(1,j))   ! V=MC
          pp = ptot/(1+0.2*ma*ma)**3.5   ! isentropic
          tt = ttot/(1+0.2*ma*ma)        ! isentropic
          rr = pp/(rgas*tt)

c         uu = u(1,j)
c         vv = v(1,j)
c         pp = ptot - 0.5*rinf*(uu*uu + vv*vv) ! bernoulli
c         tt = ttot - 0.5*(uu*uu + vv*vv)/cp   ! bernoulli
c         rr = pp/(rgas*tt)

c         pp = p(1,j)
c         tt = tinf
c         uu = uinf
c         vv = vinf
c         rr = pp/(rgas*tt)

        endif

C-------right wall
        if ((i.eq.ni+1).and.(j.ge.1).and.(j.le.nj)) then
          uu = u(ni,j)
          vv = v(ni,j)
          tt = t(ni,j)

          vm = sqrt(uu*uu + vv*vv) ! Vmag
          ma = vm/sqrt(ga*rgas*tt) ! V=MC

          if (ma.ge.1.0) pp = p(ni,j)
          if (ma.lt.1.0) pp = pinf

          rr = pp/(rgas*tt)
        endif

C-------main block of inner cells (including obstacle)
        if ((i.ge.1).and.(i.le.ni).and.(j.ge.1).and.(j.le.nj)) then
          rr = r(i,j)
          uu = u(i,j)
          vv = v(i,j)
          pp = p(i,j)
          tt = t(i,j)
        endif

C-------for all cell centres
        ee = rr*cv*tt + 0.5*rr*(uu*uu + vv*vv)

        Q(i,j,1) = rr
        Q(i,j,2) = rr*uu
        Q(i,j,3) = rr*vv
        Q(i,j,4) = ee

        E(i,j,1) = rr*uu
        E(i,j,2) = rr*uu*uu + pp
        E(i,j,3) = rr*uu*vv
        E(i,j,4) = uu*(ee + pp)

        F(i,j,1) = rr*vv
        F(i,j,2) = rr*vv*uu
        F(i,j,3) = rr*vv*vv + pp
        F(i,j,4) = vv*(ee + pp)

      enddo
      enddo

       



C-----step 2 update Q-vector (full time step)
      resid = 0.0

!$acc  parallel loop collapse(2)
!$acc& default(none)
!$acc& private(src,rhs)
!$acc& reduction(+:resid)
       
      do i = 1, ni
      do j = 1, nj

      if (live(i,j).eq.1) then

        src(1) = eps*(q(i+1,j,1) - 2.0*q(i,j,1) + q(i-1,j,1))
     &         + eps*(q(i,j+1,1) - 2.0*q(i,j,1) + q(i,j-1,1))

        src(2) = eps*(q(i+1,j,2) - 2.0*q(i,j,2) + q(i-1,j,2))
     &         + eps*(q(i,j+1,2) - 2.0*q(i,j,2) + q(i,j-1,2))

        src(3) = eps*(q(i+1,j,3) - 2.0*q(i,j,3) + q(i-1,j,3))
     &         + eps*(q(i,j+1,3) - 2.0*q(i,j,3) + q(i,j-1,3))

        src(4) = eps*(q(i+1,j,4) - 2.0*q(i,j,4) + q(i-1,j,4))
     &         + eps*(q(i,j+1,4) - 2.0*q(i,j,4) + q(i,j-1,4))

        cs = sqrt(ga*rgas*t(i,j))
        vm = sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j))
        dt = cfl*dx/(vm + cs)

        do n = 1,4

c---------1st order
          rhs(n) = dt*(E(i+1,j,n) - E(i-1,j,n))/(2.0*dx)
     &           + dt*(F(i,j+1,n) - F(i,j-1,n))/(2.0*dy)
     &           - dt*src(n)/(dx*dx)

          Q(i,j,n) = Q(i,j,n) - rhs(n)

          resid(n) = resid(n) + abs(rhs(n))

        enddo

        r(i,j) =  Q(i,j,1)
        u(i,j) =  Q(i,j,2)/r(i,j)
        v(i,j) =  Q(i,j,3)/r(i,j)
        t(i,j) = (Q(i,j,4)/r(i,j)
     &         - 0.5*(u(i,j)*u(i,j) + v(i,j)*v(i,j)) )/cv
        p(i,j) = r(i,j)*rgas*t(i,j) ! (ideal gas law)

      endif
      enddo
      enddo


C-----step 3 tangential velocity cut-cells
!$acc  parallel loop collapse(2)
!$acc& default(none)
!$acc& private(src,rhs)
!$acc& reduction(+:resid)

      do i = 1, ni
      do j = 1, nj

      if (live(i,j).eq.0) then  ! obstacle

        ccln = 0   ! reset ccln
        psum = 0.0 ! reset psum
        rsum = 0.0 ! reset rsum
        tsum = 0.0 ! reset tsum
        usum = 0.0 ! reset usum
        vsum = 0.0 ! reset vsum

        do ii = i-1,i+1
        do jj = j-1,j+1

          if ( live(ii,jj).eq.1 ) then ! average from live cells
            ccln = ccln + 1
            rsum = rsum + r(ii,jj)
            psum = psum + p(ii,jj)
            tsum = tsum + t(ii,jj)
            usum = usum + u(ii,jj)
            vsum = vsum + v(ii,jj)
          endif

        enddo
        enddo

        if (ccln.gt.0) then ! ccln > 0
                    
          r(i,j) = rsum/real(ccln)  !  set dens at cut-cell
          p(i,j) = psum/real(ccln)  !  set pres at cut-cell
          t(i,j) = tsum/real(ccln)  !  set temp at cut-cell

          uold = usum/real(ccln)
          vold = vsum/real(ccln)

          xx = real(i-iloc)
          yy = real(j-jloc)
          rr = sqrt(xx*xx + yy*yy)

          nrm1 = xx/rr ! normal
          nrm2 = yy/rr ! normal

          dotp = nrm1*uold + nrm2*vold
          unrm = nrm1*dotp ! norm vel
          vnrm = nrm2*dotp ! norm vel

          unew1 = uold - unrm
          vnew1 = vold - vnrm

          unew2 = uold + unrm
          vnew2 = vold + vnrm

          diff1 = SQRT( unew1**2.0 + vnew1**2.0 )
          diff2 = SQRT( unew2**2.0 + vnew2**2.0 )

          IF (diff1.LT.diff2) THEN
            unew = unew1
            vnew = vnew1
          ELSE
            unew = unew2
            vnew = vnew2
          ENDIF

          mag1 = SQRT( uold**2.0 + vold**2.0 )
          mag2 = SQRT( unew**2.0 + vnew**2.0 )
          tweak = (mag1 + 0.001)/(mag2 + 0.001)

          u(i,j) = unew*tweak
          v(i,j) = vnew*tweak

        endif ! ccln gt 0
        endif ! obstacle

      enddo
      enddo


C-----step 4 smoothing
      if (.false.) then
!$acc  parallel loop collapse(2)
!$acc& default(none)
!$acc& private(src,rhs)
!$acc& reduction(+:resid)

      do i = 1, ni
      do j = 1, nj

        if (live(i,j).eq.1) then

        im1 = max(1,  i-1)
        jm1 = max(1,  j-1)
        ip1 = min(ni, i+1)
        jp1 = min(nj, j+1)

        r(i,j) = urf*r(i,j)
     &         + (1-urf)*(r(im1,j)+r(ip1,j)+r(i,jm1)+r(i,jp1))/4.0
        u(i,j) = urf*u(i,j)
     &         + (1-urf)*(u(im1,j)+u(ip1,j)+u(i,jm1)+u(i,jp1))/4.0
        v(i,j) = urf*v(i,j)
     &         + (1-urf)*(v(im1,j)+v(ip1,j)+v(i,jm1)+v(i,jp1))/4.0
        p(i,j) = urf*p(i,j)
     &         + (1-urf)*(p(im1,j)+p(ip1,j)+p(i,jm1)+p(i,jp1))/4.0
        t(i,j) = urf*t(i,j)
     &         + (1-urf)*(t(im1,j)+t(ip1,j)+t(i,jm1)+t(i,jp1))/4.0
        endif

      enddo
      enddo
      endif


C-----save the initial residuals
      if ((sa.eq.1).and.(it.le.100)) then
        resid1 = resid
CCC     write(6,*)it,(resid(n),n = 1,4)
      endif

C-----write to screen and file
      if (mod(it,nout).eq.0) then
        write(6, '(2I8,4F12.6)') sa,it,log10(resid/resid1)
        write(10,'(2I8,4F12.6)') sa,it,log10(resid/resid1)
      endif

      enddo ! its


C=============
!$acc end data
C=============


C-----write binary vtk file
      write(unit=string, fmt='(I3.3)') sa
      outfile = './results/lax-wendroff-'//string//'.vtk'

      OPEN(unit=20, file=outfile, form='unformatted',
     &  access='stream',status='replace',convert="big_endian")

      write(str1(1:10),'(i10)') ni+1
      write(str2(1:10),'(i10)') nj+1
      write(str3(1:10),'(i10)') nk+1
      write(str4(1:10),'(i10)') ni*nj*nk

      write(20)'# vtk DataFile Version 3.0'//lf
      write(20)'vtk output'//lf
      write(20)'BINARY'//lf
      write(20)'DATASET RECTILINEAR_GRID'//lf
      write(20)'DIMENSIONS '//str1//str2//str3//lf
      write(20)'X_COORDINATES '//str1//' float'//lf
      write(20)(xn(i),i=1,ni+1)
      write(20)'Y_COORDINATES '//str2//' float'//lf
      write(20)(yn(j),j=1,nj+1)
      write(20)'Z_COORDINATES '//str3//' float'//lf
      write(20)(zn(k),k=1,nk+1)
      write(20)'CELL_DATA '//str4//lf
C-----live
      write(20)'SCALARS live int'//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)((live(i,j),i=1,ni),j=1,nj)
C-----dens
      write(20)'SCALARS dens float  '//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)((r(i,j),i=1,ni),j=1,nj)
C-----temp
      write(20)'SCALARS temp float  '//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)((t(i,j),i=1,ni),j=1,nj)
C-----pres
      write(20)'SCALARS pres float  '//lf
      write(20)'LOOKUP_TABLE default'//lf
      write(20)((p(i,j),i=1,ni),j=1,nj)
C-----velocity
      write(20)'VECTORS vel float'//lf
      write(20)((u(i,j),v(i,j),0.0,i=1,ni),j=1,nj)
      close(20)



C-----write ascii vtk file
      if (.false.) then
      write(unit=string, fmt='(I3.3)') sa
      outfile = './results/lax-wendroff-'//string//'.vtk'
      OPEN(UNIT=20, FILE=outfile)
      WRITE(20,10)'# vtk DataFile Version 2.0'
      WRITE(20,10)'# sample rectilinear grid'
      WRITE(20,10)'ASCII'
      WRITE(20,10)'DATASET RECTILINEAR_GRID'
      WRITE(20,20)'DIMENSIONS ',ni+1,nj+1,2
      WRITE(20,30)'X_COORDINATES ',ni+1,' float'
      WRITE(20,*) (xn(i),i=1,ni+1)
      WRITE(20,30)'Y_COORDINATES ',nj+1,' float'
      WRITE(20,*) (yn(j),j=1,nj+1)
      WRITE(20,30)'Z_COORDINATES ',2,' float'
      WRITE(20,*) '0.0',dz
      WRITE(20,40)'CELL_DATA ',ni*nj
C-----scalar 
      WRITE(20,10)'SCALARS live float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((live(i,j),i=1,ni),j=1,nj)
C-----scalar 
      WRITE(20,10)'SCALARS dens float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((r(i,j),i=1,ni),j=1,nj)
C-----scalar 
      WRITE(20,10)'SCALARS pres float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((p(i,j),i=1,ni),j=1,nj)
C-----scalar 
      WRITE(20,10)'SCALARS temp float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((t(i,j),i=1,ni),j=1,nj)
C-----vector
      WRITE(20,10)'VECTORS velocity float'
      WRITE(20,*)((u(i,j),v(i,j),0.0,i=1,ni),j=1,nj)
      CLOSE(20)
      endif ! true/false
C-----end write ascii vtk file


      enddo ! nsav


      close(10) ! residuals
C-----iterations

   10 FORMAT(A)
   20 FORMAT(A,3I4)
   30 FORMAT(A,I4,A)
   40 FORMAT(A,I9)

      END 






