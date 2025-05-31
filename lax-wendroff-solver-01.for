
      program main
      implicit none

      integer, parameter :: ni=200,nj=200
      integer live(0:ni+1,0:nj+1)

      real xn(ni+1),yn(nj+1)
      real r(ni,nj),u(ni,nj),v(ni,nj)
      real p(ni,nj),t(ni,nj)
      real Q(0:ni+1,0:nj+1,4)
      real E(0:ni+1,0:nj+1,4)
      real F(0:ni+1,0:nj+1,4)
      real Qf(4),  sb(4)
      real Ef(ni+1,nj,4)
      real Ff(ni,nj+1,4)
      real rhs(4),resid(4),resid1(4)

      integer i,j,n,it,sa,nits,nsav,gc
      integer ccln,ii,jj,nout,iloc,jloc

      real dx,dy,dz,dt,dtdx,dtdy
      real uu,vv,rr,pp,tt,ee,eps
      real rinf,uinf,vinf,pinf,tinf
      real rgas,cp,cv,kth,ga,ma,cs
      real ptot,ttot,einf,sigma

      real unrm,vnrm,unew,vnew,tweak
      real unew1,vnew1,unew2,vnew2
      real usum,vsum,rsum,tsum,psum
      real nrm1,nrm2,mag1,mag2,dotp
      real diff1,diff2,uold,vold,xx,yy

      character(31) outfile
      character(3) string

C-----specify constants
      dx = 0.1      ! spacing
      dy = 0.1      ! spacing
      dz = 0.1      ! spacing
      dt = 0.00001  ! timestep
      eps = 1.0     ! smoothing (0.05)

C-----air properties
      rgas = 287.1
      cp   = 1004.5
      cv   = 717.5
      ga   = 1.4
      kth  = 0.026 

      rinf = 1.225 
      uinf = 10.0
      vinf = 0.0
      pinf = 101325.0
      tinf = 288.0
      einf = rinf*cv*tinf + 0.5*rinf*(uinf*uinf + vinf*vinf)

      ptot = 101386.25  ! pstatic + hrv2
      ttot = 288.04877  ! tstatic + hv2/cp

      nsav = 40
      nits = 1000
      nout = 100

      dtdx = dt/dx
      dtdy = dt/dy

     
C-----specify live/dead regions (if needed)
      live = 1   ! everything live
      iloc = 100 ! blockage centre
      jloc = 100 ! blockage centre
       
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

      do i = 2, ni+1
        xn(i) = xn(i-1) + dx
      enddo

      do j = 2, nj+1
        yn(j) = yn(j-1) + dy
      enddo


C-----initialise arrays with height pulse
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

      do it = 1,nits
      
C-----step1 - vectors at cell centres
      do i = 0,ni+1
      do j = 0,nj+1

C-------bottom wall
        if ((j.eq.0).and.(i.ge.1).and.(i.le.ni)) then
          rr = r(i,1)
          uu = u(i,1)
          vv = v(i,1)
          pp = p(i,1)
          tt = t(i,1)
        endif

C-------top wall
        if ((j.eq.nj+1).and.(i.ge.1).and.(i.le.ni)) then
          rr = r(i,nj)
          uu = u(i,nj)
          vv = v(i,nj)
          pp = p(i,nj)
          tt = t(i,nj)
        endif

C-------left wall
        if ((i.eq.0).and.(j.ge.1).and.(j.le.nj)) then
          uu = u(1,j)
          vv = vinf
          pp = ptot - 0.5*rinf*(uu*uu + vv*vv)   ! bernoulli
          tt = ttot - 0.5*(uu*uu + vv*vv)/cp   ! bernoulli
          rr = pp/(rgas*tt)

c           pp = p(1,j)
c           tt = tinf
c           uu = uinf
c           vv = vinf
c           rr = pp/(rgas*tt)
        endif

C-------right wall
        if ((i.eq.ni+1).and.(j.ge.1).and.(j.le.nj)) then
          pp = pinf
          tt = t(ni,j)
          uu = u(ni,j)
          vv = v(ni,j)
          rr = pp/(rgas*tt)
        endif

C-------main block of inner cells
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

       

C-----step3 i-face vectors (half step)
      do j = 1, nj   ! j-cells
      do i = 1, ni+1 ! i-faces

      if ((live(i,j).eq.1).or.(live(i-1,j).eq.1)) then

        do n = 1,4
          Qf(n) = 0.5*(Q(i,j,n) + Q(i-1,j,n))
     &     - 0.5*dtdx*(E(i,j,n) - E(i-1,j,n))
        enddo

        rr =  Qf(1)
        uu =  Qf(2)/rr
        vv =  Qf(3)/rr
        tt = (Qf(4)/rr - 0.5*(uu*uu + vv*vv))/cv

        pp = rr*rgas*tt ! (ideal gas law)
        ee = rr*cv*tt + 0.5*rr*(uu*uu + vv*vv)

        Ef(i,j,1) = rr*uu
        Ef(i,j,2) = rr*uu*uu + pp
        Ef(i,j,3) = rr*uu*vv
        Ef(i,j,4) = uu*(ee + pp)

      endif

      enddo
      enddo


C-----step4 j-face vectors (half step)
      do i = 1, ni   ! i-cells 
      do j = 1, nj+1 ! j-faces 1,nj+1

      if ((live(i,j).eq.1).or.(live(i,j-1).eq.1)) then

        do n = 1,4
          Qf(n) = 0.5*(Q(i,j,n) + Q(i,j-1,n))
     &     - 0.5*dtdy*(F(i,j,n) - F(i,j-1,n))
        enddo

        rr =  Qf(1)
        uu =  Qf(2)/rr
        vv =  Qf(3)/rr
        tt = (Qf(4)/rr - 0.5*(uu*uu + vv*vv))/cv

        pp = rr*rgas*tt ! (ideal gas law)
        ee = rr*cv*tt + 0.5*rr*(uu*uu + vv*vv)

        Ff(i,j,1) = rr*vv
        Ff(i,j,2) = rr*vv*uu
        Ff(i,j,3) = rr*vv*vv + pp
        Ff(i,j,4) = vv*(ee + pp)

      endif

      enddo
      enddo


C-----step5 Q-vector (full step)
      resid = 0.0
       
      do i = 1, ni
      do j = 1, nj

      if (live(i,j).eq.1) then


        sb = 0.0

C-------sponge

        sigma = 0.0
c        if (i.le.40)  sigma = 1.0*(40-i)/40.0
c        if (i.ge.160) sigma = 1.0*(i-160)/40.0

          rr = r(i,j)
          uu = u(i,j)
          vv = v(i,j)
          pp = p(i,j)
          tt = t(i,j)
          ee = rr*cv*tt + 0.5*rr*(uu*uu + vv*vv)

        sb(1) = eps*(q(i+1,j,1) - 2.0*q(i,j,1) + q(i-1,j,1))
     &        + eps*(q(i,j+1,1) - 2.0*q(i,j,1) + q(i,j-1,1))
     &        - sigma*(rr - rinf)

        sb(2) = eps*(q(i+1,j,2) - 2.0*q(i,j,2) + q(i-1,j,2))
     &        + eps*(q(i,j+1,2) - 2.0*q(i,j,2) + q(i,j-1,2))
     &        - sigma*(rr*uu - rinf*uinf)

        sb(3) = eps*(q(i+1,j,3) - 2.0*q(i,j,3) + q(i-1,j,3))
     &        + eps*(q(i,j+1,3) - 2.0*q(i,j,3) + q(i,j-1,3))
     &        - sigma*(rr*vv - rinf*vinf)


        sb(4) = eps*(q(i+1,j,4) - 2.0*q(i,j,4) + q(i-1,j,4))
     &        + eps*(q(i,j+1,4) - 2.0*q(i,j,4) + q(i,j-1,4))
     &        - sigma*(ee - einf)

        do n = 1,4

c         1st order
          rhs(n) = dtdx*(E(i+1,j,n) - E(i-1,j,n))/2.0
     &           + dtdy*(F(i,j+1,n) - F(i,j-1,n))/2.0
     &           - dt*sb(n)/(dx*dx)

c         2nd order
c          rhs(n) = dtdx*(Ef(i+1,j,n) - Ef(i,j,n))
c     &           + dtdy*(Ff(i,j+1,n) - Ff(i,j,n))
c     &           - dt*sb(n)/(dx*dx)

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

          nrm1 = xx/rr ! surf normal
          nrm2 = yy/rr ! surf normal

          dotp = nrm1*uold + nrm2*vold
          unrm = nrm1*dotp ! surf norm vel
          vnrm = nrm2*dotp ! surf norm vel

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

C-----write vtk file
      write(unit=string, fmt='(I3.3)') sa
      outfile = './results/shallow-water-'//string//'.vtk'
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
C-----end write vtk file

      enddo ! nsav


      close(10) ! residuals
C-----iterations

   10 FORMAT(A)
   20 FORMAT(A,3I4)
   30 FORMAT(A,I4,A)
   40 FORMAT(A,I9)

      END 






