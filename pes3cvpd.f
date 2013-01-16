c      subroutine mysrf (r1,r2,theta,phi,v)
c      real*8  r1,r2,theta,phi,v
c      write(6,*) 'This is still a dymmy surface. Edit surfaces/mysurf.f'
c      stop
c      end
c
C--------------------------------------------------------------------------
C Place here your own potential energy routine (and remove the lines above).
C E.g.:
C     subroutine mysrf(r1,r2,theta,phi,v)
C        .... FORTRAN TEXT ....
C     return
C     end
C
C NOTE: You also have to modify the file source/opfuncs/usersrf.F
C       Subroutine uvpoint   and   subroutine usersurfinfo. 
C       Please edit the text accordingly.
C
C       You may name the coordinates as you like, but make sure that their 
C       ordering is appropriate. The variable v contains the energy value 
C       on return, i.e. v = V(r1,r2,theta,phi) 
C--------------------------------------------------------------------------

        subroutine pes3cvpd (dr1,dr2,dr4,r3,angle1,angle2,dihedral,v)

        implicit none

        double precision v, cart_in(5,3)
        double precision dr1,dr2,dr4,r3(3)
        double precision angle1,angle2,dihedral
        integer i, j
        integer icallnum
        save icallnum
        data icallnum /0/

        call coordtransfvpd
     +        (dr1,dr2,dr4,r3,angle1,angle2,dihedral,cart_in)

        if(icallnum .eq. 0)then
          call prepot3()
          icallnum=1
        endif

        call calcpot3(v,cart_in)

c NOTE the global minimum is at zero 

c dani TEST
c        if(v.lt.0.d0) v = 2.48d0 / 27.21138386d0 
        if(v.lt.0.d0) then     ! dani
          write(6,*) 'v= ', v  ! dani
          write(6,'(a,9(1x,f8.3))') 
     +      'coord=', dr1,dr2,dr4,r3,angle1,angle2,dihedral ! dani
c          v = 2.48d0 / 27.21138386d0
        end if ! dani

        return

        end  

        subroutine prepot3()

        implicit double precision (a-h,o-z)
        implicit integer (i-n)

        double precision V, cart_in(5,3), v0, coef(0:2066)

        double precision dc0(0:4,0:4), dw0(0:4,0:4)

        integer i,j,k,i1,j1,k1

        common/NCOE/m,mr
        common/h3o2coef/dc0,dw0,coef

        m=2022 ; mr=15

        open(20,file='h3o2.pes3C.ifchcm.dat',status='old')

        read(20,*)
        read(20,*)dc0
        read(20,*)
        read(20,*)dw0
        read(20,*)
        read(20,*)
        read(20,*)(coef(i1),i1=0,m+3*mr-1)
!        write(*,*)(coef(i1),i1=0,m+3*mr-1)
        close(20)

        return
        end ! prepot3

!*******************************************************************


        subroutine calcpot3(V,cart_in)

! this subroutine calculate the H3O2- ion potential
! input cart_in(5,3) contains cartesian coor of H O H O H
! we need to convert cart_in into cart1, the atom order in cart1 is H H H O O 

        implicit none

        integer m,mr
        double precision dc0(0:4,0:4),dw0(0:4,0:4),coef(0:2066)

        common/NCOE/m,mr
        common/h3o2coef/dc0,dw0,coef

        double precision V, cart_in(5,3), vec(0:2066)

        double precision cart0(5,3),cart1(5,3)
        integer i,j,k,l,i1,j1,k1,l1,i2,j2,k2,l2

        double precision xnuc(0:2,0:4)

        cart0=cart_in
        cart1(1,:)=cart_in(1,:)
        cart1(4,:)=cart_in(2,:)
        cart1(2,:)=cart_in(3,:)
        cart1(5,:)=cart_in(4,:)
        cart1(3,:)=cart_in(5,:)

        xnuc=transpose(cart1)

        call getvec3 (m, mr, xnuc, dc0, dw0, vec)
        V = dot_product(coef,vec)

! set C1 min as potential zero-point --- PES-3C
        V=(V+69.33697232493064d0)/219474.63067d0

! set C2 linear SP as potential zero-point --- PES-3C
!       V=(V-4.986933138719327d0)/219474.63067d0

!        V=V*219474.63067d0

        return
        end ! calcpot3
!********************************************************
        
!**********************************
        subroutine getd03 (nk, r0, dc0, dw0, d0)
        implicit none
        integer nk
        double precision r0(0:nk-1,0:nk-1), dc0(0:nk-1,0:nk-1), 
     $     dw0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1)
        integer i, j
        do i = 0, nk-1
         d0(i,i) = 0
         do j = i+1, nk-1
!          d0(i,j)=(dlog(r0(i,j))-dc0(i,j))/dw0(i,j)
          d0(i,j) = (dexp(-r0(i,j)/3.d0)-dc0(i,j))/dw0(i,j)
          d0(j,i) = d0(i,j)
         enddo
        enddo
        return
        end ! getd03
        
        subroutine getfit3 (ms, mr, dc0, dw0, coef, xn, f0, gf0, vec)
        implicit none
        integer nk, ms, mr
        parameter (nk=5)
        double precision dc0(0:nk-1,0:nk-1), dw0(0:nk-1,0:nk-1), 
     $coef(0:ms+3*mr-1),xn(0:2,0:nk-1),f0,gf0(0:2,0:nk-1),
     $vec(0:ms+3*mr-1)
        double precision  dd
        parameter (dd=1.0e-6)
        integer i, j
        double precision xn1(0:2,0:nk-1), t0, t1
!----------------------------------------------------------------------
        do i = 0, 2
         do j = 0, nk-1
          xn1 = xn ; xn1(i,j) = xn(i,j)-dd
          call getvec3 (ms, mr, xn1, dc0, dw0, vec)
          t0 = dot_product(coef,vec)
          xn1 = xn ; xn1(i,j) = xn(i,j)+dd
          call getvec3 (ms, mr, xn1, dc0, dw0, vec)
          t1 = dot_product(coef,vec)
! Note that gf0 will be the negative gradient
          gf0(i,j) = -(t1-t0)/(2*dd)
         enddo
        enddo
! f0 is last, so vec will contain a sensible return value
        call getvec3 (ms, mr, xn, dc0, dw0, vec)
        f0 = dot_product(coef,vec)
        return
        end !getvec3
        
        subroutine getr03 (nk, xn, r0)
        implicit none
        integer nk
        double precision xn(0:2,0:nk-1), r0(0:nk-1,0:nk-1)
        integer i, j
        do i = 0, nk-1
         r0(i,i) = 0
         do j = i+1, nk-1
          r0(i,j) = dsqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+ 
     $        (xn(2,j)-xn(2,i))**2)
          r0(j,i) = r0(i,j)
         enddo
        enddo
        return
        end ! getr03
        
        
        subroutine getrvec3 (m, r, vec)
        implicit none
! version for C3H2
        integer nk, m
        parameter (nk=5)
        double precision r(0:nk-1,0:nk-1), vec(0:m-1)
        integer i, j
        double precision x(0:2), r1(0:nk-1,0:nk-1), t0, t1
!-----------------------------------------------------------------------
! Test for compatibility
        if (.not.(m.eq.1.or.m.eq.4)) then
         stop 'getrvec - wrong dimension'
        endif
! Computation
        x = 0
        do i = 0, nk-1
         do j = 0, nk-1
          if (i.eq.j) then
           r1(i,j) = 0
          else
           r1(i,j) = dexp(-r(i,j))/r(i,j)
          endif
         enddo
        enddo
! CH distance
        x(0) = sum(r1(0:2,3:4))/6
! CC distance
        t0 = 0
        do i = 0, 2
         do j = i+1, 2
          t0 = t0+r1(i,j)/3
         enddo
        enddo
        x(1) = t0
! HH distance
        x(2) = r1(3,4)
! set vec
        vec(0) = 1
        if (4.le.m) then
         vec(1:3) = x
        endif
        return
        end !getrvec3
        
        
        subroutine getscale3 (n, xnuc, dc0, dw0)
        implicit none
! Version for C3H2
        integer n, nk
        parameter (nk=5)
        double precision xnuc(0:2,0:nk-1,0:n-1), dc0(0:nk-1,0:nk-1),
     $ dw0(0:nk-1,0:nk-1)
        integer ip, i, j
        double precision d0(0:nk-1,0:nk-1), dt0(0:nk-1,0:nk-1), 
     $r0(0:nk-1,0:nk-1),dc1(0:nk-1,0:nk-1), dw1(0:nk-1,0:nk-1), t0, t1
! evaluate dc0
        dt0 = 0 ; dc1 = 0 ; dw1 = 1
        do ip = 0, n-1
         call getr03 (nk, xnuc(0:2,0:nk-1,ip), r0)
         call getd03 (nk, r0, dc1, dw1, d0)
         dt0 = dt0+d0/n
        enddo
        do i = 0, nk-1
         dc0(i,i) = 0
        enddo
        t0 = 0
        do i = 0, 2
         do j = i+1, 2
          t0 = t0+dt0(i,j)/3
         enddo
        enddo
        do i = 0, 2
         do j = i+1, 2
          dc0(i,j) = t0
          dc0(j,i) = t0
         enddo
        enddo
        dc0(3,4) = dt0(3,4)
        dc0(4,3) = dt0(3,4)
        t0 = sum(dt0(0:2,3:4))/6
        dc0(0:2,3:4) = t0
        dc0(3:4,0:2) = t0
! evaluate dw0
        dt0 = 0
        do ip = 0, n-1
         call getr03 (nk, xnuc(0:2,0:nk-1,ip), r0)
         call getd03 (nk, r0, dc0, dw1, d0)
         dt0 = dt0+d0**2/n
        enddo
        t0 = 0
        do i = 0, 2
         t0 = t0+dt0(i,i)/3
        enddo
        do i = 0, 2
         dw0(i,i) = dsqrt(t0)
        enddo
        do i = 3, 4
         t0 = t0+dt0(i,i)/2
        enddo
        do i = 3, 4
         dw0(i,i) = dsqrt(t0)
        enddo
        t0 = 0
        do i = 0, 2
         do j = i+1, 2
          t0 = t0+dt0(i,j)/3
         enddo
        enddo
        do i = 0, 2
         do j = i+1, 2
          dw0(i,j) = dsqrt(t0)
          dw0(j,i) = dsqrt(t0)
         enddo
        enddo
        dw0(3,4) = dsqrt(dt0(3,4))
        dw0(4,3) = dsqrt(dt0(3,4))
        t0 = sum(dt0(0:2,3:4))/6
        dw0(0:2,3:4) = dsqrt(t0)
        dw0(3:4,0:2) = dsqrt(t0)
        return
        end !getscale3
        
        
        subroutine getsvec3 (m, d, vec)
        implicit none
        integer nk, m
        parameter (nk=5)
        double precision  d(0:nk-1,0:nk-1), vec(0:m-1)
! Version for C3H2
        integer k0, k1, k2, k3, k4, k5, k6, k7, k8
        Parameter (k0=1, k1=k0+3, k2=k1+11, k3=k2+33,k4=k3+91,k5=k4+225,
     $  k6=k5+525, k7=k6+1133, k8=k7+2321)
!! According to the Molien Series:
!! k6=k5+525, k7=k6+1133, k8=k7+2321, k9=k8+4511
        integer i, j, k, i0, i1, i2, j0, j1, j2
        double precision  d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1)
        double precision d4(0:nk-1,0:nk-1), d5(0:nk-1,0:nk-1)
        double precision d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)
        double precision  x(0:2), y(0:3), z(0:1), w(0:0), 
     $ ys(0:0), zs(0:5),us(0:7),vs(0:5),ws(0:11),ws7(0:13),ws8(0:8), t0
! There are (3,4,2,0,0,1) primaries at degrees (1..6) respectively.
! There are (1,6,8,6,12,14,9,8,5,2) secondaries at degrees (2..11).
        double precision  her2, her3, her4, her5, her6, her7
        her2(t0) = (4*t0**2-2)/dsqrt(dble(8*2))
        her3(t0) = (8*t0**2-12)*t0/dsqrt(dble(16*6))
        her4(t0) = ((16*t0**2-48)*t0**2+12)/dsqrt(dble(32*24))
        her5(t0) = ((32*t0**2-160)*t0**2+120)*t0/dsqrt(dble(64*120))
        her6(t0) = (((64*t0**2-480)*t0**2+720)*t0**2-120)/dsqrt(dble(128
     $*720))
        her7(t0) = (((128*t0**2-1344)*t0**2+3360)*t0**2-1680)*t0/dsqrt(d
     $ble(256*5040))
!-----------------------------------------------------------------------
! version for C3H2.
!-----------------------------------------------------------------------
! Test for compatibility
        if (.not.(m.eq.k0.or.m.eq.k1.or.m.eq.k2.or.m.eq.k3.or.m.eq.k4.or
     $.m.eq.k5.or.m.eq.k6.or.m.eq.k7.or.m.eq.k8)) then
         stop 'getsvec - wrong dimension'
        endif
! auxiliary distances
        do i = 0, nk-1
         do j = i+1, nk-1
          d2(i,j) = her2(d(i,j))
          d2(j,i) = d2(i,j)
          d3(i,j) = her3(d(i,j))
          d3(j,i) = d3(i,j)
          d4(i,j) = her4(d(i,j))
          d4(j,i) = d4(i,j)
          d5(i,j) = her5(d(i,j))
          d5(j,i) = d5(i,j)
          d6(i,j) = her6(d(i,j))
          d6(j,i) = d6(i,j)
          d7(i,j) = her7(d(i,j))
          d7(j,i) = d7(i,j)
         enddo
        enddo
! Primary Invariants
        x = 0 ; y = 0 ; z = 0 ; w = 0
        x(0) = (d(0,1)+d(0,2)+d(1,2))/3
        x(1) = sum(d(0:2,3:4))/6
        x(2) = d(3,4)
        y(0) = (d2(0,1)+d2(0,2)+d2(1,2))/3
        y(1) = sum(d2(0:2,3:4))/6
        z(0) = (d3(0,1)+d3(0,2)+d3(1,2))/3
        w(0) = sum(d6(0:2,3:4))/6
        do i0 = 0, 2
         y(2) = y(2)+her2(sum(d(i0,3:4))/2)/3
         z(1) = z(1)+her3(sum(d(i0,3:4))/2)/3
        enddo
        do j0 = 3, 4
         y(3) = y(3)+her2(sum(d(0:2,j0))/3)/2
        enddo
! Secondary Invariants
!! us(0:7), vs(0:5), ws(0:11), ws7(0:13), ws8(0:8)
        ys = 0 ; zs = 0 ; us = 0 ; vs = 0 ; ws = 0 ; ws7 = 0 ; ws8 = 0
        do i0 = 0, 2
         zs(0) = zs(0)+(d3(i0,3)+d3(i0,4))/6
         us(0) = us(0)+(d4(i0,3)+d4(i0,4))/6
         us(1) = us(1)+(d3(i0,3)*d(i0,4)+d(i0,3)*d3(i0,4))/6
         vs(0) = vs(0)+(d5(i0,3)+d5(i0,4))/6
         do i1 = 0, 2
         if (i1.ne.i0) then
          zs(1) = zs(1)+d(i0,i1)*d(i0,3)*d(i0,4)/3
          do j0 = 3, 4
           ys(0) = ys(0)+d(i0,i1)*d(i0,j0)/6
           zs(2) = zs(2)+d2(i0,i1)*d(i0,j0)/6
           zs(3) = zs(3)+d(i0,i1)*d2(i0,j0)/6
           zs(4) = zs(4)+d(i0,i1)*d(i0,j0)*d(i1,j0)/6
           zs(5) = zs(5)+d2(i0,j0)*d(i1,j0)/6
           us(2) = us(2)+d(i0,i1)*d3(i0,j0)/6
           us(3) = us(3)+d2(i0,i1)*d2(i0,j0)/6
           us(4) = us(4)+d3(i0,j0)*d(i1,j0)/6
           us(5) = us(5)+d(i0,i1)*d2(i0,j0)*d(i1,j0)/6
           us(6) = us(6)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/6
           vs(1) = vs(1)+d(i0,i1)*d4(i0,j0)/6
          enddo
         endif
         enddo
        enddo
        us(7) = her2(ys(0))
        vs(2) = ys(0)*zs(0)
        vs(3) = ys(0)*zs(3)
        vs(4) = ys(0)*zs(4)
        vs(5) = ys(0)*zs(5)
        ws(0) = ys(0)*us(0)
        ws(1) = ys(0)*us(1)
        ws(2) = ys(0)*us(3)
        ws(3) = ys(0)*us(4)
        ws(4) = ys(0)*us(6)
        ws(5) = her2(zs(0))
        ws(6) = zs(0)*zs(2)
        ws(7) = zs(0)*zs(3)
        ws(8) = zs(1)*zs(2)
        ws(9) = zs(2)*zs(5)
        ws(10) = her2(zs(3))
        ws(11) = zs(3)*zs(5)
        ws7(0) = ys(0)*vs(0)
        ws7(1) = ys(0)*vs(1)
        ws7(2) = zs(0)*us(0)
        ws7(3) = zs(0)*us(2)
        ws7(4) = zs(0)*us(3)
        ws7(5) = zs(1)*us(3)
        ws7(6) = zs(1)*us(6)
        ws7(7) = zs(3)*us(0)
        ws7(8) = zs(3)*us(2)
        ws7(9) = zs(4)*us(2)
        ws7(10) = zs(4)*us(3)
        ws7(11) = zs(5)*us(0)
        ws7(12) = zs(5)*us(2)
        ws7(13) = zs(5)*us(3)
        ws8(0) = zs(0)*vs(0)
        ws8(1) = zs(0)*vs(1)
        ws8(2) = zs(1)*vs(1)
        ws8(3) = zs(2)*vs(1)
        ws8(4) = zs(3)*vs(0)
        ws8(5) = zs(3)*vs(1)
        ws8(6) = zs(4)*vs(1)
        ws8(7) = zs(5)*vs(1)
        ws8(8) = us(0)*us(3)
! Compute vec(0:m-1)
        vec = 0
! constant term
        vec(0) = 1
! first degree terms
        if (k1.le.m) then
         vec(k0:k0+2) = x(0:2)
        endif
! second degree terms
        if (k2.le.m) then
         vec(k1) = her2(x(0))
         vec(k1+1:k1+2) = x(0)*x(1:2)
         vec(k1+3) = her2(x(1))
         vec(k1+4) = x(1)*x(2)
         vec(k1+5) = her2(x(2))
         vec(k1+6:k1+9) = y(0:3)
         vec(k1+10) = ys(0)
        endif
! third degree terms
        if (k3.le.m) then
         vec(k2) = her3(x(0))
         vec(k2+1:k2+2) = her2(x(0))*x(1:2)
         vec(k2+3:k2+10) = x(0)*vec(k1+3:k1+10)
         vec(k2+11) = her3(x(1))
         vec(k2+12) = her2(x(1))*x(2)
         vec(k2+13:k2+18) = x(1)*vec(k1+5:k1+10)
         vec(k2+19) = her3(x(2))
         vec(k2+20:k2+24) = x(2)*vec(k1+6:k1+10)
         vec(k2+25:k2+26) = z(0:1)
         vec(k2+27:k2+32) = zs(0:5)
        endif
! fourth degree terms
        if (k4.le.m) then
         vec(k3) = her4(x(0))
         vec(k3+1:k3+2) = her3(x(0))*x(1:2)
         vec(k3+3:k3+10) = her2(x(0))*vec(k1+3:k1+10)
         vec(k3+11:k3+32) = x(0)*vec(k2+11:k2+32)
         vec(k3+33) = her4(x(1))
         vec(k3+34) = her3(x(1))*x(2)
         vec(k3+35:k3+40) = her2(x(1))*vec(k1+5:k1+10)
         vec(k3+41:k3+54) = x(1)*vec(k2+19:k2+32)
         vec(k3+55) = her4(x(2))
         vec(k3+56:k3+60) = her2(x(2))*vec(k1+6:k1+10)
         vec(k3+61:k3+68) = x(2)*vec(k2+25:k2+32)
         k = k3+69
         do i = 0, 3
          vec(k) = her2(y(i))
          k = k+1
          vec(k:k+3-i) = y(i)*vec(k1+7+i:k1+10)
          k = k+4-i
         enddo
         if (k.ne.k3+83) then
          stop 'getsvec count failure (4)'
         endif
         vec(k3+83:k3+90) = us(0:7)
        endif
! fifth degree terms
        if (k5.le.m) then
         vec(k4) = her5(x(0))
         vec(k4+1:k4+2) = her4(x(0))*x(1:2)
         vec(k4+3:k4+10) = her3(x(0))*vec(k1+3:k1+10)
         vec(k4+11:k4+32) = her2(x(0))*vec(k2+11:k2+32)
         vec(k4+33:k4+90) = x(0)*vec(k3+33:k3+90)
         vec(k4+91) = her5(x(1))
         vec(k4+92) = her4(x(1))*x(2)
         vec(k4+93:k4+98) = her3(x(1))*vec(k1+5:k1+10)
         vec(k4+99:k4+112) = her2(x(1))*vec(k2+19:k2+32)
         vec(k4+113:k4+148) = x(1)*vec(k3+55:k3+90)
         vec(k4+149) = her5(x(2))
         vec(k4+150:k4+154) = her3(x(2))*vec(k1+6:k1+10)
         vec(k4+155:k4+162) = her2(x(2))*vec(k2+25:k2+32)
         vec(k4+163:k4+184) = x(2)*vec(k3+69:k3+90)
         k = k4+185
         do i = 0, 3
          vec(k:k+7) = y(i)*vec(k2+25:k2+32)
          k = k+8
         enddo
         do i = 0, 1
          vec(k) = z(i)*vec(k1+10)
          k = k+1
         enddo
         if (k.ne.k4+219) then
          stop 'getsvec count failure (5)'
         endif
         vec(k4+219:k4+224) = vs(0:5)
        endif
! sixth degree terms
        if (k6.le.m) then
         vec(k5) = her6(x(0))
         vec(k5+1:k5+2) = her5(x(0))*x(1:2)
         vec(k5+3:k5+10) = her4(x(0))*vec(k1+3:k1+10)
         vec(k5+11:k5+32) = her3(x(0))*vec(k2+11:k2+32)
         vec(k5+33:k5+90) = her2(x(0))*vec(k3+33:k3+90)
         vec(k5+91:k5+224) = x(0)*vec(k4+91:k4+224)
         vec(k5+225) = her6(x(1))
         vec(k5+226) = her5(x(1))*x(2)
         vec(k5+227:k5+232) = her4(x(1))*vec(k1+5:k1+10)
         vec(k5+233:k5+246) = her3(x(1))*vec(k2+19:k2+32)
         vec(k5+247:k5+282) = her2(x(1))*vec(k3+55:k3+90)
         vec(k5+283:k5+358) = x(1)*vec(k4+149:k4+224)
         vec(k5+359) = her6(x(2))
         vec(k5+360:k5+364) = her4(x(2))*vec(k1+6:k1+10)
         vec(k5+365:k5+372) = her3(x(2))*vec(k2+25:k2+32)
         vec(k5+373:k5+394) = her2(x(2))*vec(k3+69:k3+90)
         vec(k5+395:k5+434) = x(2)*vec(k4+185:k4+224)
         vec(k5+435) = her3(y(0))
         vec(k5+436:k5+439) = her2(y(0))*vec(k1+7:k1+10)
         vec(k5+440:k5+456) = y(0)*vec(k3+74:k3+90)
         vec(k5+457) = her3(y(1))
         vec(k5+458:k5+460) = her2(y(1))*vec(k1+8:k1+10)
         vec(k5+461:k5+473) = y(1)*vec(k3+78:k3+90)
         vec(k5+474) = her3(y(2))
         vec(k5+475:k5+476) = her2(y(2))*vec(k1+9:k1+10)
         vec(k5+477:k5+486) = y(2)*vec(k3+81:k3+90)
         vec(k5+487) = her3(y(3))
         vec(k5+488) = her2(y(3))*vec(k1+10)
         vec(k5+489:k5+496) = y(3)*vec(k3+83:k3+90)
         vec(k5+497) = her2(z(0))
         vec(k5+498:k5+504) = z(0)*vec(k2+26:k2+32)
         vec(k5+505) = her2(z(1))
         vec(k5+506:k5+511) = z(1)*vec(k2+27:k2+32)
         vec(k5+512) = w(0)
         vec(k5+513:k5+524) = ws(0:11)
        endif
! seventh degree terms
        if (k7.le.m) then
         vec(k6) = her7(x(0))
         vec(k6+1:k6+2) = her6(x(0))*x(1:2)
         vec(k6+3:k6+10) = her5(x(0))*vec(k1+3:k1+10)
         vec(k6+11:k6+32) = her4(x(0))*vec(k2+11:k2+32)
         vec(k6+33:k6+90) = her3(x(0))*vec(k3+33:k3+90)
         vec(k6+91:k6+224) = her2(x(0))*vec(k4+91:k4+224)
         vec(k6+225:k6+524) = x(0)*vec(k5+225:k5+524)
         vec(k6+525) = her7(x(1))
         vec(k6+526) = her6(x(1))*x(2)
         vec(k6+527:k6+532) = her5(x(1))*vec(k1+5:k1+10)
         vec(k6+533:k6+546) = her4(x(1))*vec(k2+19:k2+32)
         vec(k6+547:k6+582) = her3(x(1))*vec(k3+55:k3+90)
         vec(k6+583:k6+658) = her2(x(1))*vec(k4+149:k4+224)
         vec(k6+659:k6+824) = x(1)*vec(k5+359:k5+524)
         vec(k6+825) = her7(x(2))
         vec(k6+826:k6+830) = her5(x(2))*vec(k1+6:k1+10)
         vec(k6+831:k6+838) = her4(x(2))*vec(k2+25:k2+32)
         vec(k6+839:k6+860) = her3(x(2))*vec(k3+69:k3+90)
         vec(k6+861:k6+900) = her2(x(2))*vec(k4+185:k4+224)
         vec(k6+901:k6+990) = x(2)*vec(k5+435:k5+524)
         vec(k6+991:k6+998) = her2(y(0))*vec(k2+25:k2+32)
         vec(k6+999:k6+1030) = y(0)*vec(k4+193:k4+224)
         vec(k6+1031:k6+1038) = her2(y(1))*vec(k2+25:k2+32)
         vec(k6+1039:k6+1062) = y(1)*vec(k4+201:k4+224)
         vec(k6+1063:k6+1070) = her2(y(2))*vec(k2+25:k2+32)
         vec(k6+1071:k6+1086) = y(2)*vec(k4+209:k4+224)
         vec(k6+1087:k6+1094) = her2(y(3))*vec(k2+25:k2+32)
         vec(k6+1095:k6+1102) = y(3)*vec(k4+217:k4+224)
         vec(k6+1103:k6+1110) = z(0)*vec(k3+83:k3+90)
         vec(k6+1111:k6+1118) = z(1)*vec(k3+83:k3+90)
         vec(k6+1119:k6+1132) = ws7(0:13)
        endif
        return
        end ! getsvec3
        
        
        subroutine getvec3 (ms, mr, xn, dc0, dw0, vec)
        implicit none
! version for C3H2
        integer nk, ms, mr
        parameter (nk=5)
        double precision  xn(0:2,0:nk-1), dc0(0:nk-1,0:nk-1)
        double precision dw0(0:nk-1,0:nk-1), vec(0:ms+3*mr-1)
        integer k, l
        double precision  rvec(0:3),d0(0:nk-1,0:nk-1),r0(0:nk-1,0:nk-1)
        call getr03 (nk, xn, r0)
        call getd03 (nk, r0, dc0, dw0, d0)
        call getsvec3 (ms, d0, vec(0:ms-1))
        call getrvec3 (4, r0, rvec)
        do l = 0, mr-1
         do k = 0, 2
          vec(ms+3*l+k) = rvec(k+1)*vec(l)
         enddo
        enddo
        return
        end  ! getvec3



c---------------------------------------------------
c---------------------------------------------------
c  Subroutine of coord transf fo H2O3-
c---------------------------------------------------
c---------------------------------------------------

        subroutine coordtransfvpd
     +            (dr1,dr2,dr4,r3,angle1,angle2,dihedral,cart)

        implicit none

c----------------------------------------------
c   Parameters 
c    (1 amu / 1822.88839 me)
c----------------------------------------------
        
        real*8 amu_in_me, bohr_in_ang, dd
        real*8 Hamu, Oamu, massH, massO
        
        parameter (amu_in_me = 5.4857990945d-4)  
        parameter (Hamu = 1.007825032d0)
        parameter (Oamu = 15.99491462d0)
        parameter (massH = Hamu / amu_in_me)
        parameter (massO = Oamu / amu_in_me)
        parameter (bohr_in_ang = 1.d0 / 0.52917720859d0)
        parameter (dd=1.746d0)

c----------------------------------------------
c   END of Parameters
c----------------------------------------------

        real*8 cart(5,3), cart1(3), cart2(3)
        real*8 oh_com1(3), oh_com2(3)
        real*8 rmat(3,3), tt1(3), tt2(3)
        real*8 dr1, dr2, dr4, r3(3)
        real*8 angle1, angle2, dihedral
        real*8 ppi, mass(5), zpp
        integer i, j, k
        logical lprintcoord, loo

c----------------------------------------------
c   Useful definitions:
c----------------------------------------------

        ppi= acos(-1.d0)

! if loo, dr4 is the distance between both oxygens 

c        loo = .true.
        loo = .false.

c----------------------------------------------
c initialize coordinates
c   (distances in bohrs)
c   (angles in radians)
c----------------------------------------------

        do j = 1, 3
          do i = 1, 5
            cart(i,j) = 0.d0
          end do 
        end do 

        do i = 1, 3
          oh_com1(i) = 0.d0
          oh_com2(i) = 0.d0
          cart1(i)   = 0.d0
          cart2(i)   = 0.d0
        end do 

c----------------------------------------------
c assign masses
c----------------------------------------------

        mass(1) = massH
        mass(2) = massO
        mass(3) = massH
        mass(4) = massO
        mass(5) = massH

c----------------------------------------------
c place bridging H(1)
c----------------------------------------------

        do i = 1, 3
          cart(1,i) = r3(i)
        enddo

c----------------------------------------------
c   Place vector O2H3 along z axis
c    (H towards negative direction z-axis)
c----------------------------------------------

c Oxygen at dr1 from origin

       cart(2,3) = -dr1  

c Hydrogen at the origin

       cart(3,3) = 0.d0 

c Compute its centre of masses (oh_com1)

        do i = 1, 3
          cart1(i) = cart(2,i)
          cart2(i) = cart(3,i)
        enddo

       if(.not.loo)
     +    call com_oh3(oh_com1,cart1,cart2,mass(2),mass(3))

       if(loo) then

c Shift O2 to origin and define H3 wrt O2

        do i = 1, 3
          cart(3,i) = cart2(i) - cart(2,i) ! H wrt O
          cart(2,i) = 0.d0                 ! O to origin
          cart1(i)  = 0.d0
          cart2(i)  = 0.d0
        enddo

       else

c Shift COM1 to origin

          do i = 1, 3
            cart(2,i) = cart1(i) - oh_com1(i)
            cart(3,i) = cart2(i) - oh_com1(i)
            cart1(i)  = 0.d0
            cart2(i)  = 0.d0
          enddo

       end if

c   Rotate O2H3 around y-axis by angle1
c   (Recall rotation counterclockwise)

       do i = 1, 3
         do j = 1, 3
           rmat(i,j) = 0.d0
         end do 
       end do 

       call rotmaty3(angle1, rmat)

       do i = 1, 3
         tt1(i) = 0.d0
         tt2(i) = 0.d0
         do j = 1, 3
           tt1(i) = tt1(i) + rmat(i,j) * cart(2,j)
           tt2(i) = tt2(i) + rmat(i,j) * cart(3,j)
         end do
       end do

       do i = 1, 3
         cart(2,i) = tt1(i)
         cart(3,i) = tt2(i)
       end do

c Displace from origin

       cart(2,3) = cart(2,3) - dr4 / 2.d0
       cart(3,3) = cart(3,3) - dr4 / 2.d0

c----------------------------------------------
c   Place vector O4H5 along z axis
c----------------------------------------------

c Oxygen at dr2 from origin

       cart(4,3) = -dr2  

c Hydrogen at the origin

       cart(5,3) = 0.d0 

c Compute its centre of masses (oh_com2)

       do i = 1, 3
         cart1(i) = cart(4,i)
         cart2(i) = cart(5,i)
       enddo

       if(.not.loo)
     +   call com_oh3(oh_com2,cart1,cart2,mass(4),mass(5))

       do i = 1, 3
         cart(4,i) = cart1(i)
         cart(5,i) = cart2(i)
       enddo

       if(loo) then

c Shift O4 to origin and define H5 wrt to O

         do i = 1, 3
           cart(5,i) = cart2(i) - cart(4,i) ! H wrt O
           cart(4,i) = 0.d0                 ! O to origin
           cart1(i)  = 0.d0
           cart2(i)  = 0.d0
         enddo

       else

c Shift COM2 to origin

         do i = 1, 3
           cart(4,i) = cart1(i) - oh_com2(i)
           cart(5,i) = cart2(i) - oh_com2(i)
           cart1(i)  = 0.d0
           cart2(i)  = 0.d0
         end do

       end if

c   Rotate O4H5 around y-axis by angle2

       do i = 1, 3
         do j = 1, 3
          rmat(i,j) = 0.d0
         end do
       end do

       call rotmaty3(angle2, rmat)

       do i = 1, 3
         tt1(i) = 0.d0
         tt2(i) = 0.d0
         do j = 1, 3
           tt1(i) = tt1(i) + rmat(i,j) * cart(4,j)
           tt2(i) = tt2(i) + rmat(i,j) * cart(5,j)
         end do
       end do

       do i = 1, 3
         cart(4,i) = tt1(i)
         cart(5,i) = tt2(i)
       end do

c Displace from origin

       cart(4,3)   = cart(4,3) + dr4 / 2.d0
       cart(5,3)   = cart(5,3) + dr4 / 2.d0

c   Rotate O4H5 to set 0. dihedral wrt to O2H3
c
c       do i = 1, 3
c         do j = 1, 3
c          rmat(i,j) = 0.d0
c         end do
c       end do
c
c       call rotmatz3(ppi, rmat)
c
c       do i = 1, 3
c         tt1(i) = 0.d0
c         tt2(i) = 0.d0
c         do j = 1, 3
c           tt1(i) = tt1(i) + rmat(i,j) * cart(4,j)
c           tt2(i) = tt2(i) + rmat(i,j) * cart(5,j)
c         end do
c       end do
c
c       do i = 1, 3
c         cart(4,i) = tt1(i)
c         cart(5,i) = tt2(i)
c       end do
c
c   Rotate O4H5 around z-axis by dihedral

       dihedral = - dihedral

       do i = 1, 3
         do j = 1, 3
          rmat(i,j) = 0.d0
         end do
       end do

       call rotmatz3(dihedral, rmat)

       do i = 1, 3
         tt1(i) = 0.d0
         tt2(i) = 0.d0
         do j = 1, 3
           tt1(i) = tt1(i) + rmat(i,j) * cart(4,j)
           tt2(i) = tt2(i) + rmat(i,j) * cart(5,j)
         end do
       end do

       do i = 1, 3
         cart(4,i) = tt1(i)
         cart(5,i) = tt2(i)
       end do

       lprintcoord = .false.
c       lprintcoord = .true.

c       if(lcutalong.or.lspenergy) then

       if(lprintcoord) then

c  Centre of masses O2H3

         do i = 1, 3
           cart1(i) = cart(2,i)
           cart2(i) = cart(3,i)
         enddo
 
         call com_oh3(oh_com1,cart1,cart2,mass(2),mass(3))

c  Centre of masses O2H3

         do i = 1, 3
           cart1(i) = cart(4,i)
           cart2(i) = cart(5,i)
         enddo
  
         if(.not.loo)
     +     call com_oh3(oh_com2,cart1,cart2,mass(4),mass(5))

         zpp = r3(3) / (dr4 - 2.d0*dd) 

         write(6,'(/a)')'-------------------------------------'
         write(6,'(/a)')'Cartesian coordinates '
         write(6,'(a)') ' X correspond to:' 
         write(6,'(a)') ' (1) O2H2 Centre of Masses'
         write(6,'(a)') ' (2) O2H3 Centre of Masses'
         write(6,'(a)') ' (2) O4H5 Centre of Masses'
         write(6,'(/a)') '(Dist= AA, Angles= Degrees)'
         write(6,'(a)')'-------------------------------------'

         if(loo) then
           write(6,'(1x,i2)') 6
         else
           write(6,'(1x,i2)') 8
         end if

         write(6,*)
         write(6,'(a,3(1x,f12.6))') ' X' , 0.d0, 0.d0, 0.d0

         if(.not.loo) then
           write(6,'(a,3(1x,f12.6))') 
     +     ' X' , (oh_com1(i) / bohr_in_ang, i=1,3)
           write(6,'(a,3(1x,f12.6))') 
     +     ' X' , (oh_com2(i) / bohr_in_ang, i=1,3)
         end if

         write(6,'(a,3(1x,f12.6))')
     +   ' H' , (cart(1,i) / bohr_in_ang, i=1,3)
         write(6,'(a,3(1x,f12.6))') 
     +   ' O' , (cart(2,i) / bohr_in_ang, i=1,3)
         write(6,'(a,3(1x,f12.6))') 
     +   ' H' , (cart(3,i) / bohr_in_ang, i=1,3)
         write(6,'(a,3(1x,f12.6))') 
     +   ' O' , (cart(4,i) / bohr_in_ang, i=1,3)
         write(6,'(a,3(1x,f12.6))') 
     +   ' H' , (cart(5,i) / bohr_in_ang, i=1,3)

       end if

       return

       end subroutine

c----------------------------------------------
c  Centre of masses (OH)
c----------------------------------------------

        subroutine com_oh3(com, cart1, cart2, mass1, mass2)

          implicit none

          real*8   mass1, mass2
          real*8   cart1(3), cart2(3)
          real*8   com(3), xcom(3)
          real*8   totmass
          integer  i, j

          totmass = 0.d0
          xcom    = 0.d0

          totmass = mass1 + mass2

          xcom(1) = mass1 * cart1(1) + mass2 * cart2(1)
          xcom(2) = mass1 * cart1(2) + mass2 * cart2(2)
          xcom(3) = mass1 * cart1(3) + mass2 * cart2(3)

          xcom = xcom / totmass

          com = xcom

        end subroutine com_oh3

c----------------------------------------------
c  Rotation matrices are generated here
c----------------------------------------------

      subroutine rotmaty3(ang, rmat)

        implicit none

        real*8  rmat(3,3)
        real*8  ang

        rmat(1,1) =  cos(ang)
        rmat(1,3) =  sin(ang)

        rmat(2,2) = 1.d0

        rmat(3,1) = -sin(ang)
        rmat(3,3) =  cos(ang)

        return

      end subroutine 

      subroutine rotmatz3(ang, rmat)

        implicit none

        real*8  rmat(3,3)
        real*8  ang

        rmat(1,1) =  cos(ang)
        rmat(1,2) =  sin(ang)

        rmat(2,1) = -sin(ang)
        rmat(2,2) =  cos(ang)

        rmat(3,3) = 1.d0

        return

      end subroutine 



