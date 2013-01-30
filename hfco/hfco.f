      Block Data
      Implicit Real*8(a-h,o-z)
      Logical Init
      Common /Initialize/ xJaco(6),Init
      Data Init/.true./
      End

      subroutine hfco(r1,r2,r3,o1,o2,phi,v)
c==================================================
c     r1 distance C-H   r2 distance C-F  r3 distance C=O 
c     o1 angle H-C-O    o2 angle F-C-O 
c     phi angle diedre H hors du plan
c     Atome 1 : C , 2 : O , 3 : H , 4 : F 
c==================================================
      Implicit Real*8(a-h,o-z)
      Real*8 x(4,3)
      Logical Init
      Common /Initialize/ xJaco(6),Init
      Data seuil/0.0001/
      Data aumass/1822.802/,autoeV/27.21183/,cmtoau/4.55633e-6/
     &     ,autoag/0.529d0/,zero/0d0/,rad/0.01745329252/
      
      xJaco(1)=r1
      xJaco(2)=r2
      xJaco(3)=r3
      xJaco(4)=o1
      xJaco(5)=o2
      xJaco(6)=phi


c     coordonnees cartesiennes des 4 atomes a partir
c     de valence

      x(1,1)=zero
      x(1,2)=zero
      x(1,3)=zero
      x(2,1)=zero
      x(2,2)=zero
      x(2,3)=r3
      x(3,1)=r1*sin(acos(o1))*cos(phi)
      x(3,2)=r1*sin(phi)
      x(3,3)=-r1*o1*cos(phi)
      x(4,1)=r2*sin(acos(o2))
      x(4,2)=zero
      x(4,3)=r2*o2

c Appel des subroutines de Kato
      if(Init) Then
         call init_potential
         Init=.false.
      Endif
      call calc_potential(x,v)
c Shift de la surface
      v = v + 213.4575280d0

      return
      end

c  Subroutines de Kato 


c     initialization of HFCO potential
c==================================================
      subroutine init_potential
c==================================================
      implicit real*8(a-h,o-z)
      call set_pot_index
      call set_pot_param
      return
      end


c     NOTE: (sub.)set_pot_index has to be called beforehand.
c==================================================
      subroutine set_pot_param
c==================================================

c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3)
c     ----------------------------------------
      common /mass/ zmc,zmo,zmh,zmf,zmco,zmhf,zmtot
      common /mass2/ zm1,zm2,zm3
      common /energy_pot/ v110,vts0,v220
c     ----------------------------------------
      common /v11coef_pot/ fcmeq(6,6),v11_coef(Ncoefmax)
      common /vtscoef_pot/ fcmts(6,6),vts_coef(Ncoefmax)
      common /index_pot/ ind(6,Ncoefmax),Ncoef
c     ----------------------------------------
      common /v11geom_pot/ xeq(4,3),bceq(6),zjceq(6)
      common /vtsgeom_pot/ xts(4,3),bcts(6),zjcts(6)
      common /s1geom_pot/ xs1(4,3),bcs1(6),zjcs1(6)
c     ----------------------------------------
      pi = dacos(-1.0d0)

c     write(2,1100)
 1100 format(//20x,40('#')/30x,'Kato Potential Surface '/20x,40('#'))

c
c     Mass. (atomic unit)
c
      zmc = 12.0d0 * amu
      zmo = 16.0d0 * amu
      zmh = 1.0d0  * amu
      zmf = 19.0d0 * amu

      zmco  = zmc + zmo
      zmhf  = zmh + zmf
      zmtot = zmco + zmhf

      zm1 = zmc  * zmo   / zmco 
      zm2 = zmco * zmhf  / zmtot
      zm3 = zmh  * zmf   / zmhf

c
c     Energy.
c
      v110 = -213.457457236198d0
      vts0 = -213.375843733324d0
      v220 = -213.4565070210d0

c
c     Geometry.
c
      ! S0 equilibrium.
      !--------------------
      bceq(1)=     1.1839304d0 / bohr
      bceq(2)=     1.0853751d0 / bohr
      bceq(3)=     1.3483149d0 / bohr
      bceq(4)=   109.0623746d0 * pi/180.0d0
      bceq(5)=   122.9633241d0 * pi/180.0d0
      bceq(6)=   0.0d0
      call btoc(bceq,xdum)
      call ctoj(xdum,zjceq)

      
      ! S0 transition state.
      !--------------------
      bcts(1) =   1.1402939d0 / bohr
      bcts(2) =   1.1259456d0 / bohr
      bcts(3) =   1.8428087d0 / bohr
      bcts(4) =  48.7644395d0 * pi/180.0d0
      bcts(5) = 122.2232022d0 * pi/180.0d0
      bcts(6) =    .0000000d0 * pi/180.0d0
      call btoc(bcts,xdum)
      call ctoj(xdum,zjcts)


      ! S1 equilibrium.
      !--------------------
      bcs1(1) =    1.3843896654d0 / bohr
      bcs1(2) =    1.0769979077d0 / bohr
      bcs1(3) =    1.3386149085d0 / bohr
      bcs1(4) =  114.4671755332d0 * pi/180.0d0
      bcs1(5) =  109.2908406289d0 * pi/180.0d0
      bcs1(6) =   50.5589711437d0 * pi/180.0d0
      call btoc(bcs1,xdum)
      call ctoj(xdum,zjcs1)

c     ------------------------------
c     potential coefficients 
c     ------------------------------
      call set_pot_coef

c     ------------------------------
c     parameters for v22
c     ------------------------------
      call set_v22_param

c     ------------------------------
c     parameters for normal mode
c     ------------------------------
      call set_normal_mode

      return
      end


c     Set polynomial index.
c     Direct Taylor-expansion indices are adopted.
c==================================================
      subroutine set_pot_index
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      integer i(10)
      common /index_pot/ ind(6,Ncoefmax),Ncoef

c      header = 'set polynomial index (direct ordering)'
c      call make_header


      Ncoef = 0

c     3rd order.
      !--------------------
      norder = 3
      do 300 i1 =  1,6
      do 300 i2 = i1,6
      do 300 i3 = i2,6
         i(1) = i1
         i(2) = i2
         i(3) = i3

         idum = 0
         do k1 = 1,norder
            if( i(k1).eq.6 ) idum = idum + 1
         enddo
         ires = mod(idum, 2)

         if(ires.eq.0) then
            Ncoef = Ncoef + 1
            do icord = 1,6
               idum = 0
               do k1 = 1,norder
                  if( i(k1).eq.icord ) idum = idum + 1
               enddo
               ind(icord, Ncoef) = idum
            enddo
         endif

 300  continue !---<


c     4th order.
      !--------------------
      norder = 4
      do 400 i1 =  1,6
      do 400 i2 = i1,6
      do 400 i3 = i2,6
      do 400 i4 = i3,6
         i(1) = i1
         i(2) = i2
         i(3) = i3
         i(4) = i4


         idum = 0
         do k1 = 1,norder
            if( i(k1).eq.6 ) idum = idum + 1
         enddo
         ires = mod(idum, 2)

         if(ires.eq.0) then
            Ncoef = Ncoef + 1
            do icord = 1,6
               idum = 0
               do k1 = 1,norder
                  if( i(k1).eq.icord ) idum = idum + 1
               enddo
               ind(icord, Ncoef) = idum
            enddo
         endif

 400  continue !---<


c     5th order.
      !--------------------
      norder = 5
      do 500 i1 =  1,6
      do 500 i2 = i1,6
      do 500 i3 = i2,6
      do 500 i4 = i3,6
      do 500 i5 = i4,6
         i(1) = i1
         i(2) = i2
         i(3) = i3
         i(4) = i4
         i(5) = i5


         idum = 0
         do k1 = 1,norder
            if( i(k1).eq.6 ) idum = idum + 1
         enddo
         ires = mod(idum, 2)

         if(ires.eq.0) then
            Ncoef = Ncoef + 1
            do icord = 1,6
               idum = 0
               do k1 = 1,norder
                  if( i(k1).eq.icord ) idum = idum + 1
               enddo
               ind(icord, Ncoef) = idum
            enddo
         endif

 500  continue !---<


c     6th order.
      !--------------------
      norder = 6
      do 600 i1 =  1,6
      do 600 i2 = i1,6
      do 600 i3 = i2,6
      do 600 i4 = i3,6
      do 600 i5 = i4,6
      do 600 i6 = i5,6
         i(1) = i1
         i(2) = i2
         i(3) = i3
         i(4) = i4
         i(5) = i5
         i(6) = i6

         idum = 0
         do k1 = 1,norder
            if( i(k1).eq.6 ) idum = idum + 1
         enddo
         ires = mod(idum, 2)

         if(ires.eq.0) then
            Ncoef = Ncoef + 1
            do icord = 1,6
               idum = 0
               do k1 = 1,norder
                  if( i(k1).eq.icord ) idum = idum + 1
               enddo
               ind(icord, Ncoef) = idum
            enddo
         endif

 600  continue !---<

c      write(6,*) '   * Ncoef = ',Ncoef
      return
      end


c     calculate total potential
c     * xdum : Cartesian(a.u.)
c     * vfit : a.u.
c==================================================
      subroutine calc_potential(xdum,vfit)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
      common /switch_pot/ sw_v11, sw_vts, sw_v22
c     ----------------------------------------

c
c     [ version 2 ]
c
      call calc_v11( xdum, v11 )
      call calc_vts( xdum, vts )
      call calc_v22( xdum, v22 )

      call correct_v11( xdum, v11 )
      call correct_pot( v11 )
      call correct_pot( vts )
      call correct_pot( v22 )

      call calc_switch2( xdum, sw_vts )
      call calc_vouter( v11, v22, vouter )
      vfit = (one-sw_vts)*vouter + sw_vts*vts

      return
      end

      
c     switching function for <iversion>=2
c==================================================
      subroutine calc_switch2( xdum, sw_vts )
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
      common /vtsgeom_pot/ xts(4,3),bcts(6),Djcts(6)
c     ----------------------------------------

c     [ switching between vts and others ]
c     --------------------
      Dnorm_lower =  0.7d0
      Dnorm_upper =  1.0d0
c     --------------------

      call calc_norm(xdum,Dnorm)
      sum = 0.5d0*( Dnorm_upper + Dnorm_lower )
      dif = 0.5d0*( Dnorm_upper - Dnorm_lower )
      arg = 3.0d0*( Dnorm-sum )/dif
      sw_vts = 0.5d0*( one-dtanh(arg) )

      return
      end


c     calculate norm from TS
c     * Jacobi coordinate
c==================================================
      subroutine calc_norm(xdum,Dnorm)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),Djc(6),q(6)
      real*8 base_plus(6),base_minus(6)
c     ----------------------------------------
      common /vtsgeom_pot/ xts(4,3),bcts(6),Djcts(6)
c     ----------------------------------------
      pi = dacos(-1.0d0)

c
c     [ version 2 ]
c

c     ------------------------------
      base_plus (1) = 0.15d0 / bohr
      base_minus(1) = 0.15d0 / bohr
c     --------------------
      base_plus (2) = 0.50d0 / bohr
      base_minus(2) = 0.40d0 / bohr
c     --------------------
      base_plus (3) = 0.50d0 / bohr
      base_minus(3) = 0.30d0 / bohr
c     --------------------
      base_plus (4) = 40.0d0 * pi/180.0d0
      base_minus(4) = 25.0d0 * pi/180.0d0
c     --------------------
      base_plus (5) = 20.0d0 * pi/180.0d0
      base_minus(5) = 25.0d0 * pi/180.0d0
c     --------------------
      base_plus (6) = 50.0d0 * pi/180.0d0
      base_minus(6) = 50.0d0 * pi/180.0d0
c     ------------------------------
      npow = 6
c     ------------------------------

      call ctoj(xdum,Djc)
      do i1 = 1,6
         q(i1) = Djc(i1) - Djcts(i1)
      enddo
      dum2 = q(2)
      dum3 = q(3)
      q(2) = (dum2+dum3)/dsqrt(two)
      q(3) = (dum2-dum3)/dsqrt(two)

      Dnorm = zero
      do i1 = 1,6
         if( q(i1).gt.zero ) then
            dum = q(i1)/base_plus(i1)
         else
            dum = q(i1)/base_minus(i1)
         endif
         Dnorm = Dnorm + dum**npow
      enddo

      Dnorm = Dnorm**(one/dble(npow))
      return
      end


c     calculate distance from EQ (equilibrium geometry)
c==================================================
      subroutine calc_EQnorm(xdum,Dnorm)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),bc(6)
      real*8 base_plus(6),base_minus(6)
c     ----------------------------------------
      common /v11geom_pot/ xeq(4,3),bceq(6),Djceq(6)
c     ----------------------------------------
      pi = dacos(-1.0d0)

c     ------------------------------
      base_plus (1) = 0.37d0 / bohr
      base_minus(1) = 0.20d0 / bohr
c     --------------------
      base_plus (2) = 0.78d0 / bohr
      base_minus(2) = 0.32d0 / bohr
c     --------------------
      base_plus (3) = 0.8d0 / bohr
      base_minus(3) = 0.3d0 / bohr
c     --------------------
      base_plus (4) = 66.0d0 * pi/180.0d0
      base_minus(4) = 66.0d0 * pi/180.0d0
c     --------------------
      base_plus (5) = 47.0d0 * pi/180.0d0
      base_minus(5) = 40.0d0 * pi/180.0d0
c     --------------------
      base_plus (6) = 90.0d0 * pi/180.0d0
      base_minus(6) = 90.0d0 * pi/180.0d0
c     ------------------------------
      npow = 6
c     ------------------------------

      call ctob(xdum,bc)
      do i1 = 1,6
         bc(i1) = bc(i1) - bceq(i1)
      enddo

      Dnorm = zero
      do i1 = 1,6
         if( bc(i1).gt.zero ) then
            dum = bc(i1)/base_plus(i1)
         else
            dum = bc(i1)/base_minus(i1)
         endif
         Dnorm = Dnorm + dum**npow
      enddo
      Dnorm = Dnorm**(one/dble(npow))

      return
      end


c     give ceiling correction to v11 using EQ-norm
c==================================================
      subroutine correct_v11( xdum, v11 )
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
c     ----------------------------------------

c     --------------------
      Dnorm_lower = 1.0d0
      Dnorm_upper = 1.15d0
      E_ceiling   = vts0 + 30.0d0/kcal
c     --------------------

      call calc_EQnorm(xdum,Dnorm)

      sum = 0.5d0*(Dnorm_upper + Dnorm_lower)
      dif = 0.5d0*(Dnorm_upper - Dnorm_lower)
      arg = 3.0d0*(Dnorm-sum)/dif
      switch = 0.5d0*(one-dtanh(arg))

      v11 = switch*v11 + (one-switch)*E_ceiling
      return
      end


c     give floor and ceiling correction to energy
c==================================================
      subroutine correct_pot( vtmp )
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
c     ----------------------------------------

c     --------------------
      E_lower1   = v110 - 10.0d0/kcal
      E_lower2   = v110 -  7.0d0/kcal
      E_upper1   = vts0 + 27.0d0/kcal
      E_upper2   = vts0 + 30.0d0/kcal
c     --------------------

      if( vtmp.lt.E_lower2 ) then
         sum = 0.5d0*( E_lower2 + E_lower1 )
         dif = 0.5d0*( E_lower2 - E_lower1 )
         arg = 3.0d0*(vtmp-sum)/dif
         switch = 0.5d0*( one+dtanh(arg) )
         vtmp = switch*vtmp + (one-switch)*E_lower1
         return
      endif
      
      if( vtmp.gt.E_upper1 ) then
         sum = 0.5d0*( E_upper2 + E_upper1 )
         dif = 0.5d0*( E_upper2 - E_upper1 )
         arg = 3.0d0*(vtmp-sum)/dif
         switch = 0.5d0*( one-dtanh(arg) )
         vtmp = switch*vtmp + (one-switch)*E_upper2
         return
      endif

      return
      end


c     calculate outer potential by linking v11 and
c     v22 with the use of EVB-like method
c==================================================
      subroutine calc_vouter( v11, v22, vouter )
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header

c     --------------------
      v12 = 2.0d0 / kcal
c     --------------------
      sum = 0.5d0*( v11+v22 )
      dif = 0.5d0*( v11-v22 )
      vouter = sum - dsqrt( dif**2 + v12**2 )

      return
      end


c     local potential for interaction region
c==================================================
      subroutine calc_v11(xdum,v11)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),bc(6),u(6,0:Nordermax)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
      common /v11geom_pot/ xeq(4,3),bceq(6),zjceq(6)
      common /v11coef_pot/ fcmeq(6,6),v11_coef(Ncoefmax)
      common /index_pot/ ind(6,Ncoefmax),Ncoef
c     ----------------------------------------
      
      call ctob(xdum,bc)

      do i1 = 1,6
         u(i1,0) = one
         u(i1,1) = bc(i1) - bceq(i1)
         do iorder = 2,Nordermax
            u(i1,iorder) = u(i1,iorder-1) * u(i1,1)
         enddo
      enddo

c     quadratic term
      dum = zero
      do i1 = 1,6
         do i2 = 1,6
            dum = dum + fcmeq(i1,i2) * (bc(i1)-bceq(i1))
     $           * (bc(i2)-bceq(i2))
         enddo
      enddo
      v11   = v110 + 0.5d0*dum

c     higher order term
      do icoef = 1,Ncoef
         dum = one
         do i1 = 1,6
            dum = dum * u(i1,ind(i1,icoef))
         enddo
         v11 = v11 + v11_coef(icoef)*dum
      enddo

      return
      end


c     local potential for transition-state region
c==================================================
      subroutine calc_vts(xdum,vts)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),bc(6),u(6,0:Nordermax)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
      common /vtsgeom_pot/ xts(4,3),bcts(6),zjcts(6)
      common /vtscoef_pot/ fcmts(6,6),vts_coef(Ncoefmax)
      common /index_pot/ ind(6,Ncoefmax),Ncoef
c     ----------------------------------------
      
      call ctob(xdum,bc)

      do i1 = 1,6
         u(i1,0) = one
         u(i1,1) = bc(i1) - bcts(i1)
         do iorder = 2,Nordermax
            u(i1,iorder) = u(i1,iorder-1) * u(i1,1)
         enddo
      enddo

c     quadratic term
      dum = zero
      do i1 = 1,6
         do i2 = 1,6
            dum = dum + fcmts(i1,i2) * (bc(i1)-bcts(i1))
     $           * (bc(i2)-bcts(i2))
         enddo
      enddo
      vts   = vts0 + 0.5d0*dum

c     higher order term
      do icoef = 1,Ncoef
         dum = one
         do i1 = 1,6
            dum = dum * u(i1,ind(i1,icoef))
         enddo
         vts = vts + vts_coef(icoef)*dum
      enddo

      return
      end


c==================================================
      subroutine set_v22_param
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
c     ----------------------------------------
      common /diatomics_pot/ rco0,rhf0,De_co,De_hf,a_co,a_hf
      common /diatomics2_pot/ b_co(10),b_hf(10)
      common /diatomics3_pot/ wfreq_co,wfreq_hf
c     ----------------------------------------
      common /pvrep_potv22/ pvrep(4,10),n_vrep
      common /order_potv22/ n_Amp
      common /pAmp_potv22/ pAmp_co(4,10),pAmp_hf(4,10)
c     ----------------------------------------
      common /CTsite_v22/ rct_H0,rct_F0,rct_HF0,dct_HF0
      common /CTsite2_v22/ damp_tot,damp_site(4),damp_co,damp_hf
      common /CTsite3_v22/ r_ch0,thH0,r_cf0,thF0
      common /CTsite4_v22/ dy_H0
c     ----------------------------------------
      pi = dacos(-1.0d0)

c      header = 'set v22 parameters'
c      call make_header

c
c     CT-site parameters.
c

      ! Initial guess.
      !--------------------
c      r_ch0 = 1.126d0 / bohr
c      thH0 = 10.0d0 * pi/180.0d0
c      r_cf0 = 1.843d0 / bohr
c      thF0 = 60.0d0 * pi/180.0d0
c      rct_H0 =  .7101735853115629d0 / bohr
c      rct_F0 =  3.08587076847724d0 / bohr
c      rct_HF0 = 1.1d0 / bohr
c      dct_HF0 = 0.2d0 / bohr
c      damp_tot = 0.5d0
c      damp_site(1) = 0.45d0/3.0d0
c      damp_site(2) = 0.45d0/3.0d0
c      damp_hf = 0.45d0/3.0d0
c      dy_H0 = 1.0d0


      ! Gaussian-type CT-site.
      ! CT-site(H,F) is used, but F-site
      ! has almost no effect (it is unity
      ! almost everywhere).
      !--------------------
      rct_H0 =      2.05410733d0
      rct_F0 =     13.73818140d0
      damp_tot =       .44581541d0
      r_ch0 =      1.74972194d0
      thH0 =       .80504821d0
      r_cf0 =      2.29269875d0
      thF0 =      2.66005678d0
      rct_HF0 =      1.84037102d0
      dct_HF0 =       .20261232d0


      ! * (970315) znorm>1.3, v_upper=vts+20.0kcal
      ! * CT-site(F) is not used in evaluating
      !   damping strength.
      ! * <dy_H0> is introduced.
      !--------------------
      rct_H0 =      2.06091142d0
      damp_tot =       .42966525d0
      r_ch0 =      1.64598242d0
      thH0 =       .74780062d0
      rct_HF0 =      1.90455802d0
      dct_HF0 =       .26369038d0
      dy_H0 =       .66235880d0


c
c     Diatomics potential.
c

      rco0 = 1.1283d0 / bohr
      rhf0 = 0.9168d0 / bohr

      ! Parameters for Morse potential.
      De_co =  .4124910598767912d0  
      De_hf =  .22510798955892d0

      a_co =  1.22119123854612d0
      a_hf =  1.17794137850757d0

      ! Parameters for extended Rydberg potential.
      b_co( 1) =  2.31354004629352d0
      b_co( 2) =  1.18472254353607d0
      b_co( 3) =  1.22056768077309d0
      b_co( 4) =  .8301270321159142d0
      b_co( 5) =  -1.93044383367946d0

      b_hf( 1) =  1.66555396537341d0
      b_hf( 2) =  9.718156860755167d-03
      b_hf( 3) =  .2721590699435786d0
      b_hf( 4) =  -.353308946538763d0
      b_hf( 5) =  .100594754195106d0

      wfreq_co = 2169.81358d0 * hartree
      wfreq_hf = 4138.32d0    * hartree

c
c     Expansion order.
c

      ! site-site repulsion functions
      n_vrep = 2

      ! angle modulation functions
      n_Amp = 3

c
c     v(rep) between each site.
c

      ! Old parameters. (970313)
      !--------------------

      ! CH
      pvrep(1,1) =         6.1059484771d0 * bohr
      pvrep(1,2) =        12.4098110850d0
      pvrep(1,3) =        17.6143595209d0 * bohr
      pvrep(1,4) =         3.2272243904d0 * bohr*bohr

      ! CF
      pvrep(2,1) =         3.7190927745d0 * bohr
      pvrep(2,2) =         9.9321952199d0
      pvrep(2,3) =        28.7268619975d0 * bohr
      pvrep(2,4) =          .8399432042d0 * bohr*bohr

      ! OH
      pvrep(3,1) =         6.0550232979d0 * bohr
      pvrep(3,2) =         4.3316477225d0
      pvrep(3,3) =        18.9570554616d0 * bohr
      pvrep(3,4) =         9.7772319986d0 * bohr*bohr

      ! OF
      pvrep(4,1) =         3.5636370755d0 * bohr
      pvrep(4,2) =        27.6329685104d0
      pvrep(4,3) =        27.1318304162d0 * bohr
      pvrep(4,4) =       -13.6946562928d0 * bohr*bohr


      ! (970314)
      !--------------------
      pvrep(1,1)=      2.82987021d0
      pvrep(1,2)=     13.74465005d0
      pvrep(2,1)=      1.59543280d0
      pvrep(2,2)=     18.14836000d0
      pvrep(3,1)=      2.63098230d0
      pvrep(3,2)=      7.66606206d0
      pvrep(4,1)=      1.73298550d0
      pvrep(4,2)=     17.62397863d0


c
c     Angle modulation.
c
      ! (970314)
      !--------------------

      pAmp_co(1,1)=       .22369783d0
      pAmp_co(1,2)=       .05310495d0
      pAmp_co(1,3)=      -.57093464d0
      pAmp_co(1,4)=       .15576494d0

      pAmp_co(2,1)=      -.63557990d0
      pAmp_co(2,2)=       .36516197d0
      pAmp_co(2,3)=       .26159898d0
      pAmp_co(2,4)=       .01658816d0

      pAmp_co(3,1)=       .25322070d0
      pAmp_co(3,2)=      -.26802399d0
      pAmp_co(3,3)=      -.16409634d0
      pAmp_co(3,4)=       .04933701d0

      pAmp_co(4,1)=      -.04168087d0
      pAmp_co(4,2)=      -.00323396d0
      pAmp_co(4,3)=      -.04307188d0
      pAmp_co(4,4)=       .11446998d0

      pAmp_hf(1,1)=      -.79397896d0
      pAmp_hf(1,2)=      1.12801931d0
      pAmp_hf(1,3)=      -.06432629d0
      pAmp_hf(1,4)=      -.21583338d0

      pAmp_hf(2,1)=       .12077888d0
      pAmp_hf(2,2)=      -.02030453d0
      pAmp_hf(2,3)=      -.03862673d0
      pAmp_hf(2,4)=      -.00310208d0

      pAmp_hf(3,1)=      -.29021047d0
      pAmp_hf(3,2)=       .69769590d0
      pAmp_hf(3,3)=      -.05091486d0
      pAmp_hf(3,4)=      -.17679996d0

      pAmp_hf(4,1)=       .26486282d0
      pAmp_hf(4,2)=      -.13288284d0
      pAmp_hf(4,3)=      -.11266260d0
      pAmp_hf(4,4)=       .00849603d0

      return
      end


c==================================================
      subroutine calc_v22(xdum,v22)
c      subroutine calc_v22_CTsite(xdum,v22)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),sd(6)
c     ----------------------------------------
      common /energy_pot/ v110,vts0,v220
c     ----------------------------------------
      common /CTsite_v22/ rct_H0,rct_F0, rct_HF0,dct_HF0
      common /CTsite2_v22/ damp_tot,damp_site(4),damp_co,damp_hf
c     ----------------------------------------
      common /vparts_potv22/ vrep_site(4),vco,vhf
c     ----------------------------------------

      call calc_v22_nct(xdum,v22nct)
      call calc_CT_length(xdum,rct_H,rct_F)
      call ctos(xdum,sd)

      dum1 = (rct_H / rct_H0)**2
      dum2 = (rct_F / rct_F0)**2
      dum3 = (sd(6)-rct_HF0) / dct_HF0

c      ratio = dexp(-dum1)*dexp(-dum2)
c     $     * 0.5d0*(one + dtanh(3.0d0*dum3))

      ratio = dexp(-dum1)
     $     * 0.5d0*(one + dtanh(3.0d0*dum3))


      v_CHF = (one - damp_tot*ratio) * (vrep_site(1) + 
     $     vrep_site(2) + vhf)
      v_OHF = vrep_site(3) + vrep_site(4) + vco

      v22 = v220 + v_OHF + v_CHF

      return
      end


c==================================================
      subroutine calc_CT_length(xdum,rct_H,rct_F)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      parameter (error_tol = 1.0d-6, eps=1.0d-4)
      Integer stat
      real*8 xdum(4,3)
      real*8 v_ch(3),v_cf(3),v_oc(3),zn_H(3),zn_F(3)
      real*8 x_H(3),x_H0(3),x_F(3),x_F0(3)
c     real*8 xdum_(4,3),sd(6),sd_(6)
c     ----------------------------------------
      common /CTsite3_v22/ r_ch0,thH0, r_cf0,thF0
      common /CTsite4_v22/ dy_H0
c     ----------------------------------------
      pi = dacos(-1.0d0)

      r_ch = zero
      r_cf = zero
      r_oc = zero
      do i1 = 1,3
         v_ch(i1) = xdum(3,i1) - xdum(1,i1)
         v_cf(i1) = xdum(4,i1) - xdum(1,i1)
         v_oc(i1) = xdum(1,i1) - xdum(2,i1)
         r_ch = r_ch + v_ch(i1)*v_ch(i1)
         r_cf = r_cf + v_cf(i1)*v_cf(i1)
         r_oc = r_oc + v_oc(i1)*v_oc(i1)
      enddo
      r_ch = dsqrt(r_ch)
      r_cf = dsqrt(r_cf)
      r_oc = dsqrt(r_oc)
      do i1 = 1,3
         v_ch(i1) = v_ch(i1) / r_ch
         v_cf(i1) = v_cf(i1) / r_cf
         v_oc(i1) = v_oc(i1) / r_oc
      enddo

      cos_thF = zero
      cos_thH = zero
      do i1 = 1,3
         cos_thF = cos_thF + v_oc(i1)*v_cf(i1)
         cos_thH = cos_thH + v_oc(i1)*v_ch(i1)
      enddo
      stat=979
      thF = dacos2( cos_thF ,stat)
      if(stat.eq.0) goto 999
      stat=982
      thH = dacos2( cos_thH ,stat)
      if(stat.eq.0) goto 999

      zn_F(1) =  v_cf(2)*v_oc(3) - v_cf(3)*v_oc(2)
      zn_F(2) =  v_cf(3)*v_oc(1) - v_cf(1)*v_oc(3)
      zn_F(3) =  v_cf(1)*v_oc(2) - v_cf(2)*v_oc(1)

      zn_H(1) =  v_ch(2)*v_oc(3) - v_ch(3)*v_oc(2)
      zn_H(2) =  v_ch(3)*v_oc(1) - v_ch(1)*v_oc(3)
      zn_H(3) =  v_ch(1)*v_oc(2) - v_ch(2)*v_oc(1)

      rn_H = zero
      rn_F = zero
      do i1 = 1,3
         rn_H = rn_H + zn_H(i1)**2
         rn_F = rn_F + zn_F(i1)**2
      enddo
      rn_H = dsqrt(rn_H)
      rn_F = dsqrt(rn_F)

      if((rn_H.lt.eps).or.(rn_F.lt.eps)) then
         phai = zero
      else
         do i1 = 1,3
            zn_H(i1) = zn_H(i1) / rn_H
            zn_F(i1) = zn_F(i1) / rn_F
         enddo
         cos_phai = zero
         do i1 = 1,3
            cos_phai = cos_phai + zn_F(i1)*zn_H(i1)
         enddo
         stat=1014
         phai = dacos2( cos_phai ,stat)
         if(stat.eq.0) goto 999
      endif

c     Location of H and F.
      x_F(1) = r_cf * dsin(thF)
      x_F(2) = zero
      x_F(3) = r_cf * dcos(thF)

      x_H(1) = r_ch * dsin(thH) * dcos(phai)
      x_H(2) = r_ch * dsin(thH) * dsin(phai)
      x_H(3) = r_ch * dcos(thH)


c     Conversion check.
c      goto 1000 ! SKIP--->>
c      xdum_(1,1) = zero
c      xdum_(1,2) = zero
c      xdum_(1,3) = zero
c
c      xdum_(2,1) = zero
c      xdum_(2,2) = zero
c      xdum_(2,3) = -r_oc
c
c      do i1 = 1,3
c         xdum_(3,i1) = x_H(i1)
c         xdum_(4,i1) = x_F(i1)
c      enddo
c
c      call ctos(xdum ,sd )
c      call ctos(xdum_,sd_)
c      error = zero
c      do i1 = 1,6
c         error = error + dabs( sd(i1)-sd_(i1) )
c      enddo
c      if(error.gt.error_tol) then
c         write(6,*) 'Conversion failed.'
c         write(6,*) 'error = ',error
c         do i1 = 1,6
c            write(6,*) 'sd(',i1,') = ',sd(i1)
c            write(6,*) 'sd_(',i1,') = ',sd_(i1)
c         enddo
c         pause
c      endif
c 1000 continue !---<<


c     Location of CT site.
      x_F0(1) = r_cf0 * dsin(thF0)
      x_F0(2) = zero
      x_F0(3) = r_cf0 * dcos(thF0)

      x_H0(1) = r_ch0 * dsin(thH0)
      x_H0(2) = zero
      x_H0(3) = r_ch0 * dcos(thH0)

      rct_F = zero
      do i1 = 1,3
         rct_F = rct_F + (x_F(i1)-x_F0(i1))**2
      enddo
      rct_F = dsqrt(rct_F)

      rct_H = (x_H(1)-x_H0(1))**2 + (x_H(3)-x_H0(3))**2
     $     + ( (x_H(2)-x_H0(2))/dy_H0 )**2
      rct_H = dsqrt(rct_H)

      return
  999 Write (6,*) 'r_ch=',r_ch
      Write (6,*) 'r_cf=',r_cf
      Write (6,*) 'r_oc=',r_oc
      stop
      end


c     V22 : No Charge Transfer effect included
c     * NOTE: v220 not included.
c==================================================
      subroutine calc_v22_nct(xdum,v22nct)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),sd(6)
      real*8 cos_co(4,0:10),cos_hf(4,0:10)
c     ----------------------------------------
      common /pvrep_potv22/ pvrep(4,10),n_vrep
      common /order_potv22/ n_Amp
      common /pAmp_potv22/ pAmp_co(4,10),pAmp_hf(4,10)
c     ----------------------------------------
      common /vparts_potv22/ vrep_site(4),vco,vhf
c     ----------------------------------------

c     Site-site distances and angles.
      !--------------------
      call ctos(xdum,sd)
      call get_site_angle(xdum, cos_co, cos_hf, n_Amp)


c     Site-site exhange-repulsion energy.
      !--------------------
      do ibond = 1,4
         r = sd( 1 + ibond )
         vrep_site(ibond) = pvrep(ibond,2)
     $        * dexp(-pvrep(ibond,1)*r)
      enddo


      vrep_tot = zero
      do ibond = 1,4

         dum_co = pAmp_co(ibond,1)
         dum_hf = pAmp_hf(ibond,1)

         do n = 1,n_Amp
            dum_co = dum_co 
     $           + pAmp_co(ibond,n+1)*cos_co(ibond,n)
            dum_hf = dum_hf
     $           + pAmp_hf(ibond,n+1)*cos_hf(ibond,n)
         enddo
         Amp_co = dexp(dum_co)
         Amp_hf = dexp(dum_hf)

         vrep_site(ibond) = Amp_co*Amp_hf * vrep_site(ibond)
         vrep_tot = vrep_tot + vrep_site(ibond)

      enddo


c     diatomics potential.
      !--------------------
      call calc_diatomics_pot(sd(1),sd(6),vco,vhf)

c     total potential
      !--------------------
      v22nct = vrep_tot + vco + vhf
      return
      end


c==================================================
      subroutine get_site_angle(xdum,cos_co,cos_hf,norder)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 xdum(4,3),cos_co(4,0:10),cos_hf(4,0:10)
      external bond_angle_cos

      do ibond = 1,4
         cos_co(ibond,0) = one
         cos_hf(ibond,0) = one
      enddo

      cos_co(1,1) = bond_angle_cos(xdum,3,1,2)
      cos_co(2,1) = bond_angle_cos(xdum,4,1,2)
      cos_co(3,1) = bond_angle_cos(xdum,3,2,1)
      cos_co(4,1) = bond_angle_cos(xdum,4,2,1)

      cos_hf(1,1) = bond_angle_cos(xdum,1,3,4)
      cos_hf(2,1) = bond_angle_cos(xdum,1,4,3)
      cos_hf(3,1) = bond_angle_cos(xdum,2,3,4)
      cos_hf(4,1) = bond_angle_cos(xdum,2,4,3)

      do ibond = 1,4
         do iorder = 2,norder
            cos_co(ibond,iorder) = 
     $           two*cos_co(ibond,1)*cos_co(ibond,iorder-1)
     $           - cos_co(ibond,iorder-2)
            cos_hf(ibond,iorder) = 
     $           two*cos_hf(ibond,1)*cos_hf(ibond,iorder-1)
     $           - cos_hf(ibond,iorder-2)
         enddo
      enddo

      return
      end


c     _diatomics potential 
c     _RHF/MP2 (target energy < 70.0 kcal/mol)
c==================================================
      subroutine calc_diatomics_pot(rco,rhf,vco,vhf)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
c     ----------------------------------------
      common /diatomics_pot/ rco0,rhf0,De_co,De_hf,a_co,a_hf
      common /diatomics2_pot/ b_co(10),b_hf(10)
c     ----------------------------------------

      drco = rco - rco0
      drhf = rhf - rhf0

c     extended Rydberg potential
      sum_co = one + b_co(1)*drco + b_co(2)*(drco**2)
     $     + b_co(3)*(drco**3) + b_co(4)*(drco**4) 
     $     + b_co(5)*(drco**5)

      sum_hf = one + b_hf(1)*drhf + b_hf(2)*(drhf**2)
     $     + b_hf(3)*(drhf**3) + b_hf(4)*(drhf**4) 
     $     + b_hf(5)*(drhf**5)

      vco = De_co - De_co*sum_co *dexp(-b_co(1)*drco)
      vhf = De_hf - De_hf*sum_hf *dexp(-b_hf(1)*drhf)

      return
      end


c     _diatomics potential 
c     _Morse (exp.)
c     _NOTE: input -> a.u. (trajectory calculation
c     still uses angstrom as unit.)
c==================================================
      subroutine calc_diatomics_pot2(rco,rhf,vco,vhf)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      common /diatomics_pot/ rco0,rhf0,De_co,De_hf,a_co,a_hf

      drco = rco - rco0
      drhf = rhf - rhf0
      dum = one - dexp(-a_co*drco)
      vco = De_co * dum*dum
      dum = one - dexp(-a_hf*drhf)
      vhf = De_hf * dum*dum 

      return
      end


c     set coefficients in v11 and vts
c==================================================
      subroutine set_pot_coef
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
c     ----------------------------------------
      common /v11coef_pot/ fcmeq(6,6),v11_coef(Ncoefmax)
      common /vtscoef_pot/ fcmts(6,6),vts_coef(Ncoefmax)
c     ----------------------------------------

c
c     FCM at S0:equilibrium
c
      fcmeq(1,1) =          .9397423948d0
      fcmeq(1,2) =          .0162793929d0
      fcmeq(1,3) =          .0840682729d0
      fcmeq(1,4) =         -.0548419065d0
      fcmeq(1,5) =          .0076820560d0
      fcmeq(1,6) =          .0000000348d0

      fcmeq(2,1) =          .0162793929d0
      fcmeq(2,2) =          .3485437426d0
      fcmeq(2,3) =          .0143461237d0
      fcmeq(2,4) =          .0034382189d0
      fcmeq(2,5) =         -.0174870099d0
      fcmeq(2,6) =         -.0000000020d0

      fcmeq(3,1) =          .0840682729d0
      fcmeq(3,2) =          .0143461237d0
      fcmeq(3,3) =          .3516980522d0
      fcmeq(3,4) =          .0697067228d0
      fcmeq(3,5) =          .0673367283d0
      fcmeq(3,6) =         -.0000000126d0

      fcmeq(4,1) =         -.0548419065d0
      fcmeq(4,2) =          .0034382189d0
      fcmeq(4,3) =          .0697067228d0
      fcmeq(4,4) =          .2940189054d0
      fcmeq(4,5) =          .1328827854d0
      fcmeq(4,6) =         -.0000000004d0

      fcmeq(5,1) =          .0076820560d0
      fcmeq(5,2) =         -.0174870099d0
      fcmeq(5,3) =          .0673367283d0
      fcmeq(5,4) =          .1328827854d0
      fcmeq(5,5) =          .4621674610d0
      fcmeq(5,6) =         -.0000000083d0

      fcmeq(6,1) =          .0000000348d0
      fcmeq(6,2) =         -.0000000020d0
      fcmeq(6,3) =         -.0000000126d0
      fcmeq(6,4) =         -.0000000004d0
      fcmeq(6,5) =         -.0000000083d0
      fcmeq(6,6) =          .0882721498d0
c
c     FCM at S0:transition-state
c
      fcmts(1,1) =         1.1452843955d0
      fcmts(1,2) =         -.0127422527d0
      fcmts(1,3) =          .0385564229d0
      fcmts(1,4) =         -.0181903455d0
      fcmts(1,5) =          .0107883319d0
      fcmts(1,6) =          .0000009942d0

      fcmts(2,1) =         -.0127422527d0
      fcmts(2,2) =          .2354460736d0
      fcmts(2,3) =          .0076748773d0
      fcmts(2,4) =          .1886463228d0
      fcmts(2,5) =          .0363382165d0
      fcmts(2,6) =         -.0000000705d0

      fcmts(3,1) =          .0385564229d0
      fcmts(3,2) =          .0076748773d0
      fcmts(3,3) =          .1260692975d0
      fcmts(3,4) =          .1494575355d0
      fcmts(3,5) =         -.0082880551d0
      fcmts(3,6) =          .0000000643d0

      fcmts(4,1) =         -.0181903455d0
      fcmts(4,2) =          .1886463228d0
      fcmts(4,3) =          .1494575355d0
      fcmts(4,4) =         -.1904797241d0
      fcmts(4,5) =         -.0581915197d0
      fcmts(4,6) =          .0000002968d0

      fcmts(5,1) =          .0107883319d0
      fcmts(5,2) =          .0363382165d0
      fcmts(5,3) =         -.0082880551d0
      fcmts(5,4) =         -.0581915197d0
      fcmts(5,5) =          .1504878796d0
      fcmts(5,6) =         -.0000002089d0

      fcmts(6,1) =          .0000009942d0
      fcmts(6,2) =         -.0000000705d0
      fcmts(6,3) =          .0000000643d0
      fcmts(6,4) =          .0000002968d0
      fcmts(6,5) =         -.0000002089d0
      fcmts(6,6) =          .0519242677d0

c
c     higher-order coefficients in v11
c

      v11_coef(  1) =          -.559882661205916D+00
      v11_coef(  2) =          -.232906796305885D-01
      v11_coef(  3) =          -.726521776088361D-01
      v11_coef(  4) =           .313730987579553D-01
      v11_coef(  5) =           .686715011911285D-02
      v11_coef(  6) =           .114076992862326D-01
      v11_coef(  7) =          -.108824169397772D-01
      v11_coef(  8) =           .298069559251759D-01
      v11_coef(  9) =           .434367747907133D-01
      v11_coef( 10) =          -.312111362424967D-01
      v11_coef( 11) =           .720472605841483D-02
      v11_coef( 12) =          -.100870572959878D+00
      v11_coef( 13) =          -.481260379113531D-01
      v11_coef( 14) =          -.106670258996291D+00
      v11_coef( 15) =          -.169800514127854D+00
      v11_coef( 16) =          -.391127300200841D-01
      v11_coef( 17) =          -.180292705907534D+00
      v11_coef( 18) =           .934735407586513D-02
      v11_coef( 19) =           .853973553150172D-02
      v11_coef( 20) =           .330599655488930D-02
      v11_coef( 21) =          -.197989453524872D-01
      v11_coef( 22) =          -.482693622221962D-01
      v11_coef( 23) =           .679196788248472D-02
      v11_coef( 24) =          -.276838007638719D-01
      v11_coef( 25) =          -.171408301879854D-01
      v11_coef( 26) =          -.140175524069521D-01
      v11_coef( 27) =          -.103201057267949D-01
      v11_coef( 28) =          -.212666239926390D+00
      v11_coef( 29) =          -.209306076274773D-01
      v11_coef( 30) =          -.403387665689342D-01
      v11_coef( 31) =          -.739297836478008D-01
      v11_coef( 32) =           .126952596049335D-01
      v11_coef( 33) =          -.130447370339530D+00
      v11_coef( 34) =          -.478388119946187D-02
      v11_coef( 35) =          -.226778908811394D-01
      v11_coef( 36) =           .815842740591878D-01
      v11_coef( 37) =           .901248879025723D-01
      v11_coef( 38) =           .160621262072408D-01
      v11_coef( 39) =          -.506466268263633D-01
      v11_coef( 40) =           .249253681523608D-01
      v11_coef( 41) =           .409908370592899D+00
      v11_coef( 42) =           .128056784769732D-01
      v11_coef( 43) =           .256565898269457D-01
      v11_coef( 44) =           .227426661812664D-01
      v11_coef( 45) =          -.544504196295375D-02
      v11_coef( 46) =           .157843463740714D-01
      v11_coef( 47) =          -.217008445074869D-02
      v11_coef( 48) =          -.263477443120749D-01
      v11_coef( 49) =          -.550322941001859D-01
      v11_coef( 50) =           .364548174280287D-01
      v11_coef( 51) =          -.500304513894019D-01
      v11_coef( 52) =           .680087850725646D-01
      v11_coef( 53) =           .609562565004424D-02
      v11_coef( 54) =           .908300809847826D-02
      v11_coef( 55) =           .296086287465174D-01
      v11_coef( 56) =           .101061005655809D-01
      v11_coef( 57) =          -.528724291715540D-02
      v11_coef( 58) =           .686236105198694D-02
      v11_coef( 59) =          -.192316663422767D-02
      v11_coef( 60) =           .159729201549044D-01
      v11_coef( 61) =          -.818061246600772D-02
      v11_coef( 62) =           .757082536358028D-02
      v11_coef( 63) =          -.388003879411136D-01
      v11_coef( 64) =          -.504642510081322D-02
      v11_coef( 65) =           .270355044420737D-01
      v11_coef( 66) =           .330554844026565D-01
      v11_coef( 67) =           .146491484614300D-01
      v11_coef( 68) =          -.549107029170821D-02
      v11_coef( 69) =           .561211126907189D-02
      v11_coef( 70) =           .816949455240398D-01
      v11_coef( 71) =           .265625690032769D-01
      v11_coef( 72) =           .198235544541230D-01
      v11_coef( 73) =           .158063826707978D+00
      v11_coef( 74) =           .518157984651908D-02
      v11_coef( 75) =          -.194786540079733D-01
      v11_coef( 76) =          -.918653097809285D-01
      v11_coef( 77) =          -.142513622206340D+00
      v11_coef( 78) =          -.198276412339180D-01
      v11_coef( 79) =           .449831038459943D-02
      v11_coef( 80) =          -.137451251111253D-01
      v11_coef( 81) =           .120212288659951D+00
      v11_coef( 82) =          -.240333496778197D-02
      v11_coef( 83) =           .617282938057679D-02
      v11_coef( 84) =           .111815650309869D-01
      v11_coef( 85) =           .403253733675198D-02
      v11_coef( 86) =          -.341557760650541D-02
      v11_coef( 87) =           .997505023667787D-02
      v11_coef( 88) =          -.222875511106359D-02
      v11_coef( 89) =           .458373368709093D-02
      v11_coef( 90) =          -.140450500433884D-01
      v11_coef( 91) =           .212409775964780D-02
      v11_coef( 92) =          -.131577872850431D-02
      v11_coef( 93) =           .272198012803195D-01
      v11_coef( 94) =          -.272294337737091D-01
      v11_coef( 95) =           .366354618775654D-01
      v11_coef( 96) =          -.117844571927281D-01
      v11_coef( 97) =           .176420868524685D-01
      v11_coef( 98) =          -.145586491575891D-02
      v11_coef( 99) =          -.137806966081022D-01
      v11_coef(100) =          -.102766964225267D-01
      v11_coef(101) =           .148504272988635D-02
      v11_coef(102) =           .317577545633464D-02
      v11_coef(103) =          -.187110000955316D-01
      v11_coef(104) =          -.477804657332538D-02
      v11_coef(105) =           .182379995999720D+00
      v11_coef(106) =          -.420280009561954D-03
      v11_coef(107) =           .190021742660192D-01
      v11_coef(108) =           .147981801898131D-01
      v11_coef(109) =          -.171372769431311D-01
      v11_coef(110) =           .151978388031368D-01
      v11_coef(111) =          -.123265202337139D-04
      v11_coef(112) =           .404090594398497D-01
      v11_coef(113) =          -.462640595738796D-01
      v11_coef(114) =          -.385580949158572D-01
      v11_coef(115) =          -.417649420202371D-02
      v11_coef(116) =           .870475125649155D-01
      v11_coef(117) =          -.108752298274461D-01
      v11_coef(118) =          -.128223636490120D-02
      v11_coef(119) =          -.457998518405063D-01
      v11_coef(120) =           .336936921979798D-01
      v11_coef(121) =          -.961683811080050D-02
      v11_coef(122) =           .551504811495002D-01
      v11_coef(123) =           .215859461497132D-01
      v11_coef(124) =           .920352550846166D-01
      v11_coef(125) =          -.957317617906504D-02
      v11_coef(126) =           .289035416819562D-02
      v11_coef(127) =          -.295865726406382D+00
      v11_coef(128) =           .878113845650774D-01
      v11_coef(129) =           .540210774363566D-01
      v11_coef(130) =           .153826294458007D-01
      v11_coef(131) =          -.700218486345063D-01
      v11_coef(132) =           .350863908787738D-01
      v11_coef(133) =          -.312857034318616D-01
      v11_coef(134) =           .201358795851120D-01
      v11_coef(135) =           .627787978347838D-01
      v11_coef(136) =           .125083358082766D-01
      v11_coef(137) =          -.844610498185968D-01
      v11_coef(138) =          -.516336081649091D-01
      v11_coef(139) =          -.436701383143835D-01
      v11_coef(140) =           .235928458887405D-01
      v11_coef(141) =           .266872392004131D-01
      v11_coef(142) =           .123346064370267D-01
      v11_coef(143) =           .528746089573426D-02
      v11_coef(144) =          -.663377639988611D-02
      v11_coef(145) =          -.426933518657627D-01
      v11_coef(146) =          -.220353857001529D-01
      v11_coef(147) =          -.208494477148412D-01
      v11_coef(148) =           .138280234271630D-01
      v11_coef(149) =          -.312909873230440D-01
      v11_coef(150) =          -.368016404543060D-01
      v11_coef(151) =          -.118628337764509D+00
      v11_coef(152) =          -.262802278161427D-01
      v11_coef(153) =          -.387766918423711D-02
      v11_coef(154) =           .463661313501080D-01
      v11_coef(155) =          -.880165608043143D-01
      v11_coef(156) =          -.139361634940262D+00
      v11_coef(157) =          -.284769824694008D-01
      v11_coef(158) =          -.739643145666842D-01
      v11_coef(159) =          -.653200252942953D-01
      v11_coef(160) =           .243440516685545D-02
      v11_coef(161) =           .513221810043133D-01
      v11_coef(162) =           .144931387373819D+00
      v11_coef(163) =          -.414660140763216D-02
      v11_coef(164) =          -.131806071058625D-01
      v11_coef(165) =          -.438478473707597D-01
      v11_coef(166) =          -.281234783065008D-01
      v11_coef(167) =          -.221662377928734D-01
      v11_coef(168) =           .635328853437169D-02
      v11_coef(169) =          -.266217515168855D-01
      v11_coef(170) =          -.174281850139875D-01
      v11_coef(171) =          -.173759753684801D-01
      v11_coef(172) =          -.188857178137826D-01
      v11_coef(173) =          -.176615091372149D-01
      v11_coef(174) =           .222088328451502D-01
      v11_coef(175) =           .112669618460430D-01
      v11_coef(176) =          -.169402164212500D-01
      v11_coef(177) =          -.392269529718434D-02
      v11_coef(178) =          -.529704805065017D-02
      v11_coef(179) =          -.171571019350840D-01
      v11_coef(180) =          -.643877969416979D-02
      v11_coef(181) =          -.307165073351887D-01
      v11_coef(182) =          -.412223428357761D-01
      v11_coef(183) =          -.469500599059302D-01
      v11_coef(184) =          -.370962913147788D-02
      v11_coef(185) =           .510175755587685D-03
      v11_coef(186) =           .504553110409467D-01
      v11_coef(187) =           .326862399680935D-01
      v11_coef(188) =          -.760964308162337D-03
      v11_coef(189) =          -.190589075769478D-01
      v11_coef(190) =           .136456940387877D-01
      v11_coef(191) =           .139993589881527D-01
      v11_coef(192) =          -.104395359701188D-01
      v11_coef(193) =          -.385599660709391D-01
      v11_coef(194) =           .876750109927230D-02
      v11_coef(195) =          -.532056308871629D-02
      v11_coef(196) =          -.830006925898215D-01
      v11_coef(197) =           .121028594313734D-01
      v11_coef(198) =           .246905531377026D-02
      v11_coef(199) =           .555625880857731D-01
      v11_coef(200) =          -.194673279843153D-01
      v11_coef(201) =           .413774497507078D-02
      v11_coef(202) =          -.143749125706364D+00
      v11_coef(203) =           .603662152187368D-02
      v11_coef(204) =          -.118323036346506D-01
      v11_coef(205) =          -.602340231329006D-01
      v11_coef(206) =          -.971152679542949D-01
      v11_coef(207) =           .166068490415198D-01
      v11_coef(208) =          -.471174074859633D-01
      v11_coef(209) =          -.184774648803902D-01
      v11_coef(210) =          -.407193652535901D-01
      v11_coef(211) =           .154755617357465D-01
      v11_coef(212) =          -.195528352787702D-03
      v11_coef(213) =          -.507742189848796D-01
      v11_coef(214) =          -.325262208702358D-01
      v11_coef(215) =          -.262041181138888D-01
      v11_coef(216) =          -.385881004843293D-01
      v11_coef(217) =          -.655533428432900D-02
      v11_coef(218) =          -.176874841254003D-02
      v11_coef(219) =           .549885753617307D-02
      v11_coef(220) =           .408910494123857D-03
      v11_coef(221) =          -.206538117436390D-01
      v11_coef(222) =          -.709960956709433D-02
      v11_coef(223) =          -.110103226630496D-02
      v11_coef(224) =          -.104079945981196D-01
      v11_coef(225) =          -.234020056741342D-01
      v11_coef(226) =          -.272906394593809D-01
      v11_coef(227) =          -.152207796210721D-01
      v11_coef(228) =          -.221688400625713D-01
      v11_coef(229) =          -.112949772282996D-02
      v11_coef(230) =          -.367394809270406D-02
      v11_coef(231) =           .108824144073965D-01
      v11_coef(232) =           .916315724639774D-02
      v11_coef(233) =          -.938136897374941D-02
      v11_coef(234) =          -.711237917886920D-03
      v11_coef(235) =          -.493433457446326D-02
      v11_coef(236) =          -.534119316475695D-02
      v11_coef(237) =           .148900368835776D-01
      v11_coef(238) =           .829275163930797D-02
      v11_coef(239) =           .112576826444461D-01
      v11_coef(240) =          -.287038159817244D-01
      v11_coef(241) =           .313191286737516D-02
      v11_coef(242) =           .129405504479530D-01
      v11_coef(243) =          -.539799107715994D-02
      v11_coef(244) =          -.511950640758362D-02
      v11_coef(245) =           .216620212464002D-01
      v11_coef(246) =           .138425149344594D-01
      v11_coef(247) =           .260237787600616D-03
      v11_coef(248) =           .342791569626061D-01
      v11_coef(249) =           .246914837479954D-02
      v11_coef(250) =          -.782797095872977D-02
      v11_coef(251) =          -.450493345324109D-01
      v11_coef(252) =          -.486452144876258D-01
      v11_coef(253) =          -.249167254206512D-02
      v11_coef(254) =          -.248948094558414D-01
      v11_coef(255) =          -.169448384072330D-01
      v11_coef(256) =          -.934619257643325D-02
      v11_coef(257) =          -.539259305748807D-02
      v11_coef(258) =          -.165168744072807D-02
      v11_coef(259) =          -.109735494149397D+00
      v11_coef(260) =          -.382344632214011D-01
      v11_coef(261) =          -.327571890897241D-01
      v11_coef(262) =          -.119811340993719D-01
      v11_coef(263) =          -.356296835698013D-01
      v11_coef(264) =          -.240704022691441D-01
      v11_coef(265) =           .264102908027330D-03
      v11_coef(266) =          -.264135425108261D-01
      v11_coef(267) =          -.213819729781410D-01
      v11_coef(268) =          -.704519061886446D-01
      v11_coef(269) =          -.889105060032099D-02
      v11_coef(270) =          -.526549567858301D-01
      v11_coef(271) =          -.180857295216679D-01
      v11_coef(272) =          -.423069360138198D-01
      v11_coef(273) =          -.340279963353305D-01
      v11_coef(274) =          -.652486584349325D-01
      v11_coef(275) =          -.845820488984991D-02
      v11_coef(276) =          -.596792213495462D-01
      v11_coef(277) =          -.171876842246426D-01
      v11_coef(278) =          -.101775347311093D+00
      v11_coef(279) =          -.972798746997494D-02
      v11_coef(280) =          -.357598532283934D-03
      v11_coef(281) =           .857280893926254D-02
      v11_coef(282) =           .405553074454651D-01
      v11_coef(283) =           .494041107674973D-01
      v11_coef(284) =          -.207247119551366D-01
      v11_coef(285) =           .778159438434529D-01
      v11_coef(286) =          -.635318130044935D-02
      v11_coef(287) =           .425897167543368D-01
      v11_coef(288) =           .964255760398373D-02
      v11_coef(289) =          -.194757739140621D-02
      v11_coef(290) =          -.578410131810600D-01
      v11_coef(291) =          -.256870112753302D-01
      v11_coef(292) =          -.244874169186960D-03
      v11_coef(293) =           .338367556744255D+00
      v11_coef(294) =          -.166413326939488D+00
      v11_coef(295) =          -.901039448645037D-01
      v11_coef(296) =          -.309651729399113D+00
      v11_coef(297) =           .446517738036249D-01
      v11_coef(298) =          -.160598372575501D+00
      v11_coef(299) =           .137263023618034D+00
      v11_coef(300) =           .316831999711686D-02
      v11_coef(301) =           .193176415542934D-01
      v11_coef(302) =          -.303592555315285D-01
      v11_coef(303) =           .195012788159949D+00
      v11_coef(304) =           .346411793361092D-01
      v11_coef(305) =           .153955315548244D+00
      v11_coef(306) =           .393253913927521D-01
      v11_coef(307) =          -.111042771692550D-01
      v11_coef(308) =          -.333784030593727D-01
      v11_coef(309) =           .395024503003954D-02
      v11_coef(310) =          -.485838680988090D-01
      v11_coef(311) =           .478716500853634D-01
      v11_coef(312) =          -.938431011057889D-01
      v11_coef(313) =           .101468129850207D+00
      v11_coef(314) =          -.354514004893759D-01
      v11_coef(315) =           .287157991789547D+00
      v11_coef(316) =           .108579458829806D+00
      v11_coef(317) =           .120482223349789D+00
      v11_coef(318) =          -.266515216803023D-01
      v11_coef(319) =           .323943464187057D-01
      v11_coef(320) =          -.117541546081103D-01
      v11_coef(321) =           .602043780450081D-01
      v11_coef(322) =           .629888869842922D-01
      v11_coef(323) =          -.816754375236510D-01
      v11_coef(324) =           .127649138639020D+00
      v11_coef(325) =           .809579086918817D-02
      v11_coef(326) =          -.293706685658037D-01
      v11_coef(327) =          -.665473576725596D-01
      v11_coef(328) =          -.354491505247614D-01
      v11_coef(329) =           .224286603307050D+00
      v11_coef(330) =           .384923879662238D-01
      v11_coef(331) =           .205996843746790D+00
      v11_coef(332) =           .653403524796266D-01
      v11_coef(333) =          -.913459151883995D-02
      v11_coef(334) =           .239608654023763D-01
      v11_coef(335) =           .208844779940361D-01
      v11_coef(336) =           .617909556087610D-01
      v11_coef(337) =          -.821410817615651D-02
      v11_coef(338) =           .429251003474216D-01
      v11_coef(339) =          -.273654495775435D-02
      v11_coef(340) =          -.118184922038127D+00
      v11_coef(341) =          -.313498794992026D-01
      v11_coef(342) =          -.392241197417683D-01
      v11_coef(343) =          -.247643957224468D-01
      v11_coef(344) =           .562072830978347D-01
      v11_coef(345) =           .797804594199157D-01
      v11_coef(346) =           .153813073819638D+00
      v11_coef(347) =           .549399148312443D-01
      v11_coef(348) =           .109393695236805D+00
      v11_coef(349) =          -.630236935202841D-01
      v11_coef(350) =          -.103176460838620D-02
      v11_coef(351) =          -.654766656715895D-01
      v11_coef(352) =          -.308994259042248D+00
      v11_coef(353) =          -.152565551565462D+00
      v11_coef(354) =          -.466955250868466D-02
      v11_coef(355) =           .572512151257380D-01
      v11_coef(356) =           .266490129358826D-01
      v11_coef(357) =          -.461499555200563D-01
      v11_coef(358) =           .122896879618257D+00
      v11_coef(359) =           .708972773043692D-01
      v11_coef(360) =           .603193833190886D-01
      v11_coef(361) =           .376928573348212D-01
      v11_coef(362) =           .114650253993825D+00
      v11_coef(363) =          -.259612484342596D-01
      v11_coef(364) =           .959807325018328D-01
      v11_coef(365) =          -.131715858840476D-01
      v11_coef(366) =           .150071303849336D+00
      v11_coef(367) =           .820129766318929D-02
      v11_coef(368) =           .348884292513990D-01
      v11_coef(369) =           .148329175391929D-01
      v11_coef(370) =           .928700669405673D-01
      v11_coef(371) =           .276187616614044D+00
      v11_coef(372) =           .172880469922167D+00
      v11_coef(373) =          -.220332427673396D-01
      v11_coef(374) =           .523964924968657D-02
      v11_coef(375) =          -.247558941403721D-01
      v11_coef(376) =          -.112428032472341D+00
      v11_coef(377) =          -.219152491038873D-01
      v11_coef(378) =           .295934723853665D-02
      v11_coef(379) =           .160467783877390D-01
      v11_coef(380) =          -.780466918232561D-02
      v11_coef(381) =           .210349140002074D-01
      v11_coef(382) =          -.465119575001727D-02
      v11_coef(383) =           .232143676785465D-01
      v11_coef(384) =           .554897775069779D-02
      v11_coef(385) =           .530072213621465D-01
      v11_coef(386) =           .484358244330546D-02
      v11_coef(387) =           .220066291925206D-01
      v11_coef(388) =           .281390089970629D-01
      v11_coef(389) =           .585477011499957D-02
      v11_coef(390) =           .259417559068442D-02
      v11_coef(391) =           .459359976940168D-02
      v11_coef(392) =          -.272886663152264D-01
      v11_coef(393) =          -.338198661627150D-01
      v11_coef(394) =          -.105345114163524D-01
      v11_coef(395) =           .281410304768763D-01
      v11_coef(396) =           .756327518590527D-02
      v11_coef(397) =           .478904535016618D-01
      v11_coef(398) =           .165290098375131D-01
      v11_coef(399) =           .616043874468421D-01
      v11_coef(400) =           .352292714142499D-02
      v11_coef(401) =           .429489830450771D-01
      v11_coef(402) =          -.157778381421182D-01
      v11_coef(403) =           .134161522012849D-01
      v11_coef(404) =           .344467191514942D-01
      v11_coef(405) =           .555535038590945D-01
      v11_coef(406) =           .753685435242634D-01
      v11_coef(407) =           .482140604740348D-01
      v11_coef(408) =           .418949676728744D-01
      v11_coef(409) =           .441071518111609D-02
      v11_coef(410) =           .172299944498003D-01
      v11_coef(411) =           .241861300814859D-01
      v11_coef(412) =           .106031385122276D+00
      v11_coef(413) =           .228267354010760D-03
      v11_coef(414) =           .601072806061447D-01
      v11_coef(415) =          -.216429573918303D-01
      v11_coef(416) =           .199021213666722D-01
      v11_coef(417) =           .114464184868705D+00
      v11_coef(418) =           .206997302873351D+00
      v11_coef(419) =          -.836349144036638D-03
      v11_coef(420) =           .820161330592509D-01
      v11_coef(421) =           .973126943865046D-02
      v11_coef(422) =           .476042660762714D-01
      v11_coef(423) =           .120164757668479D-01
      v11_coef(424) =           .165873606279228D-02
      v11_coef(425) =          -.111171232005366D-01
      v11_coef(426) =          -.802852546501819D-02
      v11_coef(427) =           .166289951913007D-02
      v11_coef(428) =          -.285470042814435D-01
      v11_coef(429) =           .128339125017237D-01
      v11_coef(430) =           .272169704735406D-01
      v11_coef(431) =          -.104963908322166D-01
      v11_coef(432) =          -.549933384384853D-02
      v11_coef(433) =           .213331636222677D-01
      v11_coef(434) =           .131567423741525D+00
      v11_coef(435) =           .200450869195918D-01
      v11_coef(436) =           .999622624965039D-01
      v11_coef(437) =           .136143558317656D-01
      v11_coef(438) =          -.134917109684706D-01
      v11_coef(439) =           .593887438573648D-01
      v11_coef(440) =           .127432519851403D+00
      v11_coef(441) =           .334116182981812D-02
      v11_coef(442) =           .157420108768887D+00
      v11_coef(443) =           .511552892998489D-02
      v11_coef(444) =           .148684888176530D+00
      v11_coef(445) =           .749606816645737D-02
      v11_coef(446) =          -.115672936092500D-02
      v11_coef(447) =          -.142595195631329D-01
      v11_coef(448) =          -.433421257202979D-01
      v11_coef(449) =           .566882758651695D-01
      v11_coef(450) =           .160448135224549D-01
      v11_coef(451) =           .117358061166737D+00
      v11_coef(452) =          -.658394329604236D-02
      v11_coef(453) =           .162561526600649D+00
      v11_coef(454) =           .329603856418617D-01
      v11_coef(455) =           .216716141944740D-02
      v11_coef(456) =           .142288259262482D+00
      v11_coef(457) =           .281221121858149D-01
      v11_coef(458) =           .464437154468969D-02
      v11_coef(459) =           .959584196858190D-02
      v11_coef(460) =           .234062695738124D-01
      v11_coef(461) =           .150063993705476D-01
      v11_coef(462) =           .230721596753827D-01
      v11_coef(463) =           .464558546272196D-02
      v11_coef(464) =           .124055406376626D-01
      v11_coef(465) =          -.141477455069323D-01
      v11_coef(466) =          -.436111814978379D-02
      v11_coef(467) =           .349728503700614D-02
      v11_coef(468) =           .791118969122994D-02
      v11_coef(469) =          -.528613182388442D-03
      v11_coef(470) =           .174863455209975D-01
      v11_coef(471) =           .144936304498201D-01
      v11_coef(472) =           .347542311189408D-01
      v11_coef(473) =           .233158696495471D-01
      v11_coef(474) =           .389436823154680D-01
      v11_coef(475) =           .189287515764432D-01
      v11_coef(476) =           .655144710987667D-02
      v11_coef(477) =          -.328125061566843D-02
      v11_coef(478) =          -.878980687217015D-03
      v11_coef(479) =           .934425229715450D-02
      v11_coef(480) =          -.487146398773719D-02
      v11_coef(481) =           .212416429002233D-01
      v11_coef(482) =          -.735915986139348D-03
      v11_coef(483) =          -.155936854683210D-02
      v11_coef(484) =           .701988134922257D-02
      v11_coef(485) =          -.318350375446837D-02
      v11_coef(486) =          -.624055472991454D-02
      v11_coef(487) =          -.323374583657686D-01
      v11_coef(488) =          -.242879622021703D-01
      v11_coef(489) =           .154625726607483D-02
      v11_coef(490) =          -.452844930541334D-02
      v11_coef(491) =          -.267232386349589D-01
      v11_coef(492) =          -.320318422328518D-01
      v11_coef(493) =           .122942004331669D-01
      v11_coef(494) =          -.476055755759628D-01
      v11_coef(495) =           .108289056674010D-01
      v11_coef(496) =          -.804331810263514D-03
      v11_coef(497) =           .304643500084620D-01
      v11_coef(498) =          -.230569640163026D-01
      v11_coef(499) =           .374505574180547D-02
      v11_coef(500) =          -.547169377028415D-01
      v11_coef(501) =           .161432436980168D-01
      v11_coef(502) =           .120737710018527D-01
      v11_coef(503) =           .118124837634381D-02
      v11_coef(504) =          -.177154399110840D-02
      v11_coef(505) =          -.780893658288623D-02
      v11_coef(506) =          -.127953686291192D-01
      v11_coef(507) =           .761088237570405D-02
      v11_coef(508) =           .131870719366509D-01
      v11_coef(509) =           .224070752720717D-01
      v11_coef(510) =           .101422592959260D-01
      v11_coef(511) =           .327599015240514D-02
      v11_coef(512) =           .154139380767396D-01
      v11_coef(513) =          -.519485320153888D-03
      v11_coef(514) =           .188201833127421D-01
      v11_coef(515) =          -.804499788904267D-02
      v11_coef(516) =          -.453000398177471D-02
      v11_coef(517) =           .382480528017517D-02
      v11_coef(518) =           .118896523384318D-02
      v11_coef(519) =          -.121078621920566D-01
      v11_coef(520) =          -.178211593572931D-02
      v11_coef(521) =           .115438639613402D-01
      v11_coef(522) =          -.126865974807260D-01
      v11_coef(523) =           .285413072540745D-02
      v11_coef(524) =          -.101761043749581D+00
      v11_coef(525) =           .203896781514576D-02
      v11_coef(526) =           .232370976864494D-02
      v11_coef(527) =           .261267345353668D-01
      v11_coef(528) =          -.207717571877320D-01
      v11_coef(529) =          -.672260569513158D-01
      v11_coef(530) =           .151941109849023D-01
      v11_coef(531) =          -.555253139116409D-01
      v11_coef(532) =           .302754119611336D-02
      v11_coef(533) =          -.578170924715681D-01
      v11_coef(534) =          -.233260339854286D-01
      v11_coef(535) =          -.388146731485132D-02
      v11_coef(536) =          -.189395778263176D-01
      v11_coef(537) =           .264983014097932D-01
      v11_coef(538) =          -.262490346073127D-02
      v11_coef(539) =           .312553631988926D-01
      v11_coef(540) =           .269910212845928D-01
      v11_coef(541) =           .189708136437921D-01
      v11_coef(542) =           .139658543551476D-01
      v11_coef(543) =           .266230080586718D-01
      v11_coef(544) =           .370651959833806D-01
      v11_coef(545) =          -.104617467076160D-02
      v11_coef(546) =           .173933239950098D-01
      v11_coef(547) =           .269886703382371D-01
      v11_coef(548) =           .641231588435077D-01
      v11_coef(549) =           .507865670616374D-02
      v11_coef(550) =           .102981705899407D-01
      v11_coef(551) =           .143744093883426D-01
      v11_coef(552) =           .205592012375713D-01
      v11_coef(553) =           .175234031387142D-01
      v11_coef(554) =           .102906719568718D-01
      v11_coef(555) =           .307722673680579D-02
      v11_coef(556) =           .350181813104456D-01
      v11_coef(557) =           .162553928731895D-01
      v11_coef(558) =           .576805770770454D-01
      v11_coef(559) =           .531372577315336D-02
      v11_coef(560) =          -.948215499536081D-03
      v11_coef(561) =          -.838357853061022D-02
      v11_coef(562) =          -.306470712165796D-01
      v11_coef(563) =          -.116989237323511D-01
      v11_coef(564) =           .403266646357264D-02
      v11_coef(565) =          -.138175063359312D-01
      v11_coef(566) =           .545013435060785D-02
      v11_coef(567) =           .250305633306568D-01
      v11_coef(568) =          -.767214396839990D-03
      v11_coef(569) =           .691628191823850D-03
      v11_coef(570) =           .657906294834191D-01
      v11_coef(571) =           .141638964147098D-01
      v11_coef(572) =          -.485692292570995D-03
      v11_coef(573) =          -.203710799310868D-02
      v11_coef(574) =           .315106907228132D-01
      v11_coef(575) =           .415129722358783D-01
      v11_coef(576) =          -.671535857547793D-02
      v11_coef(577) =           .513659879381438D-01
      v11_coef(578) =          -.911933529102954D-02
      v11_coef(579) =           .900511267096139D-02
      v11_coef(580) =           .217163050643882D-01
      v11_coef(581) =          -.407234866501102D-02
      v11_coef(582) =          -.194290835091664D-01
      v11_coef(583) =           .933935456287878D-02
      v11_coef(584) =           .206134811943064D-02
      v11_coef(585) =           .200733842054957D-01
      v11_coef(586) =          -.167218404406672D-01
      v11_coef(587) =          -.341305970395033D-02
      v11_coef(588) =           .113926905559910D-03

c
c     higher-order coeffecients in vts
c

      vts_coef(  1) =          -.700785857851808D+00
      vts_coef(  2) =          -.256410285760228D-01
      vts_coef(  3) =          -.464875430300624D-02
      vts_coef(  4) =           .193583501881592D-01
      vts_coef(  5) =           .274411526626973D-01
      vts_coef(  6) =           .908830977229824D-02
      vts_coef(  7) =          -.897664872868108D-02
      vts_coef(  8) =          -.503343744960494D-02
      vts_coef(  9) =          -.476378197823153D-01
      vts_coef( 10) =          -.224539451836736D-01
      vts_coef( 11) =          -.358807308382561D-01
      vts_coef( 12) =           .152458195506508D-01
      vts_coef( 13) =          -.278905197114855D-01
      vts_coef( 14) =          -.469875475056177D-01
      vts_coef( 15) =          -.397920868095802D-01
      vts_coef( 16) =          -.251959736095409D-01
      vts_coef( 17) =          -.143100558169463D+00
      vts_coef( 18) =           .191605419654552D-01
      vts_coef( 19) =           .836728728893659D-01
      vts_coef( 20) =          -.174750788525364D-01
      vts_coef( 21) =          -.267280175251251D-01
      vts_coef( 22) =          -.244614200848969D-01
      vts_coef( 23) =          -.270812546901320D-02
      vts_coef( 24) =          -.237366169106116D+00
      vts_coef( 25) =          -.308096013694611D-01
      vts_coef( 26) =          -.245484623437817D-01
      vts_coef( 27) =          -.198919075642687D-01
      vts_coef( 28) =          -.325578025900052D-01
      vts_coef( 29) =          -.136932478746946D+00
      vts_coef( 30) =           .645235750348750D-02
      vts_coef( 31) =          -.622090630881867D-01
      vts_coef( 32) =           .103249519001122D+00
      vts_coef( 33) =          -.672731673428250D-01
      vts_coef( 34) =          -.264865263230371D-02
      vts_coef( 35) =          -.221913997573494D-01
      vts_coef( 36) =          -.446276646707346D-01
      vts_coef( 37) =           .304986172870667D-01
      vts_coef( 38) =           .249051423327123D-01
      vts_coef( 39) =          -.102285474208043D-01
      vts_coef( 40) =          -.369849128355745D-02
      vts_coef( 41) =           .598095850417904D+00
      vts_coef( 42) =          -.301448041205120D+00
      vts_coef( 43) =           .201752762062047D+00
      vts_coef( 44) =          -.114150570266662D+00
      vts_coef( 45) =           .679932478355907D-01
      vts_coef( 46) =           .239480075827712D+00
      vts_coef( 47) =          -.384246971092749D-01
      vts_coef( 48) =           .684722434783497D+00
      vts_coef( 49) =          -.259309254423934D+00
      vts_coef( 50) =           .693662130840478D-01
      vts_coef( 51) =          -.140901940188162D+00
      vts_coef( 52) =           .302206442617509D+00
      vts_coef( 53) =           .177477173706272D+00
      vts_coef( 54) =           .422104555719124D+00
      vts_coef( 55) =           .376176517867204D-01
      vts_coef( 56) =           .105983031802273D-01
      vts_coef( 57) =           .121467609278021D-01
      vts_coef( 58) =           .138810932498329D-01
      vts_coef( 59) =          -.484831966608905D-01
      vts_coef( 60) =           .871929283050537D-01
      vts_coef( 61) =          -.106034570217459D-01
      vts_coef( 62) =           .738518718475931D-01
      vts_coef( 63) =          -.117763981146575D-01
      vts_coef( 64) =           .125499014424839D+00
      vts_coef( 65) =          -.777934343401038D-01
      vts_coef( 66) =           .167523957513042D+00
      vts_coef( 67) =           .323619860594206D-01
      vts_coef( 68) =           .593014352159515D-02
      vts_coef( 69) =          -.483713668603926D-01
      vts_coef( 70) =          -.876713357351299D-01
      vts_coef( 71) =          -.287678365425564D+00
      vts_coef( 72) =          -.680912784470543D-01
      vts_coef( 73) =          -.113500734213726D+00
      vts_coef( 74) =          -.194002856014124D-01
      vts_coef( 75) =          -.301509981293087D-01
      vts_coef( 76) =          -.254002498955981D+00
      vts_coef( 77) =           .394428317271503D-01
      vts_coef( 78) =          -.506060767817123D-01
      vts_coef( 79) =          -.755631143039891D-01
      vts_coef( 80) =           .515897309609449D-01
      vts_coef( 81) =           .122399746430911D+00
      vts_coef( 82) =          -.899570255934039D-01
      vts_coef( 83) =          -.100963595504469D+00
      vts_coef( 84) =          -.129264966130576D-01
      vts_coef( 85) =           .608285225140576D-01
      vts_coef( 86) =          -.625739015467774D-01
      vts_coef( 87) =           .410192886492970D-01
      vts_coef( 88) =          -.344151317468881D+00
      vts_coef( 89) =           .287723033811313D-01
      vts_coef( 90) =           .316016585391232D-01
      vts_coef( 91) =           .910854166421132D-02
      vts_coef( 92) =          -.117142856240357D-03
      vts_coef( 93) =           .358916273508908D-01
      vts_coef( 94) =          -.600813226659610D-02
      vts_coef( 95) =           .343613194785494D+00
      vts_coef( 96) =           .292559055658030D-02
      vts_coef( 97) =          -.366377535468716D-02
      vts_coef( 98) =          -.824946183598161D-02
      vts_coef( 99) =           .501808120877796D-01
      vts_coef(100) =          -.152142475725208D+00
      vts_coef(101) =          -.454476514007403D-02
      vts_coef(102) =          -.139791870209256D-01
      vts_coef(103) =          -.212832433923531D-01
      vts_coef(104) =           .192377540419800D-02
      vts_coef(105) =           .149187775322191D-01
      vts_coef(106) =           .853367343917652D-01
      vts_coef(107) =          -.341039576765497D-02
      vts_coef(108) =           .855387261971914D-01
      vts_coef(109) =           .189135805239344D-01
      vts_coef(110) =           .372471064368073D-01
      vts_coef(111) =          -.485841239026914D-02
      vts_coef(112) =           .810642390945937D-02
      vts_coef(113) =           .386821012940341D-01
      vts_coef(114) =           .100107181123604D+00
      vts_coef(115) =          -.235127022446562D-01
      vts_coef(116) =           .539961798768331D-02
      vts_coef(117) =          -.208280152521792D-01
      vts_coef(118) =           .452997355781822D+00
      vts_coef(119) =           .970719054795656D-01
      vts_coef(120) =           .993952553519890D-01
      vts_coef(121) =          -.353062019297381D-01
      vts_coef(122) =           .589452590934055D-01
      vts_coef(123) =          -.407431138593605D-02
      vts_coef(124) =           .901797818819266D-02
      vts_coef(125) =          -.115806763374825D-01
      vts_coef(126) =          -.230607923264092D-02
      vts_coef(127) =          -.177088443425026D-01
      vts_coef(128) =           .112408723686948D+01
      vts_coef(129) =          -.681913683299385D+00
      vts_coef(130) =           .230108797458810D-01
      vts_coef(131) =           .179045282836246D+00
      vts_coef(132) =          -.700624529323003D+00
      vts_coef(133) =           .145583856261140D+01
      vts_coef(134) =          -.105491892474285D+01
      vts_coef(135) =          -.221162866435194D+00
      vts_coef(136) =          -.705579342344599D+00
      vts_coef(137) =           .105816324575670D+01
      vts_coef(138) =          -.780317650921450D-01
      vts_coef(139) =          -.355998686245382D+00
      vts_coef(140) =           .101257760934130D+00
      vts_coef(141) =          -.852620088864061D+00
      vts_coef(142) =           .753847019455853D-01
      vts_coef(143) =          -.378278020832507D+00
      vts_coef(144) =           .301460334173006D+00
      vts_coef(145) =          -.614630731374315D+00
      vts_coef(146) =           .727402840695987D+00
      vts_coef(147) =          -.538814632468429D+00
      vts_coef(148) =          -.351866909264936D+00
      vts_coef(149) =          -.155528763706048D+01
      vts_coef(150) =           .392212037647259D+00
      vts_coef(151) =           .160409146815567D+00
      vts_coef(152) =          -.244689968664319D+00
      vts_coef(153) =          -.120548887923292D+00
      vts_coef(154) =           .183817504508807D+00
      vts_coef(155) =          -.532993701085256D-01
      vts_coef(156) =           .557305724461756D+00
      vts_coef(157) =          -.108381786612190D+01
      vts_coef(158) =          -.177470572534647D+01
      vts_coef(159) =           .300673881896437D+00
      vts_coef(160) =          -.553740707579137D-01
      vts_coef(161) =          -.268303874871718D+00
      vts_coef(162) =          -.951550526815552D+00
      vts_coef(163) =          -.818053214333986D+00
      vts_coef(164) =          -.241392724671540D+00
      vts_coef(165) =           .296383932113084D+00
      vts_coef(166) =          -.129612814291386D+00
      vts_coef(167) =          -.347260855676481D-01
      vts_coef(168) =           .890908462690230D-01
      vts_coef(169) =           .975922876298838D-01
      vts_coef(170) =           .110881652264478D-01
      vts_coef(171) =          -.398507371711818D-02
      vts_coef(172) =           .315391711302042D+00
      vts_coef(173) =           .677455135970386D-01
      vts_coef(174) =          -.663046910967574D+00
      vts_coef(175) =           .144958401341242D+01
      vts_coef(176) =          -.332454305381836D+00
      vts_coef(177) =           .382806730891755D-01
      vts_coef(178) =          -.202128332208664D+00
      vts_coef(179) =          -.110198693187174D+01
      vts_coef(180) =           .582807606944795D-01
      vts_coef(181) =          -.673146292026748D+00
      vts_coef(182) =          -.446983707015308D+00
      vts_coef(183) =           .353362866944260D+00
      vts_coef(184) =          -.784505738560501D-01
      vts_coef(185) =          -.960663576603207D+00
      vts_coef(186) =           .145000572815117D+01
      vts_coef(187) =          -.754556760374984D-01
      vts_coef(188) =           .148515624952196D+00
      vts_coef(189) =           .155775053252601D+00
      vts_coef(190) =           .328145479511563D-02
      vts_coef(191) =           .620324803028645D-01
      vts_coef(192) =           .173574747366009D+00
      vts_coef(193) =           .103082449427952D+00
      vts_coef(194) =          -.422841760210180D+00
      vts_coef(195) =           .273127383055519D+00
      vts_coef(196) =           .250313105968807D-01
      vts_coef(197) =          -.809780228989513D-03
      vts_coef(198) =          -.433263440490701D+00
      vts_coef(199) =          -.824429001537848D-01
      vts_coef(200) =          -.188721141028543D-01
      vts_coef(201) =          -.726061919130893D-01
      vts_coef(202) =          -.567867407024245D-01
      vts_coef(203) =          -.102693580254160D+00
      vts_coef(204) =          -.143734487725581D+00
      vts_coef(205) =          -.320643539532896D+00
      vts_coef(206) =           .102568875353014D+01
      vts_coef(207) =           .137146917094955D+00
      vts_coef(208) =          -.450776803211225D+00
      vts_coef(209) =           .684633231491316D-01
      vts_coef(210) =           .620209569271602D-01
      vts_coef(211) =          -.102398087702097D+00
      vts_coef(212) =           .452335078891675D-02
      vts_coef(213) =          -.664885695517196D-01
      vts_coef(214) =           .109277607055143D+00
      vts_coef(215) =           .432876541872606D-01
      vts_coef(216) =           .250742675463093D-01
      vts_coef(217) =          -.140285522608181D+00
      vts_coef(218) =           .596739174469895D-02
      vts_coef(219) =          -.983548858902155D-01
      vts_coef(220) =           .325891928514221D+00
      vts_coef(221) =          -.328434140000514D-01
      vts_coef(222) =          -.324361048842937D-01
      vts_coef(223) =          -.566635845269292D-02
      vts_coef(224) =           .940821573993470D-01
      vts_coef(225) =          -.298224090046186D-01
      vts_coef(226) =           .761674153697664D-01
      vts_coef(227) =          -.659233277011709D+00
      vts_coef(228) =           .776521913251219D-01
      vts_coef(229) =          -.132729687539573D-01
      vts_coef(230) =           .293327334301865D-01
      vts_coef(231) =          -.337440438281374D+00
      vts_coef(232) =          -.571667987357650D-01
      vts_coef(233) =           .164920902665429D-01
      vts_coef(234) =           .386609636199271D-01
      vts_coef(235) =           .423938893562620D-02
      vts_coef(236) =           .155868522947074D-01
      vts_coef(237) =          -.503083377593375D-01
      vts_coef(238) =          -.325278292111508D-01
      vts_coef(239) =          -.903089974508740D-01
      vts_coef(240) =           .508201490507810D+00
      vts_coef(241) =          -.405814349802903D+00
      vts_coef(242) =          -.353568075280347D-01
      vts_coef(243) =          -.233830968769662D-01
      vts_coef(244) =           .539351121908607D+00
      vts_coef(245) =          -.224109566120736D+00
      vts_coef(246) =          -.344939158126378D+00
      vts_coef(247) =          -.582174158791842D-01
      vts_coef(248) =           .321057580283565D-01
      vts_coef(249) =          -.343302081842629D-01
      vts_coef(250) =           .428883528451083D+00
      vts_coef(251) =           .355293072741652D-01
      vts_coef(252) =          -.415555679632874D+00
      vts_coef(253) =           .494820321280011D-01
      vts_coef(254) =          -.426768901127833D-01
      vts_coef(255) =          -.351184076558852D-01
      vts_coef(256) =           .249443217907335D-01
      vts_coef(257) =           .408507694328180D-02
      vts_coef(258) =           .141189448797326D-02
      vts_coef(259) =           .656537091623796D-02
      vts_coef(260) =           .948817979849126D-02
      vts_coef(261) =           .276024149957511D-01
      vts_coef(262) =          -.248066655040594D+00
      vts_coef(263) =           .114280730331307D+00
      vts_coef(264) =           .396287521983953D-01
      vts_coef(265) =          -.625600809610784D-02
      vts_coef(266) =          -.611786556528715D+00
      vts_coef(267) =          -.164932624494550D+00
      vts_coef(268) =           .277283870692965D+00
      vts_coef(269) =          -.222982818665999D-01
      vts_coef(270) =          -.385039105636265D-01
      vts_coef(271) =           .303551288522741D-01
      vts_coef(272) =          -.765802905198138D+00
      vts_coef(273) =          -.684004651834833D+00
      vts_coef(274) =           .461258607158342D+00
      vts_coef(275) =           .156371240827039D-01
      vts_coef(276) =           .314102430180770D-01
      vts_coef(277) =           .238899374102570D-01
      vts_coef(278) =          -.583465993064014D-01
      vts_coef(279) =           .941887178577657D-02
      vts_coef(280) =          -.365730502016275D-04
      vts_coef(281) =          -.580794492108725D+00
      vts_coef(282) =           .158437507020977D+00
      vts_coef(283) =          -.257672751429824D+00
      vts_coef(284) =           .430896377914213D-01
      vts_coef(285) =           .367451057435395D+00
      vts_coef(286) =          -.909018350397069D-01
      vts_coef(287) =          -.710375183433897D-01
      vts_coef(288) =           .439954377212657D-01
      vts_coef(289) =           .237931095009891D-02
      vts_coef(290) =          -.335178788427873D-01
      vts_coef(291) =          -.139339419664156D-01
      vts_coef(292) =          -.169101953010853D-02
      vts_coef(293) =          -.136741898592394D+01
      vts_coef(294) =           .184453242241875D+01
      vts_coef(295) =           .225355097616369D+00
      vts_coef(296) =           .380578972402804D+01
      vts_coef(297) =           .114810385747253D+01
      vts_coef(298) =           .382862521257530D+00
      vts_coef(299) =           .122143786909263D+01
      vts_coef(300) =           .451462079620589D+01
      vts_coef(301) =           .282747620423474D+01
      vts_coef(302) =          -.182625426066183D+01
      vts_coef(303) =          -.303254238803865D+01
      vts_coef(304) =          -.369159596833052D+01
      vts_coef(305) =           .140282590154692D+01
      vts_coef(306) =           .160418055927872D+01
      vts_coef(307) =          -.152068809043982D+01
      vts_coef(308) =          -.463200627933537D+00
      vts_coef(309) =           .102199824228865D+01
      vts_coef(310) =          -.180717162508762D+01
      vts_coef(311) =           .197037085136738D+01
      vts_coef(312) =           .263440229810947D+00
      vts_coef(313) =           .105908541353620D+01
      vts_coef(314) =           .161761489579766D+00
      vts_coef(315) =           .230064377675563D+00
      vts_coef(316) =           .283699508286581D+01
      vts_coef(317) =           .863581935397112D+00
      vts_coef(318) =           .769393702120367D+00
      vts_coef(319) =          -.672261719229268D-01
      vts_coef(320) =          -.275311440512597D+00
      vts_coef(321) =          -.681237166775425D+00
      vts_coef(322) =          -.111755837433252D+01
      vts_coef(323) =          -.886063696840805D+00
      vts_coef(324) =          -.167084146442141D+01
      vts_coef(325) =          -.292179045332423D+01
      vts_coef(326) =           .726250577527915D+00
      vts_coef(327) =          -.161340866801356D+01
      vts_coef(328) =           .939204437357875D+00
      vts_coef(329) =          -.589737197319132D+01
      vts_coef(330) =           .102828641757048D+01
      vts_coef(331) =          -.593164207330797D+00
      vts_coef(332) =           .444133644399952D+00
      vts_coef(333) =           .851330825841588D-01
      vts_coef(334) =          -.418379823656323D+00
      vts_coef(335) =          -.145846308117020D+01
      vts_coef(336) =          -.431906976802528D+00
      vts_coef(337) =           .945549357893709D+00
      vts_coef(338) =           .188464633537255D+01
      vts_coef(339) =           .512181965141656D+00
      vts_coef(340) =          -.574849231112551D+01
      vts_coef(341) =          -.416415063157415D+00
      vts_coef(342) =          -.162331136153871D+00
      vts_coef(343) =           .533891867192677D-01
      vts_coef(344) =          -.591779323970214D+00
      vts_coef(345) =          -.526780643133376D+00
      vts_coef(346) =           .170446008141895D+00
      vts_coef(347) =           .471958941489183D+01
      vts_coef(348) =          -.247286145805387D+01
      vts_coef(349) =           .685267615970524D+00
      vts_coef(350) =           .156267211611140D+00
      vts_coef(351) =          -.645608039097460D+01
      vts_coef(352) =          -.274284303280372D+01
      vts_coef(353) =          -.447020911450267D-01
      vts_coef(354) =           .171039732980869D-01
      vts_coef(355) =          -.460180291978673D+00
      vts_coef(356) =           .912506089342619D-01
      vts_coef(357) =           .442692599345668D+00
      vts_coef(358) =           .268186925160464D+01
      vts_coef(359) =          -.271874437413266D+00
      vts_coef(360) =           .536937594122371D+01
      vts_coef(361) =           .363322720356417D+01
      vts_coef(362) =          -.141441724533210D+01
      vts_coef(363) =          -.225820192024898D+00
      vts_coef(364) =           .717644374490289D+01
      vts_coef(365) =           .757045727249801D+00
      vts_coef(366) =           .116920544167587D+01
      vts_coef(367) =          -.500620567917158D+00
      vts_coef(368) =          -.964457870380175D+00
      vts_coef(369) =          -.153901856863069D+00
      vts_coef(370) =          -.537085477256555D+00
      vts_coef(371) =          -.381330255334614D+01
      vts_coef(372) =           .534278757159880D+00
      vts_coef(373) =          -.397050329876786D+00
      vts_coef(374) =          -.106730879654165D+01
      vts_coef(375) =          -.967611999121862D+00
      vts_coef(376) =          -.224746455094079D+00
      vts_coef(377) =           .145766941441959D+00
      vts_coef(378) =           .827197060966096D-02
      vts_coef(379) =           .240972865617016D-01
      vts_coef(380) =          -.107640142713938D+00
      vts_coef(381) =           .213774516063609D-01
      vts_coef(382) =          -.261319299112538D-01
      vts_coef(383) =           .109558631846483D+00
      vts_coef(384) =          -.346202900942426D+00
      vts_coef(385) =          -.634171285646528D-01
      vts_coef(386) =           .565525852466905D+00
      vts_coef(387) =          -.699961387627779D+00
      vts_coef(388) =           .115733500204958D+00
      vts_coef(389) =          -.350767382423527D-01
      vts_coef(390) =           .696485424277745D-01
      vts_coef(391) =           .106303370967301D+01
      vts_coef(392) =           .482382272234691D-01
      vts_coef(393) =           .241154267234787D+00
      vts_coef(394) =           .389131959084386D+00
      vts_coef(395) =          -.221202690811707D+00
      vts_coef(396) =           .114212325899002D+00
      vts_coef(397) =          -.428572886467885D-01
      vts_coef(398) =           .773622782910685D+00
      vts_coef(399) =          -.534783503666799D+00
      vts_coef(400) =          -.329542616944484D-01
      vts_coef(401) =          -.240045910255696D+00
      vts_coef(402) =           .147479842528248D-01
      vts_coef(403) =          -.842816674504978D-01
      vts_coef(404) =          -.769317569080377D+00
      vts_coef(405) =           .422393454604714D-01
      vts_coef(406) =          -.114637674146957D+01
      vts_coef(407) =           .131122255641673D+01
      vts_coef(408) =          -.140185097658412D+00
      vts_coef(409) =          -.138497322306081D+00
      vts_coef(410) =          -.266307387521704D+01
      vts_coef(411) =           .188801753495950D+01
      vts_coef(412) =           .704620069386020D+00
      vts_coef(413) =          -.465786141959451D-01
      vts_coef(414) =           .478607612446819D-02
      vts_coef(415) =           .677271679609858D-02
      vts_coef(416) =          -.195955368524237D+01
      vts_coef(417) =           .187717800879650D+01
      vts_coef(418) =          -.525624790844074D+00
      vts_coef(419) =           .111487988282530D+00
      vts_coef(420) =          -.607328836768854D+00
      vts_coef(421) =          -.138196414993366D+00
      vts_coef(422) =          -.718745847845114D-01
      vts_coef(423) =           .810688177171872D-01
      vts_coef(424) =          -.793557879036296D-02
      vts_coef(425) =           .637933410650138D-01
      vts_coef(426) =           .475984842698139D+00
      vts_coef(427) =           .617816555768019D-01
      vts_coef(428) =           .108179733709134D+01
      vts_coef(429) =           .117636343518963D+00
      vts_coef(430) =           .408030137946453D+00
      vts_coef(431) =           .820451020244380D-01
      vts_coef(432) =           .107691988121377D+01
      vts_coef(433) =           .296364566755865D+01
      vts_coef(434) =          -.508836466221659D-01
      vts_coef(435) =           .909649085402351D-01
      vts_coef(436) =           .556723665598395D+00
      vts_coef(437) =          -.196020265090125D-01
      vts_coef(438) =           .119460715279219D+00
      vts_coef(439) =           .604912553325024D+01
      vts_coef(440) =           .229483222865539D-01
      vts_coef(441) =           .282477983331233D+00
      vts_coef(442) =           .941360442652403D+00
      vts_coef(443) =          -.414838598448886D+00
      vts_coef(444) =           .169132321741675D+00
      vts_coef(445) =           .807403390603875D-01
      vts_coef(446) =          -.643955840737119D-02
      vts_coef(447) =           .456678592382777D-01
      vts_coef(448) =           .998022916917933D+00
      vts_coef(449) =           .414565745007650D+01
      vts_coef(450) =           .591282910882827D+00
      vts_coef(451) =          -.153712317794784D+01
      vts_coef(452) =          -.288996634006850D+00
      vts_coef(453) =           .162102526273802D+00
      vts_coef(454) =          -.102952017960468D+00
      vts_coef(455) =          -.168545036816072D-01
      vts_coef(456) =           .479079139607619D-01
      vts_coef(457) =          -.117881826553221D-01
      vts_coef(458) =          -.416255423600099D-01
      vts_coef(459) =           .160106584753755D-01
      vts_coef(460) =          -.480491392363025D-01
      vts_coef(461) =          -.314657844136920D-01
      vts_coef(462) =          -.131327676920993D-01
      vts_coef(463) =           .980061934265521D-01
      vts_coef(464) =           .852356647136213D-01
      vts_coef(465) =           .594689832105523D-01
      vts_coef(466) =          -.153455611637951D+00
      vts_coef(467) =           .101729299230930D-02
      vts_coef(468) =           .104743766301494D-01
      vts_coef(469) =           .373937430904842D-03
      vts_coef(470) =          -.859727713486856D-01
      vts_coef(471) =           .130008397943503D+00
      vts_coef(472) =          -.752676595246769D-01
      vts_coef(473) =           .119905132229140D+01
      vts_coef(474) =           .284814978894916D-01
      vts_coef(475) =          -.132858749565338D-01
      vts_coef(476) =          -.166303555760728D-01
      vts_coef(477) =           .103946476954948D+01
      vts_coef(478) =           .141604682180868D+00
      vts_coef(479) =          -.884238190406692D-01
      vts_coef(480) =          -.488148864348646D-01
      vts_coef(481) =           .103563330240406D-01
      vts_coef(482) =          -.140873374383606D-01
      vts_coef(483) =           .278144219103063D-01
      vts_coef(484) =          -.293950092158801D+00
      vts_coef(485) =           .118480665017130D+00
      vts_coef(486) =          -.130815338324550D+01
      vts_coef(487) =           .411347087042135D+00
      vts_coef(488) =           .562651942078483D-01
      vts_coef(489) =           .421598989789659D-01
      vts_coef(490) =           .709743915563441D-01
      vts_coef(491) =           .897720081729636D+00
      vts_coef(492) =           .307499315819176D+00
      vts_coef(493) =           .143719935951829D+00
      vts_coef(494) =          -.375821796749104D-01
      vts_coef(495) =           .329519933048003D-01
      vts_coef(496) =           .314057654924715D+01
      vts_coef(497) =           .684277968174551D+00
      vts_coef(498) =           .182444390414411D+00
      vts_coef(499) =          -.170643684482936D-01
      vts_coef(500) =           .452692064271201D-02
      vts_coef(501) =           .536340899824772D-02
      vts_coef(502) =          -.297740801349461D-01
      vts_coef(503) =          -.141717990232128D-01
      vts_coef(504) =           .744276612141756D-03
      vts_coef(505) =           .216211660136100D-01
      vts_coef(506) =           .253016081868742D+00
      vts_coef(507) =          -.205788234974489D-01
      vts_coef(508) =           .534551355788518D+00
      vts_coef(509) =          -.102232465160955D+00
      vts_coef(510) =           .124039482059658D-01
      vts_coef(511) =          -.506271731123118D-01
      vts_coef(512) =          -.134728398482233D+01
      vts_coef(513) =          -.853816313292838D-01
      vts_coef(514) =          -.199058049068401D+00
      vts_coef(515) =          -.158011228922090D+00
      vts_coef(516) =           .155171519523586D-01
      vts_coef(517) =          -.403499999239535D-01
      vts_coef(518) =          -.604354122929476D+01
      vts_coef(519) =           .154520919093154D+01
      vts_coef(520) =          -.909800423178109D+00
      vts_coef(521) =          -.924402173058308D-02
      vts_coef(522) =           .189198567387995D-01
      vts_coef(523) =          -.139968208643395D+00
      vts_coef(524) =           .459564285706454D-01
      vts_coef(525) =           .336620409815506D-01
      vts_coef(526) =          -.896862362509914D-02
      vts_coef(527) =          -.377007708873109D+01
      vts_coef(528) =           .221699604470947D+01
      vts_coef(529) =          -.607137549610912D+00
      vts_coef(530) =           .490999441926281D-01
      vts_coef(531) =          -.282024065765973D+00
      vts_coef(532) =          -.105784479643000D+00
      vts_coef(533) =           .840188984079502D-01
      vts_coef(534) =          -.700595459640855D-01
      vts_coef(535) =          -.759087862371096D-02
      vts_coef(536) =           .203104315696550D-01
      vts_coef(537) =          -.314511698641439D-02
      vts_coef(538) =          -.162035176397892D-02
      vts_coef(539) =          -.160704303194935D-01
      vts_coef(540) =          -.112309794120046D+00
      vts_coef(541) =          -.189731466924688D-01
      vts_coef(542) =          -.184558803659893D+00
      vts_coef(543) =          -.955188865450939D-01
      vts_coef(544) =          -.233169953764872D-01
      vts_coef(545) =           .357493463531960D-01
      vts_coef(546) =           .284657300562287D+00
      vts_coef(547) =          -.286222409069410D+00
      vts_coef(548) =          -.125626343148688D+00
      vts_coef(549) =           .900871521466028D-01
      vts_coef(550) =           .875857941420052D-01
      vts_coef(551) =           .644494856871849D-01
      vts_coef(552) =           .196781997102561D+01
      vts_coef(553) =          -.126348397488717D+01
      vts_coef(554) =           .206241539308708D+00
      vts_coef(555) =          -.516453085396876D-02
      vts_coef(556) =          -.199941913584919D+00
      vts_coef(557) =           .221069946544982D+00
      vts_coef(558) =           .798572449846910D-01
      vts_coef(559) =           .814571342412546D-01
      vts_coef(560) =           .137050199269369D-01
      vts_coef(561) =           .222520418717798D+01
      vts_coef(562) =          -.992495018338478D+00
      vts_coef(563) =          -.296158119762579D+00
      vts_coef(564) =           .165856919234275D+00
      vts_coef(565) =          -.387054089349473D-01
      vts_coef(566) =           .180265039955004D+00
      vts_coef(567) =          -.275859449506677D+00
      vts_coef(568) =           .228186935281694D+00
      vts_coef(569) =           .170053709522876D-01
      vts_coef(570) =           .768360845437401D-01
      vts_coef(571) =           .908432859908446D-01
      vts_coef(572) =           .630154824806790D-02
      vts_coef(573) =          -.626618048891549D+00
      vts_coef(574) =           .131925951800416D+01
      vts_coef(575) =          -.152707784533978D+01
      vts_coef(576) =           .208910840768494D+00
      vts_coef(577) =           .639366514620213D+00
      vts_coef(578) =          -.110393327322501D+00
      vts_coef(579) =          -.423991941469817D+00
      vts_coef(580) =           .170677344875360D+00
      vts_coef(581) =           .133323183039123D-01
      vts_coef(582) =           .548196841451506D-01
      vts_coef(583) =           .334267795367497D-01
      vts_coef(584) =          -.115976689240370D-01
      vts_coef(585) =          -.165504046310203D-02
      vts_coef(586) =           .219314187335149D-01
      vts_coef(587) =          -.832795360650561D-03
      vts_coef(588) =           .441070909465990D-03

      return
      end


c     * normal frequency
c     * Cartesian coordinate
c     * normal mode (two different normalization)
c==================================================
      subroutine set_normal_mode
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      real*8 Dmass(4,3)
c     ----------------------------------------
      common /mass/ zmc,zmo,zmh,zmf,zmco,zmhf,zmtot
      common /v11geom_pot/ xeq(4,3),bceq(6),zjceq(6)
      common /vtsgeom_pot/ xts(4,3),bcts(6),zjcts(6)
      common /s1geom_pot/  xs1(4,3),bcs1(6),zjcs1(6)
c     ----------------------------------------
      common /s0eq_nmode/ w_s0eq(6),xnm_s0eq(6,4,3),
     $     xnm_s0eq_au(6,4,3)
      common /s0ts_nmode/ w_s0ts(6),xnm_s0ts(6,4,3),
     $     xnm_s0ts_au(6,4,3)
      common /s1eq_nmode/ w_s1eq(6),xnm_s1eq(6,4,3),
     $     xnm_s1eq_au(6,4,3)
c     ----------------------------------------

      do i1 = 1,3
         Dmass(1,i1) = zmc
         Dmass(2,i1) = zmo
         Dmass(3,i1) = zmh
         Dmass(4,i1) = zmf
      enddo

c
c     --------------------
c     S0 equilibrium
c     --------------------
c
      ! normal frequency
      w_s0eq(1) =   663.50d0 * hartree
      w_s0eq(2) =  1036.93d0 * hartree
      w_s0eq(3) =  1071.60d0 * hartree
      w_s0eq(4) =  1404.61d0 * hartree
      w_s0eq(5) =  1846.98d0 * hartree
      w_s0eq(6) =  3145.59d0 * hartree


      ! Cartesian coordinate of EQ on S0.
      ! Normal-mode analysis were
      ! performed at this geometry. <xnm_s0eq> should
      ! be used together with <xts>.
      !--------------------
      call btoc( bceq, xeq )


      ! normal mode ; mass-unweighted and normalized.
      !--------------------
      xnm_s0eq(1,1,1) =      -.43038091d0
      xnm_s0eq(1,1,2) =      -.00000001d0
      xnm_s0eq(1,1,3) =       .10028029d0
      xnm_s0eq(1,2,1) =      -.05216262d0
      xnm_s0eq(1,2,2) =      -.00000001d0
      xnm_s0eq(1,2,3) =      -.56989868d0
      xnm_s0eq(1,3,1) =      -.38780344d0
      xnm_s0eq(1,3,2) =       .00000015d0
      xnm_s0eq(1,3,3) =       .22347460d0
      xnm_s0eq(1,4,1) =       .33616504d0
      xnm_s0eq(1,4,2) =       .00000000d0
      xnm_s0eq(1,4,3) =       .40481027d0
      xnm_s0eq(2,1,1) =       .00000007d0
      xnm_s0eq(2,1,2) =      -.22577111d0
      xnm_s0eq(2,1,3) =      -.00000023d0
      xnm_s0eq(2,2,1) =      -.00000003d0
      xnm_s0eq(2,2,2) =       .06269932d0
      xnm_s0eq(2,2,3) =       .00000003d0
      xnm_s0eq(2,3,1) =       .00000010d0
      xnm_s0eq(2,3,2) =       .97139130d0
      xnm_s0eq(2,3,3) =      -.00000016d0
      xnm_s0eq(2,4,1) =      -.00000003d0
      xnm_s0eq(2,4,2) =       .03866700d0
      xnm_s0eq(2,4,3) =       .00000013d0
      xnm_s0eq(3,1,1) =      -.27477188d0
      xnm_s0eq(3,1,2) =      -.00000042d0
      xnm_s0eq(3,1,3) =       .61513506d0
      xnm_s0eq(3,2,1) =       .12588763d0
      xnm_s0eq(3,2,2) =       .00000012d0
      xnm_s0eq(3,2,3) =      -.02468446d0
      xnm_s0eq(3,3,1) =      -.29485041d0
      xnm_s0eq(3,3,2) =       .00000183d0
      xnm_s0eq(3,3,3) =       .52855267d0
      xnm_s0eq(3,4,1) =       .08304984d0
      xnm_s0eq(3,4,2) =       .00000007d0
      xnm_s0eq(3,4,3) =      -.39553797d0
      xnm_s0eq(4,1,1) =       .01007065d0
      xnm_s0eq(4,1,2) =       .00000000d0
      xnm_s0eq(4,1,3) =      -.05334690d0
      xnm_s0eq(4,2,1) =      -.05391130d0
      xnm_s0eq(4,2,2) =       .00000000d0
      xnm_s0eq(4,2,3) =       .00057288d0
      xnm_s0eq(4,3,1) =       .35873118d0
      xnm_s0eq(4,3,2) =      -.00000001d0
      xnm_s0eq(4,3,3) =       .92994838d0
      xnm_s0eq(4,4,1) =       .02015758d0
      xnm_s0eq(4,4,2) =       .00000000d0
      xnm_s0eq(4,4,3) =      -.01573399d0
      xnm_s0eq(5,1,1) =      -.54174233d0
      xnm_s0eq(5,1,2) =       .00000006d0
      xnm_s0eq(5,1,3) =      -.38309862d0
      xnm_s0eq(5,2,1) =       .41640610d0
      xnm_s0eq(5,2,2) =      -.00000002d0
      xnm_s0eq(5,2,3) =       .25964931d0
      xnm_s0eq(5,3,1) =      -.43885235d0
      xnm_s0eq(5,3,2) =      -.00000025d0
      xnm_s0eq(5,3,3) =       .35512427d0
      xnm_s0eq(5,4,1) =       .01459354d0
      xnm_s0eq(5,4,2) =      -.00000001d0
      xnm_s0eq(5,4,3) =       .00461473d0
      xnm_s0eq(6,1,1) =       .08713525d0
      xnm_s0eq(6,1,2) =       .00000000d0
      xnm_s0eq(6,1,3) =      -.02126037d0
      xnm_s0eq(6,2,1) =      -.00664678d0
      xnm_s0eq(6,2,2) =       .00000000d0
      xnm_s0eq(6,2,3) =      -.00398493d0
      xnm_s0eq(6,3,1) =      -.93961588d0
      xnm_s0eq(6,3,2) =       .00000001d0
      xnm_s0eq(6,3,3) =       .33017727d0
      xnm_s0eq(6,4,1) =       .00001815d0
      xnm_s0eq(6,4,2) =       .00000000d0
      xnm_s0eq(6,4,3) =      -.00059440d0

      ! normal mode: This can be used directly with
      ! normal coordinates {Q} in atomic unit.
      !--------------------
      do imode = 1,6

         Dnorm = zero
         do i1 = 1,4
            do i2 = 1,3
               Dnorm = Dnorm + Dmass(i1,i2)
     $              * xnm_s0eq( imode,i1,i2 )**2
            enddo
         enddo
         Dnorm = dsqrt(Dnorm)
         do i1 = 1,4
            do i2 = 1,3
               xnm_s0eq_au( imode,i1,i2 ) = 
     $              xnm_s0eq( imode,i1,i2 ) / Dnorm
            enddo
         enddo

      enddo

c
c     --------------------
c     S0 transition state
c     --------------------
c
      ! normal frequency
      w_s0ts(1) =  -1528.23d0 * hartree
      w_s0ts(2) =    346.30d0 * hartree
      w_s0ts(3) =    805.23d0 * hartree
      w_s0ts(4) =    834.86d0 * hartree
      w_s0ts(5) =   1972.70d0 * hartree
      w_s0ts(6) =   2799.26d0 * hartree


      ! Cartesian coordinate of TS (a.u.)
      ! FCM calculation and normal-mode analysis were
      ! performed at this geometry. <xnm_s0ts> should
      ! be used together with <xts>.
      !--------------------
      xts(1,1) =    .0210217441d0
      xts(1,2) =    .0000000000d0     
      xts(1,3) =    .0221637224d0

      xts(2,1) =   -.0150921703d0
      xts(2,2) =    .0000000000d0 
      xts(2,3) =  -2.1323770544d0

      xts(3,1) =   -.2770156181d0
      xts(3,2) =    .0000000000d0 
      xts(3,3) =   2.1289159150d0

      xts(4,1) =  -2.8934790390d0
      xts(4,2) =    .0000000000d0 
      xts(4,3) =   1.9281604367d0


      ! normal mode: mass-unweighted and normalized.
      !--------------------
      xnm_s0ts(1,1,1) =       .11555157d0
      xnm_s0ts(1,1,2) =      -.00000017d0
      xnm_s0ts(1,1,3) =      -.04875082d0
      xnm_s0ts(1,2,1) =      -.01015227d0
      xnm_s0ts(1,2,2) =       .00000036d0
      xnm_s0ts(1,2,3) =      -.01660497d0
      xnm_s0ts(1,3,1) =      -.98492120d0
      xnm_s0ts(1,3,2) =      -.00000474d0
      xnm_s0ts(1,3,3) =       .11020617d0
      xnm_s0ts(1,4,1) =      -.01258916d0
      xnm_s0ts(1,4,2) =       .00000028d0
      xnm_s0ts(1,4,3) =       .03897111d0
      xnm_s0ts(2,1,1) =       .31922627d0
      xnm_s0ts(2,1,2) =       .00000969d0
      xnm_s0ts(2,1,3) =       .30677628d0
      xnm_s0ts(2,2,1) =      -.51933705d0
      xnm_s0ts(2,2,2) =      -.00000759d0
      xnm_s0ts(2,2,3) =       .32146502d0
      xnm_s0ts(2,3,1) =      -.11607523d0
      xnm_s0ts(2,3,2) =       .00000870d0
      xnm_s0ts(2,3,3) =       .35445281d0
      xnm_s0ts(2,4,1) =       .24185269d0
      xnm_s0ts(2,4,2) =      -.00000433d0
      xnm_s0ts(2,4,3) =      -.48304160d0
      xnm_s0ts(3,1,1) =      -.46227577d0
      xnm_s0ts(3,1,2) =       .00000500d0
      xnm_s0ts(3,1,3) =       .05545649d0
      xnm_s0ts(3,2,1) =       .09150805d0
      xnm_s0ts(3,2,2) =      -.00000191d0
      xnm_s0ts(3,2,3) =       .04362505d0
      xnm_s0ts(3,3,1) =      -.81285639d0
      xnm_s0ts(3,3,2) =      -.00002954d0
      xnm_s0ts(3,3,3) =       .19760398d0
      xnm_s0ts(3,4,1) =       .25771222d0
      xnm_s0ts(3,4,2) =       .00000034d0
      xnm_s0ts(3,4,3) =      -.08216076d0
      xnm_s0ts(4,1,1) =      -.00000308d0
      xnm_s0ts(4,1,2) =      -.14519589d0
      xnm_s0ts(4,1,3) =      -.00000001d0
      xnm_s0ts(4,2,1) =       .00000008d0
      xnm_s0ts(4,2,2) =       .05416968d0
      xnm_s0ts(4,2,3) =       .00000097d0
      xnm_s0ts(4,3,1) =      -.00000976d0
      xnm_s0ts(4,3,2) =       .98790132d0
      xnm_s0ts(4,3,3) =       .00000277d0
      xnm_s0ts(4,4,1) =       .00000240d0
      xnm_s0ts(4,4,2) =      -.00589730d0
      xnm_s0ts(4,4,3) =      -.00000095d0
      xnm_s0ts(5,1,1) =      -.02152836d0
      xnm_s0ts(5,1,2) =       .00000011d0
      xnm_s0ts(5,1,3) =      -.36427534d0
      xnm_s0ts(5,2,1) =       .00444570d0
      xnm_s0ts(5,2,2) =       .00000010d0
      xnm_s0ts(5,2,3) =       .33045251d0
      xnm_s0ts(5,3,1) =      -.00546543d0
      xnm_s0ts(5,3,2) =      -.00000407d0
      xnm_s0ts(5,3,3) =      -.87033503d0
      xnm_s0ts(5,4,1) =       .01013269d0
      xnm_s0ts(5,4,2) =       .00000005d0
      xnm_s0ts(5,4,3) =      -.00240322d0
      xnm_s0ts(6,1,1) =       .00592821d0
      xnm_s0ts(6,1,2) =      -.00000004d0
      xnm_s0ts(6,1,3) =       .13753716d0
      xnm_s0ts(6,2,1) =      -.00415499d0
      xnm_s0ts(6,2,2) =      -.00000008d0
      xnm_s0ts(6,2,3) =      -.04687659d0
      xnm_s0ts(6,3,1) =      -.16824498d0
      xnm_s0ts(6,3,2) =       .00000171d0
      xnm_s0ts(6,3,3) =      -.97490399d0
      xnm_s0ts(6,4,1) =       .00861192d0
      xnm_s0ts(6,4,2) =      -.00000009d0
      xnm_s0ts(6,4,3) =       .00392148d0


      ! normal mode: This can be used directly with
      ! normal coordinates {Q} in atomic unit.
      !--------------------
      do imode = 1,6

         Dnorm = zero
         do i1 = 1,4
            do i2 = 1,3
               Dnorm = Dnorm + Dmass(i1,i2)
     $              * xnm_s0ts( imode,i1,i2 )**2
            enddo
         enddo
         Dnorm = dsqrt(Dnorm)
         do i1 = 1,4
            do i2 = 1,3
               xnm_s0ts_au( imode,i1,i2 ) = 
     $              xnm_s0ts( imode,i1,i2 ) / Dnorm
            enddo
         enddo

      enddo


c
c     --------------------
c     S1 equilibrium
c     --------------------

      ! normal frequency
      w_s1eq( 1) =   514.44d0 * hartree
      w_s1eq( 2) =  1102.75d0 * hartree
      w_s1eq( 3) =  1158.87d0 * hartree
      w_s1eq( 4) =  1184.87d0 * hartree
      w_s1eq( 5) =  1423.77d0 * hartree
      w_s1eq( 6) =  3282.11d0 * hartree


      ! Cartesian coordinate of EQ:S1.
      ! FCM calculation and normal-mode analysis were
      ! performed at this geometry. <xnm_s1eq> should
      ! be used together with <xs1>.
      !------------------------------
      xs1(1,1) =    .2229542731d0
      xs1(1,2) =   -.4250403107d0
      xs1(1,3) =    .0888733270d0

      xs1(2,1) =  -1.9677217163d0
      xs1(2,2) =    .1344387329d0
      xs1(2,3) =  -1.2271797232d0

      xs1(3,1) =   1.9116094145d0
      xs1(3,2) =    .3571100404d0
      xs1(3,3) =   -.7350517681d0

      xs1(4,1) =   -.1053939356d0
      xs1(4,2) =    .1334915374d0
      xs1(4,3) =   2.5341111838d0


      ! normal mode: mass unweighted and normalized
      !--------------------
      xnm_s1eq(1,1,1) =      -.36304550d0
      xnm_s1eq(1,1,2) =       .15504910d0
      xnm_s1eq(1,1,3) =       .18912515d0
      xnm_s1eq(1,2,1) =      -.12059195d0
      xnm_s1eq(1,2,2) =      -.05639273d0
      xnm_s1eq(1,2,3) =      -.60820968d0
      xnm_s1eq(1,3,1) =      -.36349346d0
      xnm_s1eq(1,3,2) =       .02919630d0
      xnm_s1eq(1,3,3) =       .11415626d0
      xnm_s1eq(1,4,1) =       .34997424d0
      xnm_s1eq(1,4,2) =      -.05197378d0
      xnm_s1eq(1,4,3) =       .38672088d0
      xnm_s1eq(2,1,1) =      -.22042801d0
      xnm_s1eq(2,1,2) =      -.07288350d0
      xnm_s1eq(2,1,3) =      -.00429780d0
      xnm_s1eq(2,2,1) =       .17128092d0
      xnm_s1eq(2,2,2) =      -.00483622d0
      xnm_s1eq(2,2,3) =       .05846537d0
      xnm_s1eq(2,3,1) =      -.60254398d0
      xnm_s1eq(2,3,2) =       .73788502d0
      xnm_s1eq(2,3,3) =      -.05501002d0
      xnm_s1eq(2,4,1) =       .02669397d0
      xnm_s1eq(2,4,2) =       .01126823d0
      xnm_s1eq(2,4,3) =      -.04362434d0
      xnm_s1eq(3,1,1) =       .18995742d0
      xnm_s1eq(3,1,2) =      -.25022229d0
      xnm_s1eq(3,1,3) =      -.02265129d0
      xnm_s1eq(3,2,1) =      -.12532150d0
      xnm_s1eq(3,2,2) =       .07067341d0
      xnm_s1eq(3,2,3) =      -.07092955d0
      xnm_s1eq(3,3,1) =      -.17072595d0
      xnm_s1eq(3,3,2) =       .88177807d0
      xnm_s1eq(3,3,3) =       .24884099d0
      xnm_s1eq(3,4,1) =      -.00545363d0
      xnm_s1eq(3,4,2) =       .05211131d0
      xnm_s1eq(3,4,3) =       .06093933d0
      xnm_s1eq(4,1,1) =       .00595776d0
      xnm_s1eq(4,1,2) =       .03434730d0
      xnm_s1eq(4,1,3) =       .46575209d0
      xnm_s1eq(4,2,1) =      -.08502889d0
      xnm_s1eq(4,2,2) =       .02025763d0
      xnm_s1eq(4,2,3) =      -.04278376d0
      xnm_s1eq(4,3,1) =       .07382045d0
      xnm_s1eq(4,3,2) =       .21991775d0
      xnm_s1eq(4,3,3) =       .78877464d0
      xnm_s1eq(4,4,1) =       .06395519d0
      xnm_s1eq(4,4,2) =      -.05032671d0
      xnm_s1eq(4,4,3) =      -.29964524d0
      xnm_s1eq(5,1,1) =       .07855718d0
      xnm_s1eq(5,1,2) =       .00013797d0
      xnm_s1eq(5,1,3) =       .12157393d0
      xnm_s1eq(5,2,1) =      -.00654111d0
      xnm_s1eq(5,2,2) =       .01853761d0
      xnm_s1eq(5,2,3) =      -.02430197d0
      xnm_s1eq(5,3,1) =      -.39354373d0
      xnm_s1eq(5,3,2) =      -.03113387d0
      xnm_s1eq(5,3,3) =      -.90631469d0
      xnm_s1eq(5,4,1) =      -.02339393d0
      xnm_s1eq(5,4,2) =      -.01405913d0
      xnm_s1eq(5,4,3) =      -.00861795d0
      xnm_s1eq(6,1,1) =       .07050519d0
      xnm_s1eq(6,1,2) =       .02869892d0
      xnm_s1eq(6,1,3) =      -.03321795d0
      xnm_s1eq(6,2,1) =       .00010543d0
      xnm_s1eq(6,2,2) =       .00091275d0
      xnm_s1eq(6,2,3) =       .00029737d0
      xnm_s1eq(6,3,1) =      -.83641784d0
      xnm_s1eq(6,3,2) =      -.36959333d0
      xnm_s1eq(6,3,3) =       .39611401d0
      xnm_s1eq(6,4,1) =      -.00059639d0
      xnm_s1eq(6,4,2) =       .00055801d0
      xnm_s1eq(6,4,3) =      -.00011876d0

      ! normal mode: This can be used directly with
      ! normal coordinates {Q} in atomic unit.
      !--------------------
      do imode = 1,6

         Dnorm = zero
         do i1 = 1,4
            do i2 = 1,3
               Dnorm = Dnorm + Dmass(i1,i2)
     $              * xnm_s1eq( imode,i1,i2 )**2
            enddo
         enddo
         Dnorm = dsqrt(Dnorm)
         do i1 = 1,4
            do i2 = 1,3
               xnm_s1eq_au( imode,i1,i2 ) = 
     $              xnm_s1eq( imode,i1,i2 ) / Dnorm
            enddo
         enddo

      enddo
      
      return
      end


c     get bond-angle in radian.
c     * use (func.)dacos2
c============================================================
      real*8 function bond_angle(xdum,iatom1,iatom2,iatom3)
c============================================================
      implicit real*8(a-h,o-z)
      parameter (zero=0.0d0)
      Integer stat
      real*8 xdum(4,3)
      real*8 x1(3),x2(3),x3(3),u12(3),u23(3)

      do i1 = 1,3
         x1(i1) = xdum(iatom1,i1)
         x2(i1) = xdum(iatom2,i1)
         x3(i1) = xdum(iatom3,i1)
      enddo

      r12 = zero
      r23 = zero
      do i1 = 1,3
         u12(i1) = x1(i1) - x2(i1)
         u23(i1) = x2(i1) - x3(i1)
         r12 = r12 + u12(i1)*u12(i1)
         r23 = r23 + u23(i1)*u23(i1)
      enddo

      r12 = dsqrt(r12)
      r23 = dsqrt(r23)

      do i1 = 1,3
         u12(i1) = u12(i1)/r12
         u23(i1) = u23(i1)/r23
      enddo

      dum = zero
      do i1 = 1,3
         dum = dum + u12(i1)*u23(i1)
      enddo
      stat=2987
      bond_angle = dacos2(dum,stat)
      Return
      end


c     get bond-angle in cosine.
c============================================================
      real*8 function bond_angle_cos(xdum,iatom1,iatom2,iatom3)
c============================================================
      implicit real*8(a-h,o-z)
      parameter (zero=0.0d0)
      real*8 xdum(4,3)
      real*8 x1(3),x2(3),x3(3),u12(3),u23(3)

      do i1 = 1,3
         x1(i1) = xdum(iatom1,i1)
         x2(i1) = xdum(iatom2,i1)
         x3(i1) = xdum(iatom3,i1)
      enddo

      r12 = zero
      r23 = zero
      do i1 = 1,3
         u12(i1) = x1(i1) - x2(i1)
         u23(i1) = x2(i1) - x3(i1)
         r12 = r12 + u12(i1)*u12(i1)
         r23 = r23 + u23(i1)*u23(i1)
      enddo

      r12 = dsqrt(r12)
      r23 = dsqrt(r23)

      do i1 = 1,3
         u12(i1) = u12(i1)/r12
         u23(i1) = u23(i1)/r23
      enddo

      dum = zero
      do i1 = 1,3
         dum = dum + u12(i1)*u23(i1)
      enddo
      bond_angle_cos = dum
      end



c
c     Jacobi(radian) to Cartesian
c     * set /mass/ before calling this routine.
c==================================================
      subroutine jtoc(zjc,xdum)
c==================================================
c============================================================
c
c     -- pot.inc -- (970924)
c     include file for HFCO potential program
c
c============================================================
      implicit real*8(a-h,o-z)
c     ----------------------------------------
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      parameter (bohr=0.52917706d0)
      real*8 kcal
      parameter (kcal = 627.5075d0)
      parameter (amu = 1836.1516d0)
      parameter (hartree = 4.556037d-6)
c     ----------------------------------------
      parameter (Ncoefmax = 600)
c     ----------------------------------------
      ! order of polynomial expansion (v11,vts)
      parameter (Nordermax = 8)
c     ----------------------------------------
c      character header*50
c      common /header_com/ header
      dimension xdum(4,3),zjc(6),xj(3,3)
      common /mass/ zmc,zmo,zmh,zmf,zmco,zmhf,zmtot

c     jacobi coordinates
      r1   = zjc(1)
      r2   = zjc(2)
      r3   = zjc(3)
      th1  = zjc(4)
      th2  = zjc(5)
      phai = zjc(6)

      s1 = dsin(th1)
      c1 = dcos(th1)
      s2 = dsin(th2)
      c2 = dcos(th2)
      sp = dsin(phai)
      cp = dcos(phai)

      xj(1,1) = r1 * s1 
      xj(1,2) = zero
      xj(1,3) = r1 * c1

      xj(2,1) = zero
      xj(2,2) = zero
      xj(2,3) = r2

      xj(3,1) = r3 * s2 * cp
      xj(3,2) = r3 * s2 * sp
      xj(3,3) = r3 * c2

      do i1 = 1,3
         xdum(1,i1) = - (zmhf/zmtot) * xj(2,i1)
     |                + (zmo/zmco) * xj(1,i1)
         xdum(2,i1) = - (zmhf/zmtot) * xj(2,i1)
     |                - (zmc/zmco) * xj(1,i1)
         xdum(3,i1) =   (zmco/zmtot) * xj(2,i1)
     |                + (zmf/zmhf) * xj(3,i1)
         xdum(4,i1) =   (zmco/zmtot) * xj(2,i1)
     |                - (zmh/zmhf) * xj(3,i1)
      enddo

      return
      end


c
c     Cartesian to Jacobi(radian)
c     * set /mass/ before calling this routine.
c==================================================
      subroutine ctoj(x,zjc)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0 )
      dimension x(4,3),xj(3,3),zjc(6)
      common /mass/ zmc,zmo,zmh,zmf,zmco,zmhf,zmtot
      pi = dacos(-1.0d0)

c     --------------------
      epsa = 1.0d-8
c     --------------------


      do i1 = 1,3
         xj(1,i1) = x(1,i1) - x(2,i1)
         xj(2,i1) =   ( zmh*x(3,i1) + zmf*x(4,i1) ) / zmhf
     |              - ( zmc*x(1,i1) + zmo*x(2,i1) ) / zmco
         xj(3,i1) = x(3,i1) - x(4,i1)
      enddo

c     r1,r2,r3
      r1 = xj(1,1)*xj(1,1) + xj(1,2)*xj(1,2) + xj(1,3)*xj(1,3)
      r2 = xj(2,1)*xj(2,1) + xj(2,2)*xj(2,2) + xj(2,3)*xj(2,3)
      r3 = xj(3,1)*xj(3,1) + xj(3,2)*xj(3,2) + xj(3,3)*xj(3,3)
      r1 = dsqrt( r1 )
      r2 = dsqrt( r2 )
      r3 = dsqrt( r3 )

c     inner product
      p12 = xj(1,1)*xj(2,1) + xj(1,2)*xj(2,2) + xj(1,3)*xj(2,3)
      p23 = xj(2,1)*xj(3,1) + xj(2,2)*xj(3,2) + xj(2,3)*xj(3,3)
      p13 = xj(1,1)*xj(3,1) + xj(1,2)*xj(3,2) + xj(1,3)*xj(3,3)

c     th1,th2
      th1 = dacos2( p12 / (r1*r2) ,3114)
      th2 = dacos2( p23 / (r2*r3) ,3115)

      dum1 = pi - th1
      dum2 = pi - th2

      if( (th1.lt.epsa).or.(dum1.lt.epsa)
     |     .or.(th2.lt.epsa).or.(dum2.lt.epsa)) then
         phai = zero
         goto 1000
      endif


c     phai
*     ablolute value
      dum1 = r2*r2*p13 - p23*p12
      dum2 = r1*r1*r2*r2 - p12*p12
      dum3 = r2*r2*r3*r3 - p23*p23
      dum4 = dsqrt( dum2*dum3 )
      phai_abs = dacos2( dum1 / dum4 ,3133)

*     sign
      dum =   xj(2,1) * ( xj(1,2)*xj(3,3) - xj(1,3)*xj(3,2) )
     |      + xj(2,2) * ( xj(1,3)*xj(3,1) - xj(1,1)*xj(3,3) )
     |      + xj(2,3) * ( xj(1,1)*xj(3,2) - xj(1,2)*xj(3,1) )
      if(dum.ge.zero) then
         phai = phai_abs
      else
         phai = -phai_abs
      endif


 1000 continue
      zjc(1) = r1
      zjc(2) = r2
      zjc(3) = r3
      zjc(4) = th1
      zjc(5) = th2
      zjc(6) = phai

      return
      end


c
c     Cartesian to site-site distance
c==================================================
      subroutine ctos(x,sd)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0 )
      dimension x(4,3),sd(6)

      do i = 1,6
         sd(i) = zero
      enddo
c     
c     sd(1) ... rco
c       (2) ... rch
c       (3) ... rcf
c       (4) ... roh
c       (5) ... rof
c       (5) ... rhf
c
      do i=1,3
         sd(1) = sd(1) + (x(1,i)-x(2,i))**2
         sd(2) = sd(2) + (x(1,i)-x(3,i))**2
         sd(3) = sd(3) + (x(1,i)-x(4,i))**2
         sd(4) = sd(4) + (x(2,i)-x(3,i))**2
         sd(5) = sd(5) + (x(2,i)-x(4,i))**2
         sd(6) = sd(6) + (x(3,i)-x(4,i))**2
      enddo
      do i = 1,6
         sd(i) = dsqrt( sd(i) )
      enddo
      return
      end


c
c     Cartesian to bond(radian,pole-CF)
c     * tau [-pi:pi]
c==================================================
      subroutine ctob(xdum,bc)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0 )
      Integer stat
      dimension xdum(4,3),bc(6)
      dimension dx(4,3)
      Logical Init
      Common /Initialize/ xJaco(6),Init
      pi = dacos(-1.0d0)

c     --------------------
      epsa = 1.0d-8
c     --------------------

c
c     rco,rch,rcf
c
      do i1 = 2,4
         do i2 = 1,3
            dx(i1,i2) = xdum(i1,i2) - xdum(1,i2)
         enddo
      enddo

      r2 = dx(2,1)*dx(2,1) + dx(2,2)*dx(2,2) + dx(2,3)*dx(2,3)
      r3 = dx(3,1)*dx(3,1) + dx(3,2)*dx(3,2) + dx(3,3)*dx(3,3)
      r4 = dx(4,1)*dx(4,1) + dx(4,2)*dx(4,2) + dx(4,3)*dx(4,3)
      r2 = dsqrt( r2 )
      r3 = dsqrt( r3 )
      r4 = dsqrt( r4 )

c
c     th1,th2
c

*     inner product
      p23 = dx(2,1)*dx(3,1) + dx(2,2)*dx(3,2) + dx(2,3)*dx(3,3)
      p34 = dx(3,1)*dx(4,1) + dx(3,2)*dx(4,2) + dx(3,3)*dx(4,3)
      p24 = dx(2,1)*dx(4,1) + dx(2,2)*dx(4,2) + dx(2,3)*dx(4,3)

      stat=3245
      th1 = dacos2( p34 / (r3*r4) ,stat)
      if(stat.eq.0) goto 999
      th2 = dacos2( p24 / (r2*r4) ,3235)

      dum1 = pi - th1
      dum2 = pi - th2
      if( (th1.lt.epsa).or.(dum1.lt.epsa)
     |     .or.(th2.lt.epsa).or.(dum2.lt.epsa)) then
         tau = zero
         goto 1000
      endif

c
c     tau
c

*     absolute value
      dum1 = p24*p34 - r4*r4*p23
      dum2 = r3*r3*r4*r4 - p34*p34
      dum3 = r2*r2*r4*r4 - p24*p24
      dum4 = dsqrt( dum2*dum3 )
      
      tau_abs = dacos2( dum1 / dum4 ,3255)

*     sign
      dum =   dx(4,1) * (dx(3,2)*dx(2,3) - dx(3,3)*dx(2,2)) 
     $      + dx(4,2) * (dx(3,3)*dx(2,1) - dx(3,1)*dx(2,3)) 
     $      + dx(4,3) * (dx(3,1)*dx(2,2) - dx(3,2)*dx(2,1)) 
      
      if(dum.ge.zero) then
         tau =  tau_abs
      else
         tau = -tau_abs
      endif

 1000 continue
      bc(1) = r2
      bc(2) = r3
      bc(3) = r4
      bc(4) = th1
      bc(5) = th2
      bc(6) = tau

      return
  999 Write (6,1120) xJaco
 1120 format('Jacobi :',f10.3/(8x,f10.3))
      Write (6,1100) ((xdum(i1,i2),i2=1,3),i1=1,4)
 1100 format('xyz-coord :',3f10.3/(11x,3f10.3))
      Write (6,1110) r2,r3,r4
 1110 format('r[i],i=2..4 :',3f10.3)
      Stop
      end

c
c     Cartesian to bond(radian,pole-CF)
c     * tau [-pi:pi]
c==================================================
      subroutine ctob_old(x,bc)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0 )
      dimension x(4,3),bc(6)
      dimension r1(3),r2(3),r3(3),zn1(3),zn2(3)
      pi = dacos(-1.0d0)

c     ----------------------------------------
      eps = 1.0d-8
c     ----------------------------------------
c
c     get rco,rch,rcf
c     * bond vector R1,R2,R3
c
      do i=1,3
         r1(i) = x(2,i)-x(1,i) ! co
         r2(i) = x(3,i)-x(1,i) ! ch
         r3(i) = x(4,i)-x(1,i) ! cf
      enddo

      rco = zero
      rch = zero
      rcf = zero
      do i=1,3
         rco = rco + r1(i)**2
         rch = rch + r2(i)**2
         rcf = rcf + r3(i)**2
      enddo
      rco = dsqrt(rco)
      rch = dsqrt(rch)
      rcf = dsqrt(rcf)

*     length of vector { rco / rch / rcf } is zero.
      if((rco.lt.eps).or.(rch.lt.eps).or.(rcf.lt.eps)) then
         ierr = 1
         goto 30000
      endif
c
c     get s1,s2
c     * S1 = dacos(R2*R3/|R2*R3|)
c     * S2 = dacos(R1*R3/|R1*R3|)
c
      dum1 = zero
      dum2 = zero
      do i = 1,3
         dum1 = dum1 + r2(i)*r3(i)
         dum2 = dum2 + r1(i)*r3(i)
      enddo
      dum1 = dum1/(rch*rcf)
      dum2 = dum2/(rco*rcf)
      s1 = dacos2(dum1,3334)
      s2 = dacos2(dum2,3335)
c
c     get tau
c     * normal vector N1,N2
c
      zn1(1) = r3(2)*r2(3) - r3(3)*r2(2)
      zn1(2) = r3(3)*r2(1) - r3(1)*r2(3)
      zn1(3) = r3(1)*r2(2) - r3(2)*r2(1) 

      zn2(1) = r1(2)*r3(3) - r1(3)*r3(2)
      zn2(2) = r1(3)*r3(1) - r1(1)*r3(3)
      zn2(3) = r1(1)*r3(2) - r1(2)*r3(1) 

      rzn1 = zero
      rzn2 = zero
      do i = 1,3
         rzn1 = rzn1 + zn1(i)**2
         rzn2 = rzn2 + zn2(i)**2
      enddo
      rzn1 = dsqrt(rzn1)
      rzn2 = dsqrt(rzn2)

c     [(FCH) or (FCO) is linear : tau = 0 ]
      if((rzn1.lt.eps).or.(rzn2.lt.eps)) then
         tau = zero
         goto 10000
      endif

      do i = 1,3
         zn1(i) = zn1(i)/rzn1
         zn2(i) = zn2(i)/rzn2
      enddo
c
c     tau0 = dacos ( N1 * N2 )
c
      dum = zero
      do i = 1,3
         dum = dum + zn1(i)*zn2(i)
      enddo
      tau0 = dacos2(dum,3374)
c
c     tau <- tau0
c
c     [planar]
      dum1 = dabs(tau0)
      dum2 = dabs(tau0-pi)
      if((dum1.lt.eps).or.(dum2.lt.eps)) then
         tau = tau0
         goto 10000
      endif

c     [non-planar]
      dum = zero
      dum = dum + r3(1)*( zn1(2)*zn2(3) - zn1(3)*zn2(2) )
      dum = dum + r3(2)*( zn1(3)*zn2(1) - zn1(1)*zn2(3) )
      dum = dum + r3(3)*( zn1(1)*zn2(2) - zn1(2)*zn2(1) )
      if(dum.lt.-eps) tau =  tau0
      if(dum.gt. eps) tau = -tau0

10000 continue
      bc(1) = rco
      bc(2) = rch
      bc(3) = rcf
      bc(4) = s1
      bc(5) = s2
      bc(6) = tau
      return

c     ----------------------------------------
c     singular structure
c     ----------------------------------------
30000 continue
      write(6,*) 'ctob: singular structure appeared.'
      if(ierr.eq.1) then
         write(6,*) 'distance between CO or CH or CF is nearly zero.'
      endif
      do i1 = 1,4
         do i2 = 1,3
            write(6,*) 'x(',i1,',',i2,')=',x(i1,i2)
         enddo
      enddo
      stop
      end



c
c     bond(pole-CF,radian) to Cartesian
c==================================================
      subroutine btoc(bc,x)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0 )
      dimension bc(6),x(4,3)

      rco = bc(1)
      rch = bc(2)
      rcf = bc(3)
      s1  = bc(4)
      s2  = bc(5)
      tau = bc(6)
      
c     C
      x(1,1) = zero
      x(1,2) = zero
      x(1,3) = zero

c     O
      x(2,1) = -rco*dsin(s2)
      x(2,2) = zero
      x(2,3) =  rco*dcos(s2)

c     H
      x(3,1) = rch * dsin(s1) * dcos(tau)
      x(3,2) = rch * dsin(s1) * dsin(tau)
      x(3,3) = rch * dcos(s1)

c     F
      x(4,1) = zero
      x(4,2) = zero
      x(4,3) = rcf
      
      return
      end


c     dacos with boundary check
c==================================================
      real*8 function dacos2(x,stat)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( zero=0.0d0, one=1.0d0 )
      parameter ( eps=1.0d-4 )
      Integer stat
      pi = dacos(-1.0d0)

      if((x.le.one).and.(x.ge.-one)) then
         dacos2 = dacos(x)
         return
      endif

      if((x.gt.one).and.(x.le.(one+eps))) then
         dacos2 = zero
         return
      endif

      if((x.lt.-one).and.(x.ge.-(one+eps))) then
         dacos2 = pi
         return
      endif

c     error case
      write(6,'("dacos2[",i4,"] : argument is out of [-1:1]")') stat
      write(6,'("        eps = 1.0d-4                 ")')
      write(6,'("        argument = ",f40.20)') x
c     stat=0
      Return
      end

c
c     change unit
c     * mode = 1  ... a.u.  -> angs.
c     * mode = -1 ... angs. -> a.u.
c==================================================
      subroutine unit_change(q1,q2,mode)
c==================================================
      implicit real*8(a-h,o-z)
      parameter ( bohr=0.52917706d0 )
      real*8 q1(6),q2(6),qdum(6)
      pi = dacos(-1.0d0)

      do i1 = 1,6
         qdum(i1) = q1(i1)
      enddo

      if(mode.eq.1) then
         qdum(1) = qdum(1) * bohr
         qdum(2) = qdum(2) * bohr
         qdum(3) = qdum(3) * bohr
         qdum(4) = qdum(4) * 180.0d0/pi
         qdum(5) = qdum(5) * 180.0d0/pi
         qdum(6) = qdum(6) * 180.0d0/pi
      else
         qdum(1) = qdum(1) / bohr
         qdum(2) = qdum(2) / bohr
         qdum(3) = qdum(3) / bohr
         qdum(4) = qdum(4) * pi/180.0d0
         qdum(5) = qdum(5) * pi/180.0d0
         qdum(6) = qdum(6) * pi/180.0d0
      endif

      do i1 = 1,6
         q2(i1) = qdum(i1)
      enddo

      return
      end


c==================================================
      real*8 function escale(v)
c==================================================
      implicit real*8(a-h,o-z)
      common /energy_pot/ v110,vts0,v220
      escale = (v - vts0) * 627.5075d0
      end









C--------------------------------------------------------------------------
C Place here your own potential energy routine (and remove the lines above).
C E.g.:
C     subroutine mysrf(r1,r2,theta,phi,v)
C        .... FORTRAN TEXT ....
C     return
C     end
C
C NOTE: You also have to modify the file $MCTDH_DIR/source/opfuncs/funcsrf.F
C       Subroutine defsrf, line 497-507    # line numbers change when a
C       Subroutine vpoint, line 1517-1521  # new surface is added.
C                                          # search for 'mysurf'.
C
C       You may name the coordinates as you like, but make sure that their 
C       ordering is appropriate. The variable v contains the energy value 
C       on return, i.e. v = V(r1,r2,theta,phi) 
C--------------------------------------------------------------------------

