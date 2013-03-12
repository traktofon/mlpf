!=======================================================================
module units_m
!=======================================================================

   use base_m
   use tokenize_m
   use map_str2dbl_m
   implicit none
   private

   public :: &
      parse_value, parse_energy, parse_mass, parse_length, parse_momentum, parse_time, &
      init_units


   type :: unitclass_t
      type(map_str2dbl_t) :: map
      character(len=c1)   :: desc
   end type unitclass_t

   type(unitclass_t),save :: &
      units_energy, units_mass, units_length, units_momentum, units_time, units

   real(dbl),parameter :: ufs = 41.34137333656d0
   real(dbl),parameter :: uev = 27.21138386d0
   real(dbl),parameter :: invcm = 2.1947463137d+5
   real(dbl),parameter :: kcal = 6.27503d+2
   real(dbl),parameter :: kjoule = 2.6255d+3
   real(dbl),parameter :: kelvin = 3.15777d+5
   real(dbl),parameter :: attojoule = 1.602177d-1*uev
   real(dbl),parameter :: atomicmassunit = 1822.88848325d0
   real(dbl),parameter :: protonmass = 1836.15267247d0
   real(dbl),parameter :: hydrogenmass = 1837.15264409d0
   real(dbl),parameter :: deuteriummass = 3671.482934845d0
   real(dbl),parameter :: angstroem = 0.52917720859d0
   real(dbl),parameter :: picometer = 0.52917720859d+2
   real(dbl),parameter :: nanometer = 0.52917720859d-1
   real(dbl),parameter :: debye = 0.39343d0

   contains

!=======================================================================

   subroutine init_units
      call init_units_energy
      call init_units_mass
      call init_units_length
      call init_units_momentum
      call init_units_time
      units%desc = "any"
      call add_units_from(units_energy)
      call add_units_from(units_mass)
      call add_units_from(units_length)
      call add_units_from(units_momentum)
      call add_units_from(units_time)
      contains
      subroutine add_units_from(u)
         type(unitclass_t),intent(inout) :: u
         character(len=c1)               :: key
         real(dbl)                       :: val
         call map_iter_start(u%map)
         do while (map_iter_next(u%map,key,val))
            call map_put(units%map,key,val)
         enddo
      end subroutine add_units_from
   end subroutine init_units


   subroutine init_units_energy
      units_energy%desc = "energy"
      call put("au",       1.d0)
      call put("mH",       1.d-3)
      call put("eV",       1.d0/uev)
      call put("meV",      1.d-3/uev)
      call put("eV-1",     uev)
      call put("cm-1",     1.d0/invcm)
      call put("kcal/mol", 1.d0/kcal)
      call put("kJ/mol",   1.d0/kjoule)
      call put("kelvin",   1.d0/kelvin)
      call put("aJ",       1.d0/attojoule)
      call put("debye",    debye)
      contains
      subroutine put(label,cval)
         character(len=*),intent(in) :: label
         real(dbl),intent(in)        :: cval
         call map_put(units_energy%map, label, cval)
      end subroutine put
   end subroutine init_units_energy


   subroutine init_units_mass
      units_mass%desc = "mass"
      call put("au",     1.d0)
      call put("AMU",    atomicmassunit)
      call put("H-mass", hydrogenmass)
      call put("D-mass", deuteriummass)
      call put("p-mass", protonmass)
      contains
      subroutine put(label,cval)
         character(len=*),intent(in) :: label
         real(dbl),intent(in)        :: cval
         call map_put(units_mass%map, label, cval)
      end subroutine put
   end subroutine init_units_mass


   subroutine init_units_length
      units_length%desc = "length"
      call put("au",    1.d0)
      call put("Angst", 1.d0/angstroem)
      call put("pm",    1.d0/picometer)
      call put("nm",    1.d0/nanometer)
      call put("PI",    PI)
      call put("deg",   PI/180.d0)
      contains
      subroutine put(label,cval)
         character(len=*),intent(in) :: label
         real(dbl),intent(in)        :: cval
         call map_put(units_length%map, label, cval)
      end subroutine put
   end subroutine init_units_length


   subroutine init_units_momentum
      units_momentum%desc = "momentum"
      call put("au",      1.d0)
      call put("Angst-1", angstroem)
      call put("pm-1",    picometer)
      call put("nm-1",    nanometer)
      contains
      subroutine put(label,cval)
         character(len=*),intent(in) :: label
         real(dbl),intent(in)        :: cval
         call map_put(units_momentum%map, label, cval)
      end subroutine put
   end subroutine init_units_momentum


   subroutine init_units_time
      units_time%desc = "time"
      call put("au",   1.d0)
      call put("fs",   ufs)
      call put("ps",   1.d3*ufs)
      call put("fs-1", 1.d0/ufs)
      contains
      subroutine put(label,cval)
         character(len=*),intent(in) :: label
         real(dbl),intent(in)        :: cval
         call map_put(units_time%map, label, cval)
      end subroutine put
   end subroutine init_units_time


!=======================================================================


   function parse_value(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units,dflt)
      else
         val = parse_with_unit(t,units)
      endif
   end function parse_value
   
   
   function parse_energy(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units_energy,dflt)
      else
         val = parse_with_unit(t,units_energy)
      endif
   end function parse_energy


   function parse_mass(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units_mass,dflt)
      else
         val = parse_with_unit(t,units_mass)
      endif
   end function parse_mass


   function parse_length(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units_length,dflt)
      else
         val = parse_with_unit(t,units_length)
      endif
   end function parse_length


   function parse_momentum(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units_momentum,dflt)
      else
         val = parse_with_unit(t,units_momentum)
      endif
   end function parse_momentum


   function parse_time(t,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      if (present(dflt)) then
         val = parse_with_unit(t,units_time,dflt)
      else
         val = parse_with_unit(t,units_time)
      endif
   end function parse_time


   function parse_with_unit(t,u,dflt) result(val)
      type(tokenizer_t),intent(inout) :: t
      type(unitclass_t),intent(in)    :: u
      real(dbl),intent(in),optional   :: dflt
      real(dbl)                       :: val
      character(len=maxtoklen) :: token
      integer :: ierr
      real(dbl) :: ufac
      token = t%get()
      read(token,*,iostat=ierr) val
      if (ierr/=0) then
         if (present(dflt)) then
            val = dflt
            return
         else
            call t%error("expected real number")
         endif
      endif
      call t%gofwd
      token = t%get()
      if (token /= ",") return
      call t%gofwd
      token = t%get()
      if (.not.map_get(u%map,token,ufac)) &
         call t%error("expected "//trim(u%desc)//" unit")
      val = val*ufac
      call t%gofwd
   end function parse_with_unit

end module units_m
