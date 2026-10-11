! Reference driver for the production DISORT kernel at any default-REAL kind.
!
! Built by build_reference.sh with gfortran in double (-fdefault-real-8) and
! quadruple precision (-freal-4-real-16 -freal-8-real-16).  Every real input
! is read as an exact integer mantissa M and binary exponent E (value =
! M * 2**E), so float32 or float64 inputs are represented exactly at any
! precision.  Outputs are written with 36 significant digits.
!
! Input file (list-directed):
!   nlyr nmom ntau numu nphi nstr plank
!   fbeam umu0 phi0 fisot accur wvnmlo wvnmhi btemp ttemp temis albedo   (M E pairs)
!   dtauc(nlyr) ssalb(nlyr) pmom(0:nmom,nlyr) temper(0:nlyr) utau(ntau) umu(numu) phi(nphi)  (M E pairs)
!
! Usage: disort_ref_driver INPUT OUTPUT

program disort_ref_driver
   implicit none
   integer :: nlyr, nmom, ntau, numu, nphi, nstr, iplank, status, lc, k
   integer :: ibcnd, maxcly, maxulv, maxumu, maxphi, maxmom
   logical :: plank, prnt(5), usrtau, usrang, lamber, onlyfl
   real :: fbeam, umu0, phi0, fisot, accur, wvnmlo, wvnmhi, btemp, ttemp, temis, albedo
   real, allocatable :: dtauc(:), ssalb(:), pmom(:,:), temper(:), utau(:), umu(:), phi(:)
   real, allocatable :: rfldir(:), rfldn(:), flup(:), dfdt(:), uavg(:), uu(:,:,:), albmed(:), trnmed(:)
   integer, external :: disort
   character(len=1024) :: infile, outfile
   real :: x

   call get_command_argument(1, infile)
   call get_command_argument(2, outfile)
   open(10, file=trim(infile), status='old', action='read')
   read(10, *) nlyr, nmom, ntau, numu, nphi, nstr, iplank
   plank = iplank /= 0
   fbeam = exact(10); umu0 = exact(10); phi0 = exact(10); fisot = exact(10); accur = exact(10)
   wvnmlo = exact(10); wvnmhi = exact(10); btemp = exact(10); ttemp = exact(10); temis = exact(10)
   albedo = exact(10)
   allocate(dtauc(nlyr), ssalb(nlyr), pmom(0:nmom, nlyr), temper(0:nlyr), utau(ntau), umu(numu), phi(nphi))
   do lc = 1, nlyr
      dtauc(lc) = exact(10)
   end do
   do lc = 1, nlyr
      ssalb(lc) = exact(10)
   end do
   do lc = 1, nlyr
      do k = 0, nmom
         pmom(k, lc) = exact(10)
      end do
   end do
   do lc = 0, nlyr
      temper(lc) = exact(10)
   end do
   do k = 1, ntau
      utau(k) = exact(10)
   end do
   do k = 1, numu
      umu(k) = exact(10)
   end do
   do k = 1, nphi
      phi(k) = exact(10)
   end do
   close(10)

   allocate(rfldir(ntau), rfldn(ntau), flup(ntau), dfdt(ntau), uavg(ntau), uu(numu, ntau, nphi), &
            albmed(numu), trnmed(numu))
   ! DISORT's first-call self-test (SLFTST) temporarily overwrites its scalar
   ! arguments, so every argument must be a variable, never a literal.
   prnt = .false.
   usrtau = .true.; usrang = .true.; lamber = .true.; onlyfl = .false.; ibcnd = 0
   maxcly = nlyr; maxulv = ntau; maxumu = numu; maxphi = nphi; maxmom = nmom
   status = disort(nlyr, dtauc, ssalb, nmom, pmom, temper, wvnmlo, wvnmhi, usrtau, ntau, utau, nstr, &
                   usrang, numu, umu, nphi, phi, ibcnd, fbeam, umu0, phi0, fisot, lamber, albedo, btemp, &
                   ttemp, temis, plank, onlyfl, accur, prnt, maxcly, maxulv, maxumu, maxphi, maxmom, rfldir, &
                   rfldn, flup, dfdt, uavg, uu, albmed, trnmed)

   open(11, file=trim(outfile), status='replace', action='write')
   write(11, '(a, i0, a, i0, a, i0)') '# status ', status, ' real_bytes ', storage_size(x) / 8, &
      ' digits ', digits(x)
   call put('rfldir', rfldir); call put('rfldn', rfldn); call put('flup', flup)
   call put('dfdt', dfdt); call put('uavg', uavg); call put('uu', reshape(uu, [size(uu)]))
   call put('albmed', albmed); call put('trnmed', trnmed)
   close(11)

contains

   real function exact(unit)
      integer, intent(in) :: unit
      integer(kind=8) :: m
      integer :: e
      read(unit, *) m, e
      exact = scale(real(m), e)
   end function exact

   subroutine put(name, values)
      character(len=*), intent(in) :: name
      real, intent(in) :: values(:)
      integer :: i
      write(11, '(a, 1x, i0)') name, size(values)
      do i = 1, size(values)
         write(11, '(es45.35e4)') values(i)
      end do
   end subroutine put

end program disort_ref_driver
