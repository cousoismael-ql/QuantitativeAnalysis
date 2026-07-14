program ssfm_nlse_normalizada
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'

  integer, parameter :: dp = C_DOUBLE

  real(dp), parameter :: pi = acos(-1.0_dp)

  integer :: nt
  integer :: nxi
  integer :: save_every

  real(dp) :: taumin
  real(dp) :: taumax
  real(dp) :: ximax
  real(dp) :: Nsol

  complex(C_DOUBLE_COMPLEX), parameter :: Iunit = cmplx(0.0_dp, 1.0_dp, kind=C_DOUBLE_COMPLEX)

  real(dp) :: dtau, dxi, window
  real(dp) :: tau, xi
  integer :: i, n

  complex(C_DOUBLE_COMPLEX), allocatable :: u(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: uhat(:)
  real(dp), allocatable :: omega(:)
  real(dp), allocatable :: I_spec(:)

  type(C_PTR) :: plan_forward, plan_backward

  integer :: j, idx, kshift
  real(dp) :: omega_s, Imax, I_dB
  real(dp) :: norm_fft

  integer :: ios, unit_params
  character(len=256) :: param_file

  namelist /parametros/ nt, nxi, save_every, taumin, taumax, ximax, Nsol

  ! ============================================================
  ! Leer parámetros desde archivo
  ! ============================================================

  nt = 1024
  nxi = 5000
  save_every = 50

  taumin = -40.0_dp
  taumax =  40.0_dp

  ximax = 5.0_dp
  Nsol = 3.0_dp

  param_file = "params.txt"

  if (command_argument_count() >= 1) then
     call get_command_argument(1, param_file)
  end if

  open(newunit=unit_params, file=trim(param_file), status="old", action="read", iostat=ios)

  if (ios /= 0) then
     print *, "Error: no se pudo abrir el archivo de parámetros: ", trim(param_file)
     stop
  end if

  read(unit_params, nml=parametros, iostat=ios)

  if (ios /= 0) then
     print *, "Error: no se pudo leer correctamente el archivo de parámetros."
     stop
  end if

  close(unit_params)

  if (nt <= 0) stop "Error: nt debe ser positivo."
  if (nxi <= 0) stop "Error: nxi debe ser positivo."
  if (save_every <= 0) stop "Error: save_every debe ser positivo."
  if (taumax <= taumin) stop "Error: taumax debe ser mayor que taumin."
  if (ximax <= 0.0_dp) stop "Error: ximax debe ser positivo."
  if (Nsol <= 0.0_dp) stop "Error: Nsol debe ser positivo."

  dtau = (taumax - taumin) / real(nt, dp)
  dxi = ximax / real(nxi, dp)
  window = taumax - taumin

  allocate(u(0:nt-1))
  allocate(uhat(0:nt-1))
  allocate(omega(0:nt-1))
  allocate(I_spec(0:nt-1))

  ! Frecuencias angulares
  do i = 0, nt-1
     if (i <= nt/2) then
        omega(i) = 2.0_dp*pi*real(i, dp)/window
     else
        omega(i) = 2.0_dp*pi*real(i-nt, dp)/window
     end if
  end do

  ! Condición inicial:
  ! u(0,tau) = N sech(tau)
  do i = 0, nt-1
     tau = taumin + real(i, dp)*dtau
     u(i) = Nsol * sech(tau)
  end do

  ! Crear planes FFTW
  plan_forward = fftw_plan_dft_1d(nt, u, uhat, FFTW_FORWARD, FFTW_ESTIMATE)
  plan_backward = fftw_plan_dft_1d(nt, uhat, u, FFTW_BACKWARD, FFTW_ESTIMATE)

  open(unit=10, file="nlse_normalizada.dat", status="replace", action="write")

  call save_field(10, 0.0_dp, u, nt, dtau, taumin)

   do n = 1, nxi

      call ssfm_step_trapezoid(u, uhat, nt, omega, dxi, plan_forward, plan_backward)

      if (mod(n, save_every) == 0 .or. n == nxi) then
         xi = real(n, dp)*dxi
         call save_field(10, xi, u, nt, dtau, taumin)
      end if

   end do

  close(10)

  ! ============================================================
  ! Espectro de salida
  ! ============================================================

   call fftw_execute_dft(plan_forward, u, uhat)

   norm_fft = dtau / sqrt(2.0_dp*pi)

   do j = 0, nt-1
      I_spec(j) = abs(norm_fft * uhat(j))**2
   end do

   Imax = maxval(I_spec)

   if (Imax > 0.0_dp) then
      I_spec = I_spec / Imax
   end if

   open(unit=20, file='espectro_salida.dat', status='replace', action='write')

   write(20,*) '# omega_normalizada   Intensidad_normalizada   Intensidad_dB'

   do j = 0, nt-1

      kshift = j - nt/2
      idx = modulo(kshift, nt)

      omega_s = 2.0_dp*pi*real(kshift, dp)/(real(nt, dp)*dtau)

      I_dB = 10.0_dp*log10(I_spec(idx) + 1.0e-30_dp)

      write(20,'(3ES22.12)') omega_s, I_spec(idx), I_dB

   end do

   close(20)

   call fftw_destroy_plan(plan_forward)
   call fftw_destroy_plan(plan_backward)

   print *, "Simulación terminada."
   print *, "Archivo generado: nlse_normalizada.dat"
   print *, "Archivo generado: espectro_salida.dat"
   print *, "Columnas campo: xi, tau, Re(u), Im(u), |u|^2"
   print *, "Columnas espectro: omega, I_norm, I_dB"

contains

  function sech(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: y

    y = 1.0_dp / cosh(x)
  end function sech


  subroutine linear_half_step(u, uhat, nt, omega, dxi, plan_forward, plan_backward)
    integer, intent(in) :: nt
    real(dp), intent(in) :: omega(0:nt-1), dxi
    complex(C_DOUBLE_COMPLEX), intent(inout) :: u(0:nt-1)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: uhat(0:nt-1)
    type(C_PTR), intent(in) :: plan_forward, plan_backward

    real(dp) :: phase
    integer :: k

    ! Transformada directa: u(tau) -> uhat(omega)
    call fftw_execute_dft(plan_forward, u, uhat)

    ! Medio paso lineal en Fourier
    do k = 0, nt-1
       phase = -omega(k)**2 * dxi / 4.0_dp
       uhat(k) = uhat(k) * exp(cmplx(0.0_dp, phase, kind=C_DOUBLE_COMPLEX))
    end do

    ! Transformada inversa: uhat(omega) -> u(tau)
    call fftw_execute_dft(plan_backward, uhat, u)

    ! FFTW no normaliza la transformada inversa
    u = u / real(nt, dp)

  end subroutine linear_half_step

  subroutine ssfm_step_trapezoid(u, uhat, nt, omega, dxi, plan_forward, plan_backward)
      integer, intent(in) :: nt
      real(dp), intent(in) :: omega(0:nt-1), dxi
      complex(C_DOUBLE_COMPLEX), intent(inout) :: u(0:nt-1)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: uhat(0:nt-1)
      type(C_PTR), intent(in) :: plan_forward, plan_backward

      integer, parameter :: max_iter = 4
      real(dp), parameter :: tol = 1.0e-8_dp

      complex(C_DOUBLE_COMPLEX) :: u_old(0:nt-1)
      complex(C_DOUBLE_COMPLEX) :: u_half(0:nt-1)
      complex(C_DOUBLE_COMPLEX) :: u_guess(0:nt-1)
      complex(C_DOUBLE_COMPLEX) :: u_new(0:nt-1)

      integer :: i, iter
      real(dp) :: phase
      real(dp) :: err, num, denom

      ! Campo actual: A(z)
      u_old = u

      ! Primer medio paso lineal:
      ! A_L = exp(hD/2) A(z)
      u_half = u_old
      call linear_half_step(u_half, uhat, nt, omega, dxi, plan_forward, plan_backward)

      ! ------------------------------------------------------------
      ! Predictor:
      ! Se inicia con N(A_sig^(0)) = N(A_actual)
      ! Por tanto, se usa |u_old|^2 para la primera fase no lineal.
      ! ------------------------------------------------------------
      u_guess = u_half

      do i = 0, nt-1
         phase = dxi * abs(u_old(i))**2
         u_guess(i) = u_guess(i) * exp(Iunit * phase)
      end do

      ! Segundo medio paso lineal para obtener la primera estimación
      ! de A(z + h)
      call linear_half_step(u_guess, uhat, nt, omega, dxi, plan_forward, plan_backward)

      ! ------------------------------------------------------------
      ! Corrector:
      ! integral N dz ≈ h/2 [N(A(z)) + N(A(z+h))]
      ! ------------------------------------------------------------
      do iter = 1, max_iter

         u_new = u_half

         do i = 0, nt-1
            phase = 0.5_dp * dxi * ( abs(u_old(i))**2 + abs(u_guess(i))**2 )
            u_new(i) = u_new(i) * exp(Iunit * phase)
         end do

         ! Segundo medio paso lineal
         call linear_half_step(u_new, uhat, nt, omega, dxi, plan_forward, plan_backward)

         ! Error relativo tipo norma L2:
         ! sqrt(sum |u_new - u_guess|^2) / sqrt(sum |u_guess|^2)
         num = sqrt(sum(abs(u_new - u_guess)**2))
         denom = sqrt(sum(abs(u_guess)**2))

         if (denom < tiny(denom)) denom = 1.0_dp

         err = num / denom

         if (err < tol) exit

         u_guess = u_new

      end do

      ! Campo actualizado: A(z + h)
      u = u_new

      end subroutine ssfm_step_trapezoid


  subroutine save_field(unit, xi, u, nt, dtau, taumin)
    integer, intent(in) :: unit, nt
    real(dp), intent(in) :: xi, dtau, taumin
    complex(C_DOUBLE_COMPLEX), intent(in) :: u(0:nt-1)

    integer :: i
    real(dp) :: tau

    do i = 0, nt-1
       tau = taumin + real(i, dp)*dtau
       write(unit,'(5ES22.12)') xi, tau, real(u(i)), aimag(u(i)), abs(u(i))**2
    end do

    write(unit,*)
  end subroutine save_field

end program ssfm_nlse_normalizada
