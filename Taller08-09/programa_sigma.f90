!!!!!!!!!! USAR MODULO Y CALCULAR (I) EN t[0,1], Y OBTENER σ !!!!!!!!!!

program programa_sigma
    use modulo_romberg
    implicit none

    real(dp), parameter :: tol = 1.0e-8_dp                                      !definir la tolerancia pedida
    integer, parameter  :: max_k = 25                                           !definir el numero maximo de niveles
    real(dp), parameter :: h_val  = 6.62607015e-34_dp                           !constante de Planck [J*s]
    real(dp), parameter :: c_val  = 299792458.0_dp                              !velocidad de la luz [m/s]
    real(dp), parameter :: kb_val = 1.380649e-23_dp                             !constante de Boltzmann [J/K]

    real(dp) :: integral_t, integral_x, sigma_t, sigma_x, factor_pre            !definir las variables que usaremos
    real(dp), allocatable :: S(:,:)                                             !Arreglo de maximo de niveles
    integer :: k_fin_t, neval_t, k_fin_x, neval_x, u_out                        !definir las variables para los archivos de salida

    allocate(S(0:max_k, 0:max_k))                                               !Reservar memoria para el arreglo

    !! Factor multiplicativo para obtener sigma a partir de I !!
    factor_pre = (2.0_dp * acos(-1.0_dp) * kb_val**4) / (h_val**3 * c_val**2)    !definicion del factor de precalculo para no saturar el codigo despues


    !! Calculo con cambio de variable en t in [0, 1] !!
    call integrar_romberg(0.0_dp, 1.0_dp, tol, max_k, 3, S, k_fin_t, neval_t)          !calcular la integral con cambio de variable
    integral_t = S(k_fin_t, K_fin_t)                                                   !obtener la integral final (en t)
    sigma_t = factor_pre * integral_t                                                !obtener el valor de sigma con la ecuación


    !! Calculo sin cambio de variable en x in [0, 100] !!
    call integrar_romberg(0.0_dp, 100.0_dp, tol, max_k, 4, S, k_fin_x, neval_x)        !calcular la integral sin cambio de variable
    integral_x = S(k_fin_x, K_fin_x)                                                   !obtener la integral final (en x)
    sigma_x = factor_pre * integral_x                                                  !obtener el valor de sigma con la ecuación

    !! Escribimos el archivo resultado.sigma.dat !!
    open(newunit=u_out, file='resultado_sigma.dat', status='replace', action='write')

    write(u_out, '(i4, i10, 2es24.16)') k_fin_t, neval_t, integral_t, sigma_t                !fila 1 (cambio de varaible)
    write(u_out, '(i4, i10, 2es24.16)') k_fin_x, neval_x, integral_x, sigma_x                !fila 2 (sin cambio de variable)

    close(u_out)
    deallocate(S)                                                                            !Liberar memoria

    print *, "La parte C fue ejecutada de forma correcta"
end program programa_sigma
