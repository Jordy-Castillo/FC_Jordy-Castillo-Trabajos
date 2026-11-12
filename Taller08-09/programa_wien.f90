!!!!!!!!!! ESCRIBIR LOS 3 ARCHIVOS DE SALIDA, LA CONSTANTE WEIN Y X5 Y X3 !!!!!!!!!!

program programa_wien
    use modulo_ceros, only: dp, biseccion, newton, secante
    implicit none

    real(dp), parameter :: tol = 1.0e-10_dp                                            !definimos la tolerancia epsilon dada
    integer, parameter  :: max_iter = 100                                              !maximo de iteraciones
    real(dp), parameter :: h_val  = 6.62607015e-34_dp                                  !constante de Planck [J*s]
    real(dp), parameter :: c_val  = 299792458.0_dp                                     !velocidad de la luz [m/s]
    real(dp), parameter :: kb_val = 1.380649e-23_dp                                    !constante de Boltzmann [J/K]

    real(dp) :: x5, x3, b_wien                                                         !definimos las variables que usaremos
    real(dp) :: raiz_res                                                               !raiz resultante
    integer :: it_res, u_res, u_it, u_cte                                              !variables para los archivos de salida

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GENERAR ARCHIVO RESULTADO.WIEN.DAT E ITERACIONES_WIEN.DAT !!!

!! Las iteraciones de n=5 las registra el propio modulo (u_hist), asi la raiz salen del mismo codigo !!

    open(newunit=u_res, file='resultado_wien.dat',   status='replace', action='write')      !archivo de resultados
    open(newunit=u_it,  file='iteraciones_wien.dat', status='replace', action='write')      !archivo de iteraciones

    call biseccion(5.0_dp, 1.0_dp, 10.0_dp, tol, max_iter, raiz_res, it_res, u_hist=u_it)   !biseccion n=5 (registra iteraciones)
    write(u_res, '(i2, i2, es24.16, i4)') 5, 1, raiz_res, it_res                            !se escribe en el archivo de salida

    call newton(5.0_dp, 10.0_dp, tol, max_iter, x5, it_res, u_hist=u_it)                    !Newton n=5 (registra iteraciones)
    write(u_res, '(i2, i2, es24.16, i4)') 5, 2, x5, it_res                                  !se escribe en el archivo de salida

    call secante(5.0_dp, 10.0_dp, 9.0_dp, tol, max_iter, raiz_res, it_res, u_hist=u_it)     !secante n=5 (registra iteraciones)
    write(u_res, '(i2, i2, es24.16, i4)') 5, 3, raiz_res, it_res                            !se escribe en el archivo de salida

    call biseccion(3.0_dp, 1.0_dp, 10.0_dp, tol, max_iter, raiz_res, it_res)                !se ejecuta la bisección para n=3
    write(u_res, '(i2, i2, es24.16, i4)') 3, 1, raiz_res, it_res                            !se escribe en el archivo de salida

    call newton(3.0_dp, 10.0_dp, tol, max_iter, x3, it_res)                                 !se ejecuta el Newton-Raphson para n=3
    write(u_res, '(i2, i2, es24.16, i4)') 3, 2, x3, it_res                                  !se escribe en el archivo de salida

    call secante(3.0_dp, 10.0_dp, 9.0_dp, tol, max_iter, raiz_res, it_res)                  !se ejecuta la secante para n=3
    write(u_res, '(i2, i2, es24.16, i4)') 3, 3, raiz_res, it_res                            !se escribe en el archivo de salida

    !!!! PUNTO B.3 !!!!
    call newton(5.0_dp, 0.5_dp, tol, max_iter, raiz_res, it_res)                            !se ejecuta el Newton-Raphson para n=5 con x0 = 0.5
    write(u_res, '(i2, i2, es24.16, i4)') 5, 4, raiz_res, it_res                            !se escribe en el archivo de salida

    close(u_res)
    close(u_it)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! GENERAR ARCHIVIO CONSTANTE_WIEN.DAT !!!
    b_wien = (h_val * c_val) / (kb_val * x5)                                                !se calcula la constante de Wien mediante su ecuacion

    open(newunit=u_cte, file='constante_wien.dat', status='replace', action='write')
    write(u_cte, '(3es24.16)') x5, x3, b_wien
    close(u_cte)

    print *, "La parte B ha sido ejecutada de forma correcta"
end program programa_wien
