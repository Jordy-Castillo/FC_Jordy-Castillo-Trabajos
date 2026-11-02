
import numpy as np
import matplotlib.pyplot as plt

### PARAMETROS DEL OSCILADOR ###

gamma = 0.5                   #coeficiente de amortiguamiento (s^-1)
w = 2*np.pi                   #frecuencia angular (rad/s)
w0_sq = w**2 + gamma**2

### DEFINICION DE FUNCIONES DE REFERENCIA DE x, v, a SEGUN LAS ECUACIONES ###

def x_ex(t):
  return np.exp(-gamma*t) * np.cos(w*t)

def v_ex(t):
  return -np.exp(-gamma*t) * (gamma*np.cos(w*t) + w*np.sin(w*t))

def a_ex(t):
  return np.exp(-gamma*t) * ((gamma**2 - w**2) * np.cos(w*t) + 2*gamma*w*np.sin(w*t))

### FUNCION PARA CALCULAR DERIVADAS Y ERRORES SEGUN EL dt ###

def procesar_malla(dt):
  N = int(round(3.0 / dt))              #calcular numero de pasos en 3 segundos, con redondeo al entero mas cercano
  t = np.linspace(0.0, 3.0, N+1)        #arreglo de tiempo N+1 puntos para 0.01s
  x = x_ex(t)                           #posición exacta del arreglo

  ## VELOCIDADES CON 4 FORMULAS DE C.2 ##
  v_adelantada = (x[1:] - x[:-1]) / dt                                          #rebanadas segun la tabla del codigo de teoria dado
  v_atrasada   = (x[1:] - x[:-1]) / dt
  v_central    = (x[2:] - x[:-2]) / (dt *2.0)
  v_cinco_puntos = (-x[4:] + 8.0*x[3:-1] - 8.0*x[1:-3] + x[:-4]) / (12.0*dt)

  ## ACELERACION CON 1 FORMULA DE C.2 ##
  a_central = (x[2:] - 2.0*x[1:-1] + x[:-2]) / (dt**2)

  ## ERRORES MAXIMOS ABSOLUTOS ##
  error_v_adelantada   = np.max(np.abs(v_adelantada   - v_ex(t[:-1])))        #formula del error absoluto
  error_v_atrasada     = np.max(np.abs(v_atrasada     - v_ex(t[1:])))
  error_v_central      = np.max(np.abs(v_central      - v_ex(t[1:-1])))
  error_v_cinco_puntos = np.max(np.abs(v_cinco_puntos - v_ex(t[2:-2])))
  error_a_central      = np.max(np.abs(a_central      - a_ex(t[1:-1])))

  error_fila = [dt, error_v_adelantada, error_v_atrasada, error_v_central, error_v_cinco_puntos, error_a_central]
  return t, x, v_adelantada, v_atrasada, v_central, v_cinco_puntos, a_central, error_fila

### PROCESO PARA dt = 0.01s de C.3  ###

dt_1 = 0.01
t1, x1, v_ad1, v_at1, v_cent1, v_5p1, a_cent1, error1 = procesar_malla(dt_1)

## guardar en derivadas_oscilador.dat con los 299 puntos
t_int = t1[1:-1]
x_int = x1[1:-1]
datos_salida1 = np.column_stack((t_int, x_int, v_cent1, a_cent1))
np.savetxt("derivadas_oscilador.dat", datos_salida1, fmt="%25.15e")

## generamos la grafica comparativa de derivadas_oscilador.png ##

plt.figure(figsize=(10, 6))
plt.subplot(2,1,1)
plt.plot(t_int, v_ex(t_int), '-', label="v exacta", color="red", linewidth=1.5)
plt.plot(t_int, v_cent1, '--', label="v central numerica.", color="blue", linewidth=1.2)
plt.ylabel("Velocidad [m/s]")
plt.title("Velocidad y Aceleración del Oscilador Amortiguado (dt = 0.01 s)")
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(2, 1, 2)
plt.plot(t_int, a_ex(t_int), '-', label="a exacta", color="orange", linewidth=1.5)
plt.plot(t_int, a_cent1, '--', label="a central numerica.", color="purple", linewidth=1.2)
plt.xlabel("Tiempo [s]")
plt.ylabel("Aceleración [m/s^2]")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("derivadas_oscilador.png", dpi=250)
plt.close()

### PROCESO PARA dt = 0.005s de C.4  ###

dt_2 = 0.005
t2, x2, v_ad2, v_at2, v_cent2, v_5p2, a_cent2, error2 = procesar_malla(dt_2)
matriz_errores = np.array([error1, error2])
np.savetxt("errores_derivadas.dat", matriz_errores, fmt="%25.15e")

print("Parte C ejecutada y archivos generados")
