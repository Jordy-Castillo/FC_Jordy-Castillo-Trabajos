########### Módulo de análisis ###########

import numpy as np
import matplotlib.pyplot as plt

#Leeremos el archivo de Fortran valle.dat para tener como salidas los arreglos de
#los numeros de masa A, numeros atomicos de maximo alcance Z, y energiar de enlace EB 
#y asi poder utilizarlos en este lenguaj de python

def leer_valle(archivo):                            
  """
  Lee el archivo de datos valle.dat
  
  ENTRADA: un archivo (str)
  
  SALIDA: tupla de arreglos (a,z,EB)
          a: numeros de masa
          z: numeros atomicos
          EB: energia de enlace por nucleon

  """
  datos = np.loadtxt(archivo)                       
  a = datos[:,0]
  z = datos[:,1]
  EB = datos[:,2]
  return a, z, EB


#encontramos el nucleo con mayor enlace por nucleon, es decir, el maximo alcance es
#el valor mas negativo de E_B, como se habia definido anteriormente

def encontrar_maximo(a,z,EB):   
  """
  Encuentra el nucleo con mayor enlace por nucleon de toda la curva a construir
  
  ENTRADA: arreglos (a, z, EB)
          a: numeros de masa
          z: numeros atomicos
          EB: energia de enlace por nucleon
  
  SALIDA: tupla (a_max, z_max, EB_min) 
          son valores correspondientes al nucleo mas ligado

  """                                              
  a_min = a[0]                                       
  z_min = z[0]
  EB_min = EB[0]

  for i in range(len(EB)):           
    if EB[i] < EB_min:
      a_min = a[i]
      z_min = z[i]
      EB_min = EB[i]

  return a_min, z_min, EB_min


#escribimos el nucleo mas ligado encontrado por la funcion anterior

def escribir_maximo_python(archivo_salida="maximo_python.dat", archivo_valle="valle.dat"):
  """
  Escribe el nucleo mas ligado encontrado por python en un archivo
  
  ENTRADA: 
          a_max: numeros de masa A del minimo global
          z_max: numero atomico Z del minimo global
          EB_max: energia de enlace minima (MeV)
          archivo: ruta de destino maximo_python.dat
  
  SALIDA: ninguna, ya que genera un archivo de texto estilo fotran

  """   

  A, Z, EB = leer_valle(archivo_valle)                                
  a_max, z_max, eb_max = encontrar_maximo(A, Z, EB)                   
  
  with open(archivo_salida, "w") as f:                                #abrimos el archivo de salida en modo escritura (mismo formato usado en fortran)
    f.write(f"{int(a_max):6d}{int(z_max):6d}{eb_max:18.8f}\n")


#generamos y guardamos la grafica de B/A = -EB vs A, con los puntos experimentales dados para esta parte

def graficar_curva(archivo_valle="valle.dat", experimental_data=None, imagen_salida="curva_enlace.png"):
  """
  Genera y guarda la grafica de B/A = -EB vs A con los puntos experimentales dados
  
  ENTRADA: 
          a: numeros de masa A del modelo
          EB: energias de enlace de EB del modelo (-Mev)
          datos_exp: lista de tuplas con datos AME2020 (nombre, Z, A, B/A)
          archivo_png: ruta donde se gaurdara la imagen
  
  SALIDA: ninguna, ya que genera un archivo de texto estilo fotran

  """ 

  A, Z, EB = leer_valle(archivo_valle)                                        #Leer los datos del modelo
  B_A_modelo = -EB                                                            #Calcular la energía de enlace por nucleón

  AZUL = "#162A54"                                                            #Colores de la gráfica
  ROJO = "#B02A20"

  plt.figure(figsize=(8, 5))                                                  #Dibujar la curva teórica
  plt.plot(A, B_A_modelo, color=AZUL, label="Modelo de Gota Líquida")
 
  if experimental_data is not None:           #comprueba que existan datos experimentales
        exp_A = []                              #lista para guardar coordenada X
        exp_BA = []                             #lista para guardar corrdenada Y
        exp_labels = []                         #lista para guardar los nombres
        
        for item in experimental_data:          #Obtener los datos de cada elemento
            exp_labels.append(item[0])          #toma el nombre del elemento y lo pone en la lista de nombres
            exp_A.append(item[2])               #toma el numero de nucleones y lo pone en la lista de nucleones
            exp_BA.append(item[3])              #toma la energia de enlace por nucleon y lo pone en la lista de energias de enlace

        plt.scatter(exp_A, exp_BA, color=ROJO, label="Datos experimentales (AME2020)")   #dibujar los puntos experimentales
                                                                                         #scatter recibe dos listas (NUEVO PARA MI)

        for i in range(len(exp_A)):                                             #escribir el nombre de cada elemento junto a su punto
            y_offset= 6 if i % 2 == 0 else -12                #añadi esto para que los nombres no se sobreponieran, y asi poderlos leer
            plt.annotate(exp_labels[i], (exp_A[i], exp_BA[i]), textcoords="offset points", xytext=(6, y_offset), ha="center", fontsize=8, color=ROJO)
            #escribir el nombre de la posicion i en sus coordenadas. #annotate (NUEVO PARA MI)
    
  plt.xlabel("Número de nucleones A")                                        
  plt.ylabel("Energía de enlace por nucleón B/A [MeV]")
  plt.title("Valle de Estabilidad y Curva de Enlace Nuclear")
  plt.grid(True)
  plt.legend()
  plt.tight_layout()
  plt.savefig(imagen_salida, dpi=300)
  plt.close()
