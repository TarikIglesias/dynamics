__doc__ = """
En este módulo de Python 3 se definen funciones utiles para el manejo de los
datos de tiempos de viaje que se utilizarán en la Tarea 2 del curso GF5013.

GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

5 de Septiembre de 2017

"""
import numpy as NP


######
def readTravelTimeData(filename, headerlines = 0):
    """
    lee los datos del experimento de tomografia del archivo NombreArchivo y
    extrae la ubicacion de la fuente (xF,zF), la ubicacion del receptor (xR,
    zR) y el tiempo de viaje para cada rayo registrado (TV).
    Cada fila del archivo que se lee corresponde a un rayo del experimento.
    Las coordenadas de las fuentes y receptores estan dadas en metros y el
    tiempo de viaje esta dado en milisegundos.

    La función devuelve un diccionario con los valores indicados, 
    las llaves de dicho diccionario son xF, zF, xR, zR, TV.
    """
    # open and read the datafile
    file = open(filename, 'r')
    
    # initialize containers
    xF = [] # coordenada X de fuente
    zF = [] # coordenada Z de fuente
    xR = [] # coordenada X de receptor
    zR = [] # coordenada Z de receptor
    TV = [] # Tiempo de viaje observado
    i = 0
    headers = [] # empty list
    for line in file:
       i += 1
       if i <= headerlines:
          headers.append( line )
          continue # goes to the next iteration of the for loop.
       # if the code gets here, the line is not a header
       line = line.split()
       xF.append( float( line[0] ) )
       zF.append(- float( line[1] ) )
       xR.append( float( line[2] ) )
       zR.append(- float( line[3] ) )
       TV.append( float( line[4] ) )
    
    file.close()
    
    # assemble everything into a dictionary and convert lists into numpy arrays
    data = {}
    data['headers'] = headers
    data['xF'] = NP.array( xF )
    data['zF'] = NP.array( zF )
    data['xR'] = NP.array( xR )
    data['zR'] = NP.array( zR )
    data['TV'] = NP.array( TV )
    
    return data

    






    
    
    
