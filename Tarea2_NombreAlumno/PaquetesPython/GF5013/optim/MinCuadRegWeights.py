__doc__ = """
GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

Octubre 2017
"""
from numpy.linalg import lstsq 
import numpy as NP
from .MinCuad import MinCuadSimple

# la variable COMPLETAR es necesaria para poder importar el paquete de python con los
# codigos incompletos. Donde encuentre esta variable debe completar usted el codigo.
COMPLETAR = None 

######
def MinCuadRegWeights(G, d, Wx, H, ho, Wh, epsilon):
   """
   Resuelve el problema de mínimos cuadrados regularizados espresado con matrices de
   pesos.
   Se reescribe el problema
   Min_m ||Wx*(G*m-d)||_2^2 + epsilon^2 * ||Wh*(H*m-ho)||_2^2

   como Min_m ||Fm-D||_2^2 con matrix F y vector D adecuados. 
  
   devuelve un diccionario con las siguientes variables (igual nombre de llaves)

   m = modelo estimado
   Cm = matriz de covarianza del modelo estimado


   IMPORTANT!!! ALL VARIABLES MUST BE 2D arrays (vectors must have shape (N,1))

   """

   WxG = COMPLETAR 
   Wxd = COMPLETAR
   eWhH = COMPLETAR
   eWhho = COMPLETAR
 
   F = NP.vstack( COMPLETAR )
   D = NP.vstack( COMPLETAR )
    
   return COMPLETAR

