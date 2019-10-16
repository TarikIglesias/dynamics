__doc__ = """

GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

6 de Septiembre de 2017 

"""
import numpy as NP
from matplotlib import path
######
class RayLinear2D(object):
   """
   defines an object cappable of modeling linear rays given source and 
   receiver location. 

   Also define member functions that operate in RectangularMesh instances. 
   """
   ######
   def __init__(self, xS, yS, xR, yR):
      """
      xS = X coordinate of source
      yS = Y coordinate of source
      xR = X coordinate of receiver
      yR = Y coordinate of receiver. 
      """
      # coordinates of source
      self.xS = xS
      self.yS = yS
      # coordinates of receiver
      self.xR = xR
      self.yR = yR
      # length of the ray
      self.L = NP.sqrt( (xR - xS)**2 + (yR - yS)**2 )

   ######
   def calcMeshRayLength(self, RMesh, fractional_precision = 1E-2 ):
      """
      Estimates the length of the ray on each element of the Rectangular Mesh
      instance RMesh. 
 
      fractional_precision : needed for numerical integration of the ray along 
         the mesh. The ray is subdivided in segments of length <= dR, where 
         dR = fractional_precision * length of shorter side of mesh element.

      NOTA: Esta funcion es totalmente ineficiente en cuanto al tiempo de 
	calculo, y puede hacerce en mejor forma. Para efectos del curso se
	utiliza esta funcion ya que es la mas simple de realizar y de entender.
        Basicamente se divide el rayo en muchos puntos discretos y se cuenta la
	cantidad de puntos que cae dentro de cada elemento de la malla (RMesh).
	Luego multiplicando la cantidad de puntos por el espaciamiento entre estos 
        se obtiene la longitud del rayo en un casillero de manera aproximada.
 
      """
      # define the points that represent the ray
      a = NP.min( [ RMesh.dX, RMesh.dY ] ) # length of shorter side of element.
      dR_target = a * fractional_precision
      NumR = NP.ceil( self.L / dR_target  )
      Rx = NP.linspace(self.xS, self.xR, NumR)
      Ry = NP.linspace(self.yS, self.yR, NumR)
      dR = NP.sqrt( (Rx[1] - Rx[0])**2 + (Ry[1] - Ry[0])**2 ) # effective spacing.
      Ray = NP.array(list(zip( Rx, Ry )))
      # estimate the ray length on each element of the rectangular mesh
      Nelem = len(RMesh.ID)
      Le = NP.zeros(Nelem) 
      for k in RMesh.ID:
         # get coordinates of vertices for element k
         Vx = RMesh.vert['x'][k]
         Vy = RMesh.vert['y'][k]
         Vertices = NP.array(list(zip(Vx, Vy)))
         # identify the points of the ray that lie inside of the rectangular element.
         Rectangle = path.Path( Vertices )
         # in the following two lines, raius = 1000*eps is enough to consider points 
         # at the edges as inside!. 
         tol = 1000 * NP.finfo(float).eps
         inside = Rectangle.contains_points( Ray, radius = tol )
          
         # estimate the length of the ray inside the rectangular element. 
         NumInside = NP.sum(inside)
         Le[k] = NumInside * dR

      __comment = """

	tengo que verificar que la suma de todos los largos parciales sea igual a
	la longitud del rayo. Esto porque el metodo me da un artefacto cuando el
	rayo es horizontal y pasa justo por el limite entre dos filas de
	elementos de la discretizacion, por lo que la longitud del rayo obtenida
	de sumar los elementos parciales queda el doble. Para arreglar este
	problema simplemente imponemos que la suma de las longitudes parciales
	sea igual a la longitud original del rayo (lo cual no es la solucion mas
	elegante, para hacerlo bien habria que realizar un trazado de rayos mas
	elaborado pero escapa a las espectativas de este curso).

      """
      Lray_calculated = NP.sum(Le)
      if Lray_calculated < 1.0:
          print([ self.L, Lray_calculated, self.xS, self.yS, self.xR, self.yR])
      Le = Le * self.L / Lray_calculated 

      return [Le, dR] # dR is a measure of the precision on the length calculation
