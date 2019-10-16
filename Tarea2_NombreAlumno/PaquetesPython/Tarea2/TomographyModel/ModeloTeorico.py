__doc__ = """

GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

20 de Agosto de 2017


"""
import numpy as NP
from ..Mesh import RectangularMesh
from .RayLinear2D import RayLinear2D

# la variable COMPLETAR es necesaria para poder importar el paquete de python con los
# codigos incompletos. Donde encuentre esta variable debe completar usted el codigo.
COMPLETAR = None

######
def define_Tomography_RectangularMesh(Xmin = 0, Xmax = 100, Nx = 4, \
                                      Zmin = -200, Zmax = 0, Nz = 8):
    """
    Defines the tomography mesh as a Rectangular Mesh with limits defined by 
    Xmin, Xmax, Zmin, Zmax. The size of the discretized elements are
    dx = (xmax - xmin)/Nx and dz = (zmax - zmin)/Nz
    """
    return RectangularMesh(Xmin = Xmin, Xmax = Xmax, Nx = Nx,\
                           Ymin = Zmin, Ymax = Zmax, Ny = Nz)

######
def formar_G(xF, zF, xR, zR, TV, mesh, fractional_precision = 1E-3):
    """
    Función que crea la matriz de diseño G del problema de tomografía lineal
    - xF, zF = coordenadas x,y (2 arrays de (Nobs,1)) de las fuentes para cada rayo
               del experimento de tomografía.  
    - xR, zR = coordenadas x,y (2 arrays de (Nobs,1)) de los receptores para cada rayo
               del experimento de tomografía. 
    - TV : travel time measurement of each ray
    - Nx, Nz : number of discretized elements of the mesh on each direction. 
    - mesh is a rectangular mesh object instance. 
    - fractional_precision : needed for numerical integration of the ray along 
         the mesh. The ray is subdivided in segments of length <= dR, where 
         dR = fractional_precision * length of shorter side of mesh element.

    function returns [G, d, rays]
     - G : design matrix of the direct (forward) problem. 
     - d : data vector of the inverse problem
     - rays: a list with RayLinear2D instances (same order as input data)

    """
    # define data vector
    Nd = len(TV)
    d = NP.zeros((Nd,1)) 
    d[:,0] = TV.squeeze()

    # Each row of G is the length of the ray on each element of the mesh. 
    # the size of G is :
    #   - # of rows = Number of rays measured in the experiment = Nd = len(TV) 
    #   - # of columns = Number of elements of the mesh: 
    # The columns are sorted by the order stablished in TV (sames as source and receiver
    # coordinates). 
    Npar = COMPLETAR
    # initialize G and ray list
    G = NP.zeros(( COMPLETAR, COMPLETAR))
    rays = []
    # fill G row by row
    for i in range(0, Nd):
        # initialize ray
        ray = RayLinear2D( xF[i], zF[i], xR[i], zR[i] ) 
        rays.append(ray)
        # calculate ray length vector (see ray member functions and examples notebook)
        Lray, dRray = COMPLETAR
        # fill the i-th row of G
        G[COMPLETAR, COMPLETAR] = COMPLETAR

    return [G, d, rays]    

