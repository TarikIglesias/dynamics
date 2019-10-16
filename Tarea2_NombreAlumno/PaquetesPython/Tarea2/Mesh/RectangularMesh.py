__doc__ = """
GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

5 de Septiembre de 2017


Este módulo define la topología de la malla a usar para discretizar el problema de 
tomografía. 

"""
import numpy as NP


######
class RectangularMesh(object):

   ######
   def __init__(self, Xmin = 0, Xmax = 100, Nx = 10, Ymin = 0, Ymax = 100, Ny = 10):
      """
      Initializes a rectangular mesh that ranges at [Xmin, Xmax] and [Ymin, Ymax]
      with Nx and Ny elements in the X and Y direction respectively. 

      The RectangularMesh instantiated object has information about the coordinates
      of the center of each element,
 
      """
      self.Xmin = Xmin
      self.Xmax = Xmax
      self.Nx   = Nx
      self.Ymin = Ymin
      self.Ymax = Ymax
      self.Ny   = Ny
      self.Xedge = None
      self.Yedge = None
      self.Xc = None
      self.Yc = None
      self.dX = None
      self.dY = None
      self.ID   = [] # identifier of each rectangular mesh element 
      self.cent = {'x': [], 'y' : []} # coordinates of the center of each rectangular mesh element
      self.vert = {'x': [], 'y' : []} # coordinates of the vertices of each element.  
      # each element in self.vert['x'|'y'] list, is a tuple with 4 elements
      self.neighbors = {} # each element has a tuple with the IDs of the neighbors.
                          # keys are the ID of the mesh elements
      # indices of arrays in mesh: i for Y, j for X and k for the mesh elements 
      self.arrayIndices = {'i':[], 'j':[], 'k':[]}
      self._initializeMesh()
   ######
   def _initializeMesh(self):
      """
      Initializes mesh topology, i.e., calculates vertices, centers coordinates and 
      determine the neighbors of each rectangular element. 
      """
      # create arrays with coordinates of element borders in each direction.
      self.Xedge = NP.linspace(self.Xmin, self.Xmax, self.Nx + 1)
      self.Yedge = NP.linspace(self.Ymin, self.Ymax, self.Ny + 1)

      # calculate the coordinates of the centers of each rectangular element. 
      self.Xc = 0.5 * ( self.Xedge[1:] + self.Xedge[:-1] )
      self.Yc = 0.5 * ( self.Yedge[1:] + self.Yedge[:-1] )

      # define step of mesh in each direction
      self.dX = self.Xc[1] - self.Xc[0]
      self.dY = self.Yc[1] - self.Yc[0]
 
      # Set mesh elements connectivity tables self.neighbors
      # the mesh elements are numbered in increasing order of X coordinates and then
      # increasing order of Y coordinates. 
      # Then if i indicates the row (constant value of Y) and j indicates the 
      # column (constant value of X) the k-th mesh element number is given by
      # the following formula k = i * Nx + j where i in range(n, Ny) and
      # j in range(0, Nx) 

      for i in range(0, self.Ny):
         for j in range(0, self.Nx): 
            k = self.get_element_index(i,j)
            self.ID.append(k)
            self.arrayIndices['i'].append(i)
            self.arrayIndices['j'].append(j)
            self.arrayIndices['k'].append(k)
            # determine the number of neighbors
            NumNeighbors = 0
            self.neighbors[k] = []
            # check for neighbors above and below current mesh element k
            if i == 0:
               NumNeighbors += 1
               self.neighbors[k].append( k + self.Nx  ) # only neighbor above
            elif i == self.Ny - 1:
               NumNeighbors += 1
               self.neighbors[k].append( k - self.Nx  ) # only neighbor below
            else: # has neighbors above and below
               NumNeighbors += 1
               self.neighbors[k].append( k + self.Nx  ) # neighbor above
               NumNeighbors += 1
               self.neighbors[k].append( k - self.Nx  ) # neighbor below
            # check for neighbors to the right and left sides
            if j == 0:
               NumNeighbors += 1
               self.neighbors[k].append( k + 1  ) # only neighbor to the right
            elif j == self.Nx - 1:
               NumNeighbors += 1
               self.neighbors[k].append( k - 1  ) # only neighbor to the left
            else: # has neighbors on both sides
               NumNeighbors += 1
               self.neighbors[k].append( k + 1  ) # neighbor to the right
               NumNeighbors += 1
               self.neighbors[k].append( k - 1  ) # neighbor to the left

            self.neighbors[k] = tuple(self.neighbors[k]) # so it can not be modified

            # fill tables with center and vertex coordinates. 
            # centers of rectangular mesh elements
            self.cent['x'].append( self.Xc[j] )
            self.cent['y'].append( self.Yc[i] )
            # vertices of mesh elements
            V = {}
            V['1'] = [self.Xc[j] + self.dX/2.0, self.Yc[i] + self.dY/2.0]
            V['2'] = [self.Xc[j] - self.dX/2.0, self.Yc[i] + self.dY/2.0]
            V['3'] = [self.Xc[j] - self.dX/2.0, self.Yc[i] - self.dY/2.0]
            V['4'] = [self.Xc[j] + self.dX/2.0, self.Yc[i] - self.dY/2.0]
            Vx = NP.array( [ V[vertexID][0] for vertexID in ['1', '2', '3', '4'] ] )
            Vy = NP.array( [ V[vertexID][1] for vertexID in ['1', '2', '3', '4'] ] )
            self.vert['x'].append(Vx)
            self.vert['y'].append(Vy)
      
   ######         
   def get_element_index(self, i, j):
      """
      i indicates the row (constant value of Y) and j indicates the 
      column (constant value of X) the k-th mesh element number is given by
      the following formula k = i * Nx + j where i in range(n, Ny) and
      j in range(0, Nx) 
      """
      k = i * self.Nx + j
      return k

   ######
   def get_indices_of_element(self, k):
      """
      return the array indices of element with  ID = k 
      """
      i = self.arrayIndices['i'][k]
      j = self.arrayIndices['j'][k]
      return [i, j]

   ######
   def get_elem_vertices(self, k, shrinked = True):
      """
      Get the vertices of the elemet k.
      if shrinked == True:
         returns a matrix of 4 by 2 in which each row correspond to the x and y
         coordinates of the corresponding vertex
         [[x1, y1]
          [x2, y2]
          [x3, y3]
          [x4, y4]]
      if shrinked == False:
         returns a dictionary with the coordinates of the 4 vertex, each value is 
         a numpy array with the x and y coordinate of the corresponding vertex. 
         [ NP.array([x1, y1]), 
           NP.array([x2, y2]), 
           NP.array([x3, y3]), 
           NP.array([x4, y4]) ]
      """
      Vx = self.vert['x'][k]
      Vy = self.vert['y'][k]
      
      if shrinked:
         V = NP.zeros( (4, 2) )
         V[:,0] = Vx
         V[:,1] = Vy
          
      else:
         x1, x2, x3, x4 = Vx
         y1, y2, y3, y4 = Vy
         V = []
         V.append( NP.array([x1, y1]) )
         V.append( NP.array([x2, y2]) )
         V.append( NP.array([x3, y3]) )
         V.append( NP.array([x4, y4]) )

      return V
