__doc__ = """
GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

5 de Septiembre de 2017

Este módulo contiene funciones utilitarias que se pueden usar para operar 
sobre un objeto de Rectangular Mesh.
"""
import numpy as NP
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

######
def plotRectMeshField(F, mesh, fig, ax , colorlim = None, \
                  title ='', xlabel ='', ylabel ='',\
                  colormap = 'nipy_spectral', isColormapName = True, edgewidth = 0.2,\
                 edgecolor = 'k', colorbar = True, autoscale = True, \
                  colorbar_orientation = 'vertical', colorbar_label = 'Field',\
                 axequal = True, fontsize = 12, alpha = 0.85):
   
    # shink values of F into 1 dimension. 
    F = F.squeeze()

    # transform values of field into colors. 
    if isColormapName:
        cmap = plt.get_cmap(colormap)
    else: # is a colormap instance from plt.get_cmap()
        cmap = colormap
    Fcolors = cmap(F)
    # get info for mesh
    # assemble a list of polygons to plot (each polygon is the edge of a mesh element)
    patches = [Polygon( mesh.get_elem_vertices(k, shrinked = True), closed = True )\
               for k in mesh.ID]
    pc = PatchCollection(patches, cmap = cmap, alpha = alpha, edgecolor = edgecolor,\
                         lw = edgewidth )
    pc.set_array(NP.array(F))
    if colorlim is not None:
        pc.set_clim(colorlim)

    ax.add_collection(pc)

    if colorbar:
        cb = fig.colorbar(pc, ax = ax, \
                          orientation = colorbar_orientation )
        cb.set_label(colorbar_label, size = fontsize)

    FS = 0.5
    ax.set_xlim([mesh.Xmin - FS * mesh.dX, mesh.Xmax + FS * mesh.dX])
    ax.set_ylim([mesh.Ymin - FS * mesh.dY, mesh.Ymax + FS * mesh.dY])
    ax.set_title(title, fontsize = fontsize)
    ax.set_xlabel(xlabel, fontsize = fontsize)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    if axequal:
        plt.axis('equal')
    #fig.show()
    return [fig, ax]
######
def plotRectMeshIndices(mesh, fig, ax, element = True, title = '', linewidth = 1.0, \
		        linecolor = 'r', axequal = True, IndexFontSize = 8,\
                        xlabel = '', ylabel = '', fontsize = 12):
   """
   plot the mesh grid with indices. 
   if element == True plots the global ID (k = the number of the element)
   if element == False plots the array indices (i,j), where i is the index along the Y
                coordinate and j is the index along the X coordinates.

   """
   if element:
      # get a list of strings with the ID of each element
      str2plot = [str(mesh.ID[k]) for k in mesh.ID ]
   else:
      str2plot = [str(mesh.get_indices_of_element(k)) for k in mesh.ID]

   for k in mesh.ID:
      Xv = list(mesh.vert['x'][k])
      Xv.append(Xv[0]) # to close the rectangle
      Yv = list(mesh.vert['y'][k])
      Yv.append(Yv[0]) # to close the rectangle
      ax.plot(Xv, Yv, '-' + linecolor, linewidth = linewidth)
      ax.text(mesh.cent['x'][k], mesh.cent['y'][k], str2plot[k],\
              horizontalalignment='center', verticalalignment='center',\
              size = IndexFontSize) 
   ax.set_title(title, fontsize = fontsize)
   ax.set_xlabel(xlabel, fontsize = fontsize)
   ax.set_ylabel(ylabel, fontsize = fontsize)
   ax.spines["right"].set_visible(False)
   ax.spines["top"].set_visible(False)
   if axequal:
      plt.axis('equal')
