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
######
def calcCheckerboardSynthetic(mesh, G, shift = [0,0], scale = 33.333):
    """
    Crea un modelo sintético de velocidades del medio del tipo tableto de
       ajedrez y utilizando este modelo de velocidades calcula una prediccion
       del tiempo de viaje de las ondas para la topologia de fuentes y
       receptores descritos en malla.

       return [dSyn, mSyn]
    """
    Npar = len(mesh.ID)
    backgroundVelocity = 4
    Vel = NP.ones((Npar,1)) * backgroundVelocity
    Xo = shift[0]
    Yo = shift[1]
    
    for k in range(0, Npar):
        Xc = mesh.cent['x'][k]
        Yc = mesh.cent['y'][k]
        Vel[k] = Vel[k] + 1.0* NP.sin( NP.pi * (Xc - Xo)/scale ) *\
                               NP.sin( NP.pi * (Yc - Yo)/scale ) 

    # convert velocity into slowness
    mSyn = 1/Vel
    # compute model prediction
    dSyn = G.dot(mSyn) 

    return [dSyn, mSyn] 
    
