# this defines Wh, H, and ho for Tikhonov regularization of order 0, 1. and 2. 
# Wh is the identity matrix, 
# ho is a vector of zeros, 
# H is either the identity, the gradient, or the laplacian operator (finite differences approximations).

import numpy as NP

# la variable COMPLETAR es necesaria para poder importar el paquete de python con los
# codigos incompletos. Donde encuentre esta variable debe completar usted el codigo.
COMPLETAR = None


def order_2(RMesh):
    """
    computes the laplacian of the Rectangular Mesh RMesh.
    return [H, Wh, ho] for second order Tikhonov. 
    
    """
    
    # get element ID's (each mesh element has a unique ID)
    IDs = RMesh.ID.copy()
    Npar = len(IDs)
    # the H operator is proportional to the laplacian (finite difference approximation).
    H = COMPLETAR
    # in the diagonal goes the number of neighbors of an element and the off diagonal
    # elements in a row have 0, except for the location of the neighbors that have -1.
    for elemID in IDs:
        elemNeighbors = RMesh.neighbors[elemID] 
        H[elemID, elemID] = COMPLETAR  #diagonal element.
        for neighborID in elemNeighbors:
            H[COMPLETAR, COMPLETAR] = COMPLETAR

            
    # caculate Wh and ho for Tikhonov of second order
    Wh = COMPLETAR
    ho = COMPLETAR
    # return the laplacian matrix H
    return [H, Wh, ho]
