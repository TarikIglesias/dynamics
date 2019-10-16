__doc__ = """
GF5013 - Metodos Inversos Aplicados a la Geofisica
Primavera 2017
Prof. Francisco Hernán Ortega Culaciati
ortega.francisco@u.uchile.cl
Departamento de Geofísica - FCFM
Universidad de Chile

02 de Octubre de 2017

"""
from numpy.linalg import lstsq 
import numpy as NP
######
def MinCuadSimple(G, d):
    """
    Calcula la solución de G*m = d por el método de mínimos cuadrados simples. 
    Devuelve un diccionario con m y Cm. 
    Para el cálculo de la matriz de covarianza del modelo, se asume que los 
    errores del ajuste (i.e., de la diferencia d-Gm) son i.i.d. con media nula 
    y varianza unitaria.  
    """
    Ndata, Npar = G.shape
    Cm = NP.linalg.lstsq( G.T.dot(G), NP.eye(Npar), rcond = -1 )[0]
    m = Cm.dot( G.T.dot(d) )
    return {'m' : m, 'Cm': Cm}


######
def MinCuadPesos(G, d, Wx):
    """
    Calcula la solución de Wx*G*m = Wx*d por el método de mínimos cuadrados (con pesos). 
    Devuelve un diccionario con m y Cm. 
    Para el cálculo de la matriz de covarianza de los parámetros estimados del modelo, 
    se asume que los errores del ajuste (i.e., de la diferencia d-Gm) siguen una 
    distribución normal multivariada con media nula y matriz de covarianza Cx.
    donde Cx = inv( Wx.T.dot(Wx) ) 
    """
    return MinCuadSimple( Wx.dot(G), Wx.dot(d) )


######
def MinCuadCov(G, d, Cx):
    """
    inv(Cx) = Wx.T.dot(Wx)
    Calcula la solución de Wx*G*m = Wx*d por el método de mínimos cuadrados (con pesos). 
    Devuelve un diccionario con m y Cm. 
    Para el cálculo de la matriz de covarianza de los parámetros estimados del modelo, 
    se asume que los errores del ajuste (i.e., de la diferencia d-Gm) siguen una 
    distribución normal multivariada con media nula y matriz de covarianza Cx. 
    """
    Ndata, Ndata2 = Cx.shape
    inv_Cx = NP.linalg.lstsq( Cx , NP.eye(Ndata) )[0]
    Wx = NP.linalg.cholesky( inv_Cx  ).T # as python calculates A.dot(A.T) = inv_Cx    
 
    return  MinCuadPesos(G, d, Wx)

