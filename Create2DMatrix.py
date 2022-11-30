#Use the following links to find a
#function creating 2D matrix rather than this function
#https://docs.scipy.org/doc/numpy/reference/generated/numpy.zeros.html
#numpy.zeros
import numpy as np
from numpy import matrix


def Create2DMatrix(m,n):
    Mat=[]   
    
    for i in range(0,m):
       Mat.append([])
    for i in range(0,m):
        for j in range(0,n):
            Mat[i].append(j)
            Mat[i][j]=0;
           
    return Mat


