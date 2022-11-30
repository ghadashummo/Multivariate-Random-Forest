# This function calculates the  Euclidean Distance of a group of elements
# using  the matrix that contains the values of the elements
# represented by the rows

import numpy as np
from numpy import matrix


def Euclidean_Distance(ED_array):
    # ED_array = pe.get_array(file_name='eulidean test.xlsx');
    r, c = np.matrix(ED_array).shape
    WA = np.matrix(ED_array).astype(np.float);
    sum = 0;
    j = 0;
    i = 0;
    for j in range(c):
        A = WA[:, j];
        avrg = np.mean(A);
        for i in range(r):
            y = avrg - float(WA[i, j]);
            y = y * y;
            sum = sum + y;
    #print(sum);
    return sum;


#EDA = pe.get_array(file_name='eulidean test.xlsx');
#Euclidean_Distance(ED_array=EDA)
