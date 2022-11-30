

import numpy as np
from numpy import matrix
from random import sample
import random
from reemovNestings import reemovNestings

def sampling( data,orig_map,x ):
    print('length of data',len(data))
    if len(data)<x:
       x=len(data)
    ddata = np.array(data).tolist();
    
    mattt=[];
    id_samp=[];
    samp_ind=[];
    indices=[];
    samp=[];
  
     #for i in range(x):
    while len(samp)<=x:
        sampple= random.choice(ddata);
        if sampple in samp: # sampling without replacement
           sampple= random.choice(ddata);
        else:   
           samp.append(sampple)
    #sampp=
    reemovNestings(samp)    
    r=len(samp)
   # print(samp,'samp len',r)
  # print('gff',sampple)
    id_samp=samp.sort();
   # print('id-samp',id_samp)
    sam=np.array(samp).tolist()
    #print(sam)
    for w in sam:
        #print('w       ',w)
        ind=ddata.index(w)
      
        indices.append(ind)
    #print('yaho zato',indices)    
    QS = np.matrix(orig_map);
    a, b = np.matrix(orig_map).shape
   # print('dimensions:',a,b)
   # print('indices',indices)
    mattt=QS[:, indices];
  #  print('samp:      ',samp,'matt:     ',mattt)
   
    return samp,mattt



