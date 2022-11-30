import csv
import numpy as np
import os
from numpy import matrix
import pandas as pd
from pandas import DataFrame
from Euclidean_Distance import Euclidean_Distance
from sampling import sampling
from ProcessingInputArrays import ProcessingInputArrays
import pyexcel as pe  
import multiprocessing
import time
from random import randint
import threading
import multiprocessing
import numpy as np
from multiprocessing.managers import BaseManager
from multiprocessing import Process, Queue


#os.chdir('/lustrehome/ghadashommo/GhadaMRF/miRNA-DEexpression/')

wholeArray = pe.get_array(file_name='NewDEMirMapMat-Transpose.xlsx');
trainY = pe.get_array(file_name='AllMirListExpression.xlsx');
wholeFeatures = pe.get_array(file_name='NewDE-Targets.xlsx');
mrname = pe.get_array(file_name='AllMirList.xlsx');


XY = np.matrix(wholeArray);
a,b=np.array(XY).shape
print(a,b)
MiRR = np.array(wholeFeatures).tolist();
print(len(MiRR))
Exp = np.matrix(trainY);
mrna_name=np.array(mrname)
print(len(mrna_name))
node_size=10;
counter=[0]

#os.chdir('/lustrehome/ghadashommo/GhadaMRF/miRNA-DEBestTree/miRNA-DEBestTreeNode3/')     

    
class MRF:
     
    def __init__(self):
        pass

    
    def FirstSplitTree(self,XY,MiRR,Exp,mrna_name):
        counter[0]+=1
        uu,vv=np.matrix(XY).shape
       #  print('Remaining Mirnas',len(MiRR),'uu vv',uu,vv)
        numcov=85;
        NewPrevMat=[];
        NewNameVectors=[];
        LeftNewMat=[];
        RightNewMat=[];
        RTMR=[]
        LFMR=[]
        BRK=[];
        Mp=[];
        Dt=[];
        l_idx=[];
        r_idx=[];
        BESTmir=' ';
        BESTTree=[];
        w1=0;
        w2=0;
        EDTree=[];
        MaxEd=[];
        ll_idx=[];
        rr_idx=[];
        l_idx=[];
        r_idx=[];
        BESTmir=' ';
        BESTTree=[];
        LeftMRNA=[];
        RightMRNA=[];
        left_mir=[];
        right_mir=[];
        S=0
        s1=0
        s2=0
        s3=0
       
        Dt,Mp=sampling(MiRR,XY,numcov)
        #print('Mp',Mp)
        v1,v2=np.matrix(Mp).shape
       # print(' w1,w2', v1,v2)
        for xx in range(v2):
                    
            l_idx=[];
            r_idx=[];
            for yy in range(v1):
                if (Mp[yy,xx] == 0):
                  # print('Print Mp',Mp[1,1])
                   l_idx.append(yy)
                else:
                   r_idx.append(yy)
          
           
            Left=Exp[l_idx,:]
          #  print('Left',l_idx)
            
            Right=Exp[r_idx,:]
           # print('Right',r_idx)
            if (len(l_idx) and len(r_idx))>=node_size :
            #print('Riiiiiiiiiight',Right)
               s1=Euclidean_Distance(Exp)
               s2=Euclidean_Distance(Left)
               s3=Euclidean_Distance(Right)
               S=s1-s2-s3
            else:
                S=0
               
            EDTree.append(S)
           
            
            max_ed=max(EDTree) 
            max_ed_idx=EDTree.index(max_ed) #index from the sample
               #max_idx=max_ed_idx
            
                
            if xx is max_ed_idx:
               ll_idx=l_idx;
               rr_idx=r_idx;
               BESTTree=Mp[:,max_ed_idx]
               BESTmir=Dt[max_ed_idx]
               LeftMRNA=mrna_name[ll_idx];
               RightMRNA=mrna_name[rr_idx];
               LFMR=np.array(LeftMRNA).tolist()
               RTMR=np.array(RightMRNA).tolist()
               BRK.append('BREAK')
               
      
        NewPrevMat,NewNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA=ProcessingInputArrays(XY,MiRR,Exp,mrna_name,ll_idx,rr_idx,BESTmir)
        return  BESTmir,NewNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
    def FurtherSplit(self,XY,MiRR,Exp,mrna_name):
        counter[0]+=1
        uu,vv=np.matrix(XY).shape
       #  print('Remaining Mirnas',len(MiRR),'uu vv',uu,vv)
        numcov=85;
        NewPrevMat=[];
        NewNameVectors=[];
        LeftNewMat=[];
        RightNewMat=[];
        RTMR=[]
        LFMR=[]
        BRK=[];
        Mp=[];
        Dt=[];
        l_idx=[];
        r_idx=[];
        BESTmir=' ';
        BESTTree=[];
        w1=0;
        w2=0;
        EDTree=[];
        MaxEd=[];
        ll_idx=[];
        rr_idx=[];
        l_idx=[];
        r_idx=[];
        BESTmir=' ';
        BESTTree=[];
        LeftMRNA=[];
        RightMRNA=[];
        left_mir=[];
        right_mir=[];
        S=0
        s1=0
        s2=0
        s3=0
       
        Dt,Mp=sampling(MiRR,XY,numcov)
        #print('Mp',Mp)
        v1,v2=np.matrix(Mp).shape
       # print(' w1,w2', v1,v2)
        for xx in range(v2):
                    
            l_idx=[];
            r_idx=[];
            for yy in range(v1):
                if (Mp[yy,xx] == 0):
                  # print('Print Mp',Mp[1,1])
                   l_idx.append(yy)
                else:
                   r_idx.append(yy)
          
           
            Left=Exp[l_idx,:]
          #  print('Left',l_idx)
            
            Right=Exp[r_idx,:]
           # print('Right',r_idx)
            if (len(l_idx) and len(r_idx))>=node_size :
            #print('Riiiiiiiiiight',Right)
               s1=Euclidean_Distance(Exp)
               s2=Euclidean_Distance(Left)
               s3=Euclidean_Distance(Right)
               S=s1-s2-s3
            else:
                S=0
               
            EDTree.append(S)
           
            
            max_ed=max(EDTree) 
            max_ed_idx=EDTree.index(max_ed) #index from the sample
               #max_idx=max_ed_idx
            
                
            if xx is max_ed_idx:
               ll_idx=l_idx;
               rr_idx=r_idx;
               BESTTree=Mp[:,max_ed_idx]
               BESTmir=Dt[max_ed_idx]
               LeftMRNA=mrna_name[ll_idx];
               RightMRNA=mrna_name[rr_idx];
               LFMR=np.array(LeftMRNA).tolist()
               RTMR=np.array(RightMRNA).tolist()
               BRK.append('BREAK')
               
      
        NewPrevMat,NewNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA=ProcessingInputArrays(XY,MiRR,Exp,mrna_name,ll_idx,rr_idx,BESTmir)
      
        aa,bb=np.matrix(LeftNewMat).shape
         
        cc,dd=np.matrix(RightNewMat).shape
        print('aa  cc',aa,cc)
       
       # print('BESTmir',BESTmir)
        
                     
##        if   cc>2*node_size:
##             print('2: cc>2*node_size   ',LeftMRNA,RightMRNA) 
##            # NewNameVectors=Vector
##             RBT,RBM,RBL,RBR=MRF.BestTree(self,RightNewMat,NewNameVectors,RightNewExp,RghtMRNA)
        if  aa >2*node_size and  cc>2*node_size:
##            print('1: aa >2*node_size   ',LeftMRNA,RightMRNA)           
           # Vector=
            LBT,LBM,LBL,LBR,NewMirList=MRF.FurtherSplit(self,LeftNewMat,NewNameVectors,LeftNewExp,LftMRNA)
            RBT,RBM,RBL,RBR,NewMirList=MRF.FurtherSplit(self,RightNewMat,NewNameVectors,RightNewExp,RghtMRNA)
             
             # print(' RBT,RBM,RBL,RBR', RBT,RBM,RBL,RBR)
        elif (aa< node_size ) or (cc < node_size):
            
            #print('3: aa <node_size   ',LeftMRNA,RightMRNA)
##            if  (cc < node_size) or cc>=node_size:
##                print('4: cc<node_size   ',LeftMRNA,RightMRNA)
##            #or ((aa< node_size ) and (cc < node_size)):
                with open("%s_10-de-miRNAClst_in.csv" % item, 'a') as f:
                                writer = csv.writer(f)
##                                print('RTMR',len(RTMR))
##                                print('LFMR',len(LFMR))
                                RT=np.array(RightMRNA).tolist()
                                print('RT',len(RT))
                                LF=np.array(LeftMRNA).tolist()
                                print('LF',len(LF))
                                RTLF=sum([LF,RT],[])
                                #=itertools.chain(LF, RT)
                               
                                print('RTLF',len(RTLF))
                               # writer.writerow(LFMR)
                                writer.writerow(RTLF)
                                   
        elif (aa>=node_size ) and (cc>= node_size):
             if aa>2*node_size:
                LBT,LBM,LBL,LBR,NewMirList=MRF.FurtherSplit(self,LeftNewMat,NewNameVectors,LeftNewExp,LftMRNA)

             else:
                  with open("%s_10-de-miRNAClst_in.csv" % item, 'a') as f:
                       writer = csv.writer(f)
                       writer.writerow(LFMR)
                       
             if cc>2*node_size:
                 RBT,RBM,RBL,RBR,NewMirList=MRF.FurtherSplit(self,RightNewMat,NewNameVectors,RightNewExp,RghtMRNA)
             else:
                  with open("%s_10-de-miRNAClst_in.csv" % item, 'a') as f:
                       writer = csv.writer(f)
                       writer.writerow(RTMR)

        return BESTTree,BESTmir,LeftMRNA,RightMRNA,NewNameVectors


        file.close
PROCESSES =  multiprocessing.cpu_count() - 1
item =0;


def iteratedTree(item):
    mir=''
    NNV=[];
    LNM=[];
    RNM=[];
    LNE=[];
    RNE=[]
    LMRNA=[]
    RMRNA=[]
    mrf=MRF()
    C=XY;
    D=MiRR;
    LeBT=[]
    LeBM=[]
    LeLMRNA=[]
    LeRMRNA=[]
    LeNewVecT=[]
    

    #BESTmir,NewNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA
    mir,NNV ,LNM, RNM,LNE,RNE,LMRNA,RMRNA =mrf.FirstSplitTree(XY,MiRR,Exp,mrna_name)
    LeBT,LeBM,LeLMRNA,LeRMRNA,LeNewVecT =mrf.FurtherSplit(LNM,NNV,LNE,LMRNA)
    RiBT,RiBM,RiLMRNA,RiRMRNA,RiNewVecT =mrf.FurtherSplit(RNM,LeNewVecT,RNE,RMRNA)     
    #time.sleep(randint(2,4))
    print('Exiting worker', item)
    return "ok"

if __name__ == '__main__':
    start=time.time()

    pool = multiprocessing.Pool(processes=PROCESSES)
   # itm=input("enter order no of the first cluster:  ")
    #itm=1
    #input("enter order no of the first cluster:  ")
    #N=input("enter number of iterations:  ")
    N=101
    #print(N)
    #for item in range(int(itm),int(N)+int(itm)):
   
    for item in range(N):
    
        print('R    U   N')
        iteratedTree(item)
    pool_outputs = pool.map(iteratedTree,range(N))
    
    # pool_outputs = pool.map(iteratedTree,range(item))

    pool.close()
    pool.join()
    #iteratedTree()
    print( 'Pool:', pool_outputs)    
        
    start2=time.time()

    print('Time taken =    ',start2-start)

    
        


           



        
