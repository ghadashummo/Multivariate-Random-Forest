import numpy as np
def  ProcessingInputArrays(PrevMat,PPrevNameVectors,PrevExp,PrevMRNA,ll_idx,rr_idx,s):
                # print('s',s)
                 LeftNewMat=[];
                 RightNewMat=[];
                 NewPrevMat=[];
                 NewNameVectors=[];
                # print('ProcessingInputArrays(XY,MiRR,bbb,ccc,s):',ll_idx,rr_idx)
                 PrevNameVectors=np.array(PPrevNameVectors).tolist()
##                 if s !=' ' :
##                    print('sssss',s)
                 a=PrevNameVectors.index(s);
##                 #print('a:    ' ,a)
##                 
                 NewNameVectors=np.delete(PrevNameVectors,a);
                # NewNameVectors=PrevNameVectors;
##                 else:
##                    print('no mirna to delete',s)
##                    NewNameVectors=PrevNameVectors;
##                    return PrevMat,PPrevNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA

                     
                 NewPrevMat=np.delete(PrevMat,[a],1); #([a],1) means delete column a,([a],0):means delete row a: for row 0,for column 1
                 LeftNewExp=PrevExp[ll_idx,:]
                 RightNewExp=PrevExp[rr_idx,:]
                 LftMRNA=PrevMRNA[ll_idx,:]
                 RghtMRNA=PrevMRNA[rr_idx,:]
                 LeftNewMat=NewPrevMat[ll_idx,:];
                 RightNewMat=NewPrevMat[rr_idx,:];
                 #print('new mat')
                 #print(NewPrevMat)
                 #print('old mat ')
                 #print(PrevMat)
                 return NewPrevMat,NewNameVectors,LeftNewMat,RightNewMat,LeftNewExp,RightNewExp,LftMRNA,RghtMRNA
