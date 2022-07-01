# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:14:59 2021

@author: hdec
"""
import numpy as np
def slipPlot(Nmax,Mmax,pos,cpos):
    #plot slip surface for tunnel face
    #Spoint and Epoint refers the start and end points in rhino linePlot
    Nmax1 = Nmax+1
    NumofSeg = Nmax1*Mmax
    Spoint1 = np.zeros([NumofSeg,3])
    Epoint1 = np.zeros([NumofSeg,3])
    Spoint2 = np.zeros([(Nmax1-1)*(Mmax-1),3])
    Epoint2 = np.zeros([(Nmax1-1)*(Mmax-1),3])
    
    for i in range (0,NumofSeg):
        Spoint1[i,] = pos[i,0:3]
        if i < NumofSeg - 1:
            Epoint1[i,] = pos[i+1,0:3]
        else:
            Epoint1[i,] = pos[0,0:3]

    inc = 0
    for i in range(0,Nmax1-1):
        for j in range(0,Mmax):
            k = i*Mmax + j
            if k % Mmax == 0:
                pass
            else:
                Spoint2[inc,0] = pos[k,0]
                Spoint2[inc,1] = pos[k,1]
                Spoint2[inc,2] = pos[k,2]
                inc = inc + 1
  
    inc = 0
    for i in range(1,Nmax1):
        for j in range(0,Mmax):
            k = i*Mmax + j
            if k % Mmax == 0:
                pass
            else:
                Epoint2[inc,0] = pos[k,0]
                Epoint2[inc,1] = pos[k,1]
                Epoint2[inc,2] = pos[k,2]  
                inc = inc + 1
            
    Spoint = np.append(Spoint1,Spoint2,axis = 0)
    Epoint = np.append(Epoint1,Epoint2,axis = 0)
    np.savetxt("Spoint.csv",Spoint,fmt = '%f',delimiter = ',')       
    np.savetxt("Epoint.csv",Epoint,fmt = '%f',delimiter = ',') 
    print('Export slip surface geometry completed')
    return()
    
def sForceVector(Nmax,Mmax,ptPjt,Fs):
    #plot shear force vectors calculated in slipMain
    #Point1 & Point2 refers to shear force vector at tunnel back
    #Point3 & Point4 refers to shear force vector at tunnel belly
    Scale1,Scale2 = [0.5,0.25]
    inc1,inc2 = [0,0]
    Pointer1,Pointer2 = [np.zeros([Nmax*Mmax,3]),np.zeros([Nmax*Mmax,3])]
    Pointer3,Pointer4 = [np.zeros([Nmax*Mmax,3]),np.zeros([Nmax*Mmax,3])]

    for i in range(0,Nmax):
        for j in range(0,Mmax):
             k = i*Mmax + j
             #Pointer0[k,0],Pointer0[k,1],Pointer0[k,2] = [ptPjt[k,0],ptPjt[k,1],ptPjt[k,2]]
             if ptPjt[k,4]*ptPjt[k,7]<0:
                Pointer1[inc1,0],Pointer1[inc1,1],Pointer1[inc1,2] = [ptPjt[k,0],ptPjt[k,1],ptPjt[k,2]]
                Pointer2[inc1,:] = Pointer1[inc1,:]+Fs[k,:]*Scale1
                inc1 = inc1 + 1
             else:
                Pointer3[inc2,0],Pointer3[inc2,1],Pointer3[inc2,2] = [ptPjt[k,0],ptPjt[k,1],ptPjt[k,2]]   
                Pointer4[inc2,:] = Pointer3[inc2,:]+Fs[k,:]*Scale2
                inc2 = inc2 + 1
    
    Pointer1,Pointer2 = [Pointer1[0:inc1,],Pointer2[0:inc1,]]  
    Pointer3,Pointer4 = [Pointer3[0:inc2,],Pointer4[0:inc2,]]
    
    np.savetxt("Pointer1.csv",Pointer1,fmt = '%f',delimiter = ',') 
    np.savetxt("Pointer2.csv",Pointer2,fmt = '%f',delimiter = ',') 
    np.savetxt("Pointer3.csv",Pointer3,fmt = '%f',delimiter = ',') 
    np.savetxt("Pointer4.csv",Pointer4,fmt = '%f',delimiter = ',')  
    print('Export shear force vector completed!')
    return()

def GravityPlane(Nmax,Mmax,wedge):
    for i in range(0,Nmax):
        np.savetxt("C"+"%d"%(i)+"Point.csv",wedge[i*Mmax:(i+1)*Mmax,1:4],fmt = '%f',delimiter = ',')
    return()