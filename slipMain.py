# -*- coding: utf-8 -*-
"""
Tunnel Face Limit Equilbrium Calculation
mcan.wang@foxmail.com
"""
"""
Nmax:              Number of Tunnel Slip circle Segment
Mmax：             Number of Tunnel plane segment
Rs:                Tunnel slip circle radius at the ith section
Rnmax:             Tunnel radius
R0；               Tunnel silo radius 
theta:             Slip angle at the ith circle
alpha:             Rotation angle at the jth arch
x0,y0,z0           Coordinate of the rotation centre
xc,yc,zc           Coordinate of centre at eah circle
CircleCen          Slip rotation centre
xpos,ypos,zpos     x,y,z position in Cartesian coordinates
GridPos            slip surface point in Cartesian coordinate, [x,y,z,theta,alpha].
wedge[1],[2],[3]   Restores coordinates of the wedge, wedge[4] restores  the volume.
"""

import numpy as np
from slipInte import interFric
Nmax = 20
Mmax = 90
fric = 25
Rnmax = 4.15
R0 = Rnmax/np.tan(np.radians(55+fric/2))
# R0 = Rnmax/np.tan(np.radians(73))
x0 = 0
y0 = 0
z0 = Rnmax 
Dens = 15.7
Cdepth = 8.3*3
fricRed = interFric(fric)
fcoeff = fric/100
fcoeff = np.tan(np.radians(fric))*fricRed
# fcoeff = 0.5
kcoeff = 1-np.sin(np.radians(fric))
coh = 0
coh = coh*0.67
SDratio = 0.02
#---------------------------------------------------------------------------------------
'''
The main program contains three major modules
1. First module at 'slipGeom.py' file includs geometry initialize;
2. Second module at 'slipMechh.py' file describes mechanical behaviour;
3. Third module at 'slipArch.py' exports results for visualize.
'''
#---------------------------------------------------------------------------------------
from slipGeom import slipSurface
from slipGeom import wedgeNode
from slipGeom import wedgePos
from slipGeom import melonPos
from slipGeom import surfNorm
from slipGeom import pjtPtPos
from slipGeom import pjtPlPos

from slipMech import gravityStress
from slipMech import totalStress
from slipMech import normStress
from slipMech import slipForce

from slipArch import upperStress
'''
Upper Stress surcharged on the crest of the slip torus
q2 is assessed based on Chen's [2019] model for further reference
'''
q2 = upperStress(Rnmax,fric,Dens,Cdepth,coh,SDratio)
P0 = [0,0,-q2] # Needs to Specify the direction of the P0 #Dens = 24

gridPos   = slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[0]
CircleCen = slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[4]
mesh      = wedgeNode(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[0]
wedge     = wedgePos(Nmax,Mmax,Rnmax,R0,x0,y0,z0)
melon     = melonPos(Nmax,Mmax,wedge)
snvec     = surfNorm(Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh)  
ptPjt     = np.nan_to_num(pjtPtPos(Mmax,Nmax,wedge,CircleCen,melon,mesh,snvec)) 
Ghpl      = pjtPlPos(Mmax,Nmax,melon,ptPjt,mesh)  

#--------------Main Loop for Vertical Force Equilibrium--------------------------------------------------------
'''
Call gravityStress() to initialize equilivalent Stress at wedge gravity plane
Call normStress() and slipForce() to calculate respective normal and shear forces
Calculate unbalance force for the system under gravity and reaction loads
Apply unbalance force to the system to calculate increment of reacton loads
Loop reaction loads increment until vertical equilibrium reached
IncFactor used to multiply unbalance force to apply load incrementlly, IncFactor <= 1.

'''
nstep = 20
unBalF = np.zeros([nstep,10])
eqF = np.zeros([nstep+1,2])
eqF[0,0] = 1
GStress,GArea,GVector,lcVector,SArea = gravityStress(Mmax,Nmax,mesh,melon,Ghpl,ptPjt,wedge,Dens,R0,Rnmax,x0,y0,z0)
eqStress = GStress
IncFactor = 0.5 

'''
vertical unbalance force = Gravity load + vertical surcharge + reaction support (Fn1 + Fn2) + friction force
unBalF[istep,1] records unbalance force at Z direction (Vertical)
unBalF[istep,2] records unbalance force at X direction (Symmetric)  
unBalF[istep,2] records unbalance ratio at Z direction (Vertical) 
unbalYstress records sum of horizontal force at y direction
normLSP is the normalized limit support pressure
'''
Fs10 = np.zeros([Mmax,nstep])
ts10 = np.zeros([Mmax,nstep])
Fsmid = np.zeros([Nmax,nstep])
for istep in range(0, nstep):
    #Stress     = normStress(Nmax,Mmax,P0,eqStress,GVector,mesh,GArea)
    Stress = normStress(Nmax,Mmax,P0,eqStress,lcVector,mesh,GArea,SArea)
    #Fn1,Fn3,Fs = slipForce(Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh,GVector,Stress,ptPjt,kcoeff,fcoeff,GArea)
    Fn1,Fn3,Fs,ts,SArea1,SArea2,Fn1_Norm,Fn3_Norm,ns = slipForce(Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh,lcVector,Stress,ptPjt,kcoeff,fcoeff,GArea,coh,SArea,GStress,Dens)
    unBalF[istep,0] = istep
    unBalF[istep,1] = -wedge[:,4].sum()*Dens + P0[2]*np.pi*R0**2 + \
                      -Fn1[:,2].sum() + -Fn3[:,2].sum() + Fs[:,2].sum()  
    unBalF[istep,2] = -Fn1[:,0].sum() + -Fn3[:,0].sum() + Fs[:,0].sum()         
    unBalF[istep,3] = unBalF[istep,1] / (wedge[:,4].sum()*Dens)  
    unBalF[istep,4] = unBalF[istep,2] / (wedge[:,4].sum()*Dens)    
    unbalFactor = IncFactor*unBalF[istep,3] 
    eqF[istep+1,0] = eqF[istep,0]*(1 - unbalFactor)
    eqStress,eqGstress,incStress = totalStress(Nmax,Mmax,istep,unbalFactor,eqF,GStress)
    Fs10[:,istep] = (Fs[9*Mmax:10*Mmax,0]**2+Fs[9*Mmax:10*Mmax,1]**2+Fs[9*Mmax:10*Mmax,2]**2)**0.5
    ts10[:,istep] = (ts[9*Mmax:10*Mmax,0]**2+ts[9*Mmax:10*Mmax,1]**2+ts[9*Mmax:10*Mmax,2]**2)**0.5
    for j in range(0,Nmax):
        Fsmid[j,istep] = (Fs[int(j*Mmax+Mmax/2),0]**2+Fs[int(j*Mmax+Mmax/2),1]**2+Fs[int(j*Mmax+Mmax/2),2]**2)**0.5

yFn = Fn1[:,1].sum()+Fn3[:,1].sum()
yFs = Fs[:,1].sum()
#Find the last surcharge stress
hLoad = np.zeros([Mmax,2])
for j in range (0, Mmax):
    hLoad[j,0] = hLoad[j,0]+(Stress[0][(Nmax-1)*Mmax+j][1])*GArea[(Nmax-1)*Mmax,1]
    
a=hLoad.sum()/(np.pi*Rnmax**2)/Dens/2/Rnmax
unbalYstress1 = (-Fn1[:,1].sum() - Fn3[:,1].sum() + Fs[:,1].sum()+hLoad.sum())/(np.pi*Rnmax**2)
unbalYstress = (-Fn1[:,1].sum() - Fn3[:,1].sum() + Fs[:,1].sum())/(np.pi*Rnmax**2)
normLSP = unbalYstress/(Dens*Rnmax*2)
sigmaH = (Cdepth+Rnmax)*kcoeff*Dens
normStr = unbalYstress/sigmaH

print('NormLSP = {}'.format(normLSP))
print('UnbalRatio_Z = {}'.format(unBalF[istep,3]))
print(-yFn,yFs)
#---------------End of Main Loop for Force Equilbrium Calculation------------------------------------------------

#Export files after the mainloop
from slipExpo import slipPlot
slipPlot(Nmax,Mmax,gridPos,CircleCen)
from slipExpo import sForceVector
sForceVector(Nmax,Mmax,ptPjt,Fs)