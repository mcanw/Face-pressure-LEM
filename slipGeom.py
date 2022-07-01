# -*- coding: utf-8 -*-
"""
Created on Thu May 27 09:14:59 2021
@author: CanWang
1. Save all geometry related subFunctions in slipGeom.py file
2. 
"""
import numpy as np
def slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0):
    Tmode = 2.5
    #This function was used to generate 3D rotational surface;
    Nmax1 = Nmax+1
    CircleCen = np.zeros([Nmax1,3])
    GridPos = np.zeros([Mmax*Nmax1,5])

    xpos = np.zeros([Mmax,Nmax1])
    ypos = np.zeros([Mmax,Nmax1])
    zpos = np.zeros([Mmax,Nmax1])
    Rs = np.zeros([Nmax1,1])
    theta = np.zeros([Nmax1,1])
    alpha = np.zeros([Mmax,1])
    xc = np.zeros([Nmax1,1])
    yc = np.zeros([Nmax1,1])
    zc = np.zeros([Nmax1,1])
    ppos = np.zeros([Mmax*Nmax1,3])
    xp = np.zeros([Mmax,Nmax1])
    yp = np.zeros([Mmax,Nmax1])
    zp = np.zeros([Mmax,Nmax1])
    inc = 0

    for i in range (0,Nmax1): # i step is the one at each rotation angle
        if i < Nmax:
            theta[i,] = np.pi/2*i/Nmax
            Rs[i,] = R0+4*(Rnmax-R0)/np.pi**2*theta[i,]**2
            xc[i,] = x0
            yc[i,] = y0 - np.cos(i*np.pi/2/Nmax)*(R0+(Rnmax - R0)*(i/Nmax)**Tmode)
            zc[i,] = z0 - np.sin(i*np.pi/2/Nmax)*(R0+(Rnmax - R0)*(i/Nmax)**Tmode)
        else:
            Rs[i,] = Rnmax
            theta[i,] = np.pi/2
            xc[i,] = 0
            yc[i,] = 0
            zc[i,] = 0
        CircleCen[i,0] = xc[i,]
        CircleCen[i,1] = yc[i,]
        CircleCen[i,2] = zc[i,]
        for j in range (0,Mmax):# j step is the one at each circle plane
            alpha[j,] = 2*np.pi*j/Mmax
            xp[j,i] = xc[i,]+Rs[i,]*np.sin(alpha[j,])
            yp[j,i] = yc[i,]+Rs[i,]*np.cos(alpha[j,])
            zp[j,i] = zc[i,]       
            xpos[j,i] = x0 + np.sin(2*np.pi*j/Mmax)*(R0 + (Rnmax - R0)*(i/Nmax)**Tmode)
            ypos[j,i] = y0 - np.cos(np.pi*i/2/Nmax)*(1-np.cos(2*np.pi*j/Mmax))*(R0+(Rnmax-R0)*(i/Nmax)**Tmode)
            zpos[j,i] = z0 - np.sin(np.pi*i/2/Nmax)*(1-np.cos(2*np.pi*j/Mmax))*(R0+(Rnmax-R0)*(i/Nmax)**Tmode)
            ppos[inc,0] = xp[j,i]
            ppos[inc,1] = yp[j,i]
            ppos[inc,2] = zp[j,i]
            GridPos[inc,0] = xpos[j,i]
            GridPos[inc,1] = ypos[j,i]
            GridPos[inc,2] = zpos[j,i]
            GridPos[inc,3] = theta[i,]
            GridPos[inc,4] = alpha[j,]
            inc = inc + 1
    return(GridPos,xpos,ypos,zpos,CircleCen)
    
def wedgeNode(Nmax,Mmax,Rnmax,R0,x0,y0,z0):
    #call subFunction slipSurface
    #global meshx,meshy,meshz
    '''
    Prepare sequence of all nodes in a wedge defined by 123456.
    Rules of Nodes preparation in the mesh connection:
    1. Node 1,2,3,4 are in the primeter. node 5 and 6 are in the centre;
    2. Node 1 to 2 and node 3 to 4 are clockwise;
    3. Node 1 is above 3. node 2 is above 4;
    4. Node 5 is above 6;
    '''
    CircleCen = slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[4]
    meshx = np.zeros([Nmax*Mmax,7])  
    meshy = np.zeros([Nmax*Mmax,7])  
    meshz = np.zeros([Nmax*Mmax,7]) 
    mesh =  [meshx,meshy,meshz]
    mpos = [slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[1],\
            slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[2],\
            slipSurface(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[3]]
    mcpos = [CircleCen[:,0],CircleCen[:,1],CircleCen[:,2]]
    #global meshid
    for i in range (0,Nmax):
        for j in range(0,Mmax):  
            for k in range (0,3):
                inc = i*Mmax+j
                if j < Mmax-1:
                    mesh[k][inc,0] = inc
                    mesh[k][inc,1] = mpos[k][j,i]
                    mesh[k][inc,2] = mpos[k][j+1,i]
                    mesh[k][inc,3] = mpos[k][j,i+1]
                    mesh[k][inc,4] = mpos[k][j+1,i+1]
                    mesh[k][inc,5] = mcpos[k][i,]
                    mesh[k][inc,6] = mcpos[k][i+1,]
                
                else:
                    mesh[k][inc,0] = inc
                    mesh[k][inc,1] = mpos[k][j,i]
                    mesh[k][inc,2] = mpos[k][0,i]
                    mesh[k][inc,3] = mpos[k][j,i+1]
                    mesh[k][inc,4] = mpos[k][0,i+1]
                    mesh[k][inc,5] = mcpos[k][i,]
                    mesh[k][inc,6] = mcpos[k][i+1,]                                     
    return(mesh,mcpos)
    
def det(a,b,c): #Determinte of a 3*3 matrix
    det = a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
    return(det)
    
def volTetra(node1,node2,node3,node4,node5,node6):
    '''
    Calculate volume of a wedge;
    A wedge is divided by 3 Tetrahedrons of 1-2-3-5; 2-3-4-6; 2-3-5-6;
    Volume of each tetrahedron is calculate by determinate
    '''
    'Tetrahedron1 of 1-2-3-5'
    a1 = [node2[0]-node1[0],node2[1]-node1[1],node2[2]-node1[2]]
    b1 = [node3[0]-node1[0],node3[1]-node1[1],node3[2]-node1[2]]
    c1 = [node5[0]-node1[0],node5[1]-node1[1],node5[2]-node1[2]]
    vol1 = 1/6*abs(det(a1,b1,c1))    

    'Tetrahedron2 of 2-3-4-6'
    a2 = [node3[0]-node2[0],node3[1]-node2[1],node3[2]-node2[2]]
    b2 = [node4[0]-node2[0],node4[1]-node2[1],node4[2]-node2[2]]
    c2 = [node6[0]-node2[0],node6[1]-node2[1],node6[2]-node2[2]]
    vol2 = 1/6*abs(det(a2,b2,c2))
    
    'Tetrahedron3 of 2-3-5-6'
    a3 = [node3[0]-node2[0],node3[1]-node2[1],node3[2]-node2[2]]
    b3 = [node5[0]-node2[0],node5[1]-node2[1],node5[2]-node2[2]]
    c3 = [node6[0]-node2[0],node6[1]-node2[1],node6[2]-node2[2]]
    vol3 = 1/6*abs(det(a3,b3,c3)) 
    
    'Total volume composed of 3 Tetras'
    volume = vol1 + vol2 + vol3
    return(vol1, vol2, vol3, volume)
    
def cenTetra(node1,node2,node3,node4,node5,node6):
    #find center position of 3 tetras in the wedge
    'Tetrahedron1: centre of 1-2-3-5'
    x1 = 1/4*(node1[0]+node2[0]+node3[0]+node5[0])
    y1 = 1/4*(node1[1]+node2[1]+node3[1]+node5[1])
    z1 = 1/4*(node1[2]+node2[2]+node3[2]+node5[2])
    
    'Tetrahedron2: centre of 2-3-4-6'
    x2 = 1/4*(node2[0]+node3[0]+node4[0]+node6[0])
    y2 = 1/4*(node2[1]+node3[1]+node4[1]+node6[1])
    z2 = 1/4*(node2[2]+node3[2]+node4[2]+node6[2])
    
    'Tetrahedron3: centre of 2-3-5-6'
    x3 = 1/4*(node2[0]+node3[0]+node5[0]+node6[0])
    y3 = 1/4*(node2[1]+node3[1]+node5[1]+node6[1])
    z3 = 1/4*(node2[2]+node3[2]+node5[2]+node6[2])   
    return(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    
def wedgePos(Nmax,Mmax,Rnmax,R0,x0,y0,z0):
    #find center position of the wedge
    '''
    call subFunction wedgeNode
    call subFunction volTetra
    call subFunction cenTetra
    '''
    TVol = np.zeros([Nmax*Mmax,10])  
    TCen = np.zeros([Nmax*Mmax,10])   
    wedge = np.zeros([Nmax*Mmax,10]) #represents centre of wedge
    mesh = wedgeNode(Nmax,Mmax,Rnmax,R0,x0,y0,z0)[0]
    meshx = mesh[0]
    meshy = mesh[1]
    meshz = mesh[2]
    for i in range(0, Nmax):
        for j in range(0, Mmax):
            k = i*Mmax+j
            '''
            Tetrahedron 1 of points 1-2-3-5
            Tetrahedron 2 of points 2-3-4-6
            Tetrahedron 3 of points 2-3-5-6
        
            Volume of the Tetrahedron takes symmertric shapes
            Note if  k < Mmax/2, Node order is 1-2-3-4-5-6
                 if  k >= Mmax/2, Node order is 2-1-4-3-5-6.
                     The volume becomes symmetric.
            '''
            node1 = [meshx[k,1],meshy[k,1],meshz[k,1]]
            node2 = [meshx[k,2],meshy[k,2],meshz[k,2]]
            node3 = [meshx[k,3],meshy[k,3],meshz[k,3]]
            node4 = [meshx[k,4],meshy[k,4],meshz[k,4]]
            node5 = [meshx[k,5],meshy[k,5],meshz[k,5]]
            node6 = [meshx[k,6],meshy[k,6],meshz[k,6]]
    
            TVol[k,0] = k
        
            for m in range(0,4):
                if j < Mmax/2:
                    TVol[k,m+1] = volTetra(node1,node2,node3,node4,node5,node6)[m]
                else:
                    TVol[k,m+1] = volTetra(node2,node1,node4,node3,node5,node6)[m]
        
            TCen[k,0] = k
            
            for m in range(0,9):
                if j < Mmax/2:
                    TCen[k,m+1] = cenTetra(node1,node2,node3,node4,node5,node6)[m]
                else:
                    TCen[k,m+1] = cenTetra(node2,node1,node4,node3,node5,node6)[m]
    
            wedge[k,0] = k
            wedge[k,1] = (TVol[k,1]*TCen[k,1] + TVol[k,2]*TCen[k,4] + \
                         TVol[k,3]*TCen[k,7])/TVol[k,4]
            wedge[k,2] = (TVol[k,1]*TCen[k,2] + TVol[k,2]*TCen[k,5] + \
                         TVol[k,3]*TCen[k,8])/TVol[k,4]
            wedge[k,3] = (TVol[k,1]*TCen[k,3] + TVol[k,2]*TCen[k,6] + \
                         TVol[k,3]*TCen[k,9])/TVol[k,4]
            wedge[k,4] = TVol[k,4]
    
    return(wedge)
    
def melonPos(Nmax,Mmax,wedge):
    #find centre position of the melon
    melon = np.zeros([Nmax,5])
    for i in range (0,Nmax):
        xc = 0
        yc = 0
        zc = 0
        vol = 0
        for j in range (0,Mmax):
            vol = vol + wedge[i*Mmax+j,4]
        
        for j in range (0,Mmax):
            xc = xc + wedge[i*Mmax+j,1]*wedge[i*Mmax+j,4]/vol
            yc = yc + wedge[i*Mmax+j,2]*wedge[i*Mmax+j,4]/vol
            zc = zc + wedge[i*Mmax+j,3]*wedge[i*Mmax+j,4]/vol
        
        melon[i,0] = i
        melon[i,1] = xc
        melon[i,2] = yc
        melon[i,3] = zc
        melon[i,4] = vol
    return(melon)

"""
SubFunction tri1Surf and tri2Surf:
find tricentre and trinorm of surface triangles;
triangle1 = node2-1-3, clockwise
triangle2 = node2-3-4, clockwise
tricentre = [meshid, triangle1 area, triangle1 centre, triangle2 area, triangle2 centre]
trinorm = [meshid,triangle norm1, triangle norm2]
"""
  
def tri1Surf(head,i,Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh):
    """
    For the head and tail side by one surface of three nodes 
    cross = | i, j, k|
            |a1,a2,a3|
            |b1,b2,b3|
          = |a2*b3 - a3*b2|*i-|a1*b3-a3*b1|*j+|a1*b2-a2*b1|*k
    dot(|a1,a2,a3|,|b1,b2,b3|) = a1*b1+a2*b2+a3*b3
    vector01 = np.array([n1[0]-n0[0],n1[1]-n0[1],n1[2]-n0[2]])
    vector12 = np.array([n2[0]-n1[0],n2[1]-n1[1],n2[2]-n1[2]])
    dot.product = v01[0]*v12[0]+v01[1]*v12[1]+v01[2]*v12[2]
    """

    meshx = mesh[0]
    meshy = mesh[1]
    meshz = mesh[2]
    if head == 1:       
        n1 = np.array([meshx[i,1],meshy[i,1],meshz[i,1]])
        n2 = np.array([meshx[i,2],meshy[i,2],meshz[i,2]])
        n4 = np.array([meshx[i,4],meshy[i,4],meshz[i,4]])
        v21 = n1 - n2 # represents vector21
        v14 = n4 - n1
        dot = np.dot(v21,v14) #Triangle 1,2,4
        d21 = (v21[0]**2+v21[1]**2+v21[2]**2)**0.5
        d14 = (v14[0]**2+v14[1]**2+v14[2]**2)**0.5
        beta = np.arccos(dot/(d21*d14))
        area = 0.5*d21*d14*np.sin(beta)
        cross = np.cross(v21,v14)
        norm = cross/(cross[0]**2+cross[1]**2+cross[2]**2)**0.5
        centre = 1/3*(n1+n2+n4)
    else:
        n1 = np.array([meshx[i,1],meshy[i,1],meshz[i,1]])
        n2 = np.array([meshx[i,2],meshy[i,2],meshz[i,2]])
        n3 = np.array([meshx[i,3],meshy[i,3],meshz[i,3]])
        v21 = n1 - n2
        v13 = n3 - n1
        dot = np.dot(v21,v13)
        d21 = (v21[0]**2+v21[1]**2+v21[2]**2)**0.5
        d13 = (v13[0]**2+v13[1]**2+v13[2]**2)**0.5
        beta = np.arccos(dot/(d21*d13))
        area = 0.5*d21*d13*np.sin(beta)
        cross =  np.cross(v21,v13)
        norm = cross/(cross[0]**2+cross[1]**2+cross[2]**2)**0.5
        centre = 1/3*(n1+n2+n3)
    return(area,centre,norm)
    
def tri2Surf(i,Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh): 
    """
    For a general side formed by two surfaces of four nodes
    cross = | i, j, k|
            |a1,a2,a3|
            |b1,b2,b3|
          = |a2*b3 - a3*b2|*i-|a1*b3-a3*b1|*j+|a1*b2-a2*b1|*k
    dot(|a1,a2,a3|,|b1,b2,b3|) = a1*b1+a2*b2+a3*b3
    vector01 = np.array([n1[0]-n0[0],n1[1]-n0[1],n1[2]-n0[2]])
    vector12 = np.array([n2[0]-n1[0],n2[1]-n1[1],n2[2]-n1[2]])
    dot.product = v01[0]*v12[0]+v01[1]*v12[1]+v01[2]*v12[2]
    """
    meshx,meshy,meshz = [mesh[0],mesh[1],mesh[2]]
    n1 = np.array([meshx[i,1],meshy[i,1],meshz[i,1]])
    n2 = np.array([meshx[i,2],meshy[i,2],meshz[i,2]])
    n3 = np.array([meshx[i,3],meshy[i,3],meshz[i,3]])
    n4 = np.array([meshx[i,4],meshy[i,4],meshz[i,4]])
    v21 = n1 - n2 # represents vector21
    v13 = n3 - n1 
    v23 = n3 - n2
    v34 = n4 - n3
    dot1 = np.dot(v21,v13) # Triangle 2,1,3
    dot2 = np.dot(v23,v34) # Triangle 2,3,4
    d21 = (v21[0]**2+v21[1]**2+v21[2]**2)**0.5
    d13 = (v13[0]**2+v13[1]**2+v13[2]**2)**0.5
    d23 = (v23[0]**2+v23[1]**2+v23[2]**2)**0.5
    d34 = (v34[0]**2+v34[1]**2+v34[2]**2)**0.5
    beta1 = np.arccos(dot1/(d21*d13))
    beta2 = np.arccos(dot2/(d23*d34))

    area1 = 0.5*d21*d13*np.sin(beta1)
    area2 = 0.5*d23*d34*np.sin(beta2)
    cross1 = np.cross(v21,v13)
    cross2 = np.cross(v23,v34)
    norm1 = cross1/(cross1[0]**2+cross1[1]**2+cross1[2]**2)**0.5
    norm2 = cross2/(cross2[0]**2+cross2[1]**2+cross2[2]**2)**0.5
    centre1 = 1/3*(n2+n1+n3)
    centre2 = 1/3*(n2+n3+n4)
    return(area1,area2,centre1,centre2,norm1,norm2)
    
def surfNorm(Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh):
    '''
    To find the norm of each triangular facet
    If a triangle is at the head of the circle, degenerate to triangle 1,2,4
    If a triangle is at the tail of the circle, degenerate to triangle 1,2,3
    '''   
    # tricentre = np.zeros([Nmax*Mmax,9])
    nVec = np.zeros([Nmax*Mmax,7])
    for i in range (0,Nmax*Mmax):
        if i%Mmax == 0:
            head = 1
            matrix = tri1Surf(head,i,Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh)
            nVec[i,0] = i
            nVec[i,4] = matrix[2][0]
            nVec[i,5] = matrix[2][1]
            nVec[i,6] = matrix[2][2]
        
        elif (i+1)%Mmax == 0:
            head = 0
            matrix = tri1Surf(head,i,Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh)
            nVec[i,0] = i
            nVec[i,1] = matrix[2][0]
            nVec[i,2] = matrix[2][1]
            nVec[i,3] = matrix[2][2]
        
        else:    
            matrix = tri2Surf(i,Nmax,Mmax,Rnmax,R0,x0,y0,z0,mesh)
            nVec[i,0] = i
            nVec[i,1] = matrix[4][0]
            nVec[i,2] = matrix[4][1]
            nVec[i,3] = matrix[4][2]
            nVec[i,4] = matrix[5][0]
            nVec[i,5] = matrix[5][1]
            nVec[i,6] = matrix[5][2]
    return(nVec) 
    
def n3pts(pt1,pt2,pt3):  
    #anticlock wise direction
    v12 = pt2 - pt1
    v23 = pt2 - pt3
    cross = np.cross(v12,v23)
    norm = cross/(cross[0]**2+cross[1]**2+cross[2]**2)**0.5
    return (norm)
    
def intStPL(pt1,pt2,cp,vec): 
    #Intersect point between a plane and a line
    a,b,c = [pt1[0],pt1[1],pt1[2]] #point1
    d,e,f = [pt2[0],pt2[1],pt2[2]] #point2
    x1,y1,z1 = [cp[0],cp[1],cp[2]] #centre of plane
    n1,n2,n3 = [vec[0],vec[1],vec[2]] #norm of the plane
    t = ((n1*x1+n2*y1+n3*z1)-n1*a-n2*b-n3*c)/(n1*(d-a)+n2*(e-b)+n3*(f-c))
    x,y,z = [t*(d-a)+a, t*(e-b)+b, t*(f-c)+c]
    return(x,y,z)

def ptInTri(A,B,C,P):
    #Tell if a point is inside the triangle area on the same plane
    flag = 0
    v0 = C - A
    v1 = B - A
    v2 = P - A
    dot00 = np.dot(v0,v0)
    dot01 = np.dot(v0,v1)
    dot02 = np.dot(v0,v2)
    dot11 = np.dot(v1,v1)
    dot12 = np.dot(v1,v2)
    
    invDenom = 1/(dot00*dot11-dot01*dot01)
    u = (dot11*dot01 - dot01*dot12)*invDenom
    v = (dot00*dot12 - dot01*dot02)*invDenom
    if u>=0 and v>=0 and u+v<1:
       flag = 1    
    return(flag)
    
def triNorm(node1,node2,node3):
    #Looking at -z axis assume clockwise
    #Assume clockwise direction to have vertical component downwards
    n1 = node1
    n2 = node2
    n3 = node3
    v21 = n1 - n2 # represents vector21
    v13 = n3 - n1
    dot = np.dot(v13,v21) # Triangle 2,1,3
    d21 = (v21[0]**2+v21[1]**2+v21[2]**2)**0.5
    d13 = (v13[0]**2+v13[1]**2+v13[2]**2)**0.5
    beta = np.arccos(dot/(d21*d13))
    area = 0.5*d21*d13*np.sin(beta)
    cross = np.cross(v13,v21)
    norm = cross/(cross[0]**2+cross[1]**2+cross[2]**2)**0.5
    return(area,norm)
    
def pjtPtPos(Mmax,Nmax,wedge,CircleCen,melon,mesh,snvec):
    '''
    #Find the Gravity plane projected on the slip surface
    #wedge ID = i*Mmax+j
    Node1 is Gravity Plane intersection
    Node2 is wedge Gravity Centre
    cp is either node2 as shared by both 1-2-3 and 2-4-3
    vec is either triangle1-2-3 or 2-4-3
    Check the temporary projection scheme
    A Typical Plane for quick check is i = 6 and j = 17
    '''
    Gnvec = np.zeros([Nmax,3])
    GcInt = np.zeros([Nmax,3])
    ptPjt = np.zeros([Nmax*Mmax,10])
    meshx,meshy,meshz = [mesh[0],mesh[1],mesh[2]]
    # wedge matrix contains gravity centre and Area
    # if a = 1, pt4 in triangle 123, else in triangle 234
    # Project [0] - [2] Gravity + GeomCentre project on slip surface
    # Prorject[4] - [6] Normal of the array GeomCentre to GravityCentre
    for i in range(0,Nmax):
        m = i*Mmax + round(Mmax*0.3)
        nd1 = np.array([wedge[m,1],wedge[m,2],wedge[m,3]])
        m = i*Mmax + round(Mmax*0.6)
        nd2 = np.array([wedge[m,1],wedge[m,2],wedge[m,3]])
        m = i*Mmax + round(Mmax*0.9)
        nd3 = np.array([wedge[m,1],wedge[m,2],wedge[m,3]])
        '''
        Find Circular Centre Intersect with the Plane
        Melon is the Melon Centre Matrix
        vec is the triangle of wedge
        '''
        Gnvec[i,:] = triNorm(nd1,nd2,nd3)[1]
        nd4 = np.array([CircleCen[i,0],CircleCen[i,1],CircleCen[i,2]])
        nd5 = np.array([CircleCen[i+1,0],CircleCen[i+1,1],CircleCen[i+1,2]])
        cp1 = np.array([melon[i,1],melon[i,2],melon[i,3]])
        vec1 = np.array([Gnvec[i,0],Gnvec[i,1],Gnvec[i,2]])
        GcInt[i,0] = intStPL(nd4,nd5,cp1,vec1)[0]
        GcInt[i,1] = intStPL(nd4,nd5,cp1,vec1)[1]
        GcInt[i,2] = intStPL(nd4,nd5,cp1,vec1)[2]
        for j in range(0,Mmax):
            node1 = np.array([GcInt[i,0],GcInt[i,1],GcInt[i,2]])
            #node2 = np.array([melon[i,1],melon[i,2],melon[i,3]])
            k = i*Mmax+j
            node2 = np.array([wedge[k,1],wedge[k,2],wedge[k,3]])
            cp = np.array([meshx[k,2],meshy[k,2],meshz[k,2]])
       
            if k % Mmax == 0: #if the point is at the head of the circle
                #point 3 degenerate
                vec = np.array([snvec[k,4],snvec[k,5],snvec[k,6]])
                ptPjt[k,0] = intStPL(node1,node2,cp,vec)[0]
                ptPjt[k,1] = intStPL(node1,node2,cp,vec)[1]
                ptPjt[k,2] = intStPL(node1,node2,cp,vec)[2]
            elif (k+1) % Mmax == 0: #if this point is at the tail of the circle
                vec = np.array([snvec[k,1],snvec[k,2],snvec[k,3]])
                ptPjt[k,0] = intStPL(node1,node2,cp,vec)[0]
                ptPjt[k,1] = intStPL(node1,node2,cp,vec)[1]
                ptPjt[k,2] = intStPL(node1,node2,cp,vec)[2]     
            else:
                vec = np.array([snvec[k,1],snvec[k,2],snvec[k,3]])     
                ptPjt[k,0] = intStPL(node1,node2,cp,vec)[0]
                ptPjt[k,1] = intStPL(node1,node2,cp,vec)[1]
                ptPjt[k,2] = intStPL(node1,node2,cp,vec)[2]
                #if ptPjt is inside triangle 1-2-3
                pt1 = np.array([meshx[k,1],meshy[k,1],meshz[k,1]])
                pt2 = np.array([meshx[k,2],meshy[k,2],meshz[k,2]])
                pt3 = np.array([meshx[k,3],meshy[k,3],meshz[k,3]])
                pt4 = np.array([ptPjt[k,0],ptPjt[k,1],ptPjt[k,2]])
                flag = ptInTri(pt1,pt2,pt3,pt4)
                ptPjt[k,3] = flag 
             
                if flag  == 1:
                    pass
                else:
                    vec = np.array([snvec[k,4],snvec[k,5],snvec[k,6]])
                    ptPjt[k,0] = intStPL(node1,node2,cp,vec)[0] #Project intersection at slip surface
                    ptPjt[k,1] = intStPL(node1,node2,cp,vec)[1]
                    ptPjt[k,2] = intStPL(node1,node2,cp,vec)[2]
        
            dist = ((node1[0]-node2[0])**2+(node1[1]-node2[1])**2+(node1[2]-node2[2])**2)**0.5
            ptPjt[k,4] = (node2[0]-node1[0])/dist #vector direction from slip centre to slip surface
            ptPjt[k,5] = (node2[1]-node1[1])/dist
            ptPjt[k,6] = (node2[2]-node1[2])/dist     
        
            if node2[2] > node1[2]: # if the gravity wedge centre is higher (z2>z1) than melon geometry centre
                ptPjt[k,7] =  -ptPjt[k,4] #vector direction from higher position to lower position
                ptPjt[k,8] =  -ptPjt[k,5]
                ptPjt[k,9] =  -ptPjt[k,6]
            else:
                ptPjt[k,7] =  ptPjt[k,4]
                ptPjt[k,8] =  ptPjt[k,5]
                ptPjt[k,9] =  ptPjt[k,6]
    return(ptPjt)  
    
def pjtPlPos(Mmax,Nmax,melon,ptPjt,mesh):
#determine gravity plane perpendicular to slip surface
    k = 0
    m = 0
    Ghpl = np.zeros([Nmax*Mmax,10])
    #GVplane = np.zeros([Nmax*Mmax,10])
    meshx,meshy,meshz = [mesh[0],mesh[1],mesh[2]]
#Store nodal cooordinates at Ghpl matrix - gravity plane perpendicular to gravity stress
    for i in range(0,Nmax*Mmax):# i is the No. of the Tetrahedron
        k = k + 1
        if k > Mmax:
            k = 0
            m = m + 1
    
        Ghpl[i,0] = i
        if i%Mmax == 0 or (i+1)%Mmax == 0:
            pass
        else:
            node1 = np.array([melon[m,1],melon[m,2],melon[m,3]])
            node2 = np.array([ptPjt[i-1,0],ptPjt[i-1,1],ptPjt[i-1,2]])
            node3 = np.array([ptPjt[i,0],ptPjt[i,1],ptPjt[i,2]])
            node4 = np.array([ptPjt[i+1,0],ptPjt[i+1,1],ptPjt[i+1,2]])
            norm1 = triNorm(node1,node2,node3)[1]
            norm2 = triNorm(node1,node3,node4)[1]
            node5 = np.array([meshx[i,1],meshy[i,1],meshz[i,1]])
            node6 = np.array([meshx[i,3],meshy[i,3],meshz[i,3]])
            node7 = np.array([meshx[i,2],meshy[i,2],meshz[i,2]])
            node8 = np.array([meshx[i,4],meshy[i,4],meshz[i,4]])
        
        
            Ghpl[i,1] = intStPL(node5,node6,node3,norm1)[0]
            Ghpl[i,2] = intStPL(node5,node6,node3,norm1)[1]
            Ghpl[i,3] = intStPL(node5,node6,node3,norm1)[2]
        
            Ghpl[i,4] = intStPL(node7,node8,node3,norm2)[0]
            Ghpl[i,5] = intStPL(node7,node8,node3,norm2)[1]
            Ghpl[i,6] = intStPL(node7,node8,node3,norm2)[2]
        
            dist = ((Ghpl[i,1]-Ghpl[i,4])**2 + (Ghpl[i,2]-Ghpl[i,5])**2 + (Ghpl[i,3]-Ghpl[i,6])**2)**0.5
        
            Ghpl[i,7] = (Ghpl[i,4] - Ghpl[i,1])/dist
            Ghpl[i,8] = (Ghpl[i,5] - Ghpl[i,2])/dist
            Ghpl[i,9] = (Ghpl[i,6] - Ghpl[i,3])/dist
    return(Ghpl)     