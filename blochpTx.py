#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 19:04:31 2021

@author: qiwang
"""

import numpy as np
import sys
GAMMA = 26753

def calcrotmat(nx,ny,nz):
    '''
    calculate rotating matrix
    '''
    phi = np.sqrt(nx**2 + ny**2 + nz**2)
    rmat = np.zeros((3,3))
    if phi == 0:
        rmat = np.eye(3)
    else:
        ar = np.cos(phi/2)
        ai = -nz * np.sin(phi/2) / phi
        br = ny * np.sin(phi/2) / phi
        bi = -nx * np.sin(phi/2) / phi
        # rotation matrix
        rmat[0] = [(ar**2-ai**2-br**2+bi**2),
                   (-ar*ai*2-br*bi*2),
                   (-ar*br*2 + ai*bi*2)]
        rmat[1] = [(ar*ai*2 - br*bi*2),
                   (ar**2 - ai**2 + br**2 -bi**2),
                   (-ai*br*2 - ar*bi*2)]
        rmat[2] = [(ar*br*2 + ai*bi*2),
                   (ar*bi*2 - ai*br*2),
                   (ar**2 + ai**2 - br**2 - bi**2)]
    return rmat 

def precess(bx, by, bz,mag):

    '''
    where bx,by,bz are 1D arrays, mag be NUM*3 matrix
    '''
    # [allocation]dimensional checking expected
    b = np.sqrt(bx**2 + by**2 + bz**2)    
    if (b != 0):
        b_array = np.array([bx,by,bz]) / b
        
        
        c = np.sin(b/2)
        c = 2*c**2
        s = np.sin(b)
        k = mag[0] * b_array[0] + mag[1] * b_array[1]+ mag[2] * b_array[2]
        
        mag[0] += (bx * k)
        for i in range(3):
            mag[i] += (b_array[i] * k - mag[i]) * c + (mag[1-i] * b_array[i-1] - mag[i-1] * b_array[1-i]) * s
    return mag 
def relax(e1,e2,mag):
    mag[0] *= e2
    mag[1] *= e2   
    mag[2] = 1 + e1 * (mag[2] - 1)
    return mag
    
def blochsim(b1_complex, g, tsteps, ntime, e1, e2,
              df, dx, dy, dz, mx,
              my, mz, mode, nTx, s_complex):
    '''
    Parameters
    ----------
    b1_complex : complex 6544*1 array
        DESCRIPTION.
    g : 818*3 array
        818 time points.
    tsteps : 818*1 array
        818 time points.
    ntime : int
        818 time points.
    e1 : TYPE
        one element of the array with same dim as timepoints 818.
    e2 : TYPE
        DESCRIPTION.
    df : int
        dfreq[j]+b0[i], i in position loop, j in frequency loop.
    dx : int
        dfreq[i]+b0[i], i in position point loop.
    dy : TYPE
        dfreq[i]+b0[i], i in position point loop.
    dz : TYPE
        dfreq[i]+b0[i], i in position point loop.
    mx : TYPE
        DESCRIPTION.
    my : TYPE
        DESCRIPTION.
    mz : TYPE
        DESCRIPTION.
    mode : int
        0|1|2|3.
    nTx : int
        num of gradients = 8.
    s_complex : complex 462344*1 array
        all positional complex data from all 8 gradients.

    Returns
    -------
    None.

    '''
    FASTER_BLOCH = bool(0)
    # what're e1,e2?
    gamma_xyz = np.array([dx,dy,dz]) * GAMMA # 3*1 array, but dx,dy,dz are double at frequency point 
    mcurr0 = np.array([mx,my,mz]).transpose()           # 3*1 array
    amat = imat = np.eye(3)                 # 3*3
    decvec=bvec = np.zeros((1.3))           # 1*3
    decmat = np.zeros((3,3))                # 3*3
    
    for tcount in range(ntime):
        eff_b1r= eff_b1i= 0                 # int
        curt = tcount * nTx     # [0,818] * 8
        for i in range(nTx):    # 8
            eff_b1r += (b1_complex[curt+i]*s_complex[i]).real   # s corresponds every 8 b1_complex
            eff_b1i += (b1_complex[curt+i]*s_complex[i]).imag
            print(i)
        if not FASTER_BLOCH:# switch on FASTER_BLOCH
            rotz = (np.sum(-g[tcount,:]*gamma_xyz) - df*2*np.pi)*tsteps[tcount] # int
            rotx = (-eff_b1r*GAMMA*tsteps[tcount])
            roty = (eff_b1i*GAMMA*tsteps[tcount])   
            rotmat = calcrotmat(rotx,roty,rotz)
            
            if mode==1:
                arot = np.matmul(rotmat,amat )
                brot = np.matmul(bvec, rotmat)
            else:
                mcurr1 = np.matmul(mcurr0,rotmat)
                
            if e1.any() and e2.any():
                decvec[:,2] = 1 - e1[tcount]
                decmat[0,0] = e2[tcount]
                decmat[1,1] = e2[tcount]
                decmat[2,2] = e1[tcount]
                if mode == 1:
                    amat = np.matmul(decmat,arot ) # mat *mat
                    bvec = np.matmul(brot,decmat) # vec*mat
                    bvec = bvec + decvec
                else:
                    mcurr0 = np.matmul(mcurr1,decmat)
                    mcurr0 = mcurr0 + decvec
            else:
                if mode == 1:
                    amat = np.copy(arot)
                    bvec = np.copy(brot)
                else:
                    mcurr0 = np.copy(mcurr1)
        else:# FASTER_BLOCH
            rotz = (np.sum(-g[tcount,:]*gamma_xyz) - df*2*np.pi,)*tsteps[tcount]
            rotx = (-eff_b1r * GAMMA * tsteps[tcount])
            roty = (abs(eff_b1i * GAMMA * tsteps[tcount]))
            mcurr0 = precess(rotx, -roty, rotz,mcurr0)
            if e1 != None and e2 != None : 
                mcurr0 = relax(e1[tcount], e2[tcount])
                
        if mode == 2:
            mx,my,mz = mcurr0[:,0],mcurr0[:,1],mcurr0[:,2]
    
    if mode == 0:
        mx,my,mz = mcurr0[:,0],mcurr0[:,1],mcurr0[:,2]
    elif mode == 1:
        amat *= -1
        amat += imat
        imat = np.linalg.inv(amat)
        mvec = np.matmul(imat, bvec)
        mx, my, mz = mvec[0],mvec[1],mvec[2]   
    return mx,my,mz           
        
def blochsimfz(b1_complex,g,tsteps,ntime,t1,t2,
               dfreq,nfreq,dx,dy,dz,npos,mx,my,mz,mode,b0,nTx,s_complex):
    FASTER_BLOCH = bool(0)
    if (mode != 1 and FASTER_BLOCH):
        sys.exit("Steady state calculation not supported in BLOCH_FAST mode! Exiting\n")
    if (t1>0 and t2 >0):
        e1 = np.exp(-tsteps/t1) # array of size as tsteps
        e2 = np.exp(-tsteps/t2)
    else:
        e1 = e2 = None
    
    for fcount in range(nfreq):
        for pscount in range(npos):
            if mode == 2:
                idx = pscount * ntime + npos * fcount
            else:
                idx = pscount + npos*fcount 
                
            idxS = pscount * nTx #should idxS end at same index as s_complexn, now 8 index earlier?
            if mode == 3:
                blochsim(b1_complex,g,tsteps,ntime,e1,e2,
                         dfreq[fcount]+b0[pscount],
                         dx[pscount],dy[pscount],dz[pscount],mx[idx],my[idx],mz[idx], 1,
                         nTx,s_complex[idxS])
                blochsim(b1_complex,g,tsteps,ntime,e1,e2,
                         dfreq[fcount]+b0[pscount],
                         dx[pscount],dy[pscount],dz[pscount],mx[idx],my[idx],mz[idx],2,
                         nTx,s_complex[idxS])
            else:
                blochsim(b1_complex,g,tsteps,ntime,e1,e2,
                         dfreq[fcount]+b0[pscount],
                         dx[pscount],dy[pscount],dz[pscount],mx[idx],my[idx],mz[idx],mode,
                         nTx,s_complex[idxS])
                
            
        
        
                
                
                        