#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 22:57:21 2021

@author: qiwang
"""

# skipping all utilities functions, here're two main blochs:
import os
import sys
import numpy as np 

def rotation_matrix(nx,ny,nz):
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

def precess(bx, by, bz):
    b = sqrt(bx**2 + by**2 + bz**2)
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

def relax(e1,e2):
    mag[0] *= e2
    mag[1] *= e2
    mag[2] = 1 + e1 * (mag[2] - 1)
    return mag

def blochsim(b1_complex, xgrad, ygrad, zgrad, tsteps, ntime, e1, e2,
              df, dx, dy, dz, mx,
              my, mz, mode, nTx, s_complex):
    GAMMA = 26753.0
    gamma_xyz = [dx,dy,dz] * GAMMA
    mcurr0[0] = mx
    mcurr0[1] = my
    mcurr0[2] = mz
    amat,imat = np.eye(3),np.eye(3)
    
    for tcount in range(ntime):
        eff_b1r, eff_b1i = 0,0
        curt = tcount * nTx
        for i in range(nTx):
            eff_b1r += (b1_complex[curt+i]*s_complex[i]).real
            eff_b1i += (b1_complex[curt+i]*s_complex[i]).imag
        if FASTER_BLOCH:
            rotz = -(xgrad++ * gammadx + ygrad++ * gammady + zgrad++ * gammadz +
                     df * TWOPI) * tsteps # pointers to array here
            roty = (-eff_b1r * GAMMA * tsteps)
            rotz = (abs(eff_b1i) * GAMMA * tsteps++)
            rotmat = rotation_matrix(rotx, roty, rotz)
            if mode == 1:
                arot = np.matmul(amat, rotmat)
                brot = np.matmul(bvec, rotmat)
            else:
                mcurr1 = np.matmul(mcurr0, rotmat)
                
            if e1 != None and e2 != None:
                decvec[2] = 1 - e1
                decmat[0] = e2
                decmat[4] = e2++
                decmat[8] = e1++
                if mode == 1:
                    amat = np.matmul(arot, decmat)
                    bvec = np.matmul(brot, decmat)
                    bvec = bvec + decvec
                else:
                    mcurr0 = np.matmul(mcurr1, decmat)
                    mcurr0 = mcurr0 + decvec
            else:
                if mode == 1:
                    amat = np.copy(arot)
                    bvec = np.copy(brot)
                else:
                    mcurr0 = np.copy(mcurr1)
        else:
            rotz = -(xgrad++ * gammadx + ygrad++ * gammady + zgrad++ * gammadz +
                     df * TWOPI) * tsteps;
            rotx = (-eff_b1r * GAMMA * tsteps)
            roty = (abs(eff_b1i * GAMMA * tsteps++))
            mcurr0 = precess(rotx, -roty, rotz)
            if e1 != None and e2 != None : 
                mcurr0 = relax(e1++, e2++)

    if mode == 2:
        mx, my, mz = mcurr0[0],mcurr0[1],mcurr0[2]
        mx++,my++,mz++
    if mode == 0:
        mx, my, mz = mcurr0[0],mcurr0[1],mcurr0[2]
    elif mode == 1:
        amat *= -1
        amat += imat
        imat = np.linalg.inv(amat)
        mvec = np.matmul(imat, bvec)
        mx, my, mz = mvec[0],mvec[1],mvec[2]
        
def blochsimfz(b1_complex, xgrad, ygrad, zgrad, tsteps, ntime, t1, t2, dfreq, nfreq, dxpos, dypos,
               dzpos, npos, mx, my, mz, mode, b0, nTx, s_complex):
        if (mode != 1 and FAST_BLOCH):
            sys.exit("Steady state calculation not supported in BLOCH_FAST mode! Exiting\n")
        if (t1>0 and t2>0):
            e1 = malloc()
            e2 = malloc()
            
            e1ptr = e1
            e2ptr = e2
            tstepsptr = tsteps
            tstepptr2  = tstepsptr
            for count in range(ntime):
                e1ptr = exp(-tstepptr / t1)
                e2ptr = exp(-tstepptr2/ t2)
                e1ptr += 1
                e2ptr += 1
                tstepptr2 += 1
        else:
            e1, e2 = None,None
        for fcount in range(nfreq):
            # allocated for pragma in C
            for poscount in range(npos):
                if mode == 2:
                    idx = poscount * ntime + npos * fcount
                else:
                    idx = poscount + npos * fcount
                idxS = poscount * nTx
                if mode == 3:
                    blochsim(b1_complex, xgrad, ygrad, zgrad, tsteps, ntime, e1, e2,
                         dfreq[fcount] + b0[poscount], dxpos[poscount], dxpos[poscount],
                         dzpos[poscount], mx[idx], my[idx], mz[idx], 1, nTx,
                         s_complex[idxS])
                    blochsim(b1_complex, xgrad, ygrad, zgrad, tsteps, ntime, e1, e2,
                         dfreq[fcount] + b0[poscount], dxpos[poscount], dxpos[poscount],
                         dzpos[poscount], mx[idx], my[idx], mz[idx], 2, nTx,
                         s_complex[idxS])
                else:
                    blochsim(b1_complex, xgrad, ygrad, zgrad, tsteps, ntime, e1, e2,
                         dfreq[fcount] + b0[poscount], dxpos[poscount], dxpos[poscount],
                         dzpos[poscount], mx[idx], my[idx], mz[idx], mode, nTx,
                         s_complex[idxS])
                    
        del e1
        del e2
        

