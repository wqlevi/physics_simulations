#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 18:08:14 2021

@author: qiwang
"""
import matplotlib.pyplot as plt

# B1 visualization
t = np.linspace(0,output['b1r'].shape[0],output['b1r'].shape[0])
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.plot(t,output['b1r'])
ax2.plot(t,output['b1i'],color = 'red')
ax1.legend('Real')
ax2.legend('Imaginary')

# gradient visualization
fig,ax = plt.subplots(g.shape[1],1)
fig.subplots_adjust(hspace = 1, wspace=.001)
ax = ax.ravel()
for i in range(g.shape[1]):
    # ax[i] = plt.subplot(g.shape[1],1,i+1)
    ax[i].plot(g[:,i])
    ax[i].set_title('gradient %s'%(i+1))
    
# S visualization
sr_new = output['sr'].reshape(npos,nTx)
si_new = output['si'].reshape(npos,nTx)
t = np.linspace(0,sr_new.shape[0],sr_new.shape[0])
ax1 = plt.subplot(211)
ax2 = plt.subplot(212)
ax1.plot(t,sr_new)
ax2.plot(t,si_new)
ax1.set_title('Real')
ax2.set_title('Imaginary')

# NEW mx,my,mz visualization
fig,ax = plt.subplots(g_new.shape[1],1)
fig.subplots_adjust(hspace = 1, wspace=.001)
ax = ax.ravel()
for i in range(g_new.shape[1]):
    # ax[i] = plt.subplot(g.shape[1],1,i+1)
    ax[i].plot(g_new[:,i,:])
    ax[i].set_title('gradient %s'%(i+1))
# OLD mx,my,mz
fig,ax = plt.subplots(g_new.shape[1],1)
fig.subplots_adjust(hspace = 1, wspace=.001)
ax = ax.ravel()
for i in range(g_new.shape[1]):
    # ax[i] = plt.subplot(g.shape[1],1,i+1)
    ax[i].plot(g_old[i,:,:])
    ax[i].set_title('gradient %s'%(i+1))