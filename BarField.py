#!/usr/bin/env python2
# -*- coding: utf-8 -*-
% matplotlib qt
"""
Created on Mon Jan 25 11:46:24 2021

@author: qiwang
"""
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
 
x = np.linspace(-4,4,10)
y = np.linspace(-4,4,10)
z = np.linspace(-4,4,10)
 
x,y,z = np.meshgrid(x,y,z)
 
fig = plt.figure()
ax = fig.gca(projection = '3d')
 
def B(x,y):
    i = 1 # ampere
    mu = 1.26*10**(-6)
    mag = (mu*i)/(2*np.pi*np.sqrt(x**2+y**2))
    #mag = (mu/(2*np.pi))*(i/np.sqrt((x)**2+(y)**2))
    by = mag*(np.cos(np.arctan2(y,x)))
    bx = mag*(-np.sin(np.arctan2(y,x)))
    bz = z*0
    return bx,by,bz
 
def cylinder(r):
    phi = np.linspace(-2*np.pi,2*np.pi,100)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return x,y

bx,by,bz = B(x,y)
cx,cy = cylinder(0.2)

ax.quiver(x,y,z,bx,by,bz,color= 'b',length = 1,normalize = True)

for i in np.linspace(-4,4,800):
    ax.plot(cx,cy,i,label = 'Cylinder',color = 'r')
    
plt.title('Magnetic field simulation')
plt.xlabel('x')
plt.ylabel('y')
plt.show()