#!/usr/bin/env python2
# -*- coding: utf-32 -*-
"""
Created on Sat Apr  3 02:00:05 2021

This demo opens infile.dat and write contents:
    written in class, and methods


@author: qiwang
"""
import numpy as np
# Read file with ordered type encoding
infile = 'infile.dat'

class data_IO:
    '''
    I/O of binary file, and constant encoding format;
    Input: binary file to numpy array
    Output: numpy array to binary file
    '''
    def __init__(self):
        # customized encoding format w.r.t examplary data
        self.dt = np.dtype([('mode','u4'),('ntime','u4'),('nTx','u4'),('npos','u4'),
       ('nfreq','u4'),('properties','u4'),('none','f8'),('t1','f8'),('t2','f8'),
       ('b1r','(6544,1)f8'),('b1i','(6544,1)f8'),('g','(818,3)f8'),
       ('b0','(57793,1)f8'),('tsteps','(818,1)f8'),('dfreq','f8'),
       ('dpos','(57793,3)f8'),('mx','(57793,1)f8'),('my','(57793,1)f8'),
       ('mz','(57793,1)f8'),('sr','(462344,1)f8'),('si','(462344,1)f8')])
        
    def read_file(self,infile): 
        # [allocation]raise exception for null opening
        try:
            f = open(infile,'rb+')
        except IOError:
            print("IO Error: File doesn't exist.")

        np_data_in = np.fromfile(f,dtype = self.dt)
        return np_data_in
    
    def write_file(self,outfile,data):
        '''
        plug in data in ndarray type
        '''        
        # [aborted]try to use f.write("%d,%f"%(var1,var2)) to continuously write vars
        # attempt of data.tofile(f,format = '%s')
        data.tofile(f,format = '%s')
        f.close()
        
        

