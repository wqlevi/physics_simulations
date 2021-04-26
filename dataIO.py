#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 18:44:59 2021

This demo opens infile.dat and write contents:
    written in class, and methods
    
@author: qiwang
"""
import numpy as np
import numpy.lib.recfunctions as rfn
# Read file with ordered type encoding
#infile = 'infile.dat'

class data_IO:
    '''
    I/O of binary file, and constant encoding format;
    Input: binary file to numpy array
    Output: numpy array to binary file
    '''
    
    

    def __init__(self):
        # customized encoding format w.r.t examplary data
        self.name = [['mode','ntime','nTx','npos','nfreq','properties','none','t1','t2'],
                    ['b1r','b1i','g','b0','tsteps','dfreq','dpos','mx','my','mz','sr','si']]
        
        self.dtype = [6*['u4']+3*['f8'],
                    12*['f8']]    

    def read_file(self,infile): 
        '''
        WARNING: THIS FUNCTION BRING 2 OUTPUTS
        input: input filename
        output : ndarray in structures
        '''  
        # [allocation]raise exception for null opening
        dt_int = np.dtype([('mode','u4'),('ntime','u4'),('nTx','u4'),('npos','u4'),
       ('nfreq','u4'),('properties','u4'),('none','f8'),('t1','f8'),('t2','f8')])
        with open (infile,'rb+') as f:
            para_basic = np.fromfile(f,dtype=dt_int,count = 1)
            f.close()
            
        dt_mat = np.dtype([('b1r','f8',(para_basic['ntime'][0]*para_basic['nTx'][0])),('b1i','f8',(para_basic['ntime'][0]*para_basic['nTx'][0])),('g','f8',(para_basic['ntime'][0],3)),
       ('b0','f8',para_basic['npos'][0]),('tsteps','f8',para_basic['ntime'][0]),('dfreq','f8'),
       ('dpos','f8',(para_basic['npos'][0],3)),('mx','f8',para_basic['npos'][0]),('my','f8',para_basic['npos'][0]),
       ('mz','f8',para_basic['npos'][0]),('sr','f8',para_basic['npos'][0]*para_basic['nTx'][0]),('si','f8',para_basic['npos'][0]*para_basic['nTx'][0])])

        with open (infile,'rb+') as f:
            para_matrix = np.fromfile(f,dtype = dt_mat,offset = para_basic.itemsize)
            f.close()
        return para_basic,para_matrix
    
    def read_file2(self,infile):
        '''
        UPDATES: A MORE ORGANIZED FORMAT OF DATA
        '''        
        
        dt_int = np.dtype([(i,j) for i,j in zip(self.name[0],self.dtype[0])]) 
        with open(infile,'rb+') as f:
            para_bsc = np.fromfile(f,dtype = dt_int, count = 1)
            f.close()
              
        shape_mat = [
        para_bsc['ntime'][0]*para_bsc['nTx'][0],
        para_bsc['ntime'][0]*para_bsc['nTx'][0],
        (para_bsc['ntime'][0],3),
        para_bsc['npos'][0],
        para_bsc['ntime'][0],
        1,
        (para_bsc['npos'][0],3),
        para_bsc['npos'][0],
        para_bsc['npos'][0],
        para_bsc['npos'][0],
        para_bsc['npos'][0]*para_bsc['nTx'][0],
        para_bsc['npos'][0]*para_bsc['nTx'][0]
        ]
                
        dt_mat = np.dtype([(i,j,k) for i,j,k in zip(self.name[1],self.dtype[1],shape_mat)])
        with open (infile,'rb+') as f:
            para_mat = np.fromfile(f,dtype = dt_mat,offset = para_bsc.itemsize)
            f.close()
        output = rfn.merge_arrays([para_bsc,para_mat],flatten = True, usemask=False)
        return output
            
    
    def write_file(self,outfile,data):
        '''
        input: output filename, data(ndarray type)
        '''        
        f = open(outfile,'w')
        # [aborted]try to use f.write("%d,%f"%(var1,var2)) to continuously write vars
        # attempt of data.tofile(f,format = '%s')
        data.tofile(f,format = '%s')
        f.close()