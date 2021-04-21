#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 22:13:55 2021

This demo opens ini file

@author: qiwang
"""
# .ini read
import os
import configparser as cp

os.listdir(os.getcwd())

cp = cp.ConfigParser()
cp.read('UP_2D_FA7.ini')
sections = cp.sections()
# laod database to dict
parser_dict = {s:dict(cp.items(s)) for s in cp.sections()}

