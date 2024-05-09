# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:26:26 2024

@author: mb1536
"""

import pickle as pkl

with open('ADV&strain.pkl', 'rb') as f:
    [vector_date, current_mag, v_time, vlink] = pkl.load(f)
    
v_time = v_time[1:327862]
vlink = vlink[1:327862]

t = v_time[-1]
p = vector_date[-1]