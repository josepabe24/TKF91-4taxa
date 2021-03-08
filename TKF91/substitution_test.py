# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 09:50:04 2021

@author: josh
"""

from pyvolve import model
import numpy as np

model_typ="nucleotide"
param={'state_freqs': np.array([0.25, 0.25, 0.25, 0.25]), 
       'mu': {'AC': 1.0,'CA': 1.0,'AG': 1.0,'GA': 1.0,'AT': 1.0,'TA': 1.0,
              'CG': 1.0,'GC': 1.0,'CT': 1.0,'TC': 1.0,'GT': 1.0,'TG': 1.0}}

m=model.Model(model_typ,param)

m.extract_rate_matrix()
m.extract_parameters()

np.exp(m.extract_rate_matrix()*2)
