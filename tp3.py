#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 18:01:49 2019

@author: franco
"""

# In[]

# Imports 

import numpy as np
import matplotlib.pyplot as plt
from tp2 import  *  

# In[]

"""
Función para calcular el jacobiano (solo velocidad de traslación) en la terna 0
"""

def calc_Jac_tern0(t1,t2,t3):
    
    J = np.array([[-sin(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), cos(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*cos(t1)*cos(t2+t3)],
                  [cos(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), sin(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*sin(t1)*cos(t2+t3)],
                  [0, -(d4*sin(t2+t3) +a2*cos(t2)), -d4*sin(t2+t3)]])
    return J
# In[]

def graf_vel(xi, yi,zi,xf,yf,zf,dt,V_0_4):
    """
    xi,yi,zi -> posición inicial
    xf,yf,zf -> posición final
    dt -> paso
    V_0_4 -> Vector velocidad expresado en la terna 0 
    
    """
    
    # Variable homogénea inicial
    Ai = np.array([[1, 0, 0, xi], [0,1,0,yi] ,[0,0,1,zi],[0,0,0,1] ])
    t_vec = pos_prob_inv(Ai)
    
    