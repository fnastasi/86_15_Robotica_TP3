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

"""
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
"""
# In[]

"""
Función para calcular el jacobiano (solo velocidad de traslación) en la terna 0
"""

def calc_Jac_tern0(t1,t2,t3):
    t1 = t1*pi/180
    t2 = t2*pi/180
    t3 = t3*pi/180
    
    J = np.array([[-sin(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), cos(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*cos(t1)*cos(t2+t3)],
                  [cos(t1)* (d4*sin(t2+t3) + a2*cos(t2) +a1), sin(t1)*(d4*cos(t2+t3)-a2*sin(t2)), d4*sin(t1)*cos(t2+t3)],
                  [0, -(d4*sin(t2+t3) +a2*cos(t2)), -d4*sin(t2+t3)]])
    return J
# In[]




def graf_vel(xi, yi,zi,xf,yf,zf,dt,V_o_4):
    """
    xi,yi,zi -> posición inicial
    xf,yf,zf -> posición final
    dt -> paso
    V_0_4 -> Vector velocidad expresado en la terna 0 
    
    xi = 10; yi = -100; zi =500
    
    """
    
    # Inicialización de variales
    y = yi
    
    time = (yf-yi)/(V_o_4[1])
    cant_iter = int(time/dt)
    
    t1_vec = zeros(cant_iter)
    t2_vec = zeros(cant_iter)
    t3_vec = zeros(cant_iter)
    
      
    
    for i in range(cant_iter):
        # Calculo de matriz homogénea A en la posición del robot
        A = np.array([[1, 0, 0, xi], [0,1,0,y] ,[0,0,1,zi],[0,0,0,1] ])
        
        # Calculo los valores de tita en la posición del robot con una configuración y los guardo en los vectores correspondientes
        t1,t2,t3, t4,t5,t6 = pos_prob_inv(A,1,1,1)
        
        t1_vec[i] = t1
        t2_vec[i] = t2
        t3_vec[i] = t3
        
        # Calculo el Jacobiano 
        J =calc_Jac_tern0(t1,t2,t3)
        
        # Calculo los valores de tita punto
        t_vec_punto = linalg.solve(J,V_o_4)
        
        # Actualizó el valores de  y
        y = y + dt*V_o_4[1]
        
        
    
    time_vec = linspace(0,time,cant_iter)
    
    fig1,ax1 = plt.subplots()     
    #ax[0].plot(time_vec*1e3,t1_vec)
    #ax[0].set_xlabel('Tiempo')
    #ax[0].set_ylabel(r'$\theta_1 (°)$')
    #ax[0].grid(True)
    ax1.plot(time_vec*1e3,t1_vec)
    ax1.set_xlabel('Tiempo (ms)')
    ax1.set_ylabel(r'$\theta_1 (°)$')
    ax1.grid(True)
    
    fig2,ax2 = plt.subplots() 
    ax2.plot(time_vec*1e3,t2_vec)
    ax2.set_xlabel('Tiempo (ms)')
    ax2.set_ylabel(r'$\theta_2 (°)$')
    ax2.grid(True)

    fig3,ax3 = plt.subplots() 
    ax3.plot(time_vec*1e3,t3_vec)
    ax3.set_xlabel('Tiempo (ms)')
    ax3.set_ylabel(r'$\theta_3 (°)$')
    ax3.grid(True)
    
    
    fig1.savefig("Tp3_tita_1_xi_"+str(xi)+".png",dpi = 300)
    fig2.savefig("Tp3_tita_2_xi_"+str(xi)+".png",dpi = 300)
    fig3.savefig("Tp3_tita_3_xi_"+str(xi)+".png",dpi = 300)
    
# In[]
    
# Creación de gráficos

xi = 0.1; yi = -100; zi =500
xf = 0.1; yf = 100; zf =500   
dt = 1e-4    
V_o_4 = array([0,200,0])

graf_vel(xi, yi,zi,xf,yf,zf,dt,V_o_4)
