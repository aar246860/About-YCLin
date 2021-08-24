#!/usr/bin/env python
# coding: utf-8

# In[106]:


from anaflow import theis, ext_theis_2d
import scipy.special as sc
import numpy as np
from scipy.special import expi
from mpmath import *
mp.dps = 10; mp.prec = 60; mp.pretty = True
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings


# In[128]:


def app_ext_theis_2d(time, rad, S, TG, var, Len_scale, rate):
    xi = 1.6
    rw = 0.05
    R = rw + 1.49861*np.sqrt(TG*time/S) 
    lamda1 = -var*(rad*xi)**2/(2*(Len_scale**2+(rad*xi)**2))
    lamda2 = var*Len_scale**2/(2*(Len_scale**2+(rad*xi)**2))
    lamda1R = -var*(R*xi)**2/(2*(Len_scale**2+(R*xi)**2))
    lamda2R = var*Len_scale**2/(2*(Len_scale**2+(R*xi)**2))
    A1 = -(rate/(4*np.pi*TG))*(expi(lamda2R)-np.exp(var/2)*expi(lamda1R))
    A2 = -rate/(2*np.pi*TG)
    return A1 + 0.5*A2*(np.exp(var/2)*expi(lamda1)-expi(lamda2))


# In[134]:


#Example1: Varying Len_scale from 1 to 5#
#r value
t1=np.logspace(0,4,num=100,base=10)
#T value
TT1=[];
TT2=[];
TT3=[];
TT4=[];
TT5=[];
TTz1=[];
TTz2=[];
TTz3=[];
TTz4=[];
TTz5=[];
for n in range(0, len(t1)):
    TTz1.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TTz2.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 2, 0.0001))
    TTz3.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 3, 0.0001))
    TTz4.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 4, 0.0001))
    TTz5.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 5, 0.0001))
    TT1.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TT2.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 2, 0.0001))
    TT3.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 3, 0.0001))
    TT4.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 4, 0.0001))
    TT5.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 5, 0.0001))
#in semi-log scale
plt.plot(t1, TTz1, label=label_TG, color="C"+str(1), marker=".",markevery=3)
plt.plot(t1, TTz2, label=label_TG, color="C"+str(2), marker=".",markevery=3)
plt.plot(t1, TTz3, label=label_TG, color="C"+str(3), marker=".",markevery=3)
plt.plot(t1, TTz4, label=label_TG, color="C"+str(4), marker=".",markevery=3)
plt.plot(t1, TTz5, label=label_TG, color="C"+str(5), marker=".",markevery=3)
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
plt.plot(t1, TT4, label=label_ef, color="C"+str(4))
plt.plot(t1, TT5, label=label_ef, color="C"+str(5))
plt.xlabel("t (sec)")
plt.ylabel("s (m)")
plt.show()


# In[133]:


#Example2: Varying var from 1 to 5#
#r value
t1=np.logspace(0,4,num=100,base=10)
#T value
TT1=[];
TT2=[];
TT3=[];
TT4=[];
TT5=[];
TTz1=[];
TTz2=[];
TTz3=[];
TTz4=[];
TTz5=[];
for n in range(0, len(t1)):
    TTz1.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TTz2.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 2, 1, 0.0001))
    TTz3.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 3, 1, 0.0001))
    TTz4.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 4, 1, 0.0001))
    TTz5.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 5, 1, 0.0001))
    TT1.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TT2.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 2, 1, 0.0001))
    TT3.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 3, 1, 0.0001))
    TT4.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 4, 1, 0.0001))
    TT5.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 5, 1, 0.0001))
#in semi-log scale
plt.plot(t1, TTz1, label=label_TG, color="C"+str(1), marker=".",markevery=2)
plt.plot(t1, TTz2, label=label_TG, color="C"+str(2), marker=".",markevery=2)
plt.plot(t1, TTz3, label=label_TG, color="C"+str(3), marker=".",markevery=2)
plt.plot(t1, TTz4, label=label_TG, color="C"+str(4), marker=".",markevery=2)
plt.plot(t1, TTz5, label=label_TG, color="C"+str(5), marker=".",markevery=2)
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
plt.plot(t1, TT4, label=label_ef, color="C"+str(4))
plt.plot(t1, TT5, label=label_ef, color="C"+str(5))
plt.xlabel("t (sec)")
plt.ylabel("s (m)")
plt.show()


# In[ ]:




