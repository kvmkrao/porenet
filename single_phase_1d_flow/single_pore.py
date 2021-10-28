
import sys
import itertools
import pandas as pd
import tensorflow as tf
import numpy as np
import math
from sklearn import datasets
from sklearn import metrics
from sklearn import model_selection
from sklearn import preprocessing
from matplotlib import pyplot
from sklearn.metrics import mean_squared_error, r2_score

mu1=3.3e-3
mu2=0
L=1e-6
sig=1.6
small=1e-2 
r=5.5e-6;
x1=small 
t1=small 
t2=10 
p1=0 
p2=0; 
A=math.pi*r**2; 
pc=2*math.pi*r*(sig/A);
dp=p1-p2+pc;
dt=0.05; 
n=1000; 
x=np.zeros(n+1); 
xe=np.zeros(n+1); 
t=np.zeros(n+1);
x[1]=x1;
t[1]=t1;
xe[1]=x1;
q=0;
v=0;
for i in range(1,n):
    a=8*(mu1*x[i] + mu2*(L-x[i])) ;  
    b=math.pi*r**4;
    q=dp*b/a; v=q/A;
    x[i+1]=x[i]+v*dt;
    t[i+1]=t[i]+dt;
    xe[i+1]=math.sqrt((r*sig*t[i+1])/(2*mu1));

pyplot.figure(figsize=[10,6])
pyplot.plot(t,x, marker='o', linestyle='--',label='Approx',color="blue",linewidth=2.0)
pyplot.title('Infiltration')
pyplot.plot(t, xe, '*-',label='Exact',color="red",linewidth=2.0)
pyplot.xlabel('time ($ \mu s$)',fontsize=20,fontweight='bold')
pyplot.ylabel('Depth of Penetration ($ \mu m$)',fontsize=14,fontweight='bold')
pyplot.legend(loc='upper right', numpoints=1, ncol=2, fancybox=False, shadow=False,fontsize=16)
pyplot.tick_params(axis='both', which='major', labelsize=14)
pyplot.savefig('predict.png')
pyplot.show()
pyplot.close()


np.savetxt('array.csv', [t,x,xe], delimiter=',', fmt='%d')
#plot(t,x,t,xe)
#fprintf('Infiltration at t=%f should be %f',t(n+1),x(n+1))
#ma=[t;x;xe];
#mb=ma';
#csvwrite('validation.csv',mb);
