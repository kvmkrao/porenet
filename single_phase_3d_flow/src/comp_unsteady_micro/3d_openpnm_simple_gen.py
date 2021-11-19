''' Author: Murali Kotteda 
date : May 27, 2021       
openPNM : create a 3D network

command to install openpnm
$conda install -c conda-forge openpnm

execute 
$python 3d_openpnm_simple_gen.py 

instructions to load the vtp files in paraview can be found at 
https://openpnm.readthedocs.io/en/latest/getting_started/visualize_in_paraview.html

To visualize the pore data, we need to add some glyphs. First click on the Glyph button in the tool bar. Then, you can plot the pore data as spheres, where their size and/or their color can be mapped to some variables. In the images below spherical glyphs are assigned to the pores, where the diameter is linked to the pore diameters and the color to the concentration. Clicking on the Apply button renders these settings.

To visualize the throat data (like diameters or molar flow rates) we need to set up the following filters in Paraview. First, use the Shrink filter and set it up to 1. Then, the cell data needs to be transposed to the data point with CellDatatoPointdata filter. Then extract the surface with the filter ExtractSurface. Finally, the data may be plotted as tube by using the Tube filter. As previously for the pore data, either the Tube radius or color can be linked to throat data.

Throat data plotted with tubes radius proportional to throat diameter and tubes color linked to the throat mole rate:
'''

#import scipy as sp
#import matplotlib.pyplot as plt
import os
import warnings
import numpy as np
import openpnm as op
import scipy.stats
import math 

np.random.seed(11)
ws = op.Workspace()
ws.settings['loglevel'] = 40
warnings.filterwarnings('ignore')
nx = 10 
ny = 10 
nz = 10 

pn = op.network.Cubic(shape=[nx, ny, nz], spacing=1e-7, connectivity=6)

with open("link.dat", "w") as f:
    f.write(str(len(pn['throat.conns']))+'\n')
    np.savetxt(f, (pn['throat.conns']).astype(int), fmt='%i', delimiter="\t")

#lower = np.log(2)
lower = np.log(10)
upper = np.log(90)
#upper = np.log(100)
mu = np.log(20)
sigma = np.log(20)
N = nx*ny*nz 

samples = scipy.stats.truncnorm.rvs(
          (lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=N)

#plt.hist(np.exp(samples), bins= 50,alpha=0.3,label='histogram');
np.savetxt('node_dia.dat', np.exp(samples), delimiter=',',fmt='%8.6f')   # X is an array
#plt.show()

spacing= 100
file = open('node.dat', 'w')
file.write(str(N)+'\n')
for k in range(nz):
    for i in range(nx): 
        for j in range(ny): 
            xc = float(i+1)*spacing
            yc = float(j+1)*spacing 
            zc = float(k+1)*spacing 
            #print(i,j,k)
            #print(xc,yc,zc)
            file.write(" %d %d %d\n" % (xc, yc, zc ))

porosity = 0 
for k in range(N):
    tmp = 0.5*np.exp(samples[k])
    porosity = porosity + (4.0/3.0)*math.pi* tmp*tmp*tmp
    #porosity = porosity + math.pi* 100.0* tmp*tmp
    #print(np.exp(samples[k]))

print("porosity (nodes)")
print(porosity/(1000*1000*1000))


N1 = len(pn['throat.conns'])
file1 = open('link_dia.dat', 'w')
file1.write(str(N1)+'\n')

lower = np.log(8)
upper = np.log(40)
mu = np.log(20)
sigma = np.log(20)

samples = scipy.stats.truncnorm.rvs(
          (lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=N1)

porosity = 0 
for k in range(N1):
    file1.write(str(np.exp(samples[k]))+'\n')
    tmp = 0.5*np.exp(samples[k])
     # volume pi*r^2*h 
    porosity = porosity + math.pi* tmp*tmp*100.0 

print("porosity (throats)")
print(porosity/(nx*100*ny*100*nz*100))

