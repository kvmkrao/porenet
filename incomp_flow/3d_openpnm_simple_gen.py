''' Author: Murali Kotteda 
date : May 27, 2021       
openPNM : create a 3D network

command to install openpnm
$conda install -c conda-forge openpnm

execute 
$python 3d_openpnm_simple2.py 

instructions to load the vtp files in paraview can be found at 
https://openpnm.readthedocs.io/en/latest/getting_started/visualize_in_paraview.html

To visualize the pore data, we need to add some glyphs. First click on the Glyph button in the tool bar. Then, you can plot the pore data as spheres, where their size and/or their color can be mapped to some variables. In the images below spherical glyphs are assigned to the pores, where the diameter is linked to the pore diameters and the color to the concentration. Clicking on the Apply button renders these settings.

To visualize the throat data (like diameters or molar flow rates) we need to set up the following filters in Paraview. First, use the Shrink filter and set it up to 1. Then, the cell data needs to be transposed to the data point with CellDatatoPointdata filter. Then extract the surface with the filter ExtractSurface. Finally, the data may be plotted as tube by using the Tube filter. As previously for the pore data, either the Tube radius or color can be linked to throat data.

Throat data plotted with tubes radius proportional to throat diameter and tubes color linked to the throat mole rate:
'''

import scipy as sp
import matplotlib.pyplot as plt
import os
import warnings
import numpy as np
import openpnm as op

np.random.seed(11)
ws = op.Workspace()
ws.settings['loglevel'] = 40
warnings.filterwarnings('ignore')

pn = op.network.Cubic(shape=[11, 11, 11], spacing=1e-7, connectivity=6)
print(pn['throat.conns'][pn.pores('all')])
f3=np.savetxt("link.dat",(pn['throat.conns']).astype(int),fmt='%i')
#f2=np.savetxt("link.dat",(pn['throat.conns'][pn.pores('all')]))

h = pn.check_network_health()
print(h) 

geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

air = op.phases.Water(network=pn)
#air = op.phases.Air(network=pn)

#air['pore.temperature'] = 298.15
air['pore.viscosity']   = 1.5e-5

#phys = op.physics.Standard2D(network=pn, phase=air, geometry=geo)
phys = op.physics.Standard(network=pn, phase=air, geometry=geo)

sf = op.algorithms.StokesFlow(network=pn, phase=air)
pinlet  = 101325*5.0
poutlet = 101325*10.0
sf.set_value_BC(pores=pn.pores('top'),  values=pinlet)
sf.set_value_BC(pores=pn.pores('bottom'), values=poutlet)
sf.run()
air.update(sf.results())
K = sf.calc_effective_permeability()
print(K)

fig = plt.hist(geo['pore.diameter'], bins=25, edgecolor='k')
plt.savefig('pore_dia_hist.png')

#Q = perm.rate(pores=pn.pores('bottom'), mode='group')
#A = (im.shape[0] * im.shape[1]) * resolution**2
#L = im.shape[2] * resolution
#mu = air['pore.viscosity'].max()
#delta_P = 101325 - 0
#K = Q * L * mu / (A * delta_P)
print('The value of K is:', K/0.987e-12*1000, 'mD')


'''
print(pn['pore.coords'][pn.pores('all')])
f1=np.savetxt("node.dat",(pn['pore.coords'][pn.pores('all')]))

print(pn['throat.conns'][pn.pores('all')])
f3=np.savetxt("link.dat",(pn['throat.conns']).astype(int),fmt='%i')
#f2=np.savetxt("link.dat",(pn['throat.conns'][pn.pores('all')]))

Nx = 11 
Ny = 11 
Nz = 11 
spacing=1e-7
'''

#mu, sigma = 3., 1. # mean and standard deviation
#s = np.random.lognormal(mu, sigma, 1331)
#count, bins, ignored = plt.hist(s, 100, density=True, align='mid')

mean   = 50
stddev = 10
sigma_squared = np.log((stddev/mean)**2 + 1)
mu = np.log(mean) - 0.5*sigma_squared
sigma = np.sqrt(sigma_squared)

sample = np.random.lognormal(mu, sigma, size=1331)

print(np.mean(sample))
print(np.std(sample))
print(np.min(sample), np.max(sample))

print(len(sample)) 

#print(len(pn['throat.diameter'])) 
#pn['throat.diameter'] = sample
#f4=np.savetxt("throat_diameter",(pn['throat.diameter']))

op.io.VTK.save(network=pn, filename='test_file') 

#xc = np.linspace(0, 10, num=10)
#yc = np.linspace(0, 10, num=10)
#zc = np.linspace(0, 10, num=10)

#b = np.array([np.linspace(1,10,num=10), np.linspace(1,10,num=10), np.linspace(1,10,num=10)])
#print(b) 
#xc= np.zeros(1331)
#yc= np.zeros(1331)
#zc= np.zeros(1331)


mu, sigma = 50, 10 # mean and standard deviation
s = np.random.normal(mu, sigma, 1331)
#Verify the mean and the variance:
print('error mean', abs(mu - np.mean(s)))
print('error std',  abs(sigma - np.std(s, ddof=1)))
#Display the histogram of the samples, along with the probability density function:
count, bins, ignored = plt.hist(s, 30, density=True)
#fig=plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
            np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
            linewidth=2, color='r')
#plt.savefig('pore_dia_hist1.png')
plt.show()


# nodes information 
spacing= 1e-7
file = open('mnodes.dat', 'w')
for k in range(11):
    for i in range(11): 
        for j in range(11): 
            xc = float(i)*100 #*spacing
            yc = float(j)*100 #*spacing 
            zc = float(k)*100 #*spacing 
            #print(i,j,k)
            #print(xc,yc,zc)
            file.write(" %f %f %f\n" % (xc, yc, zc ))
            #file.write(" %d %d %d\n" % (xc, yc, zc ))
            

# diameter information 
file1 = open('mdiameter.dat', 'w') 
for i in range(1331):
    file1.write(" %f \n" % (sample[i]))

# connection information 
file2 = open('mlinks.dat', 'w') 
count=0 

'''
for k in range(10): 
    for j in range(10): 
        node1 = 
'''


for k in range(11):
    for i in range(11):
        for j in range(11):
            node  = k*11*11+ i*11 + j

            node1 = k*11*11+ i*11 + j+1    # + k direction 
            if(j!=10):
                file2.write("%d %d\n" %(node, node1))

#            node2 = k*11*11+ j*11 + i-1      # -k direction 
#            if(i!=0):
#                file2.write("%d %d\n" %(node, node2))

            node3 = k*11*11+ (i+1)*11 + j  # +j direction 
            if(i!=10):
                file2.write("%d %d\n" %(node, node3))

#            node4 = k*11*11+ (j-1)*11 + i  # -j direction 
#            if(j!=0):
#                file2.write("%d %d\n" %(node, node4))

            node5 = (k+1)*11*11+ i*11 + j  # +i direction 
            if(k!=10):
                file2.write("%d %d\n" %(node, node5))

#            node6 = (k-1)*11*11+ j*11 + i # -i direction 
#            if(k!=0):
#                file2.write("%d %d\n" %(node, node6))
                
            #print(node,node1,node2,node3, node4, node5,node6)
            #       1330 1331 1329  1341    1319   1451 1209
            #print(i,j,k)
            #print(xc,yc,zc)
            #file.write(" %0.7f %0.7f %0.7f\n" % (xc, yc, zc ))
