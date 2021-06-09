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
import openpnm as op
import numpy as np 
import matplotlib.pyplot as plt
np.random.seed(10)

op.Workspace().settings['loglevel'] = 50

Nx = 11 
Ny = 11 
Nz = 11 
spacing=1e-4
pn = op.network.Cubic(shape=[11, 11, 11],spacing=[spacing, spacing, spacing])
#regular network 
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
#geo.show_hist(['pore.diameter', 'throat.diameter', 'throat.length'])
#trim the network (1 in x, 4 in y, 7 in z) 
#op.topotools.trim(network=pn, pores=[1, 4, 7])

print(pn['pore.coords'][pn.pores('all')])
f=np.savetxt("node.dat",(pn['pore.coords'][pn.pores('all')]))

print(pn['pore.diameter'][pn.pores('all')])
f1=np.savetxt("pore_diameter",(pn['pore.diameter']))
#f1=np.savetxt("pore_diamter",(pn['pore.diameter'][pn.pores('all')]))

print(pn['pore.volume'][pn.pores('all')])
f2=np.savetxt("pore_volume",(pn['pore.volume']))
#f1=np.savetxt("pore_diamter",(pn['pore.diameter'][pn.pores('all')]))

print(pn['throat.conns'][pn.pores('all')])
f3=np.savetxt("link.dat",(pn['throat.conns']).astype(int),fmt='%i')
#f2=np.savetxt("link.dat",(pn['throat.conns'][pn.pores('all')]))

print(pn['throat.diameter'][pn.pores('all')])
f4=np.savetxt("throat_diameter",(pn['throat.diameter']))

print(pn['throat.area'][pn.pores('all')])
f5=np.savetxt("throat_area",(pn['throat.area']))

print(pn['throat.conduit_lengths.pore1'][pn.pores('all')])
f6=np.savetxt("throat.conduit_lengths.pore1",(pn['throat.conduit_lengths.pore1']))

print(pn['throat.conduit_lengths.pore2'][pn.pores('all')])
f7=np.savetxt("throat.conduit_lengths.pore2",(pn['throat.conduit_lengths.pore2']))

print(pn['throat.conduit_lengths.throat'][pn.pores('all')])
f8=np.savetxt("throat.conduit_lengths.throat",(pn['throat.conduit_lengths.throat']))

print(pn['throat.volume'][pn.pores('all')])
f8=np.savetxt("throat.volume",(pn['throat.volume']))

op.io.VTK.save(network=pn, filename='test_file')
#op.io.MAT.save(network=pn, filename='test_file')
#op.io.CSV.save(network=pn, filename='test_file')


fig = plt.figure()
plt.hist(pn.num_neighbors(pn.Ps), edgecolor='k')
fig.patch.set_facecolor('white')
#plt.show()
plt.savefig('num_neighbor.png')


#geo.add_model(propname='pore.seed',
#              model=op.models.geometry.pore_size.random,
#              seed=0)

'''
geo.add_model(propname='pore.diameter',
              model=op.models.geometry.pore_size.weibull,
              scale=0.5e-4, shape=0.8, loc=1e-6)

geo.add_model(propname='pore.volume',
              model=op.models.geometry.pore_volume.sphere)
geo.add_model(propname='throat.diameter',
              model=op.models.geometry.throat_size.from_neighbor_pores,
              mode='min')
geo.add_model(propname='throat.length',
              model=op.models.geometry.throat_length.classic)
geo.add_model(propname='throat.volume',
              model=op.models.geometry.throat_volume.cylinder)

'''
#
air = op.phases.Air(network=pn)
#phys = op.physics.Standard2D(network=pn, phase=air, geometry=geo)
phys = op.physics.Standard(network=pn, phase=air, geometry=geo)

sf = op.algorithms.StokesFlow(network=pn, phase=air)
sf.set_value_BC(pores=pn.pores('left'),  values=200000)
sf.set_value_BC(pores=pn.pores('right'), values=100000)
sf.run()
air.update(sf.results())
K = sf.calc_effective_permeability()
print(K)

print(pn.settings)
print(sf.settings)
'''
sim = sf['pore.pressure'][pn.pores('internal')]
temp_map = np.reshape(a=sim, newshape=[11,11,1])
plt.subplots(1, 1, figsize=(10, 5))
plt.imshow(temp_map, cmap=plt.cm.plasma);
plt.colorbar();
plt.show()
'''

'''
c=c=air['pore.pressure']
c2d = c.reshape(pn._shape).squeeze()
plt.imshow(np.rot90(c2d))
plt.colorbar()
plt.show()
'''

geo.show_hist(['pore.diameter', 'throat.diameter', 'throat.length'])
#plt.hist(pn['pore_coords'])
#plt.show() 

