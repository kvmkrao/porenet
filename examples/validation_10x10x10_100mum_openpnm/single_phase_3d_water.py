#!/usr/bin/env python
# coding: utf-8

# In[1]:


import openpnm as op
import openpnm.models.geometry as gmods
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from scipy.stats import lognorm
from scipy.stats import norm

# Define geometrical parameters
Lc = 1e-4
Nx, Ny, Nz = (10, 10, 10)
Pin, Pout = (2000000, 1000000)
# Generate network, geometry, phase, and physics
#pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc,connectivity=26)
#op.topotools.reduce_coordination(network=pn, z=5)

pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc,connectivity=6)
#Create a Geometry Object and Assign Geometric Properties to Pores and Throats
#The Network pn does not contain any information about pore and throat sizes at this point. The next step is to create
#a Geometry object to manage the geometrical properties.
#geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)


# test case https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/porenetwork_upscaling/README.md 
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

npr = pn.num_pores()
print(pn.num_pores())
print(pn.num_throats())
print(pn.props())
print(pn.pores('left'))
print(pn.labels())

geo['pore.old_diameter'] = geo.pop('pore.diameter')
#geo.add_model(propname='pore.diameter',
#              model=gmods.pore_size.weibull,
#              shape=0.5, loc=0, scale=1e-5)
sd = 10
mean = 20
geo.add_model(propname='pore.diameter',
              model=gmods.pore_size.normal,
              loc=mean, scale=sd)
              #loc=5e-5, scale=1e-5)

#s = 0.954
#x = np.linspace(lognorm.ppf(0.01, s),lognorm.ppf(0.99, s), npr)
#mean, var, skew, kurt = 20, 10, 0, 0 #lognorm.stats(s, moments='mvsk')
#geo['pore.diameter'] = lognorm.pdf(x, s)*0.0001
#plt.plot(x, lognorm.pdf(x, s),'r-', lw=5, alpha=0.6, label='lognorm pdf')

x = np.linspace(4, 100, npr)
#y = norm.pdf(x, 20, 10)
geo['pore.diameter'] = norm.pdf(x,20,5)*0.001
plt.plot(x, geo['pore.diameter'],'r-', lw=5, alpha=0.6, label='norm pdf')
plt.savefig('pore_diameter.png')

for i in range(npr):
    #print(geo['pore.diameter'][i])
    if geo['pore.diameter'][i] < 2e-6: 
        geo['pore.diameter'][i] = geo['pore.diameter'].mean()

geo.add_model(propname='throat.diameter',
              model=op.models.misc.from_neighbor_pores,prop='pore.diameter',mode='min')

#geo['throat.diameter'] = geo['throat.diameter']/2.0

#geo.add_model(propname='throat.endpoints',
#                model=op.models.geometry.throat_endpoints.spherical_pores)
#geo.add_model(propname='throat.area',
#                model=op.models.geometry.throat_cross_sectional_area.cylinder)
#geo.add_model(propname='pore.area',
#                model=op.models.geometry.pore_cross_sectional_area.sphere)
#geo.add_model(propname='throat.conduit_lengths',
#                model=op.models.geometry.throat_length.conduit_lengths)
        
#fig = plt.hist(geo['pore.diameter'], bins=15,density=True, edgecolor='k', alpha=0.5) 
print("mean pore diameter", geo['pore.diameter'].mean()) 
print("min  pore diameter", geo['pore.diameter'].min())
print("max  pore diameter", geo['pore.diameter'].max())

print("mean throat diameter", geo['throat.diameter'].mean()) 
print("min  throat diameter", geo['throat.diameter'].min())
print("max  throat diameter", geo['throat.diameter'].max())


#Hg = op.phases.Mercury(network=pn)
#phys = op.physics.Standard(network=pn, phase=Hg, geometry=geo)

# Create algorithm and run simulation
#mip = op.algorithms.Porosimetry(network=pn)
#mip.setup(phase=Hg)
#mip.set_inlets(pores=pn.pores(['left', 'right']))
#mip.run()
print(pn.check_network_health())
print(geo.models)


# In[2]:


# Generate phase and physics
water = op.phases.Water(network=pn)
water['pore.viscosity'] =0.001
phys = op.physics.Standard(network=pn, phase=water, geometry=geo)

# Create algorithm, set boundary conditions and run simulation
sf = op.algorithms.StokesFlow(network=pn, phase=water)

sf.set_value_BC(pores=pn.pores('left'), values=Pin)
sf.set_value_BC(pores=pn.pores('right'), values=Pout)
sf.run()
print(water)


# In[8]:


Q = sf.rate(pores=pn.pores('left'))
A = Ny*Nz*Lc**2
L = Nx*Lc
mu = water['pore.viscosity'].mean()
K = Q*mu*L/(A*(Pin-Pout))
print(K, Q, mu, L, A, Pin-Pout)
print("effective permeability",sf.calc_effective_permeability())
print(Q)

print(pn.shape.max())
print(pn.spacing.max()) 
vol_total = (pn._shape * pn._spacing).prod()
vol_pores = geo['pore.volume'].sum()
vol_throats = geo['throat.volume'].sum()
porosity = (vol_pores + vol_throats) / vol_total * 100
print('Porosity (%):', porosity)


# In[4]:


mip = op.algorithms.Porosimetry(network=pn, phase=water)
mip.set_inlets(pores=pn.pores('left'))
mip.run()
fig = mip.plot_intrusion_curve()


# In[5]:


fig = plt.figure()
plt.hist(pn.num_neighbors(pn.Ps), edgecolor='k')
fig.patch.set_facecolor('white')
print("mean co-ordination number",pn.num_neighbors(pn.Ps).mean())
print("minimum co-ordination number",pn.num_neighbors(pn.Ps).min())
print("maximum co-ordination number",pn.num_neighbors(pn.Ps).max())


# In[ ]:





# In[ ]:





# In[ ]:




