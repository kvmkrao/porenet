
import scipy as sp
import matplotlib.pyplot as plt
import os
import warnings
import numpy as np
import openpnm as op

import math
from scipy import stats

#b = []
#
#for i in range(1331):
#   a = 5. + np.random.standard_normal(100)
#   b.append(np.product(a))
#
#b = np.array(b) / np.min(b) # scale values to be positive
#count, bins, ignored = plt.hist(b, 100, density=True, align='mid')
#sigma = np.std(np.log(b))
#mu = np.mean(np.log(b))

x = np.linspace(2, 100, 1331)

sigma = 5 
mu    = 20 
#print(min(bins)) 
#print(max(bins)) 

print(sigma) 
print(mu) 

pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
       / (x * sigma * np.sqrt(2 * np.pi)))

file1 = open('mdiameter.dat', 'w')
#file1.write(str(pdf)) 
#for i in range(1331):
#    file1.write(" %s \n" % (str(pdf[i])))
#    #file1.write(" %s \n" % (str(pdf[i])))

plt.plot(x, pdf, color='r', linewidth=2)
plt.show()

mu, sigma = 55, 16  # mean and standard deviation
snormal = np.random.normal(mu, sigma, 1331)
for i in range(1331):
    file1.write(" %f \n" % (snormal[i]))
#    file1.write(" %0.9f \n" % (snormal[i]*1e-9))



#mu, sigma = 1.30, 1.5  # mean and standard deviation
#snormal = np.random.lognormal(mu, sigma, 1331)
#for i in range(1331):
#    file1.write(" %f \n" % (snormal[i]))
##    file1.write(" %0.9f \n" % (snormal[i]*1e-9))

# standard deviation of normal distribution
sigma = 0.859455801705594
# mean of normal distribution
mu = 0.418749176686875
# hopefully, total is the value where you need the cdf
total = 37

frozen_lognorm = stats.lognorm(s=sigma, scale=math.exp(mu))
#frozen_lognorm.cdf(total) # use whatever function and value you need here
#plt.plot(frogen_lognorm) 
xc = frozen_lognorm.cdf(1331) # use whatever function and value you need here
#print(frogen_lognorm) 
#plt.plot(xc) 
#plt.show()

