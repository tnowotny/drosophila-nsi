import numpy as np
import matplotlib.pyplot as plt
import sys


fname= sys.argv[1]

name= [ "pM1", "pM2", "pM3", "pM4", "pM5", "pM6", "g_M", "cmid_a", "cslope", "g_r", "err" ] 
para= np.loadtxt(fname, dtype= np.float32)
print(para.shape)

wd = 4
ht = 3
plt.figure(figsize=(20,10))
for x in range(wd):
    for y in range(ht):
        i = y*wd+x
        if i < para.shape[1]:
            plt.subplot(ht,wd,i+1)
            plt.plot(para[:,i])
            plt.plot(np.ones(para.shape[0])*para[0,i])
            plt.title(name[i])
plt.show()
