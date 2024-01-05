from droso_AB3_asym_fit import *
from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.optimize import minimize

net.store()

df = read_csv("ab3A_conc.csv")
print(df)
conc= [-3, -4, -5, -6, -7 ]

trg_sdf= []
for c in conc:
    trg_sdf.append(np.asarray(df[str(c)]))
trg_sdf= np.vstack(trg_sdf)

the_conc= 10.0**np.asarray(conc)

time= np.asarray(df["time"])
p= {
    "tin": -time[0]*second,
    "tdur": 5.0*second,
    "tout": (time[-1]-5.0)*second
}

def fit_fun(x, trg_sdf, conc, net, p, resfile, plot= False):
    reps= len(conc)
    print(x)
    if resfile is not None:
        for xi in x:
            resfile.write(f"{str(xi)} ")
    pdict= {
        'pM1': x[0],
        'pM2': x[1]*mV,
        'pM3': x[2]*mV,
        'pM4': x[3],
        'pM5': x[4]*mV,
        'pM6': x[5]*mV,
        'g_M': x[6]*msiemens,
        'cmid_a': x[7],
        'cslope': x[8],
        'g_r': x[9]*msiemens,
    }
    net.set_states({'sensillum': pdict})
    net.store()
    y= []
    if plot:
        plt.figure()
    for i in range(reps):
        net.run(p["tin"])
        net.set_states({'sensillum': {"conc_a": conc[i]}})
        net.run(p["tdur"])
        net.set_states({'sensillum': {"conc_a": 0.0}})
        net.run(p["tout"])
        sT= net.get_states()["spikeA"]["t"]
        rsdf= sdf(sT/ms,t0= 0.0, tmax= (p["tin"]+p["tdur"]+p["tout"])/ms, dt=10.0, sigma=50.0)
        if plot: 
            plt.plot(trg_sdf[i])
            plt.plot(rsdf["sdf"])
        y.append(trg_sdf[i]-rsdf["sdf"][:len(trg_sdf[i])])
        net.restore()
    y= np.hstack(y)
    if plot:
        plt.figure()
        plt.plot(y)
        plt.show()
    err= np.linalg.norm(y)
    if resfile is not None:
        resfile.write(f"{err}\n")
        resfile.flush()        
    return err

x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
x0= [ 0.04085458101661702, 23.324173009601324, 6.467984263287498, 0.001005189468563782, -32.19489404645946, -4.247262362895521, 0.000970760277054877, -3.071751510549535, 0.7533306082531139, 0.0012243353684708985, 660.4355835454747 ]
bounds= [ (0.02, 0.08), (10, 30), (5,9), (0.0005,0.002), (-40,-20), (-5,-3), (3e-4,5e-3),
          (-6.0, -2.0), (0.5, 2), (2e-4, 8e-3) ]

# For quick test of a parameter combination
fit_fun(x0, trg_sdf, the_conc, net, p, None, plot= True)
exit(1)

resfile= open("fit_trajectory.txt","w")
# For actual fitting
res_lsq= minimize(fit_fun,x0, args=(trg_sdf,the_conc,net,p,resfile),bounds=bounds,method='Nelder-Mead',options={'xatol': 1e-5})

print(res_lsq)
resfile.close()
