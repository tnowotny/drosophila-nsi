from droso_AB3_asym_fit import *
from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv

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

def fit_fun(x):
    reps= len(the_conc)
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
    for i in range(reps):
        net.run(p["tin"])
        net.set_states({'sensillum': {"conc_a": the_conc[i]}})
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
    err= np.linalg.norm(y)
    if resfile is not None:
        resfile.write(f"{err}\n")
        resfile.flush()        
    return err

if __name__ == '__main__':
    resfile= open("fit_trajectory_PSO.txt","w")
    x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
    bound_low= [ 0.03,  10,  5, 0.0007, -40, -5, 3e-4, -6.0, 0.5,  2e-4 ]
    bound_high= [ 0.06, 30, 9, 0.0015, -20, -3, 5e-3, -2.0, 2, 8e-3]

    from sko.PSO import PSO
    pso = PSO(func=fit_fun, n_dim=10, pop=1000, max_iter=150, lb=bound_low, ub=bound_high, w=0.8, c1=0.5, c2=0.5)
    pso.run()
    print('best_x is ', pso.gbest_x, 'best_y is', pso.gbest_y)
    resfile.close()
