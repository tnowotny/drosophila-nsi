from droso_AB3_asym_fit import *
from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.optimize import minimize
from brian2 import seed

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
    "tin": -time[0]*second-10*ms,
    "tdur": 5.0*second,
    "tout": (time[-1]-5.0)*second+10*ms
}

def fit_fun(x, trg_sdf, conc, net, p, resfile, plot= False, std= 50.0, repetitions= 5):
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
    sT= []
    for i in range(reps):
        sT.append([])
    for n in range(repetitions):
        for i in range(reps):
            seed(n*11)
            sensillum.v_a = El + (randn() * 5 - 5)*mV
            sensillum.v_b = El + (randn() * 5 - 5)*mV
            net.run(p["tin"])
            net.set_states({'sensillum': {"conc_a": conc[i]}})
            net.run(p["tdur"])
            net.set_states({'sensillum': {"conc_a": 0.0}})
            net.run(p["tout"])
            sT[i].append(net.get_states()["spikeA"]["t"]/ms)
            net.restore()
    for i in range(reps):
        sT[i]= np.hstack(sT[i])
        rsdf= sdf(sT[i],t0= 0.0, tmax= (p["tin"]+p["tdur"]+p["tout"])/ms, dt=10.0, sigma=std/repetitions)
        rsdf= rsdf["sdf"]/repetitions
        if plot:
            plt.plot(rsdf)
            plt.plot(trg_sdf[i])            
        y.append(trg_sdf[i]-rsdf[:len(trg_sdf[i])])
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
        

def fit_fun_v2(x, trg_sdf, conc, net, p, resfile, plot= False, repetitions= 10):
    # use repititions runs and do a histogram with 10 ms bins
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
    sT= []
    for i in range(reps):
        sT.append([])
    for n in range(repetitions):
        for i in range(reps):
            seed(n*11)
            sensillum.v_a = El + (randn() * 5 - 5)*mV
            sensillum.v_b = El + (randn() * 5 - 5)*mV
            net.run(p["tin"])
            net.set_states({'sensillum': {"conc_a": conc[i]}})
            net.run(p["tdur"])
            net.set_states({'sensillum': {"conc_a": 0.0}})
            net.run(p["tout"])
            sT[i].append(net.get_states()["spikeA"]["t"])
            net.restore()
    for i in range(reps):
        sT[i]= np.hstack(sT[i])
        sdf, edges= np.histogram(sT[i]/ms, int((p["tin"]+p["tdur"]+p["tout"])/(10*ms)), range=(0.0,(p["tin"]+p["tdur"]+p["tout"])/ms))
        sdf= sdf.astype(float)
        sdf *= 100.0/repetitions
        if plot:
            plt.plot(sdf)
            plt.plot(trg_sdf[i])
        y.append(trg_sdf[i]-sdf[:len(trg_sdf[i])])
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

if __name__ == '__main__':
    x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
    bounds= [ (0.02, 0.08), (10, 30), (5,9), (0.0005,0.002),
              (-40,-20), (-5,-3), (3e-4,5e-3), (-6.0, -2.0),
              (0.5, 2), (2e-4, 8e-3) ]

    resfile= open("fit_trajectory.txt","w")
    res_lsq= minimize(fit_fun,x0, args=(trg_sdf,the_conc,net,p,resfile),bounds=bounds,method='Nelder-Mead',options={'xatol': 1e-4})

    print(res_lsq)
    resfile.close()
