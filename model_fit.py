import numpy as np
import matplotlib.pyplot as plt
from utils import *
from droso_AB3_asym import *
from brian2modelfitting import *
from pandas import read_csv

df = read_csv("ab3A_conc.csv")

print(df)

conc= ["-3", "-4", "-5", "-6", "-7" ]

"""
# This is for initial testing only:
for c in conc:
    plt.figure()
    plt.plot(df["time"], df[c])
    sT= spikes(df["time"],df[c])
    sT= np.asarray(sT)*1000
    nsdf= sdf(sT, dt= 10, sigma= 50)
    plt.plot(nsdf["t"]/1000,nsdf["sdf"])
    plt.figure()
    plt.plot(sT, np.ones(len(sT)), "|")
   
plt.show()
"""

# extract the output spikes
out_spikes= []
for c in conc:
    sT= spikes(df["time"],df[c])
    sT= np.asarray(sT)*1000*ms
    sT= sT-sT[0]
    out_spikes.append(sT)


# make the input traces
time= np.asarray(df["time"])
steps= []
steps.append(int((-time[0])/(dt/second)))
steps.append(int(5.0/(dt/second)))
steps.append(int((time[-1]-5.0)/(dt/second)))
inp_traces= []
for c in conc:
    conc= 10**int(c)
    dd= []
    dd.append(np.zeros(steps[0]))
    dd.append(np.ones(steps[1])*conc)
    dd.append(np.zeros(steps[2]))
    inp_traces.append(np.hstack(dd))
inp_traces= np.vstack(inp_traces)

n_opt = SkoptOptimizer('ET')
metric = GammaFactor(delta=2*ms, time= inp_traces.shape[1]*dt)

init_sens= {
    "v_a": El,
    "m_a": 0.0,
    "h_a": 1.0,
    "n_a": 0.0,
    "v_b": El,
    "m_b": 0.0,
    "h_b": 1.0,
    "n_b": 0.0,
    "refr_b": 0*ms
}

fitter = SpikeFitter(model=eqs,
                     input=inp_traces,
                     input_var='conc_a',
                     output=out_spikes,
                     dt=dt,
                     n_samples=30,
                     threshold='v_a>-20*mV',
                     reset='',
                     method='milstein',
                     param_init= init_sens)

results, error = fitter.fit(n_rounds=2,
                            optimizer=n_opt,
                            metric=metric,
                            pM1=[0.01,0.1],
                            pM2=[-10*mV,40*mV],
                            pM3=[1*mV,10*mV],
                            pM4=[0.0001,0.01],
                            pM5=[-60*mV,0*mV],
                            pM6=[-10*mV,-1*mV],
                            g_M=[1e-4*msiemens,11e-4*msiemens]
                            )
