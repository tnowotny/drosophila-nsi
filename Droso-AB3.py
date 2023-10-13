from brian2 import *
import matplotlib.pyplot as plt
import numpy as np


def sdf(sT, t0, tmax, dt, sigma):
    tleft= t0-4*sigma
    tright= tmax+4*sigma
    n= int((tright-tleft)/dt)
    sdf= np.zeros(n)
    kwdt= 3*sigma
    i= 0
    x= np.arange(-kwdt,kwdt,dt)
    x= np.exp(-np.power(x,2)/(2*sigma*sigma))
    x= x/(sigma*np.sqrt(2.0*np.pi))*1000.0
    if sT is not None:
        for t in sT:
            if (t > t0 and t < tmax): 
                left= int((t-tleft-kwdt)/dt)
                right= int((t-tleft+kwdt)/dt)
                if left >= 0 and right <= n:
                    sdf[left:right]+=x
           
    return sdf
                  





# Parameters

Cm = 1*nfarad
gl = 0.00018685*msiemens

El = -63.563*mV
EK = -95*mV
ENa = 50*mV
g_na = 0.05*msiemens
g_kd = 0.01*msiemens
dt= 0.05*msecond
I_noise= 1*namp*dt**0.5  # note: not sure this is rescaling correctly for different timesteps
I_base= 0.4*namp
I_off= I_base
VT = -65*mV
pM1 = 0.1
pM2 = 0*mV
pM3 = 10*mV
pM4 = 0.001
pM5 = -30*mV
pM6 = -10*mV
g_M= 3e-5*msiemens
cmid= -3   # log concentration of mid-activation
g_r= 1e-5*msiemens
conc_a= 0
conc_b= 0

# sensillum
gl_sl = 0.001*msiemens
El_sl = 40*mV
Cm_sl = 50*nfarad

"""
v_a= np.arange(-100,100,0.1)*mV
plt.figure()
plt.plot(v_a,pM1/(1+exp((pM2-v_a)/pM3)))
plt.plot(v_a,pM4/(1+exp((pM5-v_a)/pM6)))
plt.show()
exit(1)
"""

"""
area = 20000*umetre**2
Cm = (1*ufarad*cm**-2) * area
gl = (5e-5*siemens*cm**-2) * area  #(5e-5*siemens*cm**-2) * area

El = -60*mV
EK = -90*mV
ENa = 50*mV
g_na = (100*msiemens*cm**-2) * area
g_kd = (30*msiemens*cm**-2) * area
dt= 0.01*msecond
I_noise= 5e-1*namp*dt**0.5 # (5e-3*namp*dt**0.5)  # note: not rescaling this to be valid for different timesteps
I_off= -0.035*namp
VT = -63*mV
"""


# The model
eqs = Equations('''
dv_a/dt = (gl*(El-v_a)-
         g_na*(m_a*m_a*m_a)*h_a*(v_a-ENa)-
         g_kd*(n_a*n_a*n_a*n_a)*(v_a-EK)-
         g_M*M_a*(v_a-EK)-
         g_r*s_r_a*(v_a-ENa-v_sl)+I_off+I_noise*xi_a)/Cm : volt
dm_a/dt = alpha_m_a*(1-m_a)-beta_m_a*m_a : 1
dn_a/dt = alpha_n_a*(1-n_a)-beta_n_a*n_a : 1
dh_a/dt = alpha_h_a*(1-h_a)-beta_h_a*h_a : 1
dM_a/dt = alpha_M_a*(1-M_a)-beta_M_a*M_a : 1


alpha_m_a = 0.32*(mV**-1)*4*mV/exprel((13*mV-v_a+VT)/(4*mV))/ms : Hz
beta_m_a = 0.28*(mV**-1)*5*mV/exprel((v_a-VT-40*mV)/(5*mV))/ms : Hz
alpha_h_a = 0.128*exp((17*mV-v_a+VT)/(18*mV))/ms : Hz
beta_h_a = 4./(1+exp((40*mV-v_a+VT)/(5*mV)))/ms : Hz
alpha_n_a = 0.032*(mV**-1)*5*mV/exprel((15*mV-v_a+VT)/(5*mV))/ms : Hz
beta_n_a = .5*exp((10*mV-v_a+VT)/(40*mV))/ms : Hz
alpha_M_a = pM1/(1+exp((pM2-v_a)/pM3))/ms : Hz
beta_M_a = pM4/(1+exp((pM5-v_a)/pM6))/ms : Hz
s_r_a = 1.0/(1.0+exp((cmid-log10(conc_a))/2)) : 1

dv_b/dt = (gl*(El-v_b)-
         g_na*(m_b*m_b*m_b)*h_b*(v_b-ENa)-
         g_kd*(n_b*n_b*n_b*n_b)*(v_b-EK)-
         g_M*M_b*(v_b-EK)-
         g_r*s_r_b*(v_a-ENa-v_sl)+I_off+I_noise*xi_b)/Cm : volt
dm_b/dt = alpha_m_b*(1-m_b)-beta_m_b*m_b : 1
dn_b/dt = alpha_n_b*(1-n_b)-beta_n_b*n_b : 1
dh_b/dt = alpha_h_b*(1-h_b)-beta_h_b*h_b : 1
dM_b/dt = alpha_M_b*(1-M_b)-beta_M_b*M_b : 1
drefr_b/dt = -1 : second

alpha_m_b = 0.32*(mV**-1)*4*mV/exprel((13*mV-v_b+VT)/(4*mV))/ms : Hz
beta_m_b = 0.28*(mV**-1)*5*mV/exprel((v_b-VT-40*mV)/(5*mV))/ms : Hz
alpha_h_b = 0.128*exp((17*mV-v_b+VT)/(18*mV))/ms : Hz
beta_h_b = 4./(1+exp((40*mV-v_b+VT)/(5*mV)))/ms : Hz
alpha_n_b = 0.032*(mV**-1)*5*mV/exprel((15*mV-v_b+VT)/(5*mV))/ms : Hz
beta_n_b = .5*exp((10*mV-v_b+VT)/(40*mV))/ms : Hz
alpha_M_b = pM1/(1+exp((pM2-v_b)/pM3))/ms : Hz
beta_M_b = pM4/(1+exp((pM5-v_b)/pM6))/ms : Hz
s_r_b= 1/(1+exp((cmid-log10(conc_b))/2)) : 1

dv_sl/dt = gl_sl*(El_sl-v_sl)/Cm_sl+(g_r*s_r_a*(v_a-ENa-v_sl)+g_r*s_r_b*(v_b-ENa-v_sl))/Cm : volt
''')

P = NeuronGroup(1, model=eqs, threshold='v_a>-20*mV', events={'b_spike': '(v_b>-20*mV) and (refr_b <= 0*ms)'},refractory=3*ms,
                method='milstein',dt=dt)

P.run_on_event('b_spike','refr_b= 3*ms')

# Initialization
P.v_a = 'El + (randn() * 5 - 5)*mV'
P.m_a= 0.0
P.h_a= 1.0
P.n_a= 0.0
P.v_b = 'El + (randn() * 5 - 5)*mV'
P.m_b= 0.0
P.h_b= 1.0
P.n_b= 0.0
P.refr_b= 0*ms

# Record a few traces
spikeA = SpikeMonitor(P)
spikeB = EventMonitor(P,'b_spike')

trace = StateMonitor(P, ('v_a','v_b','v_sl'), record= [0])
for i in range(20):
    #I_off= I_base + i*0.05*namp
    conc_a= 10**(-(20-i)/4)
    run(3.5 * second, report='text')
    #_off= I_base
    conc_b= 10**-1
    run(1 * second, report='text')
    conc_b= 0.0
    run(0.5 * second, report='text')
    conc_a= 0.0
    run(5 * second, report='text')
figure()
plot(trace.t/ms, trace[0].v_sl/mV,lw=0.2)
xlabel('t (ms)')
ylabel('"LFP" (mV)')
savefig('lfp',dpi=600)

#plot(trace.t/ms, trace[0].v_a/mV)
#plot(trace.t/ms, trace[0].v_b/mV)
#xlabel('t (ms)')
#ylabel('v_a, v_b (mV)')
#figure()
#plot(trace.t/ms, trace[0].M_a)
#plot(trace.t/ms, trace[0].M_b)

figure()
scatter(spikeA.t/ms,spikeA.i)
scatter(spikeB.t/ms,spikeB.i+1)

sdf_a= sdf(spikeA.t/ms,0.0,2e5,2,100)
sdf_b= sdf(spikeB.t/ms,0.0,2e5,2,100)

figure()
plot(sdf_a,lw=0.2)
plot(sdf_b,lw=0.2)
xlabel('t (ms)')
ylabel('f_a, f_b (Hz)')
savefig('sdfs',dpi=600)
show()
