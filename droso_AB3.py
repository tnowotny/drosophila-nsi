from brian2 import *
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

# Basic neuron Parameters
x=[4.22825245e-02, 2.15402562e+01, 7.98889005e+00, 1.16413590e-03, -3.15916662e+01, -4.51528381e+00, 8.44657433e-04, -3.34225086e+00, 8.77785161e-01, 9.75556066e-04]

Cm = 1*nfarad
Cm1 = 1*nfarad
gl = 0.00018685*msiemens
El = -63.563*mV
EK = -95*mV
ENa = 50*mV
g_na = 0.05*msiemens
g_kd = 0.01*msiemens
g_01 = 0.001*msiemens
tau_sr = 25*ms
dt = 0.05*msecond

I_noise= 1*namp*dt**0.5
I_base= 0.8*namp
I_off= I_base

VT = -65*mV

# Adaptation M-current parameters
# activation
pM1 = x[0]
pM2 = x[1]*mV
pM3 = x[2]*mV
# inactivation
pM4 = x[3]
pM5 = x[4]*mV
pM6 = x[5]*mV
g_M= x[6]*msiemens

# receptor current parameters
cmid_a = x[7] # log10 concentration of mid-activation
cmid_b = cmid_a * 1.2
cslope_a = x[8]  # log10 concentration slope of activation
cslope_b = cslope_a / 1
g_r_a = x[9]*msiemens

# scale capacitance (C), conductances (g), I_noise, I_base: by surface area ratios from Zhang et al 2019)
SA_ratio =  0.75  # ab3 ratio 0.73: B SA (91+16) / A SA (98+46)

Cm_b = Cm * SA_ratio
Cm1_b = Cm1 * SA_ratio
gl_b = gl * SA_ratio
g_na_b = g_na * SA_ratio
g_kd_b = g_kd * SA_ratio
g_01_b = g_01 * SA_ratio
I_noise_b = I_noise * SA_ratio
I_base_b = I_base * 1.5 # sym 2 / 0.75 assym 1.5 / 0.5 assym 1
I_off_b = I_base_b
g_M_b = g_M
g_r_b = g_r_a * SA_ratio

# odour parameters
conc_a= 0
conc_b= 0

# sensillum parameters
gl_sl = 0.001*msiemens
El_sl = -30*mV
Cm_sl = 8*nfarad

# The model
eqs = Equations('''
dv_a/dt = (gl*(El-v_a)-
         g_na*(m_a*m_a*m_a)*h_a*(v_a-ENa)-
         g_kd*(n_a*n_a*n_a*n_a)*(v_a-EK)-
         g_M*M_a*(v_a-EK)+g_01*(v1_a-v_a)+
         I_off+I_noise*xi_a)/Cm : volt
dv1_a/dt = (g_01*(v_a-v1_a)-g_r_a*s_r_a*(v1_a-ENa-v_sl))/Cm1 : volt
dm_a/dt = alpha_m_a*(1-m_a)-beta_m_a*m_a : 1
dn_a/dt = alpha_n_a*(1-n_a)-beta_n_a*n_a : 1
dh_a/dt = alpha_h_a*(1-h_a)-beta_h_a*h_a : 1
dM_a/dt = alpha_M_a*(1-M_a)-beta_M_a*M_a : 1
ds0_r_a/dt = (s_r_inf_a - s0_r_a)/tau_sr : 1
ds_r_a/dt = (s0_r_a - s_r_a)/tau_sr : 1

alpha_m_a = 0.32*(mV**-1)*4*mV/exprel((13*mV-v_a+VT)/(4*mV))/ms : Hz
beta_m_a = 0.28*(mV**-1)*5*mV/exprel((v_a-VT-40*mV)/(5*mV))/ms : Hz
alpha_h_a = 0.128*exp((17*mV-v_a+VT)/(18*mV))/ms : Hz
beta_h_a = 4./(1+exp((40*mV-v_a+VT)/(5*mV)))/ms : Hz
alpha_n_a = 0.032*(mV**-1)*5*mV/exprel((15*mV-v_a+VT)/(5*mV))/ms : Hz
beta_n_a = .5*exp((10*mV-v_a+VT)/(40*mV))/ms : Hz
alpha_M_a = pM1/(1+exp((pM2-v_a)/pM3))/ms : Hz
beta_M_a = pM4/(1+exp((pM5-v_a)/pM6))/ms : Hz
s_r_inf_a = 1.0/(1.0+exp((cmid_a-log10(conc_a))*cslope_a)) : 1

dv_b/dt = (gl_b*(El-v_b)-
         g_na_b*(m_b*m_b*m_b)*h_b*(v_b-ENa)-
         g_kd_b*(n_b*n_b*n_b*n_b)*(v_b-EK)-
         g_M_b*M_b*(v_b-EK)+g_01_b*(v1_b-v_b)+
         I_off_b+I_noise_b*xi_b)/Cm_b : volt
dv1_b/dt = (g_01_b*(v_b-v1_b)-g_r_b*s_r_b*(v1_b-ENa-v_sl))/Cm1_b : volt
dm_b/dt = alpha_m_b*(1-m_b)-beta_m_b*m_b : 1
dn_b/dt = alpha_n_b*(1-n_b)-beta_n_b*n_b : 1
dh_b/dt = alpha_h_b*(1-h_b)-beta_h_b*h_b : 1
dM_b/dt = alpha_M_b*(1-M_b)-beta_M_b*M_b : 1
ds0_r_b/dt = (s_r_inf_b - s0_r_b)/tau_sr : 1
ds_r_b/dt = (s0_r_b - s_r_b)/tau_sr : 1
drefr_b/dt = -1 : second

alpha_m_b = 0.32*(mV**-1)*4*mV/exprel((13*mV-v_b+VT)/(4*mV))/ms : Hz
beta_m_b = 0.28*(mV**-1)*5*mV/exprel((v_b-VT-40*mV)/(5*mV))/ms : Hz
alpha_h_b = 0.128*exp((17*mV-v_b+VT)/(18*mV))/ms : Hz
beta_h_b = 4./(1+exp((40*mV-v_b+VT)/(5*mV)))/ms : Hz
alpha_n_b = 0.032*(mV**-1)*5*mV/exprel((15*mV-v_b+VT)/(5*mV))/ms : Hz
beta_n_b = .5*exp((10*mV-v_b+VT)/(40*mV))/ms : Hz
alpha_M_b = pM1/(1+exp((pM2-v_b)/pM3))/ms : Hz
beta_M_b = pM4/(1+exp((pM5-v_b)/pM6))/ms : Hz
s_r_inf_b = 1.0/(1.0+exp((cmid_b-log10(conc_b))*cslope_b)) : 1

dv_sl/dt = (gl_sl*(El_sl-v_sl)+g_r_a*s_r_a*(v1_a-ENa-v_sl)+g_r_b*s_r_b*(v1_b-ENa-v_sl))/Cm_sl : volt
''')

sensillum = NeuronGroup(1, model=eqs, threshold='v_a>-20*mV', events={'b_spike': '(v_b>-20*mV) and (refr_b <= 0*ms)'},
                        refractory=3*ms, method='milstein',dt=dt)

sensillum.run_on_event('b_spike','refr_b= 3*ms')

# Initialization
sensillum.v_a = El + (randn() * 5 - 5)*mV
sensillum.v1_a = El
sensillum.m_a= 0.0
sensillum.h_a= 1.0
sensillum.n_a= 0.0
sensillum.v_b = El + (randn() * 5 - 5)*mV
sensillum.v1_b = El
sensillum.m_b= 0.0
sensillum.h_b= 1.0
sensillum.n_b= 0.0
sensillum.refr_b= 0*ms
