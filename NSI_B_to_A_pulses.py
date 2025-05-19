from droso_AB3 import *
from utils import *

# make experiment folder
expt = 'NSI_B_to_A_pulses'
trial = 'NSI_A4B4'
path = os.path.join('./', expt, trial)

# set concentrations
a_test_conc = 10**-4
b_test_conc = 10**-4
a_2H_conc = 0 # reponse overlap of A to B simulation
b_MH_conc = 0 # reponse overlap of B to A simulation

# Check whether the specified path exists or not
isExist = os.path.exists(path)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(path)
   print("The directory created")
if isExist:
    print("The experiment already exists, rename and restart")
    exit()  

# save config
shutil.copy('droso_AB3.py', path)
shutil.copy((expt + '.py'), path)

# Record spikes
spikeA = SpikeMonitor(sensillum)
spikeB = EventMonitor(sensillum,'b_spike')

# Record a few traces
trace = StateMonitor(sensillum, ('s_r_a','s_r_b','v1_a','v1_b','v_a','v_b','v_sl','M_a','M_b'), record= [0])

# set timing parameters (seconds)
offset = (0,0.025,0.05,0.1,0.15,0.45,0.95,1.65,2.55)
Astim = 5
Bstim = 0.5
clean = 5

x = 2.5+((Astim+clean)*(len(offset)+2))+((Bstim+clean)*2)
print("Experiment length = ", x,'s')

# run a little to get away from initial artefacts
run(2.5 * second, report='text')

# control data
conc_a= a_test_conc
conc_b= b_MH_conc
run(Astim * second, report='text')
conc_a= 0.0
conc_b= 0.0
run(clean * second, report='text')
conc_a= a_test_conc
conc_b= b_MH_conc
run(Bstim * second, report='text')
conc_a= 0.0
conc_b= 0.0
run(clean * second, report='text')
conc_b= b_test_conc
conc_a= a_2H_conc
run(Astim * second, report='text')
conc_b= 0.0
conc_a= 0.0
run(clean * second, report='text')
conc_b= b_test_conc
conc_a= a_2H_conc
run(Bstim * second, report='text')
conc_b= 0.0
conc_a= 0.0
run(clean * second, report='text')

# Loop over a few trials
for k in range(len(offset)):
    i = (offset[k])
    conc_a= a_test_conc
    conc_b= b_MH_conc
    print("conc_a = ", conc_a)
    print("conc_b = ", conc_b)
    run(i * second, report='text')
    conc_b= b_test_conc + b_MH_conc
    conc_a= a_test_conc + a_2H_conc
    print("conc_a_overlap = ", conc_a)
    print("conc_b_overlap = ", conc_b)
    run(Bstim * second, report='text')
    conc_b= b_MH_conc
    conc_a= a_test_conc
    run((Astim-i-Bstim) * second, report='text')
    conc_b= 0.0
    conc_a= 0.0
    run(clean * second, report='text')

# Plot intersting stuff

# Plot variablesfrom the StateMonitor
# figure()
# plot(trace.t/ms, trace[0].v_sl/mV,lw=0.2)
# xlabel('t (ms)')
# ylabel('"LFP" (mV)')
# xlim([0, (x*1000)])
# ylim([-100, 0])
# savefig(os.path.join(expt, trial, 'lfp'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, (trace[0].v_sl/mV)+30,lw=0.2)
# xlabel('t (ms)')
# ylabel('"LFP" (mV)')
# xlim([0, (x*1000)])
# ylim([-40, 0])
# savefig(os.path.join(expt, trial, 'lfp_delta'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, (trace[0].s_r_a/mV),lw=0.2)
# plot(trace.t/ms, (trace[0].s_r_b/mV),lw=0.2)
# xlabel('t (ms)')
# ylabel('"receptor activation"(mV)')
# xlim([0, (x*1000)])
# ylim([0, 510])
# savefig(os.path.join(expt, trial, 's_r_a'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, (trace[0].v1_a/mV),lw=0.2)
# plot(trace.t/ms, (trace[0].v1_b/mV),lw=0.2)
# xlabel('t (ms)')
# ylabel('"v1" (mV)')
# xlim([0, (x*1000)])
# ylim([-85, 0])
# savefig(os.path.join(expt, trial, 'v1'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, g_01*(trace[0].v1_a - trace[0].v_a),lw=0.2)
# plot(trace.t/ms, g_01*(trace[0].v1_b - trace[0].v_b),lw=0.2)
# xlabel('t (ms)')
# ylabel('"g_01*(v1-v)"')
# xlim([0, (x*1000)])
# savefig(os.path.join(expt, trial, 'backprop'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, g_r_a*trace[0].s_r_a*(trace[0].v1_a - ENa - trace[0].v_sl),lw=0.2)
# plot(trace.t/ms, g_r_b*trace[0].s_r_b*(trace[0].v1_b - ENa - trace[0].v_sl),lw=0.2)
# xlabel('t (ms)')
# ylabel('"R currents"')
# xlim([0, (x*1000)])
# savefig(os.path.join(expt, trial, 'R'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, g_M*trace[0].M_a*(trace[0].v1_a - EK),lw=0.2)
# plot(trace.t/ms, g_M_b*trace[0].M_b*(trace[0].v1_b - EK),lw=0.2)
# xlabel('t (ms)')
# ylabel('"M currents"')
# ylim([0, 0.65])
# savefig(os.path.join(expt, trial, 'M'+'_'+expt+'_'+trial),dpi=600)

# figure()
# plot(trace.t/ms, trace[0].M_a,lw=0.2)
# plot(trace.t/ms, trace[0].M_b,lw=0.2)
# xlabel('t (ms)')
# ylabel('"M currents"')
# xlim([0, (x*1000)])
# savefig(os.path.join(expt, trial, 'M_a_M_b'+'_'+expt+'_'+trial),dpi=600)

# # Plot spikes 
# figure()
# scatter(spikeA.t/ms,spikeA.i, s=0.1)
# scatter(spikeB.t/ms,spikeB.i+1, s=0.1)
# savefig(os.path.join(expt, trial, 'spikes'+'_'+expt+'_'+trial),dpi=600)

# Make and plot spike density functions
sdf_a= sdf(spikeA.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 50)
sdf_b= sdf(spikeB.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 50)
sdf_a["t"] = ((sdf_a["t"])/1000)
sdf_b["t"] = ((sdf_b["t"])/1000)

figure()
plot(sdf_a["t"],sdf_a["sdf"],lw=0.2)
plot(sdf_b["t"],sdf_b["sdf"],lw=0.2)
xlabel('t (s)')
ylabel('f_a, f_b (Hz)')
xlim([0, x])
ylim([0, 350])
savefig(os.path.join(expt, trial, 'sdfs'+'_'+expt+'_'+trial),dpi=600)

# save data
spikesA = spikeA.t/ms
spikesB = spikeB.t/ms
np.save(os.path.join(expt, trial, 'spikesA.npy'), spikesA)
np.save(os.path.join(expt, trial, 'spikesB.npy'), spikesB)

np.save(os.path.join(expt, trial, 'LFP.npy'), trace[0].v_sl/mV)
np.save(os.path.join(expt, trial, 'v1_a.npy'), trace[0].v_a/mV)
np.save(os.path.join(expt, trial, 'v1_b.npy'), trace[0].v_b/mV)
np.save(os.path.join(expt, trial, 's_r_a.npy'), trace[0].s_r_a/mV)
np.save(os.path.join(expt, trial, 's_r_b.npy'), trace[0].s_r_b/mV)

np.save(os.path.join(expt, trial, 'Mcur_A.npy'), (g_M*trace[0].M_a*(trace[0].v_a - EK)))
np.save(os.path.join(expt, trial, 'Mcur_B.npy'), (g_M_b*trace[0].M_b*(trace[0].v_b - EK)))
np.save(os.path.join(expt, trial, 'M_a.npy'), trace[0].M_a/mV)
np.save(os.path.join(expt, trial, 'M_b.npy'), trace[0].M_b/mV)

np.save(os.path.join(expt, trial, 'Rcur_A.npy'), (g_r_a*trace[0].s_r_a*(trace[0].v1_a - ENa - trace[0].v_sl)))
np.save(os.path.join(expt, trial, 'Rcur_B.npy'), (g_r_b*trace[0].s_r_b*(trace[0].v1_b - ENa - trace[0].v_sl)))

np.save(os.path.join(expt, trial, 'time.npy'), sdf_a["t"])
np.save(os.path.join(expt, trial, 'sdf_a.npy'), sdf_a["sdf"])
np.save(os.path.join(expt, trial, 'sdf_b.npy'), sdf_b["sdf"])
