from droso_AB3 import *

# make experiment folder
expt = 'conc_pulse'
trial = 'test1'
path = os.path.join('./', expt, trial)
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
# shutil.copy('droso_AB3.py', path)
# shutil.copy((expt + '.py'), path)

# Record spikes
spikeA = SpikeMonitor(sensillum)
spikeB = EventMonitor(sensillum,'b_spike')

# Record a few traces
trace = StateMonitor(sensillum, ('s_r_a','s_r_b','v_a','v_b','v_sl','M_a','M_b'), record= [0])

# set stimulation parameters
m_conc = 9  # minimum concentration 10^-
s_conc = 2  # concentration step size
trials = 4 # trials / concentration steps
pulse  = 3  # pulse length (s)
clean  = 5  # time between pulses (s)

x = 5 + ((pulse+clean)*2)*trials

print("Experiment length =", x,'s')

conc_a_data =[]
conc_b_data =[]

# run a little to get away from initial artefacts
run(2.5 * second, report='text')
# Loop over a few trials
for i in range(trials):  
    exp_a = (m_conc-(i*s_conc))
    exp_a = round(exp_a,2)
    exp_b = (m_conc-(i*s_conc))
    exp_b = round(exp_b,2)
    conc_a_data.append(exp_a)
    conc_b_data.append(exp_b)
    print("A conc =",exp_a)
    print("B conc =",exp_b)

    conc_a= 10**(-exp_a)
    run(pulse * second, report='text')
    conc_a= 0.0
    run(clean * second, report='text')
    conc_b= 10**(-exp_b)
    run(pulse * second, report='text')
    conc_b= 0.0
    run(clean * second, report='text')

print("A conc =",conc_a_data)
print("B conc =",conc_b_data)

# Plot intersting stuff

# Plot variablesfrom the StateMonitor
figure()
plot(trace.t/ms, trace[0].v_sl/mV,lw=0.2)
xlabel('t (ms)')
ylabel('"LFP" (mV)')
xlim([0, (x*1000)])
ylim([-100, 0])
savefig(os.path.join(expt, trial, 'lfp'+'_'+expt+'_'+trial),dpi=600)

figure()
plot(trace.t/ms, (trace[0].v_sl/mV)+30,lw=0.2)
xlabel('t (ms)')
ylabel('"LFP" (mV)')
xlim([0, (x*1000)])
ylim([-40, 0])
savefig(os.path.join(expt, trial, 'lfp_delta'+'_'+expt+'_'+trial),dpi=600)

figure()
plot(trace.t/ms, g_r*trace[0].s_r_a*(trace[0].v_a - ENa - trace[0].v_sl),lw=0.2)
plot(trace.t/ms, g_r_b*trace[0].s_r_b*(trace[0].v_b - ENa - trace[0].v_sl),lw=0.2)
xlabel('t (ms)')
ylabel('"R currents"')
xlim([0, (x*1000)])
##ylim([0, 1])
savefig(os.path.join(expt, trial, 'R'+'_'+expt+'_'+trial),dpi=600)

figure()
plot(trace.t/ms, trace[0].M_a,lw=0.2)
plot(trace.t/ms, trace[0].M_b,lw=0.2)
xlabel('t (ms)')
ylabel('"M currents"')
xlim([0, (x*1000)])
##ylim([0, 0.65])
savefig(os.path.join(expt, trial, 'M'+'_'+expt+'_'+trial),dpi=600)

# Plot spikes 
figure()
scatter(spikeA.t/ms,spikeA.i, s=0.1)
scatter(spikeB.t/ms,spikeB.i+1, s=0.1)
savefig(os.path.join(expt, trial, 'spikes'+'_'+expt+'_'+trial),dpi=600)

# Make and plot spike density functions
sdf_a= sdf(spikeA.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 50)
sdf_b= sdf(spikeB.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 50)
sdf_a["t"] = ((sdf_a["t"])/1000)
sdf_b["t"] = ((sdf_b["t"])/1000)
spikeA = spikeA.t/ms
spikeB = spikeB.t/ms

figure()
plot(sdf_a["t"],sdf_a["sdf"],lw=0.2)
plot(sdf_b["t"],sdf_b["sdf"],lw=0.2)
xlabel('t (s)')
ylabel('f_a, f_b (Hz)')
xlim([0, x])
#ylim([0, 300])
title(conc_a_data)
savefig(os.path.join(expt, trial, 'sdfs'+'_'+expt+'_'+trial),dpi=600)

# save data
np.save(os.path.join(expt, trial, 'time.npy'), sdf_a["t"])
np.save(os.path.join(expt, trial, 'sdf_a.npy'), sdf_a["sdf"])
np.save(os.path.join(expt, trial, 'sdf_b.npy'), sdf_b["sdf"])
np.save(os.path.join(expt, trial, 'spikeA.npy'), spikeA)
np.save(os.path.join(expt, trial, 'spikeB.npy'), spikeB)
np.save(os.path.join(expt, trial, 'LFP.npy'), trace[0].v_sl/mV)
np.save(os.path.join(expt, trial, 'M_a.npy'), trace[0].M_a/mV)
np.save(os.path.join(expt, trial, 'M_b.npy'), trace[0].M_b/mV)
np.save(os.path.join(expt, trial, 'R_a.npy'), g_r*trace[0].s_r_a*(trace[0].v_a - ENa - trace[0].v_sl)/mV)
np.save(os.path.join(expt, trial, 'R_b.npy'), g_r*trace[0].s_r_b*(trace[0].v_b - ENa - trace[0].v_sl)/mV)