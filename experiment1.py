from droso_AB3 import *

"""
v_a= np.arange(-100,100,0.1)*mV
plt.figure()
plt.plot(v_a,pM1/(1+exp((pM2-v_a)/pM3)))
plt.plot(v_a,pM4/(1+exp((pM5-v_a)/pM6)))
plt.show()
exit(1)
"""

# Record spikes
spikeA = SpikeMonitor(sensillum)
spikeB = EventMonitor(sensillum,'b_spike')

# Record a few traces
trace = StateMonitor(sensillum, ('v_a','v_b','v_sl'), record= [0])

# run a little to get away from initial artefacts
run(2.5 * second, report='text')
# Loop over a few trials
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

# Plot intersting stuff

# Plot variablesfrom the StateMonitor
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

# Plot spikes 
figure()
scatter(spikeA.t/ms,spikeA.i)
scatter(spikeB.t/ms,spikeB.i+1)

# Make and plot spike density functions
sdf_a= sdf(spikeA.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 100)
sdf_b= sdf(spikeB.t/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 100)

figure()
plot(sdf_a["t"],sdf_a["sdf"],lw=0.2)
plot(sdf_b["t"],sdf_b["sdf"],lw=0.2)
xlabel('t (ms)')
ylabel('f_a, f_b (Hz)')
savefig('sdfs',dpi=600)
show()
