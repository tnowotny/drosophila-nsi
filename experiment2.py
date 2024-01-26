from droso_AB3_asym_fit import *
from utils import sdf

"""
v_a= np.arange(-100,100,0.1)*mV
plt.figure()
plt.plot(v_a,pM1/(1+exp((pM2-v_a)/pM3)))
plt.plot(v_a,pM4/(1+exp((pM5-v_a)/pM6)))
plt.show()
exit(1)
"""

# Record a few traces
trace = StateMonitor(sensillum, ('v_a','v1_a','v_b','v_sl'), record= [0],name="traces")
net.add(trace)

# net.run a little to get away from initial artefacts
net.run(2.5 * second, report='text')
# Loop over a few trials
for i in range(20):
    #I_off= I_base + i*0.05*namp
    conc_a_ini= 10**(-(20-i)/4)
    conc_b= 10**-1
    sensillum.set_states({'conc_a': conc_a_ini})
    net.run(1 * second, report='text')
    conc_b = 0.0
    net.run(4.0 * second, report='text')
    conc_a_ini = 0.0
    sensillum.set_states({'conc_a': conc_a_ini})
    net.run(5 * second, report='text')

# Plot intersting stuff

# Plot variablesfrom the StateMonitor
figure()
d= net.get_states()["traces"]
plot(d["t"]/ms, d["v_sl"]/mV,lw=0.2)
xlabel('t (ms)')
ylabel('"LFP" (mV)')
savefig('lfp',dpi=600)

figure()
plot(d["t"]/ms, d["v_a"]/mV,lw=0.2)
plot(d["t"]/ms, d["v_b"]/mV,lw=0.2)
figure()
plot(d["t"]/ms, d["v1_a"]/mV)

# Plot spikes 
figure()
sa = net.get_states()["spikeA"]
scatter(sa["t"]/ms,sa["i"], s=0.1)
sb = net.get_states()["spikeB"]
scatter(sb["t"]/ms,sb["i"]+1, s=0.1)

# Make and plot spike density functions
sdf_a= sdf(sa["t"]/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 100)
sdf_b= sdf(sb["t"]/ms,t0= 0.0,tmax= 2e5,dt= 2,sigma= 100)

figure()
plot(sdf_a["t"],sdf_a["sdf"],lw=0.2)
plot(sdf_b["t"],sdf_b["sdf"],lw=0.2)
xlabel('t (ms)')
ylabel('f_a, f_b (Hz)')
savefig('sdfs',dpi=600)
show()
