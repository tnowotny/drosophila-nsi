import numpy as np

def sdf(sT, t0= None, tmax= None, dt= 1.0, sigma= 50.0):
    if t0 is None:
        t0= np.min(sT)
    if tmax is None:
        tmax= np.max(sT)
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
                if left > 0 and left+len(x) < n:
                    sdf[left:left+len(x)] += x
    t= np.arange(t0, tmax, dt)
    start= int(4*sigma/dt)
    stop= start+len(t)
    the_sdf= {
        "t": t,
        "sdf": sdf[start:stop]
    }
    return the_sdf


"""
poor man's method to make a spike train out of an SDF
"""

def spikes(t, sdf, t0= None, tmax= None):
    mindt= 0.001
    t= np.asarray(t)
    sdf= np.asarray(sdf)

    if t0 is None:
        t0= t[0]
    else:
        if t0 < t[0] or t0 >= t[-1]:
            print("the lower time limit t0 must be within the provided time range of the SDF!")
            exit(1)
    if tmax is None:
        tmax = t[-1]
    else:
        if tmax < t0 or tmax >= t[-1]:
            print("the upper time limit must be greater than the lower time limit and less than the upper end of the provided time range of the SDF")
        exit(1)

    sT= []
    tt= t0
    done= False
    sdf= sdf[:-1] # make by one shorter to fit the logical expression for grabbing the right value
    while not done:
        f= sdf[np.logical_and(t[:-1] <= tt, t[1:] > tt)] 
        if f == 0:
            tt += mindt
        else:
            dt= 1/f 
            sT.append(tt+dt)
            tt= tt+dt
        done= tt >= tmax
    return sT
        
