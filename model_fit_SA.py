from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.optimize import anneal
from model_fit_simple import *

if __name__ == '__main__':
    x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
    bound_low= [ 0.03,  10,  5, 0.0007, -40, -5, 3e-4, -6.0, 0.5,  2e-4 ]
    bound_high= [ 0.06, 30, 9, 0.0015, -20, -3, 5e-3, -2.0, 2, 8e-3]

    resfile= open("fit_trajectory_SA.txt","w")
    res_sa= anneal(fit_fun,x0, args=(trg_sdf,the_conc,net,p,resfile),lower=bound_low,upper=bound_high,full_output=True,maxeval=1000000)

    print(res_sa)
    resfile.close()
