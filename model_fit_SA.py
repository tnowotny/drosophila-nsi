from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.optimize import dual_annealing
from model_fit_simple import *

if __name__ == '__main__':
    x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
    bounds= [ (0.02, 0.08), (10, 30), (5,9), (0.0005,0.002),
              (-40,-20), (-5,-3), (3e-4,5e-3), (-6.0, -2.0),
              (0.5, 2), (2e-4, 8e-3) ]

    resfile= open("fit_trajectory_SA.txt","w")
    res_sa= dual_annealing(fit_fun,bounds=bounds, args=(trg_sdf,the_conc,net,p,resfile),x0=x0,maxfun=1000000)

    print(res_sa)
    resfile.close()
