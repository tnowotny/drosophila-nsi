from utils import *
import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
from model_fit_PSO import *


if __name__ == '__main__':
    resfile= open("fit_trajectory_GA.txt","w")
    x0= [0.04, 20, 7, 0.001, -30, -4, 1e-3, -4.1, 1, 1e-3 ]
    bound_low= [ 0.02,  10,  5, 0.0005, -40, -5, 3e-4, -6.0, 0.5,  2e-4 ]
    bound_high= [ 0.08, 30, 9, 0.002, -20, -3, 5e-3, -2.0, 2, 8e-3]

    from sko.GA import GA
    ga = GA(func=fit_fun, n_dim=10, size_pop=40, max_iter=150, lb=bound_low, ub=bound_high, prob_mut= 0.001)
    ga.run()
    print('best_x is ', ga.best_x, 'best_y is', ga.best_y)
    resfile.close()
