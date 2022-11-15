from LinearHillSSA import *
import numpy as np
from time import time


S = 2
N = 500
system = SSA(model_template(S=S))
x0 = np.ones(S)
tmax = 300
dt = 1
ts = np.arange(0,tmax,dt)

print(f"""
Simulation of housekeeping gene expression with {N} cells 
and {tmax} min total time, sampling each {dt} min""")
start = time()
system.measure_ensemble(ts,N,x0)
print(f"{time()-start:.3f} seconds")