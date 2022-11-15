from LinearHillSSA import *
import numpy as np
import matplotlib.pyplot as plt
import re

labels = re.sub("(\w)(\d)",r"$\1_\2$","m1 m2 m3 p1 p2 p3").split()
S = 6
N = 1
C = left_permutation_matrix(S//2)
system = SSA(model_template(
    S = S,
    km = 5e-3,
    γm = np.log(2)/120,
    kp = np.log(2)/6,
    γp = np.log(2)/60,
    v = 0.5,
    α = 40,
    n = -2,
    C = C,
))
x0 = np.ones(S)
x0[3] = 20
ts = np.arange(0,10_000,10)

fig,ax = plt.subplots(2,sharex=True,figsize=(8,6))

plt.sca(ax[0])
u = system.deterministic_solution(ts,x0,method="RK45")
p = plt.plot(ts,u.T)
plt.legend(p,labels)
plt.ylabel("deterministic solution")


plt.sca(ax[1])
#plt.gca().set_prop_cycle(None)
u = system.measure_path(ts,x0)
p = plt.step(ts,u.T)
plt.legend(p,labels)
plt.ylabel("stochastic jump simulation")
plt.xlim(ts.min(),ts.max())

plt.tight_layout()
plt.savefig("repressilator.svg")
plt.show()

