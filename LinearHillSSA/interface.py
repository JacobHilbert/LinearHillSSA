import numpy as np
import warnings

from scipy.optimize import root
from scipy.integrate import solve_ivp

from dataclasses import dataclass
from functools import partial

from .ssa import linearhillssa as ssa

def model_template(S=2,km=1,γm=1/5,kp=50,γp=1/30,C=None,v=1,α=1,n=1):
    I = np.eye(S//2)
    u = np.ones(S//2)
    if C is None:
        C = 0*I
    parameters = [km,γm,kp,γp,v,α,n]
    for i,p in enumerate(parameters):
        if np.shape(p) == ():
            parameters[i] = u*p
    km,γm,kp,γp,v,α,n = parameters
    return {
        "A" : np.block([
            [I*0,  I*0 ],
            [I*γm, I*0 ],
            [I*kp, I*0 ],
            [I*0,  I*γp]
        ]),
        "b" : np.ravel([km,u*0,u*0,u*0]),
        "T" : np.block([
            [I*0, C  ],
            [I*0, I*0],
            [I*0, I*0],
            [I*0, I*0]
        ]),
        "v" : np.ravel([v,u,u,u]),
        "α" : np.ravel([α,u,u,u]),
        "n" : np.ravel([n,u,u,u]),
        "R" : np.block([
            [+1*I,  0*I],
            [-1*I,  0*I],
            [ 0*I, +1*I],
            [ 0*I, -1*I]
        ])
    }
    
def left_permutation_matrix(N):
    """Gives the matrix that permutes the first column with the last, 
    and every other with its last"""
    return np.roll(np.eye(N),1,axis=0)

@dataclass
class SSA:
    parameters: dict
    
    def __post_init__(self):
        self.S = self.parameters["A"].shape[1]
    
    def __repr__(self):
        return f"<SSA model with {self.S} variables>"
        
    
    
    def path(self,iterations,x0=None):
        if x0 is None:
            x0 = np.ones(self.S)
        x0 = np.copy(x0)
        X = np.zeros((self.S,iterations),order="F")
        times = np.zeros(iterations)
        ssa.ssa_path(x0,X,times,*self.parameters.values())
        return times,X

    def measure_path(self,times,x0=None):
        if x0 is None:
            x0 = np.ones(self.S)
        x0 = np.copy(x0)
        X = np.zeros((self.S,len(times)),order="F")
        ssa.ssa_measure(x0,times,X,*self.parameters.values())
        return X
        
    def measure_ensemble(self,times,N,x0=None):
        if x0 is None:
            x0 = np.ones(self.S)
        if np.ndim(x0) == 1:
            x0 = np.array([np.copy(x0) for i in range(N)]).T
        x0 = np.asfortranarray(x0)
        X = np.zeros((
            self.S,
            len(times),
            N
        ),order="F")
        ssa.ssa_ensemble(x0,times,X,*self.parameters.values())
        return X
                
    def deterministic_solution(self,times,x0=None,**kwargs):
        if x0 is None:
            x0 = np.ones(self.S)
        sol = solve_ivp(
            ssa.jacobian,
            (times.min(),times.max()),
            x0,
            args=self.parameters.values(),
            t_eval=times,
            vectorized=True,
            **kwargs
        )
        if not sol.success:
            warnings.warn(
                "scipy.integrate.solve_ivp failed to converge "\
                f"with message: \n\t'{sol.message}'\n",RuntimeWarning)
        return sol.y
        
    def steady_state(self,x0=None,**kwargs):
        if x0 is None:
            x0 = np.ones(self.S)
        f = partial(self.deterministic_model,0.0)
        sol = root(f,x0,**kwargs)
        if not sol.success:
            warnings.warn(
                "scipy.optimize.root failed to converge "\
                f"with message: '{sol.message}'",RuntimeWarning)
        return sol.x
        
        
    















