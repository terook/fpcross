"""Package fpcross, module equation_demo.equation_duff: duffing problem.

This module contains the class EquationDum, which extends the base class
Equation and provides a description of the 3-dimensional dumbbell model, which
is a case of the Fokker-Planck equation. The formulation is based on the work
Dolgov, S. V., Khoromskij, B. N., and Oseledets, I. V. (2012). Fast Solution of
Parabolic Problems in the Tensor Train/Quantized Tensor Train Format with
Initial Application to the Fokker--Planck Equation. SIAM J. Sci. Comput. 34,
A3016–A3038; doi:10.1137/120864210.

"""
import numpy as np
import teneva


from ..equation import Equation


class EquationDum(Equation):
    """Autonomous form of forced duffing for the FPE (it works only for d=3)."""

    def __init__(self, d=3, e=1.E-5, is_full=False, name='Dum'):
        super().__init__(d, e, is_full, name)

        if is_full:
            raise ValueError('This equation works only for TT-format')

        self.set_grid(n=60, a=-10., b=+10.)
        self.set_grid_time(m=100, t=10.)

    def f(self, X, t):
        X = X.T
		res=X
		res[0]= X[1]
		res[1]=-self.damp*X[1]-self.stif*X[0]-self.beta*X[0]**3+self.favg+self.falt*sin(X[2])
		res[2]= self.omga
        return res.T

    def f1(self, X, t):
        X = X.T
        res=X*0
		res[1]=-self.damp
        return res.T

    def init(self):
        self.with_rs = False
        self.with_rt = False

        self.damp = +0.1
        self.stif = +1.0
		self.beta = -0.3
		self.favg = +0.2
		self.falt = +0.2
		self.omga = +1.3
