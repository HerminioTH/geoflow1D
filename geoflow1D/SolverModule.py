import numpy as np
import scipy.sparse.linalg as spla

class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
        self.residues = []

    def __call__(self, rk=None):
        if self.niter == 0:
            self.r0 = rk
            # self.r0 = 1.0
        self.niter += 1
        self.residues.append(rk/self.r0)
        if self._disp:
            print('iter %3i\trk = %s' % (self.niter, str(rk)))

class cg_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
        self.residues = []

    def __call__(self, rk=None):
        if self.niter == 0:
            self.r0 = np.linalg.norm(rk)
            # self.r0 = 1.0
        self.niter += 1
        self.residues.append(np.linalg.norm(rk)/self.r0)
        if self._disp:
            print('iter %3i\trk = %f' % (self.niter, self.residues[-1]))

class lgmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0
        self.residues = []

    def __call__(self, rk=None):
        if self.niter == 0:
            self.r0 = rk
            self.r0 = 1.0
        self.niter += 1
        self.residues.append(rk/self.r0)
        if self._disp:
            print('iter %3i' % (self.niter))
            # print('iter %3i\trk = %s' % (self.niter, str(np.linalg.norm(rk))))

class Solver(object):
    def __init__(self, tol=1e-8, maxiter=500, restart=200):
        self.tol = tol
        self.maxiter = maxiter
        self.restart = restart
        self.prec = None
        self.counter = gmres_counter(False)
        # self.counter = cg_counter(True)
        # self.counter = lgmres_counter(True)

    def setPreconditioner(self, prec):
        self.prec = prec

    def solve(self, M, b):
        if self.prec:
            # self.solution, a = spla.cg(M, b, tol=self.tol, maxiter=self.maxiter, callback=self.counter, M=self.prec)
            self.solution, a = spla.gmres(M, b, tol=self.tol, maxiter=self.maxiter, restart=self.restart, callback=self.counter, M=self.prec)
            # self.solution, a = spla.lgmres(M, b, tol=self.tol, maxiter=self.maxiter, callback=self.counter, M=self.prec)
            # print("ite: %i"%self.counter.niter)
        else:
            self.solution, a = spla.gmres(M, b, tol=self.tol, maxiter=self.maxiter, restart=self.restart, callback=self.counter)
