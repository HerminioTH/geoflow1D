
# 		     | sigma
# 		     |
# 		 +---V---+  ---
# 		 |       |   |
# 		 |       |   |
# 		 |       |   |
# 		 |       |   |
# 		 |       |   |
# 		 |       |   |  H
# 		 |		 |   |
# 		 |       |   |
# 		 |       |   |
#  x ^	 |       |   |
# 	 |	 |       |   |
# 	_|_  |_______|  _|_



import sys
sys.path.append("../../geoflow1D")
from GridModule import *
from FieldsModule import *
from LinearSystemModule import *
from GeoModule import *
from SolverModule import *
import numpy as np
from matplotlib import pyplot as plt

class SolidProps(object):
	def __init__(self, grid, M, rho):
		self.M = ScalarField(grid.getNumberOfRegions())
		self.M.setValue(grid.getRegions()[0], M)
		self.rho = ScalarField(grid.getNumberOfRegions())
		self.rho.setValue(grid.getRegions()[0], rho)


mm = 1000.

# -------------- GRID DATA ----------------------------
H = 10
nVertices = 15
nodesCoord, elemConn = createGridData(H, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
grid.buildStencil()
# -----------------------------------------------------

# -------------- PROPERTIES ----------------------------

M = 1.3e8		# Constrained modulus
rho = 2300.		# Solid density
props = SolidProps(grid, M, rho)
g = -9.81
# -----------------------------------------------------


# ------------- CREATE LINEAR SYSTEM ------------------
nDOF = 1
ls = LinearSystemCOO(grid.stencil, nDOF)
ls.initialize()
# -----------------------------------------------------

# -------------- NUMERICAL SOLUTION -------------------
AssemblyStiffnessMatrix(ls, grid, props, 0)
AssemblyGravityToVector(ls, grid, props, g, 0)
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
ls.applyDirichlet(0, 0)
sigma = -5e4
ls.applyNeumann(-1, sigma)
# -----------------------------------------------------

# ----------------- DEFINE SOLVER ---------------------
solver = Solver(tol=1e-8, maxiter=500)
solver.solve(ls.matrix, ls.rhs)
# -----------------------------------------------------

# ------------- ANALYTICAL SOLUTION -------------------
def analyticalSolution(M, stress, L, x, gravity, rho):
	x = np.array(x)
	return x*(-stress + rho*g*L)/M - rho*g*x*x/(2*M)
x_a = np.linspace(0, H, 100)
u_a = analyticalSolution(M, sigma, H, x_a, g, rho)
# -----------------------------------------------------

# -------------- PLOT SOLUTION ------------------------
x_n = [v.getCoordinate() for v in grid.getVertices()]
u_n = solver.solution
plt.plot(u_n*mm, x_n, 'o', label='Numeric')
plt.plot(u_a*mm, x_a, '-', label='Analytic')
plt.grid(True)
plt.xlabel('Displacement (mm)')
plt.ylabel('Coordinate X (m)')
plt.show()
# -----------------------------------------------------
