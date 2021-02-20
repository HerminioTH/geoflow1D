

# 		 ---   +-------+
# 		  |    |.......|  |\
#  	 |    |    |.......|  |-\
#  g |	  |    |.......|  |--\
#  	 |	  |    |.......|  |---\    --> p(x) = p(0) + rho*g*x
#  	 V	  |    |.......|  |----\  /
# 		H |    |.......|  |-----\/
# 		  |    |.......|  |------\
# 		  |    |.......|  |-------\
# 		  |    |.......|  |--------\
#  x ^	  |    |.......|  |---------\
# 	 |	  |    |.......|  |----------\
# 	_|_  _|_   |_______|  |___________\



import sys
sys.path.append("../../geoflow1D")
from GridModule import *
from FieldsModule import *
from LinearSystemModule import *
from FlowModule import *
from SolverModule import *
import numpy as np
from matplotlib import pyplot as plt

class FluidProps(object):
	def __init__(self, grid, permeability, density, viscosity):
		self.k = ScalarField(grid.getNumberOfRegions())
		self.k.setValue(grid.getRegions()[0], permeability)
		self.rho_f = density
		self.mu = viscosity


mm = 1000.

# -------------- GRID DATA ----------------------------
H = 10
nVertices = 25
nodesCoord, elemConn = createGridData(H, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
grid.buildStencil()
# -----------------------------------------------------


# -------------- PROPERTIES ----------------------------
k = 1.2e-12		# Permeability
rho = 1000.		# Fluid density
mu = 1e-3		# Fluid viscosity
props = FluidProps(grid, k, rho, mu)
g = -9.81
q = 0.0
# -----------------------------------------------------


# ------------- CREATE LINEAR SYSTEM ------------------
nDOF = 1
ls = LinearSystemCOO(grid.stencil, nDOF)
ls.initialize()
# -----------------------------------------------------


# ------------------ FLOW MODEL -----------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, props, 0)
AssemblyDarcyVelocitiesToVector(ls, grid, props, g, 0)
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
p_bar = 1e5
ls.applyDirichlet(0, p_bar)
# -----------------------------------------------------


# # ----------------- DEFINE SOLVER ---------------------
solver = Solver(tol=1e-14, maxiter=5000)
solver.solve(ls.matrix, ls.rhs)
p_n = solver.solution
x_n = [v.getCoordinate() for v in grid.getVertices()]
# # -----------------------------------------------------


# ------------- ANALYTICAL SOLUTION -------------------
def analyticalSolution(x, rho, g, p_bar):
	return rho*g*np.array(x) + p_bar
x_a = np.linspace(0, H, 100)
p_a = analyticalSolution(x_a, rho, g, p_bar)
# -----------------------------------------------------


# -------------- PLOT SOLUTION ------------------------
plt.plot(p_n, x_n, 'o', label='Numeric')
plt.plot(p_a, x_a, '-', label='Analytic')
plt.legend(loc=0, numpoints=1)
plt.grid(True)
plt.xlabel('pressure')
plt.ylabel('x')
plt.show()
# # -----------------------------------------------------
