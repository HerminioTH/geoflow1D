import geoflow1D
from geoflow1D.GridModule import *
from geoflow1D.FieldsModule import *
from geoflow1D.FlowModule import *
from geoflow1D.GeoModule import *
from geoflow1D.LinearSystemModule import *
from geoflow1D.CycleControllersModule import *
from geoflow1D.ResultsHandlerModule import *
from geoflow1D.PhysicalPropertiesModule import *
from geoflow1D.SolverModule import *
from geoflow1D.UtilsModule import *
from UndrainedConsolidation import *

from prettytable import PrettyTable


# ------------------ GRID DATA ------------------------
L = 10
L_mid = 4.
n_lower = 13
n_upper = 24
nVertices = n_lower + n_upper - 1

regionLowerElems = [e for e in range(n_lower-1)]
regionUpperElems = [e for e in range(n_lower-1, n_lower+n_upper-2)]
elemConn = [[i,i+1] for i in range(nVertices-1)]

nodesCoordLower = np.linspace(0, L_mid, n_lower)
nodesCoordUpper = np.linspace(L_mid, L, n_upper)
nodesCoord = np.concatenate((nodesCoordLower[:-1], nodesCoordUpper))

gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
gridData.addBoundary("TOP", nVertices-2, nVertices-1)
gridData.addBoundary("BOTTOM", 0, 0)
gridData.setElementsToRegion("LOWER_LAYER", regionLowerElems)
gridData.setElementsToRegion("UPPER_LAYER", regionUpperElems)
# -----------------------------------------------------

# --------------------- GRID --------------------------
grid = Grid_1D(gridData)
grid.buildStencil()
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
folder_settings = "settings//"
num_set = getJsonData(folder_settings + "numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler(timeStep, finalTime, initialTime)
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
props = PhysicalProperties(grid, "settings//")
# -----------------------------------------------------


# ---------------- INITIAL FIELDS ---------------------
ic = getJsonData(folder_settings + "IC.json")
p_old, u_old = computeUndrainedSolution(grid, props, folder_settings)
g = ic.get("Gravity")
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
nDOF = 2
ls = LinearSystemCOO(grid.stencil, nDOF)
ls.initialize()
pShift = 1
uShift = (1-pShift)
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
bound_p = getJsonData(folder_settings + "BC_p.json")
bound_u = getJsonData(folder_settings + "BC_u.json")
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, props, pShift)
AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, props, pShift)
AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, props, pShift)
ls.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls, grid, props, uShift)
AssemblyPorePressureToMatrix(ls, grid, props, uShift)
ls.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
# -----------------------------------------------------

# ------------- DEFINE PRECONDITIONER -----------------
M_LU = spla.spilu(ls.matrix, fill_factor=100.0)
preconditioner = lambda b : M_LU.solve(b)
prec = spla.LinearOperator((2*nVertices, 2*nVertices), preconditioner)
# -----------------------------------------------------

# ----------------- DEFINE SOLVER ---------------------
solver = Solver(tol=1e-15, maxiter=5000)
solver.setPreconditioner(prec)
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
folderName = "results//FIM//"
res_p = SaveResults(grid, "p.txt", folderName, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", folderName, 'Displacement', 'm')
res_u.copySettings("settings//", folderName)
# -----------------------------------------------------


# ----------------- PRETTY TABLE ----------------------
table = PrettyTable(['Time Level', '        Time', '     nIte', '        Residue'])
table.align['Time Level'] = 'l'
table.align['        Time'] = 'r'
table.align['     nIte'] = 'r'
table.align['        Residue'] = 'r'
table.hrules = 1
print(table)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
timeHandler.advanceTime()
timeLevel = 0
while timeHandler.isFinalTimeReached():
	ls.eraseVector()

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyDarcyVelocitiesToVector(ls, grid, props, g, pShift)
	AssemblyBiotAccumulationToVector(ls, grid, props, timeStep, p_old, pShift)
	AssemblyVolumetricStrainToVector(ls, grid, props, timeStep, u_old, pShift)
	ls.applyBoundaryConditionsToVector(grid, bound_p, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyGravityToVector(ls, grid, props, g, uShift)
	ls.applyBoundaryConditionsToVector(grid, bound_u, uShift)
	# -----------------------------------------------------

	solver.solve(ls.matrix, ls.rhs)
	ls.solution = solver.solution


	if pShift == 0:	p_new, u_new = ls.splitSolution(nVertices)
	else:			u_new, p_new = ls.splitSolution(nVertices)

	p_old.setField(p_new)
	u_old.setField(u_new)

	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	timeHandler.advanceTime()

	table.add_row([timeLevel, timeHandler.getCurrentTime(), solver.counter.niter, solver.counter.residues[-2]])
	print( "\n".join(table.get_string().splitlines()[-2:]) )
	timeLevel += 1

	solver.counter.niter = 0

	# print(ls.matrix.data)

import matplotlib.pyplot as plt

# plt.spy(ls.matrix.todense())
# plt.show()

res_p.close()
res_u.close()
