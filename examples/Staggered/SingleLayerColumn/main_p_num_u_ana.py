import geoflow1D
from geoflow1D.GridModule import *
from geoflow1D.FieldsModule import *
from geoflow1D.CycleControllersModule import *
from geoflow1D.ResultsHandlerModule import *
from geoflow1D.PhysicalPropertiesModule import *
from geoflow1D.SolverModule import *
from geoflow1D.UtilsModule import *

from StaggeredFlowModule import *
from StaggeredGeoModule import *
from StaggeredLinearSystemModule import *
from StaggeredResultsHandler import StaggeredSaveResults

from prettytable import PrettyTable
import matplotlib.pyplot as plt

import solgeom as sg
# from Terzaghi import Solution

# ------------------ GRID DATA ------------------------
L = 10
nVertices = 24
nElements = nVertices - 1
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
gridData.addBoundary("TOP", nVertices-2, nVertices-1)
gridData.addBoundary("BOTTOM", 0, 0)
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
p_old = ScalarField(nElements, ic.get("Initial Condition").get("Value_u"))
# for elem in grid.getElements():
# 	p_old.setValue(elem, 1e5)
u_old = ScalarField(nVertices, ic.get("Initial Condition").get("Value_p"))
g = ic.get("Gravity")
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
nDOF = 2
ls = LinearSystemDense(grid.stencil, nDOF)
ls.initialize()
pShift = 1
uShift = (1-pShift)
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
bound_p = getJsonData(folder_settings + "BC_p.json")
bound_u = getJsonData(folder_settings + "BC_u.json")
# -----------------------------------------------------

# -------------- ANALYTICAL SOLUTION ------------------
rock = json.load(open("settings//solid.json", "r"))
fluid = json.load(open("settings//fluid.json", "r"))
load = bound_u.get("TOP").get("Value")
terza = sg.Terzaghi.Solution(L, load, rock, fluid , 0.0)

# for vertex in grid.getVertices():
# 	u_old.setValue(vertex, terza.getDisplacementValuesAtPosition(vertex.getCoordinate(), [0])[0])
# for elem in grid.getElements():
# 	p_old.setValue(elem, terza.getPressureValuesAtPosition(elem.getCentroid(), [0])[0])
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls, grid, props, pShift)
AssemblyBiotAccumulationToMatrix(ls, grid, timeStep, props, pShift)
AssemblyVolumetricStrainToMatrix(ls, grid, timeStep, props, pShift)
# # -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
# AssemblyStiffnessMatrix(ls, grid, props, uShift)
# AssemblyPorePressureToMatrix(ls, grid, props, uShift)
# ls.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)

for vertex in grid.getVertices():
	ls.applyDirichletToMatrix(vertex.getIndex() + uShift*nVertices, 0.0)
# -----------------------------------------------------

# print(ls.matrix)
# plt.spy(ls.matrix, markersize=5)
# plt.show()

# # ------------- DEFINE PRECONDITIONER -----------------
# M_LU = spla.spilu(ls.matrix, fill_factor=100.0)
# preconditioner = lambda b : M_LU.solve(b)
# prec = spla.LinearOperator((2*nVertices, 2*nVertices), preconditioner)
# # -----------------------------------------------------

# # ----------------- DEFINE SOLVER ---------------------
# solver = Solver(tol=1e-12, maxiter=5000)
# # solver.setPreconditioner(prec)
# # -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
folderName = "results//p_num_u_ana//"
res_p = StaggeredSaveResults(grid, "p.txt", folderName, 'Pressure', 'Pa')
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
A_inv = np.linalg.inv(ls.matrix)
timeHandler.advanceTime()
timeLevel = 0
while timeHandler.isFinalTimeReached():
	ls.eraseVector()

	# --------------- FLUID FLOW MODEL --------------------
	AssemblyBiotAccumulationToVector(ls, grid, props, timeStep, p_old, pShift)
	AssemblyVolumetricStrainToVector(ls, grid, props, timeStep, u_old, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	# AssemblyGravityToVector(ls, grid, props, g, uShift)
	# ls.applyBoundaryConditionsToVector(grid, bound_u, uShift)
	for vertex in grid.getVertices():
		value = terza.getDisplacementValuesAtPosition(vertex.getCoordinate(), [timeHandler.getCurrentTime()])[0]
		# value = terza.getDisplacementValuesAtPosition(vertex.getCoordinate(), [0])[0]
		index = vertex.getIndex() + uShift*nVertices
		ls.applyDirichletToVector(index, value)
	# -----------------------------------------------------

	# solver.solve(ls.matrix, ls.rhs)
	# ls.solution = solver.solution

	ls.solution = A_inv.dot(ls.rhs)

	if pShift == 0:	p_new, u_new = ls.splitSolution(nVertices)
	else:			u_new, p_new = ls.splitSolution(nVertices)

	p_old.setField(p_new)
	u_old.setField(u_new)

	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	timeHandler.advanceTime()

	# table.add_row([timeLevel, timeHandler.getCurrentTime(), solver.counter.niter, solver.counter.residues[-2]])
	table.add_row([timeLevel, timeHandler.getCurrentTime(), 0, 0])
	print( "\n".join(table.get_string().splitlines()[-2:]) )
	timeLevel += 1

	# solver.counter.niter = 0

res_p.close()
res_u.close()
