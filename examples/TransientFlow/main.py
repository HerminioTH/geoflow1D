
#   p = 12__   ____________________   __ q = -1e-10
#           \ |                    | /
#            \|      p_0 = 10      |/
#             |____________________|
#
#             |------->
#                      X






import sys
sys.path.append("../../geoflow1D")

from GridModule import *
from FieldsModule import *
from FlowModule import *
from LinearSystemModule import *
from CycleControllersModule import *
from ResultsHandlerModule import *
from PhysicalPropertiesModule import *
from SolverModule import *
from UtilsModule import *

# -------------- GRID DATA ----------------------------
L = 10
nVertices = 40
nodesCoord, elemConn = createGridData(L, nVertices)
gridData = GridData()
gridData.setElementConnectivity(elemConn)
gridData.setNodeCoordinates(nodesCoord)
grid = Grid_1D(gridData)
grid.buildStencil()
# -----------------------------------------------------

# -------------- NUMERICAL SETTINGS -------------------
num_set = getJsonData("settings/numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler(timeStep, finalTime, initialTime)
# -----------------------------------------------------

# -------------- PROPERTIES ---------------------------
props = PhysicalProperties(grid, "settings//")
# -----------------------------------------------------

# -------------- INITIAL CONDITION --------------------
p_init = 10.0
p_old = ScalarField(grid.getNumberOfVertices(), p_init)
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
fileName = "p.txt"
folderName = "results//"
results = SaveResults(grid, fileName, folderName, 'Pressure', 'Pa')
results.copySettings("settings//", folderName)
# -----------------------------------------------------

# ------------- CREATE LINEAR SYSTEM ------------------
nDOF = 1
ls = LinearSystemCOO(grid.stencil, nDOF)
ls.initialize()
# -----------------------------------------------------

# -------------- NUMERICAL SOLUTION -------------------
g = 0.0
AssemblyDarcyVelocitiesToMatrix(ls, grid, props)
AssemblyDarcyVelocitiesToVector(ls, grid, props, g)
AssemblyFluidFlowAccumulationToMatrix(ls, grid, timeStep, props)
AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, props, p_old)
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
p_top = 2.0
lastVertex = nVertices-1
# ls.applyDirichlet(lastVertex, p_top)
flux = -1e-10
ls.applyNeumann(lastVertex, flux)

p_bottom = 12.0
ls.applyDirichlet(0, p_bottom)
# -----------------------------------------------------

# ----------------- DEFINE SOLVER ---------------------
solver = Solver(tol=1e-14, maxiter=500)
solver.solve(ls.matrix, ls.rhs)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():

	results.saveField(timeHandler.getCurrentTime(), p_old.getField())

	ls.eraseVector()
	AssemblyDarcyVelocitiesToVector(ls, grid, props, g)
	AssemblyFluidFlowAccumulationToVector(ls, grid, timeStep, props, p_old)
	ls.applyDirichletToVector(0, p_bottom)
	# ls.applyDirichletToVector(lastVertex, p_top)
	ls.applyNeumann(lastVertex, flux)

	solver.solve(ls.matrix, ls.rhs)
	p_old.setField(solver.solution)


	timeHandler.advanceTime()
results.close()
# -----------------------------------------------------
