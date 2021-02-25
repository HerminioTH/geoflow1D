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

from matplotlib import pyplot as plt
from prettytable import PrettyTable


# ------------- AUXILIARY FUNCTIONS -------------------
def computeDenominator(grid, vec1, vec2):
	denominator = computeNormL2(grid, vec1-vec2)
	if denominator == 0.: denominator = 1
	return denominator

def computeRelativeError(grid, vec1, vec2, denominator):
	numerator = computeNormL2(grid, vec1-vec2)
	return numerator/denominator
# -----------------------------------------------------


# ------------------ GRID DATA ------------------------
L = 10
nVertices = 25
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

# ---------------- TIME SETTINGS ----------------------
folder_settings = "settings//"
num_set = getJsonData(folder_settings + "numerical_settings.json")
initialTime = num_set.get("TransientCycle").get("InitialTime")
timeStep = num_set.get("TransientCycle").get("TimeStep")
finalTime = num_set.get("TransientCycle").get("FinalTime")
timeHandler = TimeHandler( timeStep, finalTime, initialTime )
# -----------------------------------------------------

# -------------- ITERATIVE SETTINGS -------------------
maxIte = num_set.get("IterativeCycle").get("MaximumNumberOfIterations")
maxTol = num_set.get("IterativeCycle").get("Tolerance")
iterativeController = IterativeCycleController(maxIte, maxTol)
# -----------------------------------------------------

# -------------- ITERATIVE SETTINGS -------------------
# Choose between Fixed-Strain (FSN) and Fixed-Stress (FSS)
splittingScheme = num_set.get("SplittingScheme")
# -----------------------------------------------------

# ------------------ PROPERTIES -----------------------
props = PhysicalProperties(grid, "settings//")
# -----------------------------------------------------

# -------------- UNDRAINED SOLUTION -------------------
# folder_settings = "settings//"
# p_old, u_old = computeUndrainedSolution(grid, folder_settings)
# ic = getJsonData(folder_settings + "IC.json")
# p_new = ScalarField(grid.getNumberOfVertices())
# u_new = ScalarField(grid.getNumberOfVertices())
# g = ic.get("Gravity")

ic = getJsonData(folder_settings + "IC.json")
p_old = ScalarField(nVertices, ic.get("Initial Condition").get("Value_u"))
u_old = ScalarField(nVertices, ic.get("Initial Condition").get("Value_p"))
p_new = ScalarField(grid.getNumberOfVertices())
u_new = ScalarField(grid.getNumberOfVertices())
g = ic.get("Gravity")
# -----------------------------------------------------


print(u_old.getField())


# ---------------- TUNING PARAMETER -------------------
delta = 2.0
deltas = np.repeat(delta, nVertices)
deltaField = ScalarField(nVertices)
deltaField.setField(deltas)
# -----------------------------------------------------


# ------------- CREATE LINEAR SYSTEM ------------------
ls_mass = LinearSystemCOO(grid.stencil, nDOF=1)
ls_mass.initialize()
ls_geom = LinearSystemCOO(grid.stencil, nDOF=1)
ls_geom.initialize()
# -----------------------------------------------------

# ------------- BOUNDARY CONDITIONS -------------------
bound_p = getJsonData(folder_settings + "BC_p.json")
bound_u = getJsonData(folder_settings + "BC_u.json")
# -----------------------------------------------------

# --------------- FLUID FLOW MODEL --------------------
AssemblyDarcyVelocitiesToMatrix(ls_mass, grid, props)
AssemblyBiotAccumulationToMatrix(ls_mass, grid, timeStep, props)
if splittingScheme == "FSS":
	AssemblyFixedStressAccumulationToMatrix(ls_mass, grid, timeStep, props, deltaField)
ls_mass.applyBoundaryConditionsToMatrix(grid, bound_p)
# -----------------------------------------------------

# -------------- GEOMECHANICAL MODEL ------------------
AssemblyStiffnessMatrix(ls_geom, grid, props)
AssemblyPorePressureToMatrix(ls_geom, grid, props)
ls_geom.applyBoundaryConditionsToMatrix(grid, bound_u)
# -----------------------------------------------------

# ----------------- DEFINE SOLVER ---------------------
solver = Solver(tol=1e-14, maxiter=500)
# -----------------------------------------------------

# --------------- RESULTS HANDLER ---------------------
folderName = "results//%s//"%splittingScheme
res_p = SaveResults(grid, "p.txt", folderName, 'Pressure', 'Pa')
res_u = SaveResults(grid, "u.txt", folderName, 'Displacement', 'm')
res_u.copySettings("settings//", folderName)
res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())
# -----------------------------------------------------

# ----------------- PRETTY TABLE ----------------------
table = PrettyTable(['Time Level', '   Time', '   nItes', '     Error'])
table.align['Time Level'] = 'l'
table.align['   Time'] = 'r'
table.align['   nItes'] = 'r'
table.align['     Error'] = 'r'
table.hrules = 1
print(table)
# -----------------------------------------------------

# -------------- TRANSIENT SOLUTION -------------------

timeLevel = 1
timeHandler.advanceTime()
while timeHandler.isFinalTimeReached():
	iterativeController.reset()
	error_list = []
	while iterativeController.keepCycling():

		# --------------- FLUID FLOW MODEL --------------------
		ls_mass.eraseVector()
		AssemblyBiotAccumulationToVector(ls_mass, grid, props, timeStep, p_old)
		AssemblyDarcyVelocitiesToVector(ls_mass, grid, props, g)
		if splittingScheme == "FSS":
			AssemblyFixedStressAccumulationToVector(ls_mass, grid, props, timeStep, deltaField, p_new)
		AssemblyVolumetricStrainToVector(ls_mass, grid, props, timeStep, u_old)
		AssemblyVolumetricStrainToVector(ls_mass, grid, props, -timeStep, u_new)
		ls_mass.applyBoundaryConditionsToVector(grid, bound_p)
		solver.solve(ls_mass.matrix, ls_mass.rhs)

		p_last = solver.solution
		if iterativeController.iteNumber == 0:
			den = computeDenominator(grid, p_new.getField(), p_last)
		error = computeRelativeError(grid, p_new.getField(), p_last, den)
		error_list.append(error)

		iterativeController.execute(error)

		p_new.setField(solver.solution)
		# -----------------------------------------------------

		# -------------- GEOMECHANICAL MODEL ------------------
		ls_geom.eraseVector()
		AssemblyGravityToVector(ls_geom, grid, props, g)
		AssemblyPorePressureToVector(ls_geom, grid, props, p_new)
		ls_geom.applyBoundaryConditionsToVector(grid, bound_u)
		solver.solve(ls_geom.matrix, ls_geom.rhs)
		u_new.setField(solver.solution)
		# -----------------------------------------------------

	table.add_row([timeLevel, timeHandler.getCurrentTime(), len(error_list), "%.4e"%error_list[-1]])
	print( "\n".join(table.get_string().splitlines()[-2:]) )
	timeLevel += 1

	p_old.setField(p_new.getField())
	u_old.setField(u_new.getField())

	res_p.saveField(timeHandler.getCurrentTime(), p_old.getField())
	res_u.saveField(timeHandler.getCurrentTime(), u_old.getField())

	timeHandler.advanceTime()
# -----------------------------------------------------

res_p.close()
res_u.close()
