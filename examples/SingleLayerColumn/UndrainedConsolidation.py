import sys
sys.path.append("../../geoflow1D")

from GridModule import *
from FieldsModule import *
from FlowModule import *
from GeoModule import *
from LinearSystemModule import *
from CycleControllersModule import *
from ResultsHandlerModule import *
from PhysicalPropertiesModule import *
from SolverModule import *
from UtilsModule import *


def computeUndrainedSolution(grid, folder_settings):
	n = grid.getNumberOfVertices()

	# -------------- NUMERICAL SETTINGS -------------------
	num_set = getJsonData(folder_settings + "numerical_settings.json")
	initialTime = num_set.get("TransientCycle").get("InitialTime")
	timeStep = num_set.get("TransientCycle").get("FinalTime")
	finalTime = 10*timeStep
	timeHandler = TimeHandler( timeStep, finalTime, initialTime )
	# -----------------------------------------------------

	# -------------- PROPERTIES ---------------------------
	props = PhysicalProperties(grid, "settings//")
	# -----------------------------------------------------

	# ---------------- INITIAL FIELDS ---------------------
	ic = getJsonData(folder_settings + "IC.json")
	p_old = ScalarField(n, ic.get("Initial Condition").get("Value_u"))
	u_old = ScalarField(n, ic.get("Initial Condition").get("Value_p"))
	g = ic.get("Gravity")
	# -----------------------------------------------------

	# ------------- CREATE LINEAR SYSTEM ------------------
	nDOF = 2
	ls = LinearSystemCOO(grid.stencil, nDOF)
	ls.initialize()
	pShift = 0
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
	# ls.applyBoundaryConditionsToMatrix(grid, bound_p, pShift)
	# -----------------------------------------------------

	# -------------- GEOMECHANICAL MODEL ------------------
	AssemblyStiffnessMatrix(ls, grid, props, uShift)
	AssemblyPorePressureToMatrix(ls, grid, props, uShift)
	ls.applyBoundaryConditionsToMatrix(grid, bound_u, uShift)
	# -----------------------------------------------------

	# ----------------- DEFINE SOLVER ---------------------
	solver = Solver(tol=1e-14, maxiter=500)
	# -----------------------------------------------------


	# -------------- TRANSIENT SOLUTION -------------------
	timeHandler.advanceTime()
	while timeHandler.isFinalTimeReached():
		ls.eraseVector()

		# --------------- FLUID FLOW MODEL --------------------
		AssemblyDarcyVelocitiesToVector(ls, grid, props, g, pShift)
		AssemblyBiotAccumulationToVector(ls, grid, props, timeStep, p_old, pShift)
		AssemblyVolumetricStrainToVector(ls, grid, props, timeStep, u_old, pShift)
		# ls.applyBoundaryConditionsToVector(grid, bound_p, pShift)
		# -----------------------------------------------------

		# -------------- GEOMECHANICAL MODEL ------------------
		AssemblyGravityToVector(ls, grid, props, g, uShift)
		ls.applyBoundaryConditionsToVector(grid, bound_u, uShift)
		# -----------------------------------------------------

		solver.solve(ls.matrix, ls.rhs)
		ls.solution = solver.solution

		if pShift == 0:	p_new, u_new = ls.splitSolution(n)
		else:			u_new, p_new = ls.splitSolution(n)

		p_old.setField(p_new)
		u_old.setField(u_new)

		timeHandler.advanceTime()

	return p_old, u_old
	'''
	'''
	# -----------------------------------------------------

if __name__ == "__main__":
	# ------------------ GRID DATA ------------------------
	L = 10
	nVertices = 40
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

	folder_settings = "settings//"
	p_old, u_old = computeUndrainedSolution(grid, folder_settings)
