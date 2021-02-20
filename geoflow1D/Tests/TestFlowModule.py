
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

import unittest
from GridLib import *
from FlowLib import *
from FieldsLib import *
from LinearSystemLib import *

class Test_Flow(unittest.TestCase):
	def setUp(self):
		L = 6.
		n = 3
		self.dx = L/(n-1.)
		A = 1.0
		nodesCoord, elemConn = createGridData(L, n)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		self.grid = Grid_1D(gridData)

		self.timeStep = 1.0

		k = 3.0
		self.permeability = ScalarField(self.grid.getNumberOfRegions())
		self.permeability.setValue(self.grid.getRegions()[0], k)

		alpha = 1.0
		self.biot = ScalarField(self.grid.getNumberOfRegions())
		self.biot.setValue(self.grid.getRegions()[0], alpha)

		self.viscosity = 1e-3
		self.density = 1000.
		self.gravity = -10.

		self.D = k*A/self.viscosity
		self.Dx = self.D/self.dx

		self.ls = LinearSystem(self.grid.getNumberOfVertices())

	def test_AssemblyDarcyVelocities(self):
		AssemblyDarcyVelocitiesToMatrix(self.ls, self.grid, self.viscosity, self.permeability, pShift=0)
		self.assertEqual(self.ls.getMatrixValue(0,0), self.Dx)
		self.assertEqual(self.ls.getMatrixValue(0,1), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,0), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,1), 2*self.Dx)
		self.assertEqual(self.ls.getMatrixValue(1,2), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(2,1), -self.Dx)
		self.assertEqual(self.ls.getMatrixValue(2,2), self.Dx)

		AssemblyDarcyVelocitiesToVector(self.ls, self.grid, self.viscosity, self.permeability, self.density, self.gravity, pShift=0)
		self.assertEqual(self.ls.getVectorValue(0), -self.D*self.density*self.gravity)
		self.assertEqual(self.ls.getVectorValue(1), 0.0)
		self.assertEqual(self.ls.getVectorValue(2), self.D*self.density*self.gravity)


	def test_VolumetricStrainToVector(self):
		self.ls.eraseVector()
		u_old = ScalarField(self.grid.getNumberOfVertices())
		[u_old.setValue(v, self.dx + v.getCoordinate()**2) for v in self.grid.getVertices()]
		AssemblyVolumetricStrainToVector(self.ls, self.grid, self.timeStep, self.biot, u_old)
		self.assertEqual(self.ls.getVector()[0], 4.5)
		self.assertEqual(self.ls.getVector()[1], 18.)
		self.assertEqual(self.ls.getVector()[2], 13.5)




if __name__ == '__main__':
	unittest.main()