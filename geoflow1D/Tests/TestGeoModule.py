
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

import unittest
from GridLib import *
from GeoLib import *
from FieldsLib import *
from LinearSystemLib import *

class Test_Assembly(unittest.TestCase):
	def setUp(self):
		L = 12.
		n = 7
		nodesCoord, elemConn = createGridData(L, n)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		gridData.setElementsToRegion('lower_layer', [0, 1])
		gridData.setElementsToRegion('upper_layer', [2, 3, 4, 5])
		self.grid = Grid_1D(gridData)

		valuesModulus = [1.2, 2.5]
		self.modulus = ScalarField(self.grid.getNumberOfRegions())
		[self.modulus.setValue(region, valuesModulus[region.getIndex()]) for region in self.grid.getRegions()]

		valuesDensity = [17, 13]
		self.density = ScalarField(self.grid.getNumberOfRegions())
		[self.density.setValue(region, valuesDensity[region.getIndex()]) for region in self.grid.getRegions()]

		self.gravity = -9.81

		self.ls = LinearSystem(self.grid.getNumberOfVertices())

	def test_Stiffness(self):
		AssemblyStiffnessMatrix(self.ls, self.grid, self.modulus, 0)
		A = [[-0.6,  0.6,  0.,   0.,   0.,   0.,   0.  ],
			 [ 0.6, -1.2,  0.6,  0.,   0.,   0.,   0.  ],
			 [ 0.,   0.6, -1.85, 1.25, 0.,   0.,   0.  ],
			 [ 0.,   0.,   1.25,-2.5,  1.25, 0.,   0.  ],
			 [ 0.,   0.,   0.,   1.25,-2.5,  1.25, 0.  ],
			 [ 0.,   0.,   0.,   0.,   1.25, -2.5,  1.25],
			 [ 0.,   0.,   0.,   0.,   0.,   1.25, -1.25]]
		for i in range(7):
			for j in range(7):
				self.assertAlmostEqual(self.ls.getMatrix()[i][j], A[i][j], 2)


	def test_Gravity(self):
		AssemblyGravityToVector(self.ls, self.grid, self.density, self.gravity, 0)
		self.assertListEqual(list(self.ls.getVector()), [ 166.77, 333.54, 294.3, 255.06, 255.06, 255.06, 127.53])





if __name__ == '__main__':
	unittest.main()
