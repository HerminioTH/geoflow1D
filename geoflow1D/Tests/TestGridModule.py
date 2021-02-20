import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

import unittest
from GridLib import *


class Test_Vertex(unittest.TestCase):
	def setUp(self):
		self.v1 = Vertex(11, 0.4)
		self.v1.addToVolume(0.23)
		self.v1.addToVolume(0.31)

	def test_Coordinate(self):
		self.assertEqual(self.v1.getCoordinate(), 0.4)

	def test_AddToVolume(self):
		self.assertEqual(self.v1.getVolume(), 0.54)

	def test_Index(self):
		self.assertEqual(self.v1.getIndex(), 11)

class Test_Face(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.4)
		v2 = Vertex(12, 0.5)
		self.face = Face([v1, v2], 2.3)

	def test_NeighborVertices(self):
		self.assertEqual(self.face.getBackwardVertex().getIndex(), 11)
		self.assertEqual(self.face.getForwardVertex().getIndex(), 12)

	def test_Area(self):
		self.assertEqual(self.face.getArea(), 2.3)

	def test_Coordinate(self):
		self.assertEqual(self.face.getCoordinate(), 0.45)

class Test_Element(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.4)
		v2 = Vertex(12, 0.52)
		self.e = Element([v1, v2], 3, 1, 2.3)

	def test_Index(self):
		self.assertEqual(self.e.getIndex(), 3)

	def test_ParentRegionIndex(self):
		self.assertEqual(self.e.getParentRegionIndex(), 1)

	def test_Area(self):
		self.assertEqual(self.e.getArea(), 2.3)

	def test_Length(self):
		self.assertEqual(self.e.getLength(), 0.12)

	def test_Volume(self):
		self.assertEqual(self.e.getVolume(), 2.3*0.12)

	def test_Vertices(self):
		self.assertEqual(self.e.getVertices()[0].getIndex(), 11)
		self.assertEqual(self.e.getVertices()[1].getIndex(), 12)

	def test_Face(self):
		self.assertEqual(self.e.getFace().getCoordinate(), (0.4 + 0.52)/2)


class Test_Region(unittest.TestCase):
	def setUp(self):
		v1 = Vertex(11, 0.40)
		v2 = Vertex(12, 0.52)
		v3 = Vertex(17, 0.60)
		self.regionIndex = 7
		e1 = Element([v1, v2], 3, self.regionIndex, 2.3)
		e2 = Element([v2, v3], 5, self.regionIndex, 2.3)
		self.region = Region("BODY", self.regionIndex)
		self.region.addElement(e1)
		self.region.addElement(e2)

	def test_Index(self):
		self.assertEqual(self.region.getIndex(), self.regionIndex)

	def test_Name(self):
		self.assertEqual(self.region.getName(), "BODY")

	def test_Elements(self):
		elements = self.region.getElements()
		self.assertEqual(elements[0].getIndex(), 3)
		self.assertEqual(elements[1].getIndex(), 5)
		self.assertEqual(elements[0].getParentRegionIndex(), self.regionIndex)
		self.assertEqual(elements[1].getParentRegionIndex(), self.regionIndex)



class Test_GridData(unittest.TestCase):
	def setUp(self):
		nodesCoord, elemConn = createGridData(10, 11)
		self.gd_1R = GridData()
		self.gd_1R.setElementConnectivity(elemConn)
		self.gd_1R.setNodeCoordinates(nodesCoord)
		self.gd_1R.initialize()

		self.gd_2R = GridData()
		self.gd_2R.setElementConnectivity(elemConn)
		self.gd_2R.setNodeCoordinates(nodesCoord)
		self.gd_2R.setElementsToRegion('lower_layer', [0, 1, 2, 3])
		self.gd_2R.setElementsToRegion('upper_layer', [4, 5, 6, 7, 8, 9])
		self.gd_2R.initialize()

	def test_RegionNames(self):
		self.assertEqual(self.gd_1R.regionNames[0], 'None')
		self.assertEqual(self.gd_2R.regionNames[0], 'lower_layer')
		self.assertEqual(self.gd_2R.regionNames[1], 'upper_layer')

	def test_ElementConnectivity(self):
		self.assertListEqual(self.gd_1R.elemConnectivity, [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]])
		self.assertListEqual(self.gd_2R.elemConnectivity, [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]])

	def test_RegionElements(self):
		self.assertListEqual(self.gd_1R.regionElements[0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
		self.assertListEqual(self.gd_2R.regionElements[0], [0, 1, 2, 3])
		self.assertListEqual(self.gd_2R.regionElements[1], [4, 5, 6, 7, 8, 9])

	def test_NodeCoordinates(self):
		self.assertListEqual(list(self.gd_1R.nodeCoordinates), [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10])
		self.assertListEqual(list(self.gd_2R.nodeCoordinates), [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10])


class Test_Grid(unittest.TestCase):
	def setUp(self):
		self.n, self.H = 11, 15
		nodesCoord, elemConn = createGridData(self.H, self.n)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		gridData.setElementsToRegion('lower_layer', [0, 1, 2, 3])
		gridData.setElementsToRegion('upper_layer', [4, 5, 6, 7, 8, 9])
		gridData.addBoundary("TOP", 9, self.n-1)
		gridData.addBoundary("BOTTOM", 0, 0)
		self.grid = Grid_1D(gridData)

	def test_VerticesOnRegions(self):
		vertices_region_0 = self.grid.getVerticesFromRegion(self.grid.getRegions()[0])
		self.assertEqual(vertices_region_0[0].getIndex(), 0)
		self.assertEqual(vertices_region_0[1].getIndex(), 1)
		self.assertEqual(vertices_region_0[2].getIndex(), 2)
		self.assertEqual(vertices_region_0[3].getIndex(), 3)
		self.assertEqual(vertices_region_0[4].getIndex(), 4)

		vertices_region_1 = self.grid.getVerticesFromRegion(self.grid.getRegions()[1])
		self.assertEqual(vertices_region_1[0].getIndex(), 4)
		self.assertEqual(vertices_region_1[1].getIndex(), 5)
		self.assertEqual(vertices_region_1[2].getIndex(), 6)
		self.assertEqual(vertices_region_1[3].getIndex(), 7)
		self.assertEqual(vertices_region_1[4].getIndex(), 8)
		self.assertEqual(vertices_region_1[5].getIndex(), 9)
		self.assertEqual(vertices_region_1[6].getIndex(), 10)

	def test_RegionNames(self):
		self.assertEqual(self.grid.getRegions()[0].getName(), "lower_layer")
		self.assertEqual(self.grid.getRegions()[1].getName(), "upper_layer")

	def test_Numbers(self):
		self.assertEqual(self.grid.getNumberOfVertices(), self.n)
		self.assertEqual(self.grid.getNumberOfElements(), self.n-1)
		self.assertEqual(self.grid.getNumberOfRegions(), 2)

	def test_Volumes(self):
		vertices = self.grid.getVertices()
		self.assertEqual(vertices[0].getVolume(), 0.75)
		self.assertEqual(vertices[-1].getVolume(), 0.75)
		for v in vertices[1:-1]:
			self.assertEqual(v.getVolume(), 1.5)

	def test_Regions(self):
		regions = self.grid.getRegions()
		elems_r1 = regions[0].getElements()
		elem_r1_indices = [0, 1, 2, 3]
		[self.assertEqual(elems_r1[i].getIndex(), elem_r1_indices[i]) for i in range(4)]
		elems_r2 = regions[1].getElements()
		elem_r2_indices = [4, 5, 6, 7, 8, 9]
		[self.assertEqual(elems_r2[i].getIndex(), elem_r2_indices[i]) for i in range(4)]

	def test_ElementParentRegions(self):
		for i,region in enumerate(self.grid.getRegions()):
			if i == 0:	[self.assertEqual(e.getParentRegionIndex(), 0) for e in region.getElements()]
			else:		[self.assertEqual(e.getParentRegionIndex(), 1) for e in region.getElements()]

	def test_Boundaries(self):
		boundaries = self.grid.getBoundaries()
		self.assertEqual(boundaries[0].getName(), "TOP")
		self.assertEqual(boundaries[0].getVertex().getIndex(), 10)
		self.assertEqual(boundaries[0].getElement().getIndex(), 9)

		self.assertEqual(boundaries[1].getName(), "BOTTOM")
		self.assertEqual(boundaries[1].getVertex().getIndex(), 0)
		self.assertEqual(boundaries[1].getElement().getIndex(), 0)

		self.assertEqual(self.grid.getBoundary("TOP").getName(), "TOP")
		self.assertEqual(self.grid.getBoundary("BOTTOM").getName(), "BOTTOM")





if __name__ == '__main__':
	unittest.main()
