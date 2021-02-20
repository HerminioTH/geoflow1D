import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))

from ResultsHandlerLib import *
from GridLib import *
import unittest

class Test_Results(unittest.TestCase):
	def setUp(self):
		L = 12.
		n = 4
		nodesCoord, elemConn = createGridData(L, n)
		gridData = GridData()
		gridData.setElementConnectivity(elemConn)
		gridData.setNodeCoordinates(nodesCoord)
		self.grid = Grid_1D(gridData)

		self.source = "settings//"
		self.directory = "FOLDER_ONE//FOLDER_TWO//"
		self.fileName = "Field.txt"
		self.fieldName = "luminosity"
		self.unitName = "cd"

	def test_SaveResults(self):
		res = SaveResults(self.grid, self.fileName, self.directory, self.fieldName, self.unitName)
		res.saveField(0.12, [50.2, 34.8, 14.7, 2.3])
		res.saveField(0.24, [50.2, 32.8, 13.7, 2.1])
		res.copySettings(self.source, self.directory)
		res.close()

		f = open(self.directory + self.fileName)
		lines = f.readlines()
		self.assertEqual(lines[0], "Coordinates:\n")
		self.assertEqual(lines[1], '0.0,4.0,8.0,12.0,\n')
		self.assertEqual(lines[2], 'Volumes:\n')
		self.assertEqual(lines[3], '2.0,4.0,4.0,2.0,\n')
		self.assertEqual(lines[5], ' Time(s), %s (%s) on vertices\n'%(self.fieldName, self.unitName))
		self.assertEqual(lines[6], '0.12,50.2,34.8,14.7,2.3,\n')
		self.assertEqual(lines[7], '0.24,50.2,32.8,13.7,2.1,\n')
		f.close()

	# def test_CopySettings(self):
	# 	res = SaveResults(self.grid, self.fileName, self.directory, self.fieldName, self.unitName)
	# 	# res.copySettings(self.source, self.directory)
	# 	res.close()

	def test_ReadResults(self):
		res = ReadResults(self.directory + self.fileName)
		self.assertListEqual(res.coord, [0.0, 4.0, 8.0, 12.0])
		self.assertListEqual(res.volumes, [2.0, 4.0, 4.0, 2.0])
		self.assertListEqual(res.times, [0.12, 0.24])
		self.assertListEqual(res.getSolutionAtTime(0.12), [50.2, 34.8, 14.7, 2.3])
		self.assertListEqual(res.getSolutionAtTime(0.24), [50.2, 32.8, 13.7, 2.1])

		with self.assertRaises(Exception) as context:
			res.getSolutionAtTime(0.27)
		self.assertTrue('Time 0.270000 does not exist.', context.exception)



if __name__ == '__main__':
	unittest.main()
