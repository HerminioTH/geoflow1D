import numpy as np

class ScalarField(object):
    def __init__(self, size, initialValue=0, name=None, unity=None ):
        self.__name = name
        self.__unity = unity
        self.__listOfValues = initialValue*np.ones(size)
        self.teste = 'oi'

    def setValue(self, entity, value):
        self.__listOfValues[entity.getIndex()] = value

    def setField(self, values):
        self.__listOfValues = values

    def getField(self):
        return self.__listOfValues

    def addValue(self, entity, value):
        self.__listOfValues[entity.getIndex()] += value

    def getValue(self, entity):
        return self.__listOfValues[entity.getIndex()]

    def getName(self):
        return self.__name

    def getUnity(self):
        return self.__unity

    def getValues(self):
        return self.__listOfValues

def updateField(fieldOld, fieldArray, grid):
    for i,vertex in enumerate(grid.getVertices()):
        fieldOld.setValue(vertex, fieldArray[i])



if __name__ == '__main__':
    from GridLib import *

    elemConn = np.array([[0,1],[1,2],[2,3],[3,4]])
    nodesCoord = np.array([0., 2., 4., 6., 8.])

    gridData = GridData()
    gridData.setElementConnectivity( elemConn )
    gridData.setNodeCoordinates( nodesCoord )

    g = Grid_1D( gridData )

    volumes = ScalarField( g.getNumberOfVertices() )
    width = ScalarField( g.getNumberOfElements() )
