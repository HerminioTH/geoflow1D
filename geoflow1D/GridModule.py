import numpy as np

class GridData( object ):
    def __init__( self ):
        self.elemConnectivity = []
        self.nodeCoordinates = []
        self.regionElements = []
        self.regionNames = []
        self.boundaries = []

    def setElementConnectivity(self, elemConnectivity):
        self.elemConnectivity = elemConnectivity

    def setNodeCoordinates(self, nodeCoordinates):
        self.nodeCoordinates = nodeCoordinates

    def setElementsToRegion(self, regionName, regionElements=None):
        if regionElements == None:
            self.regionElements.append(list(range(len(self.elemConnectivity))))
            self.regionNames.append(regionName)
        else:
            self.regionElements.append(regionElements)
            self.regionNames.append(regionName)

    def addBoundary(self, boundaryName, elemIndex, vertexIndex):
        bound = {"NAME": boundaryName, "ElementIndex": elemIndex, "VertexIndex": vertexIndex}
        self.boundaries.append(bound)

    def initialize(self):
        if len(self.regionElements) == 0:
            self.regionElements.append(list(range(len(self.elemConnectivity))))
            self.regionNames = ['None']


class Vertex( object ):
    def __init__(self, index, coord):
        self.__globalIndex = index
        self.__x = coord
        self.__volume = 0.0

    def getIndex(self):
        return self.__globalIndex

    def getCoordinate(self):
        return self.__x

    def addToVolume(self, value):
        self.__volume += value

    def getVolume(self):
        return self.__volume



class Face(object):
    def __init__(self, vertices, height=1.0):
        self.__bVertex = vertices[0]
        self.__fVertex = vertices[1]
        self.__faceCoord = ( vertices[0].getCoordinate() + vertices[1].getCoordinate() ) / 2.0
        self.__area = height

    def getBackwardVertex(self):
        return self.__bVertex

    def getForwardVertex(self):
        return self.__fVertex

    def getCoordinate(self):
        return self.__faceCoord

    def getArea(self):
        return self.__area


class Element(object):
    def __init__(self, vertices, index, parentRegionIndex, area=1.0):
        self.__globalIndex = index
        self.__parentRegionIndex = parentRegionIndex
        self.__area = area
        self.__vertices = vertices
        self.__buildFace()
        self.__elementLength = abs(vertices[1].getCoordinate() - vertices[0].getCoordinate())

    def __buildFace(self):
        self.__face = Face(self.__vertices)

    def getIndex(self):
        return self.__globalIndex

    def getParentRegionIndex(self):
        return self.__parentRegionIndex

    def getCentroid(self):
        x0 = self.__vertices[0].getCoordinate()
        x1 = self.__vertices[1].getCoordinate()
        return (x0 + x1)/2

    def getVertices(self):
        return self.__vertices

    def getFace(self):
        return self.__face

    def getLength(self):
        return self.__elementLength

    def getArea(self):
        return self.__area

    def getVolume(self):
        return self.__area*self.__elementLength

    def getSubVolume(self):
        return self.__area*self.__elementLength/2.


class Region(object):
    def __init__(self, name, index):
        self.__globalIndex = index
        self.__name = name
        self.__elements = []

    def addElement(self, element):
        self.__elements.append(element)

    def getIndex(self):
        return self.__globalIndex

    def getName(self):
        return self.__name

    def getElements(self):
        return self.__elements

class Boundary(object):
    def __init__(self, name, element, vertex):
        self.__name = name
        self.__element = element
        self.__vertex = vertex

    def getName(self):
        return self.__name

    def getElement(self):
        return self.__element

    def getVertex(self):
        return self.__vertex


class Grid_1D( object ):
    def __init__(self, gridData):
        gridData.initialize()
        self.__buildVertices(gridData)
        self.__buildRegions(gridData)
        self.__buildBoundaries(gridData)
        self.__computeVolumes()

    def __buildVertices(self, gridData):
        self.__vertices = []
        iVertices = []
        hash_table = np.zeros(len(gridData.nodeCoordinates))
        for iElem in gridData.elemConnectivity:
            for iVertex in iElem:
                if hash_table[iVertex] == 0:
                    iVertices.append( iVertex )
                    self.__vertices.append( Vertex( iVertex, gridData.nodeCoordinates[iVertex] ) )
                    hash_table[iVertex] = 1
        self.__numberOfVertices = len(self.__vertices)

    def __buildRegions(self, gridData):
        self.__regions = []
        self.__elements = []
        self.__numberOfRegions = len(gridData.regionNames)
        self.__numberOfElements = len(gridData.elemConnectivity)

        elementIndex = 0
        for regionIndex in range(self.__numberOfRegions):
            regionName = gridData.regionNames[regionIndex]
            region = Region(regionName, regionIndex)
            for iElem in gridData.regionElements[regionIndex]:
                v1 = gridData.elemConnectivity[iElem][0]
                v2 = gridData.elemConnectivity[iElem][1]
                e = Element( [ self.__vertices[v1], self.__vertices[v2] ], elementIndex, regionIndex )
                region.addElement(e)
                self.__elements.append(e)
                elementIndex += 1
            self.__regions.append(region)

    def __buildBoundaries(self, gridData):
        self.__boundaries = []
        for bound in gridData.boundaries:
            for e in self.__elements:
                if e.getIndex() == bound.get("ElementIndex"):
                    for v in e.getVertices():
                        if v.getIndex() == bound.get("VertexIndex"):
                            name = bound.get("NAME")
                            self.__boundaries.append(Boundary(name, e, v))


    def buildStencil(self):
        nVertices = len(self.__vertices)
        self.stencil = [[] for i in range(nVertices)]
        for element in self.__elements:
            localHandle = 0
            for v1 in element.getVertices():
                for v2 in element.getVertices()[localHandle:]:
                    if not v2.getIndex() in self.stencil[v1.getIndex()]:        self.stencil[v1.getIndex()].append(v2.getIndex())
                    if not v1.getIndex() in self.stencil[v2.getIndex()]:        self.stencil[v2.getIndex()].append(v1.getIndex())
                localHandle += 1



    def getVerticesFromRegion(self, region):
        verticesOnRegion = []
        for element in region.getElements():
            for v in element.getVertices():
                if verticesOnRegion.count(v) == 0:
                    verticesOnRegion.append(v)
        return verticesOnRegion

    def getVertices(self):
        return self.__vertices

    def getElements(self):
        return self.__elements

    def getRegions(self):
        return self.__regions

    def getBoundaries(self):
        return self.__boundaries

    def getBoundary(self, name):
        resp = False
        for bound in self.__boundaries:
            if bound.getName() == name:
                resp = True
                return bound
                break
        if resp == False:
            raise Exception("Boundary %s does not exist."%name)

    def getNumberOfElements( self ):
        return self.__numberOfElements

    def getNumberOfVertices( self ):
        return self.__numberOfVertices

    def getNumberOfRegions( self ):
        return self.__numberOfRegions

    def __computeVolumes( self ):
        for e in self.__elements:
            subVol = e.getVolume()/2.
            for v in e.getVertices():
                v.addToVolume( subVol )


def createGridData( L, numberOfNodes ):
    nodesCoord = np.linspace( 0, L, numberOfNodes )
    elemConn = [[0,0] for i in range(numberOfNodes-1)]
    # elemConn = np.zeros((numberOfNodes-1,2),dtype=int)
    for i in range( numberOfNodes-1 ):
        elemConn[i][0] = i
        elemConn[i][1] = i+1
    return nodesCoord, elemConn


if __name__ == '__main__':
##    L = 10.
##    nVertices = 6
##    nodesCoord, elemConn = createGridData( L, nVertices )
##
##    gridData = GridData()
##    gridData.setElementConnectivity( elemConn )
##    gridData.setNodeCoordinates( nodesCoord )
##    print gridData.regionElements
##
##    g = Grid_1D( gridData )
##    # # regions = g.getRegions()
##    # # print regions
##    # # for region in g.getRegions():
##    # #     print region.getName()
##    # #     for element in region.getElements():
##    # #         print element.getIndex(), element.getVertices()
##    # #     print '\n'


    # -------------- GRID DATA ----------------------------
    L_0 = 4.
    L_1 = 6.
    L = L_0 + L_1
    nVertices = 15
    nodesCoord, elemConn = createGridData(L, nVertices)
    gridData = GridData()
    gridData.setElementConnectivity(elemConn)
    gridData.setNodeCoordinates(nodesCoord)
    centroidCoord = []
    for e in elemConn:
        x_0 = gridData.nodeCoordinates[e[0]]
        x_1 = gridData.nodeCoordinates[e[1]]
        centroidCoord.append((x_0 + x_1)/2.)
    R1 = []
    R2 = []
    namesOfRegions = ['bottom', 'top']
    for e, x in enumerate(centroidCoord):
        if x <= L_0:
            R1.append(e)
        elif x > L_0:
            R2.append(e)
    gridData.setElementsToRegion(R1, namesOfRegions[0])
    gridData.setElementsToRegion(R2, namesOfRegions[1])
    # -----------------------------------------------------

    g = Grid_1D( gridData )
    regions = g.getRegions()
    print(regions)
    for region in g.getRegions():
        print(region.getName())
        for element in region.getElements():
            print(element.getIndex(), element.getVertices())
        print('\n')
