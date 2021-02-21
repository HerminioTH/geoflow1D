import numpy as np
import scipy.sparse.linalg as spla
from scipy.sparse import coo_matrix

class LinearSystemDense(object):
    def __init__(self, stencil, nDOF):
        self.stencil = stencil
        self.nDOF = nDOF
        self.nVertices = len(stencil)
        self.size = self.nDOF*self.nVertices

    def initialize(self):
        self.rhs = np.zeros(self.size)
        self.matrix = np.zeros((self.size, self.size))

    def addValueToMatrix(self, row, col, value):
        self.matrix[row][col] += value

    def setValueToMatrix(self, row, col, value):
        self.matrix[row][col] = value

    def getMatrixValue(self, row, col):
        return self.matrix[row][col]

    def addValueToVector(self, row, value):
        self.rhs[row] += value

    def setValueToVector(self, row, value):
        self.rhs[row] = value

    def getVectorValue(self, row):
        return self.rhs[row]

    def getMatrix(self):
        return self.matrix

    def getVector(self):
        return self.rhs

    def applyDirichlet(self, row, value):
        self.applyDirichletToMatrix(row, value)
        self.applyDirichletToVector(row, value)

    def applyDirichletToMatrix(self, row, value):
        for col in range(self.size):
            self.matrix[row][col] = 0.0
        self.matrix[row][row] = 1.0

    def applyDirichletToVector(self, row, value):
        self.rhs[row] = value

    def applyNeumann(self, row, value):
        self.rhs[row] += value

    def applyBoundaryConditionsToMatrix(self, grid, boundSettings, shift=0):
        n = grid.getNumberOfVertices()
        for bName in boundSettings.keys():
            bound = grid.getBoundary(bName)
            bType = boundSettings.get(bName).get("Type")
            bValue = boundSettings.get(bName).get("Value")
            if bType == "Dirichlet":
                self.applyDirichletToMatrix(bound.getVertex().getIndex() + shift*n, bValue)

    def applyBoundaryConditionsToVector(self, grid, boundSettings, shift=0):
        n = grid.getNumberOfVertices()
        for bName in boundSettings.keys():
            bound = grid.getBoundary(bName)
            bType = boundSettings.get(bName).get("Type")
            bValue = boundSettings.get(bName).get("Value")
            if bType == "Dirichlet":
                self.applyDirichletToVector(bound.getVertex().getIndex() + shift*n, bValue)
            elif bType == "Neumann":
                self.applyNeumann(bound.getVertex().getIndex() + shift*n, bValue)
            else:
                raise Exception("Boundary type %s is not supported."%bType)

    def computeResidue(self):
        self.residue = self.rhs - self.matrix.dot(self.solution)

    def getSolution(self):
        return self.solution

    def eraseVector(self):
        self.rhs = np.zeros(self.size)

    def eraseMatrix(self):
        self.matrix = np.zeros( (self.size, self.size) )

    def resetLinearSystem(self):
        self.eraseVector()
        self.eraseMatrix()

    def splitSolution(self, n):
        return self.solution[:n], self.solution[n:]



class LinearSystemCOO(LinearSystemDense):
    def initialize(self):
        self.counter = 0
        self.rhs = np.zeros(self.size)
        row, col = [], []
        self.indPtr = []
        ct = 0
        for dof_row in range(self.nDOF):
            for r, vertexStencil in enumerate(self.stencil):
                self.indPtr.append(ct)
                for dof_col in range(self.nDOF):
                    for c in vertexStencil:
                        row.append(r + dof_row*self.nVertices)
                        col.append(c + dof_col*self.nVertices)
                        ct += 1
        self.indPtr.append(ct)
        data = np.zeros(len(row))
        self.matrix = coo_matrix((data, (row, col)), shape=(self.size, self.size))

    def getIndex(self, row, col):
        pos1 = self.indPtr[row]
        pos2 = self.indPtr[row+1]
        for i, colIndex in enumerate(self.matrix.col[pos1:pos2]):
            if colIndex == col:
                return pos1+i

    def addValueToMatrix(self, row, col, value):
        index = self.getIndex(row, col)
        self.matrix.data[index] += value

    def setValueToMatrix(self, row, col, value):
        index = self.getIndex(row, col)
        self.matrix.data[index] += value

    def applyDirichletToMatrix(self, row, value):
        pos1 = self.indPtr[row]
        pos2 = self.indPtr[row+1]
        self.matrix.data[pos1:pos2] = [0 for i in range(pos1,pos2)]
        self.setValueToMatrix(row, row, 1.0)





def splitSolution(sol, n):
    return sol[:n],sol[n:]
