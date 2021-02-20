def AssemblyStiffnessMatrix(linearSystem, grid, props, uShift=0):
    for region in grid.getRegions():
        M = props.M.getValue(region)
        for e in region.getElements():
            dx = e.getLength()
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            forceOperator = [-M/dx, M/dx]
            localIndex = 0
            for v in e.getVertices():
                flux = forceOperator[localIndex]
                vIndex = v.getIndex() + uShift*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, vIndex, flux )
                linearSystem.addValueToMatrix( fIndex, vIndex, -flux )
                localIndex += 1

def AssemblyGravityToVector(linearSystem, grid, props, gravity, uShift=0):
    n = grid.getNumberOfVertices()
    for region in grid.getRegions():
        rho = props.rho.getValue(region)
        for elem in region.getElements():
            face = elem.getFace()
            bVertex = face.getBackwardVertex()
            fVertex = face.getForwardVertex()
            value = -rho*gravity*elem.getSubVolume()
            linearSystem.addValueToVector(bVertex.getIndex() + uShift*n, value)
            linearSystem.addValueToVector(fVertex.getIndex() + uShift*n, value)

def AssemblyPorePressureToMatrix(linearSystem, grid, props, uShift=0):
    for region in grid.getRegions():
        alpha = props.biot.getValue(region)
        for e in region.getElements():
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            for i,v in enumerate(e.getVertices()):
                col = v.getIndex() + (1-uShift)*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, col, -alpha/2 )
                linearSystem.addValueToMatrix( fIndex, col, +alpha/2 )

def AssemblyPorePressureToVector(linearSystem, grid, biot, pField, uShift=0):
    for region in grid.getRegions():
        alpha = biot.getValue(region)
        for e in region.getElements():
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            pBack = pField.getValue(backVertex)
            pFron = pField.getValue(forVertex)
            value = alpha/2.
            linearSystem.addValueToVector(bIndex, value*pBack)
            linearSystem.addValueToVector(bIndex, value*pFron)
            linearSystem.addValueToVector(fIndex, -value*pBack)
            linearSystem.addValueToVector(fIndex, -value*pFron)


# ------------------------- PHYSICAL INFLUENCE SCHEME - FULL -------------------------------
def AssemblyPisFullToGeoMatrix(linearSystem, grid, props, timeStep, uShift=0):
    for region in grid.getRegions():
        M = props.M.getValue(region)
        alpha = props.biot.getValue(region)
        k = props.k.getValue(region)
        Q = 1/(props.c_f*props.phi.getValue(region) + props.c_s.getValue(region)*(1 - props.phi.getValue(region)))
        for e in region.getElements():
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            dx = e.getLength()
            B = props.mu*dx*alpha/(8*k*timeStep)
            pisOperator = [-alpha*B, alpha*B]
            for localIndex, v in enumerate(e.getVertices()):
                coef = pisOperator[localIndex]
                vIndex = v.getIndex() + uShift*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, vIndex, coef )
                linearSystem.addValueToMatrix( fIndex, vIndex, -coef )

def AssemblyPisFullToGeoVector(linearSystem, grid, props, timeStep, uOldField, uShift=0):
    for region in grid.getRegions():
        M = props.M.getValue(region)
        alpha = props.biot.getValue(region)
        k = props.k.getValue(region)
        Q = 1/(props.c_f*props.phi.getValue(region) + props.c_s.getValue(region)*(1 - props.phi.getValue(region)))
        for e in region.getElements():
            f = e.getFace()
            dx = e.getLength()
            B = props.mu*dx*alpha/(8*k*timeStep)
            bVertex = f.getBackwardVertex()
            fVertex = f.getForwardVertex()
            ub = uOldField.getValue(bVertex)
            uf = uOldField.getValue(fVertex)
            linearSystem.addValueToVector(bVertex.getIndex() + uShift*grid.getNumberOfVertices(),  alpha*B*(uf - ub))
            linearSystem.addValueToVector(fVertex.getIndex() + uShift*grid.getNumberOfVertices(), -alpha*B*(uf - ub))

def AssemblyPisFullToPorePressureMatrix(linearSystem, grid, props, timeStep, uShift=0):
    for region in grid.getRegions():
        alpha = props.biot.getValue(region)
        k = props.k.getValue(region)
        Q = 1/(props.c_f*props.phi.getValue(region) + props.c_s.getValue(region)*(1 - props.phi.getValue(region)))
        for e in region.getElements():
            dx = e.getLength()
            f = e.getFace()
            bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
            A = props.mu*dx*dx/(16*k*Q*timeStep)
            for i,v in enumerate(e.getVertices()):
                col = v.getIndex() + (1-uShift)*grid.getNumberOfVertices()
                linearSystem.addValueToMatrix( bIndex, col, +alpha*A )
                linearSystem.addValueToMatrix( fIndex, col, -alpha*A )

def AssemblyPisFullToPorePressureVector(linearSystem, grid, props, timeStep, pOldField, uShift=0):
    for region in grid.getRegions():
        M = props.M.getValue(region)
        alpha = props.biot.getValue(region)
        k = props.k.getValue(region)
        Q = 1/(props.c_f*props.phi.getValue(region) + props.c_s.getValue(region)*(1 - props.phi.getValue(region)))
        for e in region.getElements():
            f = e.getFace()
            dx = e.getLength()
            A = props.mu*dx*dx/(16*k*Q*timeStep)
            bVertex = f.getBackwardVertex()
            fVertex = f.getForwardVertex()
            pb = pOldField.getValue(bVertex)
            pf = pOldField.getValue(fVertex)
            linearSystem.addValueToVector(bVertex.getIndex() + uShift*grid.getNumberOfVertices(),  alpha*A*(pf + pb))
            linearSystem.addValueToVector(fVertex.getIndex() + uShift*grid.getNumberOfVertices(), -alpha*A*(pf + pb))
# ------------------------------------------------------------------------------------------



# ---------------------------- LOOP BY ELEMENTS ----------------------------------
def AssemblyStiffnessMatrix_e(linearSystem, grid, modulus, uShift=0):
    for element in grid.getElements():
        value = modulus.getValue(element)
        dx = element.getLength()
        f = element.getFace()
        bIndex = f.getBackwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
        fIndex = f.getForwardVertex().getIndex() + uShift*grid.getNumberOfVertices()
        forceOperator = [-value/dx, value/dx]
        localIndex = 0
        for v in element.getVertices():
            flux = forceOperator[localIndex]
            vIndex = v.getIndex() + uShift*grid.getNumberOfVertices()
            linearSystem.addValueToMatrix( bIndex, vIndex, flux )
            linearSystem.addValueToMatrix( fIndex, vIndex, -flux )
            localIndex += 1

def AssemblyGravityToVector_e(linearSystem, grid, densityOnElements, gravity, uShift=0):
    for element in grid.getElements():
        rho = densityOnElements.getValue(element)
        face = element.getFace()
        bVertex = face.getBackwardVertex()
        fVertex = face.getForwardVertex()
        value = -rho*gravity*element.getSubVolume()
        linearSystem.addValueToVector(bVertex.getIndex() + uShift*grid.getNumberOfVertices(), value)
        linearSystem.addValueToVector(fVertex.getIndex() + uShift*grid.getNumberOfVertices(), value)

def AssemblyPorePressureToVector_e(linearSystem, grid, biotOnElements, pField, uShift=0):
    for element in grid.getElements():
        alpha = biotOnElements.getValue(element)
        f = element.getFace()
        backVertex = f.getBackwardVertex()
        forVertex = f.getForwardVertex()
        bIndex = backVertex.getIndex() + uShift*grid.getNumberOfVertices()
        fIndex = forVertex.getIndex() + uShift*grid.getNumberOfVertices()
        pBack = pField.getValue(backVertex)
        pFron = pField.getValue(forVertex)
        value = alpha/2.
        linearSystem.addValueToVector(bIndex, value*pBack)
        linearSystem.addValueToVector(bIndex, value*pFron)
        linearSystem.addValueToVector(fIndex, -value*pBack)
        linearSystem.addValueToVector(fIndex, -value*pFron)



# if __name__ == '__main__':
#     from GridLib import *
#     from FieldsLib import *
#     from LinearSystemLib import *

#     L_0 = 4.
#     L_1 = 6.
#     L = L_0 + L_1
#     nVertices = 10
#     nodesCoord, elemConn = createGridData( L, nVertices )

#     # -------------- GRID DATA ----------------------------
#     gridData = GridData()
#     gridData.setElementConnectivity( elemConn )
#     gridData.setNodeCoordinates( nodesCoord )
#     centroidCoord = []
#     for e in elemConn:
#         x_0 = gridData.nodeCoordinates[e[0]]
#         x_1 = gridData.nodeCoordinates[e[1]]
#         centroidCoord.append((x_0 + x_1)/2.)
#     region_1 = []
#     region_2 = []
#     namesOfRegions = ['bottom', 'top']
#     for e, x in enumerate(centroidCoord):
#         if x <= L_0:
#             region_1.append(e)
#         elif x > L_0:
#             region_2.append(e)
#     gridData.setElementsToRegion(region_1, 'lower_layer')
#     gridData.setElementsToRegion(region_2, 'upper_layer')
#     g = Grid_1D( gridData )

#     for region in g.getRegions():
#         print(region.getName())
#         for element in region.getElements():
#             vec = [element.getIndex()]
#             for v in element.getVertices():
#                 vec.append(v.getIndex())
#             print(vec)
#         print('\n')
#     # -----------------------------------------------------

#     # -------------- PROPERTIES ----------------------------
#     M = ScalarField(g.getNumberOfRegions())
#     M.setValue(g.getRegions()[0], 1000.)
#     M.setValue(g.getRegions()[1], 2000.)
#     # -----------------------------------------------------

#     # -------------- LINEAR SYSTEM ------------------------
#     ls = LinearSystem(g.getNumberOfVertices())
#     AssemblyStiffnessMatrix(ls, g, M, 0)
#     ls.applyDirichlet(0, 0)
#     ls.applyNeumann(-1, -1000)
#     print(g.getNumberOfVertices())
#     print(ls.getMatrix())
#     print(ls.getVector())
#     ls.solve()
#     print(ls.getSolution())
#     # -----------------------------------------------------
