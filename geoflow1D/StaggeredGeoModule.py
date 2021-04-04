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
