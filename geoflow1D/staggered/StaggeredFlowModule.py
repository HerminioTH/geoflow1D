# def AssemblyDarcyVelocitiesToMatrix(linearSystem, grid, props, pShift=0):
# 	for face in grid.getVertices()[1:-1]:
# 		eb = grid.getElements()[face.getIndex() - 1]
# 		ef = grid.getElements()[face.getIndex()]
# 		# kb = props.k.getValue()
# 		print(face.getIndex(), eb.getIndex(), ef.getIndex())

def AssemblyDarcyVelocitiesToMatrix(linearSystem, grid, props, pShift=0):
	for region in grid.getRegions():
		k = props.k.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			backVertex = face.getBackwardVertex()
			forVertex = face.getForwardVertex()
			bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
			fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
			diffusiveOperator = [k*A/props.mu, -k*A/props.mu]
			for i,v in enumerate(elem.getVertices()):
				flux = diffusiveOperator[i]
				vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
				linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
				linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)
