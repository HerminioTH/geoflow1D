def AssemblyDarcyVelocitiesToMatrix(linearSystem, grid, props, pShift=0):
	n = grid.getNumberOfVertices()
	for vertex in grid.getVertices()[1:-1]:
		eb = grid.getElements()[vertex.getIndex() - 1]
		ef = grid.getElements()[vertex.getIndex()]
		rb = grid.getRegions()[eb.getParentRegionIndex()]
		rf = grid.getRegions()[ef.getParentRegionIndex()]
		kb = props.k.getValue(rb)
		kf = props.k.getValue(rf)
		k = (kb*kf)/(kb + kf)
		dx = abs(eb.getCentroid() - ef.getCentroid())
		value = k/(props.mu*dx)
		linearSystem.addValueToMatrix(eb.getIndex() + pShift*n, eb.getIndex() + pShift*n, +value)
		linearSystem.addValueToMatrix(eb.getIndex() + pShift*n, ef.getIndex() + pShift*n, -value)
		linearSystem.addValueToMatrix(ef.getIndex() + pShift*n, ef.getIndex() + pShift*n, +value)
		linearSystem.addValueToMatrix(ef.getIndex() + pShift*n, eb.getIndex() + pShift*n, -value)

	last_elem = grid.getElements()[-1]
	linearSystem.addValueToMatrix(last_elem.getIndex() + pShift*n, last_elem.getIndex() + pShift*n, k/(props.mu*dx/2))


def AssemblyBiotAccumulationToMatrix(linearSystem, grid, timeStep, props, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			value = (props.cf*phi + cs*(alpha - phi))*element.getVolume()/timeStep
			linearSystem.addValueToMatrix(element.getIndex() + pShift*n, element.getIndex() + pShift*n, value)

def AssemblyBiotAccumulationToVector(linearSystem, grid, props, timeStep, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			value = (props.cf*phi + cs*(alpha - phi))*element.getVolume()/timeStep
			linearSystem.addValueToVector(element.getIndex() + pShift*n, value*p_old.getValue(element))


def AssemblyVolumetricStrainToMatrix(linearSystem, grid, timeStep, props, pShift=0):
	n = grid.getNumberOfVertices()
	for vertex in grid.getVertices()[1:-1]:
		eb = grid.getElements()[vertex.getIndex() - 1]
		ef = grid.getElements()[vertex.getIndex()]
		rb = grid.getRegions()[eb.getParentRegionIndex()]
		rf = grid.getRegions()[ef.getParentRegionIndex()]
		biot_b = props.biot.getValue(rb)
		biot_f = props.biot.getValue(rf)
		biot = (biot_b*biot_f)/(biot_b + biot_f)
		dx = abs(eb.getCentroid() - ef.getCentroid())
		value = biot/timeStep
		linearSystem.addValueToMatrix(eb.getIndex() + pShift*n, vertex.getIndex() + (1-pShift)*n, +value)
		linearSystem.addValueToMatrix(ef.getIndex() + pShift*n, vertex.getIndex() + (1-pShift)*n, -value)

	first_elem = grid.getElements()[0]
	last_elem = grid.getElements()[-1]
	first_vertex = grid.getVertices()[0]
	last_vertex = grid.getVertices()[-1]
	linearSystem.addValueToMatrix(first_elem.getIndex() + pShift*n, first_vertex.getIndex() + (1-pShift)*n, -value)
	linearSystem.addValueToMatrix(last_elem.getIndex() + pShift*n, last_vertex.getIndex() + (1-pShift)*n, +value)


def AssemblyVolumetricStrainToVector(linearSystem, grid, props, timeStep, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for vertex in grid.getVertices()[1:-1]:
		eb = grid.getElements()[vertex.getIndex() - 1]
		ef = grid.getElements()[vertex.getIndex()]
		rb = grid.getRegions()[eb.getParentRegionIndex()]
		rf = grid.getRegions()[ef.getParentRegionIndex()]
		biot_b = props.biot.getValue(rb)
		biot_f = props.biot.getValue(rf)
		biot = (biot_b*biot_f)/(biot_b + biot_f)
		dx = abs(eb.getCentroid() - ef.getCentroid())
		value = biot/timeStep
		linearSystem.addValueToVector(eb.getIndex() + pShift*n, +value*u_old.getValue(vertex))
		linearSystem.addValueToVector(ef.getIndex() + pShift*n, -value*u_old.getValue(vertex))

	first_elem = grid.getElements()[0]
	last_elem = grid.getElements()[-1]
	first_vertex = grid.getVertices()[0]
	last_vertex = grid.getVertices()[-1]
	linearSystem.addValueToVector(first_elem.getIndex() + pShift*n, -value*u_old.getValue(first_vertex))
	linearSystem.addValueToVector(last_elem.getIndex() + pShift*n, +value*u_old.getValue(last_vertex))


# def AssemblyVolumetricStrainToVector(linearSystem, grid, props, timeStep, u_old, pShift=0):
# 	n = grid.getNumberOfVertices()
# 	for region in grid.getRegions():
# 		value = props.biot.getValue(region)/(2*timeStep)
# 		for e in region.getElements():
# 			bVertex = e.getVertices()[0]
# 			fVertex = e.getVertices()[1]
# 			ub = u_old.getValue(bVertex)
# 			uf = u_old.getValue(fVertex)
# 			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(ub + uf))
# 			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(ub + uf))
# 	e = grid.getElements()[0]
# 	r = grid.getRegions()[e.getParentRegionIndex()]
# 	bVertex = e.getVertices()[0]
# 	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -props.biot.getValue(r)*u_old.getValue(bVertex)/(1*timeStep))

# 	e = grid.getElements()[-1]
# 	r = grid.getRegions()[e.getParentRegionIndex()]
# 	fVertex = e.getVertices()[1]
# 	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, props.biot.getValue(r)*u_old.getValue(fVertex)/(1*timeStep))




def AssemblyDarcyVelocitiesToVector(linearSystem, grid, props, gravity, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		k = props.k.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			value = k*A*props.rho_f*gravity/props.mu
			bIndex = face.getBackwardVertex().getIndex() + pShift*n
			fIndex = face.getForwardVertex().getIndex() + pShift*n
			linearSystem.addValueToVector(bIndex, -value)
			linearSystem.addValueToVector(fIndex,  value)
