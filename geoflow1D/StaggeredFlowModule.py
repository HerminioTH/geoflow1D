def AssemblyDarcyVelocitiesToMatrix(linearSystem, grid, props, pShift=0):
	n = grid.getNumberOfVertices()
	for vertex in grid.getVertices()[1:-1]:
		eb = grid.getElements()[vertex.getIndex() - 1]
		ef = grid.getElements()[vertex.getIndex()]
		rb = grid.getRegions()[eb.getParentRegionIndex()]
		rf = grid.getRegions()[ef.getParentRegionIndex()]
		kb = props.k.getValue(rb)
		kf = props.k.getValue(rf)
		k = (kb + kf)/(kb*kf)
		linearSystem.addValueToMatrix(eb.getIndex() + pShift*n, eb.getIndex() + pShift*n, -k/props.mu)
		linearSystem.addValueToMatrix(eb.getIndex() + pShift*n, ef.getIndex() + pShift*n, +k/props.mu)
		linearSystem.addValueToMatrix(ef.getIndex() + pShift*n, ef.getIndex() + pShift*n, +k/props.mu)
		linearSystem.addValueToMatrix(ef.getIndex() + pShift*n, eb.getIndex() + pShift*n, -k/props.mu)




	# for region in grid.getRegions():
	# 	k = props.k.getValue(region)
	# 	for elem in region.getElements():
	# 		dx = elem.getLength()
	# 		face = elem.getFace()
	# 		A = face.getArea()
	# 		backVertex = face.getBackwardVertex()
	# 		forVertex = face.getForwardVertex()
	# 		bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
	# 		fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
	# 		diffusiveOperator = [k*A/props.mu, -k*A/props.mu]
	# 		for i,v in enumerate(elem.getVertices()):
	# 			flux = diffusiveOperator[i]
	# 			vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
	# 			linearSystem.addValueToMatrix(bIndex, vIndex, +flux/dx)
	# 			linearSystem.addValueToMatrix(fIndex, vIndex, -flux/dx)


def AssemblyBiotAccumulationToMatrix(linearSystem, grid, timeStep, props, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			bIndex = element.getVertices()[0].getIndex()
			fIndex = element.getVertices()[1].getIndex()
			value = (props.cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyVolumetricStrainToMatrix(linearSystem, grid, timeStep, props, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = props.biot.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			linearSystem.addValueToMatrix(bVertex.getIndex() + pShift*n, bVertex.getIndex() + (1-pShift)*n, value)
			linearSystem.addValueToMatrix(bVertex.getIndex() + pShift*n, fVertex.getIndex() + (1-pShift)*n, value)
			linearSystem.addValueToMatrix(fVertex.getIndex() + pShift*n, bVertex.getIndex() + (1-pShift)*n, -value)
			linearSystem.addValueToMatrix(fVertex.getIndex() + pShift*n, fVertex.getIndex() + (1-pShift)*n, -value)
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	vIndex = e.getVertices()[0].getIndex()
	linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, -props.biot.getValue(r)/timeStep)

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	vIndex = e.getVertices()[1].getIndex()
	linearSystem.addValueToMatrix(vIndex + pShift*n, vIndex + (1-pShift)*n, props.biot.getValue(r)/timeStep)

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

def AssemblyBiotAccumulationToVector(linearSystem, grid, props, timeStep, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			value = (props.cf*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))

def AssemblyVolumetricStrainToVector(linearSystem, grid, props, timeStep, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = props.biot.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			ub = u_old.getValue(bVertex)
			uf = u_old.getValue(fVertex)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(ub + uf))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(ub + uf))
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	bVertex = e.getVertices()[0]
	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -props.biot.getValue(r)*u_old.getValue(bVertex)/(1*timeStep))

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	fVertex = e.getVertices()[1]
	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, props.biot.getValue(r)*u_old.getValue(fVertex)/(1*timeStep))
