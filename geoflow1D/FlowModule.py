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


def AssemblyFluidFlowAccumulationToMatrix(linearSystem, grid, timeStep, props, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			bIndex = element.getVertices()[0].getIndex()
			fIndex = element.getVertices()[1].getIndex()
			value = phi*(cs + props.cf)*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyFluidFlowAccumulationToVector(linearSystem, grid, timeStep, props, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		cs = props.cs.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			value = phi*(cs + props.cf)*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))

def AssemblyBiotAccumulationToMatrix(linearSystem, grid, props, timeStep, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.c_s.getValue(region)
		for element in region.getElements():
			bIndex = element.getVertices()[0].getIndex()
			fIndex = element.getVertices()[1].getIndex()
			value = (props.c_f*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyBiotAccumulationToVector(linearSystem, grid, props, timeStep, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		phi = props.phi.getValue(region)
		alpha = props.biot.getValue(region)
		cs = props.c_s.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			value = (props.c_f*phi + cs*(alpha - phi))*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))


def AssemblyFixedStressAccumulationToMatrix(linearSystem, grid, timeStep, biotOnRegions, deltaOnVertices, bulkModulusOnRegions, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value/deltaOnVertices.getValue(bVertex))
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value/deltaOnVertices.getValue(fVertex))

def AssemblyFixedStressAccumulationToVector(linearSystem, grid, timeStep, biotOnRegions, deltaOnVertices, bulkModulusOnRegions, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex)/deltaOnVertices.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex)/deltaOnVertices.getValue(fVertex))








def AssemblyFixedStressAccumulationToMatrix_2(linearSystem, grid, timeStep, biotOnRegions, deltaOnRegions, bulkModulusOnRegions, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		delta = deltaOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
			linearSystem.addValueToMatrix(bIndex + pShift*n, bIndex + pShift*n, value)
			linearSystem.addValueToMatrix(fIndex + pShift*n, fIndex + pShift*n, value)

def AssemblyFixedStressAccumulationToVector_2(linearSystem, grid, timeStep, biotOnRegions, deltaOnRegions, bulkModulusOnRegions, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		alpha = biotOnRegions.getValue(region)
		modulus = bulkModulusOnRegions.getValue(region)
		delta = deltaOnRegions.getValue(region)
		for element in region.getElements():
			bVertex = element.getVertices()[0]
			fVertex = element.getVertices()[1]
			bIndex = bVertex.getIndex()
			fIndex = fVertex.getIndex()
			value = (alpha*alpha/modulus)*element.getSubVolume()/timeStep/delta
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*p_old.getValue(bVertex))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, value*p_old.getValue(fVertex))


def AssemblyVolumetricStrainToMatrix(linearSystem, grid, props, timeStep, pShift=0):
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



def AssemblyVolumetricStrainToVector2(linearSystem, grid, timeStep, biotOnRegions, u_new, u_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		value = biotOnRegions.getValue(region)/(2*timeStep)
		for e in region.getElements():
			bVertex = e.getVertices()[0]
			fVertex = e.getVertices()[1]
			uob = u_old.getValue(bVertex)
			uof = u_old.getValue(fVertex)
			unb = u_new.getValue(bVertex)
			unf = u_new.getValue(fVertex)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, value*(uob + uof - unb - unf))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(uob + uof - unb - unf))
	e = grid.getElements()[0]
	r = grid.getRegions()[e.getParentRegionIndex()]
	bVertex = e.getVertices()[0]
	linearSystem.addValueToVector(bVertex.getIndex() + pShift*n, -biotOnRegions.getValue(r)*(u_old.getValue(bVertex) - u_new.getValue(bVertex))/timeStep)

	e = grid.getElements()[-1]
	r = grid.getRegions()[e.getParentRegionIndex()]
	fVertex = e.getVertices()[1]
	linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, biotOnRegions.getValue(r)*(u_old.getValue(fVertex) - u_new.getValue(fVertex))/timeStep)







# ---------------------------- PHYSICAL INFLUENCE SCHEME ----------------------------------

def AssemblyPisToMassMatrix(linearSystem, grid, props, timeStep, pShift=0):
	for region in grid.getRegions():
		M = props.M.getValue(region)
		alpha = props.biot.getValue(region)
		for elem in region.getElements():
			dx = elem.getLength()
			face = elem.getFace()
			A = face.getArea()
			backVertex = face.getBackwardVertex()
			forVertex = face.getForwardVertex()
			bIndex = backVertex.getIndex() + pShift*grid.getNumberOfVertices()
			fIndex = forVertex.getIndex() + pShift*grid.getNumberOfVertices()
			value = alpha*alpha*dx*A/(8*M*timeStep)
			diffusiveOperator = [value, -value]
			for i,v in enumerate(elem.getVertices()):
				coef = diffusiveOperator[i]
				vIndex = v.getIndex() + pShift*grid.getNumberOfVertices()
				linearSystem.addValueToMatrix(bIndex, vIndex, +coef)
				linearSystem.addValueToMatrix(fIndex, vIndex, -coef)


def AssemblyPisToMassVector(linearSystem, grid, props, timeStep, p_old, pShift=0):
	n = grid.getNumberOfVertices()
	for region in grid.getRegions():
		M = props.M.getValue(region)
		alpha = props.biot.getValue(region)
		for e in region.getElements():
			dx = e.getLength()
			face = e.getFace()
			A = face.getArea()
			bVertex = face.getBackwardVertex()
			fVertex = face.getForwardVertex()
			pb = p_old.getValue(bVertex)
			pf = p_old.getValue(fVertex)
			value = -alpha*alpha*dx*A/(8*M*timeStep)
			linearSystem.addValueToVector(bVertex.getIndex() + pShift*n,  value*(pf - pb))
			linearSystem.addValueToVector(fVertex.getIndex() + pShift*n, -value*(pf - pb))

# ------------------------------------------------------------------------------------------
