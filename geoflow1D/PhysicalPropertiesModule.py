from .UtilsModule import *
from .FieldsModule import *
import json

class PhysicalProperties(object):
	def __init__(self, grid, folderName):
		self.folderName = folderName
		self.getFluidProps()
		self.getSolidProps(grid)

	def getFluidProps(self):
		fluid = getJsonData(self.folderName + "fluid.json")
		self.rho_f = fluid.get("WATER").get("Density")
		self.mu = fluid.get("WATER").get("Viscosity")
		self.cf = fluid.get("WATER").get("Compressibility")

	def getSolidProps(self, grid):
		solid = getJsonData(self.folderName + "solid.json")
		self.k = ScalarField(grid.getNumberOfRegions())
		self.phi = ScalarField(grid.getNumberOfRegions())
		self.cs = ScalarField(grid.getNumberOfRegions())
		self.rho_s = ScalarField(grid.getNumberOfRegions())
		self.M = ScalarField(grid.getNumberOfRegions())
		self.K = ScalarField(grid.getNumberOfRegions())
		self.Q = ScalarField(grid.getNumberOfRegions())
		self.biot = ScalarField(grid.getNumberOfRegions())
		for region in grid.getRegions():
			regionName = region.getName()
			self.k.setValue(region, solid.get(regionName).get("Permeability"))
			self.phi.setValue(region, solid.get(regionName).get("Porosity"))
			self.cs.setValue(region, solid.get(regionName).get("Compressibility"))
			self.rho_s.setValue(region, solid.get(regionName).get("Density"))
			G = solid.get(regionName).get("ShearModulus")
			nu = solid.get(regionName).get("PoissonsRatio")
			CS = solid.get(regionName).get("Compressibility")
			bulk = 2*G*(1 + nu)/(3*(1 - 2*nu))
			pWave = bulk + 4*G/3.
			alpha = 1 - CS*bulk
			self.M.setValue(region, pWave)
			self.K.setValue(region, bulk)
			self.biot.setValue(region, alpha)
			self.Q.setValue(region, 1./(self.cf*self.phi.getValue(region) + (alpha - self.phi.getValue(region))*CS))
		self.rho = ScalarField(grid.getNumberOfRegions())
		for region in grid.getRegions():
			self.rho.setValue(region, self.phi.getValue(region)*self.rho_f + (1 - self.phi.getValue(region))*self.rho_s.getValue(region))


def getJsonData(data_file):
    with open(data_file, "r") as jsonFile:
        data= json.load(jsonFile)
    return data
