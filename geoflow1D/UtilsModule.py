import numpy as np
import matplotlib.pylab as plt
import json


def getJsonData(data_file):
    with open(data_file, "r") as jsonFile:
        data= json.load(jsonFile)
    return data

def saveDataToJson(fileName, data):
    with open(fileName, "w") as jsonFile:
        json.dump(data, jsonFile, indent=3)

def computeNormL2(vector, grid):
    soma = 0
    for i,vertex in enumerate(grid.getVertices()):
        soma += vertex.getVolume()*vector[i]*vector[i]
    return soma**0.5

def computeNormInf(vector):
    return vector.max()

def computeRate(error):
    return (np.log10(error[0]) - np.log10(max(1e-200, error[-1])))/len(error)


def computeNormL2OnRegions(vector, grid):
    L2 = [[] for i in range(grid.getNumberOfRegions())]
    for region in grid.getRegions():
        soma = 0
        for i,vertex in enumerate(grid.getVerticesFromRegion(region)):
            soma += vertex.getVolume()*vector[i]*vector[i]
        L2[region.getIndex()] = soma**0.5
    return np.array(L2)


def computeConsolidationCoefficient(data_fluid, data_solid, bodyName="BODY"):
    mu = data_fluid["WATER"]["Viscosity"]
    cf = data_fluid["WATER"]["Compressibility"]
    k = data_solid[bodyName]["Permeability"]
    nu = data_solid[bodyName]["PoissonsRatio"]
    phi = data_solid[bodyName]["Porosity"]
    G = data_solid[bodyName]["ShearModulus"]
    M = 2*G*(1 - nu)/(1 - 2*nu)
    Q = 1./(cf*phi)
    return k*Q*M/(mu*(M + Q))

def computeTime(props_folder, L):
    data_fluid = getJsonData(props_folder + "//fluid.json")
    data_solid = getJsonData(props_folder + "//solid.json")
    c = computeConsolidationCoefficient(data_fluid, data_solid)
    return L*L/c

def updateTimeStep(case_folder, timeStep):
    data = getJsonData(case_folder + "numerical_settings.json")
    data["TransientCycle"]["TimeStep"] = timeStep
    saveDataToJson(case_folder + "numerical_settings.json", data)

def updateFinalTime(case_folder, finalTime):
    data = getJsonData(case_folder + "numerical_settings.json")
    data["TransientCycle"]["FinalTime"] = finalTime
    saveDataToJson(case_folder + "numerical_settings.json", data)
