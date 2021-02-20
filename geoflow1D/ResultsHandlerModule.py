import os
import shutil

class SaveResults( object ):
    def __init__( self, grid, fileName, folderName=None, fieldName='NONE', units='NONE', writeGridData=True ):
        self.__grid = grid
        self.__firstTimeField = True
        self.__fieldName = fieldName
        self.__units = units
        if folderName == None:
            self.__file = open( fileName, 'w' )
        else:
            self.__ensure_dir( folderName )
            self.__file = open( folderName + fileName, 'w' )
        if writeGridData:
            self.__writeCoordinates()
            self.__writeVolumes()

    def __ensure_dir( self, folderName ):
        if not os.path.exists( folderName ):
            os.makedirs( folderName )

    def __writeCoordinates( self ):
        self.__file.write( 'Coordinates:\n' )
        line = ''
        for v in self.__grid.getVertices():
            coord = v.getCoordinate()
            line += str(coord)+','
        self.__file.write( line )
        self.__file.write('\n')

    def __writeVolumes( self ):
        self.__file.write( 'Volumes:\n' )
        line = ''
        for v in self.__grid.getVertices():
            coord = v.getVolume()
            line += str(coord)+','
        self.__file.write( line )
        self.__file.write('\n')

    def copySettings(self, source, destination):
        for file in os.listdir(source):
            self.__ensure_dir(destination)
            shutil.copyfile(source + file, destination + file)

    def saveField(self, currentTime, field):
        if self.__firstTimeField:
            self.__file.write( '\n Time(s), %s (%s) on vertices\n'%(self.__fieldName, self.__units) )
            self.__firstTimeField = False
        line = ''
        line += str(currentTime) + ','
        for p in field:
            line += str(p) + ','
        self.__file.write( line + '\n' )

    def close( self ):
        self.__file.close()



class ReadResults(object):
    def __init__(self, fileName):
        self.fileName = fileName
        self.lines = open(fileName).readlines()

        self.__readCoordinates()
        self.__readVolumes()
        self.__readField()


    def __readCoordinates(self):
        self.coord = []
        line = self.lines[1].split(',')
        for i in line[:-1]:
            self.coord.append(float(i))

    def __readVolumes(self):
        self.volumes = []
        line = self.lines[3].split(',')
        for i in line[:-1]:
            self.volumes.append(float(i))

    def __readField(self):
        self.times = []
        self.field = []
        for i in range(6, len(self.lines)):
            line = self.lines[i].split(',')
            self.times.append(float(line[0]))
            fieldAtTime = []
            for v in line[1:-1]:
                fieldAtTime.append(float(v))
            self.field.append(fieldAtTime)

    def getSolutionAtTime(self, time):
        try:
            i = self.times.index(time)
            return self.field[i]
        except:
            raise Exception('Time %f does not exist.'%time)

    def getSolutionAtTimeStep(self, timeStep):
        try:
            return self.field[timeStep]
        except:
            raise Exception('TimeStep %i does not exist.'%timeStep)
