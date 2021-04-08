import os
import shutil

class StaggeredSaveResults( object ):
    def __init__( self, grid, fileName, folderName=None, fieldName='NONE', units='NONE' ):
        self.__grid = grid
        self.__firstTimeField = True
        self.__fieldName = fieldName
        self.__units = units
        if folderName == None:
            self.__file = open( fileName, 'w' )
        else:
            self.__ensure_dir( folderName )
            self.__file = open( folderName + fileName, 'w' )
        self.__writeCoordinates()
        self.__writeVolumes()

    def __ensure_dir( self, folderName ):
        if not os.path.exists( folderName ):
            os.makedirs( folderName )

    def __writeCoordinates( self ):
        self.__file.write( 'Coordinates:\n' )
        line = ''
        for e in self.__grid.getElements():
            coord = (e.getVertices()[0].getCoordinate() + e.getVertices()[1].getCoordinate())/2
            line += str(coord)+','
        self.__file.write(line)
        self.__file.write('\n')

    def __writeVolumes( self ):
        self.__file.write( 'Volumes:\n' )
        line = ''
        for e in self.__grid.getElements():
            volume = e.getVolume()
            line += str(volume)+','
        self.__file.write(line)
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
