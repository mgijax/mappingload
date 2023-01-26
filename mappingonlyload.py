'''
#
# Purpose:
#
#	To preprocess curator created mapping file into the format needed by mappingload.py
#       in order to load:
#
#	MLD_Expts
#	MLD_Expt_Marker
#	MLD_Notes
#
#	To update Marker Chromosome and Band Assignments (optional)
#
#	MRK_Marker
#
# Assumes:
#
#	That no one else is adding Mapping or Accession IDs records to the 
#       database.
#
# Input(s):
#
#	A tab-delimited file in the format:
#		field 1: MGI Acc ID of Symbol
#		field 2: Chromosome
#		field 3: Update Marker Chromosome (yes/no)
#		field 4: Band (optional)
#		field 5: Assay Type
#		field 6: Description
#
# Input:
#       curator created file
#
# Output:
#	SQL file:
#		file of SQL commands for updating Marker chromosomes and bands
#	Diagnostics file of all input parameters and SQL commands
#	Error file
#
# Processing:
#
#	1. Verify the J: is valid.
#	    If the verification fails, report the error and stop.
#       2. Verify the createdBy is valid.
#           If the verification fails, report the error and stop.
#	3. Create the pre-processed file
#       4.  Determine the next available  Experiment Marker Key
#
# History:
#
# sc	0/24/2024
#
'''

import sys
import os
import getopt
import re
import db
import mgi_utils
import loadlib

# globals

CRT = '\n'
TAB = '\t'
PIPE = '|'

passwordFileName = os.getenv('MGD_DBPASSWORDFILE')

inputFile = ''		
outputFile = ''
logFile = ''		
sqlFile = ''

inputFileName = ''	
outputFileName = ''
logFileName = ''
sqlFileName = ''

nextMappingKey = 1000

jnum = ''
createdBy = ''

markerDict = {}

mode = os.getenv('MAPPINGMODE')

DEBUG = 0

if mode == 'preview':
    DEBUG = 1

def exit(status, message = None):
        '''
        # requires: status, the numeric exit status (integer)
        #           message (string)
        #
        # effects:
        # Print message to stderr and exits
        #
        # returns:
        #
        '''
 
        if message is not None:
                sys.stderr.write('\n' + str(message) + '\n')
 
        try:
                inputFile.close()
        except:
                pass

        db.useOneConnection()
        sys.exit(status)
 
def init():
        '''
        #
        # Initializes global file descriptors/file names, keys, lookups
        #
        '''
        global nextMappingKey, jnum, createdBy,jnum, createdBy, inputFile, outputFile
        global logFile, sqlFileName, sqlFile, markerDict

        results = db.sql(''' select nextval('mld_expt_marker_seq') as maxKey ''', 'auto')
        nextMappingKey = results[0]['maxKey']

        jnum = os.getenv('JNUM')
        createdBy = os.getenv('CREATEDBY')

        inputFileName = os.getenv('MAPPINGONLYDATAFILE') 

        outputFileName =  os.getenv('MAPPINGDATAFILE')

        logFileName = os.getenv('MAPPINGONLYDATALOG')

        sqlFileName = os.getenv('MAPPINGONLYSQLFILE')

        try:
            inputFile = open(inputFileName, 'r')
        except:
            exit(1, 'Could not open file %s\n' % inputFileName)
                
        try:
            outputFile = open(outputFileName, 'w')
        except:
            exit(1, 'Could not open file %s\n' % outputFileName)

        try:
            logFile = open(logFileName, 'a')
        except:
            exit(1, 'Could not open file %s\n' % logFileName)
               
        try:
            sqlFile = open(sqlFileName, 'w')
        except:
            exit(1, 'Could not open file %s\n' % sqlFileName)


        results = db.sql('''select a.accid, m._Marker_key
                from MRK_Marker m, MRK_Acc_View a 
                where a._Object_key = m._Marker_key 
                and m._Organism_key = 1''', 'auto')
        for r in results:
                markerKey = r['_Marker_key']
                markerID = r['accid']
                markerDict[markerID] = markerKey

        return 0

def processFile():
        '''
        # requires:
        #
        # effects:
        #	Reads input file
        #	Verifies and Processes each line in the input file
        #       Writes to intermediate output file
        #
        # returns:
        #	nothing
        #
        '''
        global nextMappingKey

        lineNum = 0

        # For each line in the input file

        for line in inputFile.readlines():
            #print('line: %s' % line)

            lineNum += 1
            # Split the line into tokens
            tokens = str.split(line, '\t')
            #print('tokens: %s' % tokens)
            try:
                markerID = tokens[0]
                chromosome = tokens[1]
                updateChr = tokens[2]
                band = tokens[3]
                assay = tokens[4]
                description = str.strip(tokens[5])
            except:
                exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

            outputFile.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextMappingKey, PIPE, markerID, PIPE, chromosome, PIPE, updateChr, PIPE, band, PIPE, assay, PIPE, description, PIPE, jnum, PIPE, createdBy, CRT))
            nextMappingKey += 1

            # update marker's chromosome...
            if markerID in markerDict:
                markerKey = markerDict[markerID]
                if updateChr == 'yes':
                        sqlFile.write('''update MRK_Marker
                                set modification_date = now(), chromosome = '%s' 
                                where _Marker_key = %s\n;\n''' % (chromosome, markerKey))

                # update cytogenetic band, if it is provided
                if band != "":
                        sqlFile.write('''update MRK_Marker 
                                set modification_date = now(), cytogeneticOffset = '%s'
                                where _Marker_key = %s\n;\n''' % (band, markerKey))
        outputFile.close()
        sqlFile.close()
        print ('DEBUG: %s' % DEBUG)
        if not DEBUG:
            cmd = 'psql -h %s -d %s -U %s -f %s -o %s.log' % (db.get_sqlServer(), db.get_sqlDatabase(), db.get_sqlUser(), sqlFileName, sqlFileName)
            print('cmd: %s' % cmd)
            os.system(cmd)
            db.commit()

        return 0

#
# Main
#

print('mappingonlyload:init()')
init()

print ('mappingonlyload:processFile()')
processFile()

exit(0)
