#!/usr/local/bin/python

'''
#
# Purpose:
#
#	To load new mapping records into Mapping structures:
#
#	MLD_Expts
#	MLD_Marker
#	MLD_Expt_Marker
#	MLD_Notes
#
#	To update Marker Chromosome and Band Assignments (optional)
#
#	MRK_Marker
#
# Assumes:
#
#	That no one else is adding Mapping or Accession IDs records to the database.
#
# Side Effects:
#
#	None
#
# Input(s):
#
#	A tab-delimited file in the format:
#		field 1: MGI Acc ID of Symbol
#		field 2: Chromosome
#		field 3: Update Marker Chromosome (yes/no)
#		field 4: Band (optional)
#		field 5: Assay Type
#
#	last line of file will be the experiment note (optional)
#
# Parameters:
#	-S = database server
#	-D = database
#	-U = user
#	-P = password file
#	-M = mode (incremental, full, preview)
#	-I = input file of mapping data
#	-R = Reference (J: in format #####)
#	-E = Experiment Type ("TEXT")
#
#	processing modes:
#		incremental - append the data to the existing Experiments (if they exist)
#			    - create the Experiments if they don't exist
#
#		full - delete the data from the existing Experiments (if they exist)
#		     - create the Experiments if they don't exist
#
#		preview - perform all record verifications but do not load the data or
#		          make any changes to the database.  used for testing or to preview
#			  the load.
#
# Output:
#
#       5 BCP files:
#
#       ACC_Accession.bcp               Accession records
#       MLD_Expts.bcp                   master Experiment records
#       MLD_Marker.bcp                  master Marker records for J:
#       MLD_Expt_Marker.bcp             Marker records for each Experiment
#       MLD_Notes.bcp             	master Experiment notes
#
#	1 SQL file:
#		file of SQL commands for updating Marker chromosomes and bands
#
#	Diagnostics file of all input parameters and SQL commands
#	Error file
#
# Processing:
#
#	1. Verify Mode.
#		if mode = incremental:  process records
#		if mode = full:  delete existing records and process
#		if mode = preview:  set "DEBUG" to True
#
#	2.  Verify the J: is valid.
#	    If the verification fails, report the error and stop.
#
#	3.  Create the master Experiment records and Accession records.
#	    If Experiment records already exist for the Reference, delete the details
#	    (but not the master experiment records themselves).
#
#	For each line in the input file:
#
#	1.  Verify the Marker Acc ID is valid.  Duplicates are reported as errors.
#	    If the verification fails, report the error and skip the record.
#
#	2.  Verify the Assay is valid.
#	    If the verification fails, report the error and skip the record.
#
#	3.  Verify the Chromosome is valid.
#	    If the verification fails, report the error and skip the record.
#
#	5.  Determine the Experiment key for the Chromosome.
#
#	5.  Create MLD_Marker record for the Marker.
#
#	6.  Create MLD_Expt_Marker record for the Marker.
#
# History:
#
# lec	01/30/2003
#	- TR 3928 (Fantom2 load)
#		. provide Assay value for each input record
#		. provide Update Chromosome? value for each input record
#
# lec	08/29/2002
#	- created for TR 4010, but other loads can use this module
#
'''

import sys
import os
import string
import getopt
import regsub
import db
import mgi_utils
import accessionlib

#globals

DEBUG = 0		# set DEBUG to false unless preview mode is selected

inputFile = ''		# file descriptor
diagFile = ''		# file descriptor
errorFile = ''		# file descriptor

exptFile = ''		# file descriptor
markerFile = ''		# file descriptor
exptMarkerFile = ''	# file descriptor
accFile = ''		# file descriptor
noteFile = ''		# file descriptor
sqlFile = ''		# file descriptor

diagFileName = ''	# file name
errorFileName = ''	# file name
passwordFileName = ''	# file name

exptFileName = ''	# file name
markerFileName = ''	# file name
exptMarkerFileName = ''	# file name
accFileName = ''	# file name
noteFileName = ''	# file name
sqlFileName = ''	# file name

mode = ''		# processing mode

markerDict = {}		# dictionary of marker accids and marker keys/symbols
chromosomeList = []	# list of valid mouse chromosome
exptDict = {}		# dictionary of chromosome/experiment key values
seqExptDict = {}	# dictionary of experiment marker sequence values

logicalDBKey = 1
mgiTypeKey = 4          # Experiment
mgiPrefix = "MGI:"
exptType = "TEXT"
bcpdelim = "|"

referenceKey = 0	# Reference Key
alleleKey = ''		# MLD_Expt_Marker._Allele_key
description = ''	# MLD_Expt_Marker.description
matrixData = 0		# MLD_Extt_Marker.matrixData

cdate = mgi_utils.date('%m/%d/%Y')	# current date

def showUsage():
	'''
	# requires:
	#
	# effects:
	# Displays the correct usage of this program and exits
	# with status of 1.
	#
	# returns:
	'''
 
	usage = 'usage: %s -S server\n' % sys.argv[0] + \
		'-D database\n' + \
		'-U user\n' + \
		'-P password file\n' + \
		'-M mode\n' + \
		'-I input file\n' + \
		'-R J# (####)\n' + \
		'-A Assay Type (ex. "sequence analysis")\n' + \
		'-E Experiment Type (ex. "TEXT", "TEXT-Physical Mapping")\n'
	exit(1, usage)
 
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
		diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
		errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
		diagFile.close()
		errorFile.close()
	except:
		pass

	db.useOneConnection()
	sys.exit(status)
 
def init():
	'''
	# requires: 
	#
	# effects: 
	# 1. Processes command line options
	# 2. Initializes local DBMS parameters
	# 3. Initializes global file descriptors/file names
	# 4. Initializes global keys
	#
	# returns:
	#
	'''
 
	global inputFile, diagFile, errorFile, errorFileName, diagFileName, passwordFileName
	global exptFile, markerFile, exptMarkerFile, accFile, noteFile, sqlFile
	global exptFileName, markerFileName, exptMarkerFileName, accFileName, noteFileName, sqlFileName
	global mode, exptType, referenceKey
 
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:M:I:R:E:')
	except:
		showUsage()
 
	#
	# Set server, database, user, passwords depending on options
	# specified by user.
	#
 
	server = ''
	database = ''
	user = ''
	password = ''
	inputFileName = ''
	jnum = ''
 
	for opt in optlist:
                if opt[0] == '-S':
                        server = opt[1]
                elif opt[0] == '-D':
                        database = opt[1]
                elif opt[0] == '-U':
                        user = opt[1]
                elif opt[0] == '-P':
			passwordFileName = opt[1]
                elif opt[0] == '-M':
                        mode = opt[1]
                elif opt[0] == '-I':
                        inputFileName = opt[1]
                elif opt[0] == '-R':
                        jnum = opt[1]
                elif opt[0] == '-E':
                        exptType = regsub.gsub('"', '', opt[1])
                else:
                        showUsage()

	# User must specify Server, Database, User and Password
	password = string.strip(open(passwordFileName, 'r').readline())
	if server == '' or \
	   database == '' or \
	   user == '' or \
	   password == '' or \
	   mode == '' or \
	   inputFileName == '' or \
	   jnum == '' or \
	   exptType == '':
		showUsage()

	# Initialize db.py DBMS parameters
	db.set_sqlLogin(user, password, server, database)
	db.useOneConnection(1)
 
	fdate = mgi_utils.date('%m%d%Y')	# current date
	head, tail = os.path.split(inputFileName) 
	diagFileName = tail + '.' + fdate + '.diagnostics'
	errorFileName = tail + '.' + fdate + '.error'
	exptFileName = tail + '.' + fdate + '.MLD_Expts.bcp'
	markerFileName = tail + '.' + fdate + '.MLD_Marker.bcp'
	exptMarkerFileName = tail + '.' + fdate + '.MLD_Expt_Marker.bcp'
	accFileName = tail + '.' + fdate + '.ACC_Accession.bcp'
	noteFileName = tail + '.' + fdate + '.MLD_Notes.bcp'
	sqlFileName = tail + '.' + fdate + '.sql'

	try:
		inputFile = open(inputFileName, 'r')
	except:
		exit(1, 'Could not open file %s\n' % inputFileName)
		
	try:
		diagFile = open(diagFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % diagFileName)
		
	try:
		errorFile = open(errorFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % errorFileName)
		
	try:
		exptFile = open(exptFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % exptFileName)
		
	try:
		markerFile = open(markerFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % markerFileName)
		
	try:
		exptMarkerFile = open(exptMarkerFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % exptMarkerFileName)
		
	try:
		accFile = open(accFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % accFileName)
		
	try:
		noteFile = open(noteFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % noteFileName)
		
	try:
		sqlFile = open(sqlFileName, 'w')
	except:
		exit(1, 'Could not open file %s\n' % sqlFileName)
		
	# Log all SQL
	db.set_sqlLogFunction(db.sqlLogAll)

	# Set Log File Descriptor
	db.set_sqlLogFD(diagFile)

	diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
	diagFile.write('Server: %s\n' % (server))
	diagFile.write('Database: %s\n' % (database))
	diagFile.write('User: %s\n' % (user))
	diagFile.write('Input File: %s\n' % (inputFileName))

	errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

	referenceKey = verifyReference(jnum)

def verifyMode():
	'''
	# requires:
	#
	# effects:
	#	Verifies the processing mode is valid.  If it is not valid,
	#	the program is aborted.
	#	Sets globals based on processing mode.
	#	Deletes data based on processing mode.
	#
	# returns:
	#	nothing
	#
	'''

	global DEBUG

	if mode == 'preview':
		DEBUG = 1
	elif mode not in ['incremental', 'full']:
		exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

def verifyAssay(assay):
	'''
	# requires:
	#	assay - string, the Assay term
	#
	# effects:
	#	verifies that the Assay exists by checking the database
	#	writes to the error file if the Assay is invalid
	#	initializes global assayKey
	#
	# returns:
	#	Assay key if found
	#
	'''

	assayKey = 0

	results = db.sql('select _Assay_Type_key from MLD_Assay_Types where description = "%s"' % (assay), 'auto')
	for r in results:
		assayKey = r['_Assay_Type_key']

	if assayKey == 0:
		exit(1, 'Invalid Assay: %s\n' % (assay))

	return(assayKey)

def verifyReference(referenceID):
	'''
	# requires:
	#	referenceID - the Accession ID of the Reference (J:)
	#
	# effects:
	#	verifies that the Reference exists by checking the database
	#	writes to the error file if the Reference is invalid
	#	initializes global referenceKey
	#
	# returns:
	#	Reference key if found
	#
	'''

	referenceKey = accessionlib.get_Object_key('J:' + referenceID, 'Reference')

	if referenceKey is None:
		exit(1, 'Invalid Reference: %s\n' % (referenceID))

	return(referenceKey)

def verifyMarker(markerID, lineNum):
	'''
	# requires:
	#	markerID - the Accession ID of the Marker
	#	lineNum - the line number of the record from the input file
	#
	# effects:
	#	verifies that:
	#		the Marker exists either in the marker dictionary or the database
	#	writes to the error file if the Marker is invalid
	#	addes the marker id and key to the marker dictionary if the Marker is valid
	#
	# returns:
	#	0 and '' if the Marker is invalid
	#	Marker Key and Marker Symbol if the Marker is valid
	#
	'''

	global markerDict

	markerKey = None

	if markerDict.has_key(markerID):
		errorFile.write('Duplicate Mouse Marker (%d) %s\n' % (lineNum, markerID))
		return(0, '')
	else:
		results = db.sql('select _Marker_key, symbol ' + \
			'from MRK_Mouse_View where mgiID = "%s"' % (markerID), 'auto')
		for r in results:
			markerKey = r['_Marker_key']
			markerSymbol = r['symbol']

		if markerKey is None:
			errorFile.write('Invalid Mouse Marker (%d) %s\n' % (lineNum, markerID))
			markerKey = 0
			markerSymbol = ''
		else:
			markerDict[markerID] = `markerKey` + ':' + markerSymbol

	return(markerKey, markerSymbol)

def verifyChromosome(chromosome, lineNum):
	'''
	# requires:
	#	chromosome - the chromosome
	#	lineNum - the line number of the record from the input file
	#
	# effects:
	#	verifies that:
	#		the Chromosome is valid
	#	writes to the error file if the Chromosome is invalid
	#
	# returns:
	#	0 if the Chromosome is invalid
	#	1 if the Chromosome is valid
	#
	'''

	global chromosomeList

	if chromosome in chromosomeList:
		return 1
	else:
		errorFile.write('Invalid Chromosome (%d) %s\n' % (lineNum, chromosome))
		return 0

def loadDictionaries():
	'''
	# requires:
	#
	# effects:
	#	loads global dictionaries/lists: chromosomeList for lookup
	#
	# returns:
	#	nothing
	'''

	global chromosomeList

	results = db.sql('select chromosome from MRK_Chromosome where _Species_key = 1 ' + \
		'and chromosome not in ("UN", "MT") order by sequenceNum', 'auto')
	for r in results:
		chromosomeList.append(r['chromosome'])

def createExperiments():
	'''
	# requires:
	#
	# effects:
	#	creates bcp entries for:
	#		Master Experiment table
	#		Accession table
	#
	# returns:
	#	nothing
	#
	'''

	global exptDict, seqExptDict

	results = db.sql('select _Expt_key, chromosome from MLD_Expts where _Refs_key = %d' % (referenceKey), 'auto')

	# experiment records exists

	if len(results) > 0:
		for r in results:
			exptDict[r['chromosome']] = r['_Expt_key']
			seqExptDict[r['_Expt_key']] = 1

		if mode == 'full':
			# delete the existing *details*.....
			db.sql('delete from MLD_Marker where _Refs_key = %d' % (referenceKey), 'auto', execute = not DEBUG)
			db.sql('delete MLD_Expt_Marker from MLD_Expt_Marker m, MLD_Expts e ' + \
				' where e._Refs_key = %d and e._Expt_key = m._Expt_key ' % (referenceKey), \
				'auto', execute = not DEBUG)

	# if no experiment records exist....create them

	else:
		exptKey = 1000
		accKey = 1000
		mgiKey = 1000
		exptTag = 1

        	results = db.sql('select maxKey = max(_Expt_key) + 1 from MLD_Expts', 'auto')
        	if results[0]['maxKey'] is None:
                	exptKey = 1000
        	else:
                	exptKey = results[0]['maxKey']

        	results = db.sql('select maxKey = max(_Accession_key) + 1 from ACC_Accession', 'auto')
        	if results[0]['maxKey'] is None:
                	accKey = 1000
        	else:
                	accKey = results[0]['maxKey']

        	results = db.sql('select maxKey = maxNumericPart + 1 from ACC_AccessionMax ' + \
			'where prefixPart = "%s"' % (mgiPrefix), 'auto')
        	mgiKey = results[0]['maxKey']
	
		for c in chromosomeList:

			bcpWrite(exptFile, [exptKey, referenceKey, exptType, exptTag, c, cdate, cdate])
			bcpWrite(accFile, [accKey, \
				mgiPrefix + str(mgiKey), \
				mgiPrefix, \
				mgiKey, \
				logicalDBKey, \
				exptKey, \
				mgiTypeKey, \
				0, 1, \
				cdate, cdate, cdate])

			exptDict[c] = exptKey
			seqExptDict[exptKey] = 1
			exptKey = exptKey + 1
			exptTag = exptTag + 1
        		accKey = accKey + 1
        		mgiKey = mgiKey + 1
		
		# Update the AccessionMax value

		db.sql('exec ACC_setMax %d' % (exptTag), None, execute = not DEBUG)

def processFile():
	'''
	# requires:
	#
	# effects:
	#	Reads input file
	#	Verifies and Processes each line in the input file
	#
	# returns:
	#	nothing
	#
	'''

	lineNum = 0

	# sequence number of marker in master marker list
	seq1 = 1

	# For each line in the input file

	for line in inputFile.readlines():

		error = 0
		note = ''
		lineNum = lineNum + 1

		# Split the line into tokens
		tokens = string.split(line[:-1], '\t')

		try:
			markerID = tokens[0]
			chromosome = tokens[1]
			updateChr = tokens[2]
			band = tokens[3]
			assay = tokens[4]
		except:
			# if it's not a valid line, assume it's the note
			note = line
#			exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

		markerKey, markerSymbol = verifyMarker(markerID, lineNum)
		assayKey = verifyAssay(assay)
		error = not verifyChromosome(chromosome, lineNum)

		# determine experiment key for this chromosome

		if not exptDict.has_key(chromosome):
			errorFile.write('Cannot Find Experiment Key For Chromosome (%d): %s\n' % (lineNum, chromosome))
			chrExptKey = 0
		else:
			chrExptKey = exptDict[chromosome]

		if markerKey == 0 or assayKey == 0 or chrExptKey == 0:
			# set error flag to true
			error = 1

		# if errors, continue to next record
		if error:
			continue

		# if no errors, process

		# add marker to master marker file
		bcpWrite(markerFile, [referenceKey, markerKey, seq1, cdate, cdate])
		seq1 = seq1 + 1

		# add marker to experiment marker file
		bcpWrite(exptMarkerFile, \
			[chrExptKey, \
			markerKey, \
			alleleKey, \
			assayKey, \
			seqExptDict[chrExptKey], \
			markerSymbol, \
			description, \
			matrixData, \
			cdate, cdate])

		# increment marker sequence number for the experiment
		seqExptDict[chrExptKey] = seqExptDict[chrExptKey] + 1

		# update marker's chromosome...
		if updateChr == 'yes':
			sqlFile.write('update MRK_Marker ' + \
				'set chromosome = "%s" ' % (chromosome) + \
				'where _Marker_key = %d\ngo\n' % (markerKey))

			sqlFile.write('update MRK_Offset ' + \
				'set offset = -1.0 ' + \
				'where _Marker_key = %d\ngo\n' % (markerKey))

		# update cytogenetic band, if it is provided
		if band != "":
			sqlFile.write('update MRK_Marker ' + \
				'set cytogeneticOffset = "%s" ' % (band) + \
				'where _Marker_key = %d\ngo\n' % (markerKey))

#	end of "for line in inputFile.readlines():"

	if len(note) > 0:

		noteSeq = 1
		
		while len(note) > 255:
			bcpWrite(noteFile, [referenceKey, noteSeq, note[:255], cdate, cdate])
			newnote = note[255:]
			note = newnote
			noteSeq = noteSeq + 1

		if len(note) > 0:
			bcpWrite(noteFile, [referenceKey, noteSeq, note, cdate, cdate])

def bcpWrite(fp, values):
	'''
	#
	# requires:
	#	fp; file pointer of bcp file
	#	values; list of values
	#
	# effects:
	#	converts each value item to a string and writes out the values
	#	to the bcpFile using the appropriate delimiter
	#
	# returns:
	#	nothing
	#
	'''

	# convert all members of values to strings
	strvalues = []
	for v in values:
		strvalues.append(str(v))

	fp.write('%s\n' % (string.join(strvalues, bcpdelim)))

def bcpFiles():
	'''
	# requires:
	#
	# effects:
	#	BCPs the data into the database
	#
	# returns:
	#	nothing
	#
	'''

	exptFile.close()
	markerFile.close()
	exptMarkerFile.close()
	accFile.close()
	noteFile.close()
	sqlFile.close()

	cmd1 = 'cat %s | bcp %s..%s in %s -c -t\"%s" -S%s -U%s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MLD_Expts', exptFileName, bcpdelim, db.get_sqlServer(), db.get_sqlUser())

	cmd2 = 'cat %s | bcp %s..%s in %s -c -t\"%s" -S%s -U%s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MLD_Marker', markerFileName, bcpdelim, db.get_sqlServer(), db.get_sqlUser())

	cmd3 = 'cat %s | bcp %s..%s in %s -c -t\"%s" -S%s -U%s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MLD_Expt_Marker', exptMarkerFileName, bcpdelim, db.get_sqlServer(), db.get_sqlUser())

	cmd4 = 'cat %s | bcp %s..%s in %s -c -t\"%s" -S%s -U%s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'ACC_Accession', accFileName, bcpdelim, db.get_sqlServer(), db.get_sqlUser())

	cmd5 = 'cat %s | bcp %s..%s in %s -c -t\"%s" -S%s -U%s' \
		% (passwordFileName, db.get_sqlDatabase(), \
	   	'MLD_Notes', noteFileName, bcpdelim, db.get_sqlServer(), db.get_sqlUser())

	cmd6 = 'cat %s | isql -S%s -D%s -U%s -i%s' \
		% (passwordFileName, db.get_sqlServer(), db.get_sqlDatabase(), db.get_sqlUser(), sqlFileName)

	diagFile.write('%s\n' % cmd1)
	diagFile.write('%s\n' % cmd2)
	diagFile.write('%s\n' % cmd3)
	diagFile.write('%s\n' % cmd4)
	diagFile.write('%s\n' % cmd5)
	diagFile.write('%s\n' % cmd6)

	if DEBUG:
		return

	os.system(cmd1)
	os.system(cmd2)
	os.system(cmd3)
	os.system(cmd4)
	os.system(cmd5)
	os.system(cmd6)
#	db.sql('dump transaction %s with truncate_only' % (db.get_sqlDatabase()), None, execute = not DEBUG)

#
# Main
#

init()
verifyMode()
loadDictionaries()
createExperiments()
processFile()
bcpFiles()
exit(0)

