#!/bin/csh -f -x

#
# MGSC V3 Assemply Mapping Load (TR 4010)
#
# Usage:
# 	mgsc.csh DBSERVER DBNAME inputfile mode
#
# Example:
#	mgsc.csh MGD mgd UNmgiID_EnsemblChr_Band incremental
#
# Purpose:
#	executes mappingload.py to update chromosome/bands for markers
#	and to create physical mapping experiment records
#

setenv DBSERVER		$1
setenv DBNAME		$2
setenv INPUTFILE	$3
setenv MODE		$4

setenv JNUM		78475
setenv EXPERIMENTTYPE	"Text-Physical Mapping"

setenv DBUTILITIESPATH		/usr/local/mgi/dbutils/mgidbutilities
setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password
setenv MAPPINGLOAD		/usr/local/mgi/dataload/mappingload

cd `dirname $0`
setenv LOG	$0.log

echo 'MGSC V3 Assembly Mapping Load' > $LOG
date >>& $LOG

# load the Annotation File
${MAPPINGLOAD}/mappingload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${INPUTFILE} -R${JNUM} -E\"${EXPERIMENTTYPE}\" >>& $LOG

date >>& $LOG

