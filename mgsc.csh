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
setenv ASSAY		"assembly"
setenv EXPERIMENTTYPE	"TEXT-Physical Mapping"

setenv DBUTILITIESPATH		/usr/local/mgi/dbutils/mgidbutilities
setenv DBUSER			mgd_dbo
setenv DBPASSWORDFILE		${DBUTILITIESPATH}/.mgd_dbo_password

cd `dirname $0`
setenv LOG	$0.log

echo 'MGSC V3 Assembly Mapping Load' > $LOG
date >>& $LOG

set loaddir = `dirname $0`

# load the Annotation File
${loaddir}/mappingload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${INPUTFILE} -R${JNUM} -A\"${ASSAY}\" -E\"${EXPERIMENTTYPE}\" >>& $LOG

date >>& $LOG

