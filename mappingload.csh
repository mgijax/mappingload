#!/bin/csh -f

#
# Wrapper script to create & load new mapping experiments
#
# Usage:  mappingload.csh
#

cd `dirname $0`

# DB schema directory; its Configuration file will set up all you need
setenv SCHEMADIR $1
source ${SCHEMADIR}/Configuration

setenv INPUTFILE	specifictoyourload
setenv MODE		preview|incremental|full
setenv JNUM		J:
setenv EXPERIMENTTYPE	"TEXT-Physical Mapping"
setenv CREATEDBY	specifictoyourload

setenv MAPPINGLOAD		/usr/local/mgi/dataload/mappingload

cd `dirname $0`
setenv LOG	$0.log

date >>& $LOG

${MAPPINGLOAD}/mappingload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${MODE} -I${INPUTFILE} -R${JNUM} -E"${EXPERIMENTTYPE}" -C${CREATEDBY} >>& $LOG

date >>& $LOG

