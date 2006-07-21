#!/bin/csh -f

#
# Wrapper script to create & load new mapping experiments
#
# Usage:  mappingload.csh configFile
#

setenv CONFIGFILE $1

source ${CONFIGFILE}

setenv MAPPINGLOAD	${DATALOAD}/mappingload/mappingload.py

cd ${MAPPINGDATADIR}

rm -rf ${MAPPINGLOG}
touch ${MAPPINGLOG}

date >> ${MAPPINGLOG}

${MAPPINGLOAD} -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${MAPPINGMODE} -I${MAPPINGDATAFILE} -R${JNUM} -E"${EXPERIMENTTYPE}" -C${CREATEDBY} >>& ${MAPPINGLOG}

date >> ${MAPPINGLOG}

