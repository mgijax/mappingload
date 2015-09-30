#!/bin/sh

#
# Wrapper script to create & load new mapping experiments
#
# Usage:  mappingload.csh configFile
#

CONFIG_FILE=$1
. ${CONFIG_FILE}

cd ${MAPPINGDATADIR}

rm -rf ${MAPPINGLOG}
touch ${MAPPINGLOG}

date >> ${MAPPINGLOG}

${MAPPINGLOAD}/mappingload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${MAPPINGMODE} -I${MAPPINGDATAFILE} -E"${EXPERIMENTTYPE}" >>& ${MAPPINGLOG}

date >> ${MAPPINGLOG}

