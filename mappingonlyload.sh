#!/bin/sh

#
# Wrapper script to create & load new mapping experiments
#
# Usage:  pp_mappingload.sh configFile
#

CONFIG_FILE=$1
. ${CONFIG_FILE}

cd ${MAPPINGDATADIR}
rm -rf ${MAPPINGONLYDATALOG}
touch ${MAPPINGONLYDATALOG}
rm -rf ${MAPPINGLOG}
touch ${MAPPINGLOG}

date >> ${MAPPINGONLYDATALOG}
${PYTHON} ${MAPPINGLOAD}/mappingonlyload.py >> ${MAPPINGONLYDATALOG}


${PYTHON} ${MAPPINGLOAD}/mappingload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${MAPPINGMODE} -I${MAPPINGDATAFILE} -E"${EXPERIMENTTYPE}" >> ${MAPPINGLOG}

date >>  ${MAPPINGONLYDATALOG}


