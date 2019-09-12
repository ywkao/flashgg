rm -rf Taggers/data
rm -rf Systematics/data
rm -rf MicroAOD/data

source MetaData/scripts/fggPopulateDataDir.sh -d Taggers/data
source MetaData/scripts/fggPopulateDataDir.sh -d Systematics/data
source MetaData/scripts/fggPopulateDataDir.sh -d MicroAOD/data
