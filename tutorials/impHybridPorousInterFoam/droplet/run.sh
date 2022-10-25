#!/bin/bash
set -e

# blockMesh
cd "${0%/*}" || exit 1

# Source tutorial run functions
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"

runApplication foamCleanTutorials

runApplication blockMesh

cp ./0/eps.org ./0/eps
cp ./0/alpha.wetting.org ./0/alpha.wetting

runApplication setFields

runApplication decomposePar

application=$(getApplication)
runParallel $application

runApplication reconstructPar
