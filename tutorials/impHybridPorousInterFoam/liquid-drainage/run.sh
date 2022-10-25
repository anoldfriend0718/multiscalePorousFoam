#!/bin/bash
set -e

# blockMesh
cd "${0%/*}" || exit 1

# Source tutorial run functions
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"


runApplication foamCleanTutorials

runApplication blockMesh

cp ./0.org/eps ./0/eps
cp ./0.org/alpha.wetting ./0/alpha.wetting

runApplication setFields

runApplication renumberMesh -overwrite

runApplication decomposePar

application=$(getApplication)
runParallel $application

runApplication reconstructPar

