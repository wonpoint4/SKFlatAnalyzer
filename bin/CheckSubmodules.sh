#!/bin/bash
submodules=$(git submodule --quiet foreach --recursive 'echo $name')
for submodule in $submodules
do
  echo "@@@@ Checking submodule at \""$submodule"\""
  cd $submodule
  echo "@@@@   Fetching.."
  git fetch --all -q
  logDiff=$(git log HEAD..remotes/min/master)
  if [ ! -z "$logDiff" ]
  then
    echo "@@@@   Printing difference"
    git log HEAD..remotes/min/master
  else
    echo "@@@@   --> This module is up-to-date"
  fi
  cd $SKFlat_WD
done
