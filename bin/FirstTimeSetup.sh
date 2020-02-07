#!/bin/bash

#### [GEScaleSyst]
#### Adding Minseok's remote
echo "@@@@ external/GEScaleSyst/"
cd external/GEScaleSyst/
echo "@@@@ Trying to add remote : git remote add min git@github.com:khaosmos93/GEScaleSyst.git"
git remote add min git@github.com:khaosmos93/GEScaleSyst.git
echo "@@@@ going back.. : cd ../../"
cd ../../
