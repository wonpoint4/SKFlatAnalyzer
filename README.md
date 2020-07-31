# SKFlatAnalyzer

## Manual

https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf

## Where to put the analyzer
TAMSA 1/2 : /data6/Users/$USER/

KISTI : ~/ # home directory

KNU :  ~/scartch/

## First time setup
```
#### When first time gie clone, use the option "--recursive" to initiate the submodules
git clone --recursive git@github.com:CMSSNU/SKFlatAnalyzer.git
cd SKFlatAnalyzer
#### add your remote repo
git remote add <nickname> git@github.com:<gitaccount>/SKFlatAnalyzer.git
git checkout <your working branch>

#### First time setup script
source bin/FirstTimeSetup.sh 
source setup.sh

#### You have to edit user info
#### First, copy the temply using the command below
cp $SKFlat_WD/python/UserInfo_template.py $SKFlat_WD/python/UserInfo_${USER}.py 
#### Then, edit $SKFlat_WD/python/UserInfo_${USER}.py

#### Now, run setup script.
#### This should be done for every new shell
source setup.sh
```
Compile
```
make clean
make
```

## Test job
```bash
SKFlat.py -a ExampleRun -i DYJets -n 50 -y 2016 &
```

## Making a new Ananlyzer
```bash
cd script/MakeCycleSkeleton/
```
Then, run
```bash
python MakeCycleSkeleton.py NewAnalyzer # <- put new analyzer name
```
It will print below lines (execute the lines) :
```bash
mv NewAnalyzer.h $SKFlat_WD/Analyzers/include/
mv NewAnalyzer.C $SKFlat_WD/Analyzers/src/
```

Then, add
```bash
#pragma link C++ class NewAnalyzer+;
```
in Analyzers/include/Analyzers_LinkDef.h

## Detailed descriptions

Look Analyzers/src/ExampleRun.C

# Tips

## Making PR

Start from the CMSSNU's master branch of CMSSNU when making pull request to prevent anoying conflicts.

