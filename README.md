# SKFlatAnalyzer

## Manual

https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf

## Where to put the analyzer
TAMSA 1/2 : /data6/Users/$USER/

KISTI : ~/ # home directory

KNU :  ~/scartch/

## First time setup
```
source setup.sh # This should be done for every new shell
#### Setup user info
#### Need to be done only once
cp $SKFlat_WD/python/UserInfo_template.py $SKFlat_WD/python/UserInfo_${USER}.py 
# edit $SKFlat_WD/python/UserInfo_${USER}.py
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

# How To update your forked repo to follow the original repository

To get your forked repo from this one
```
push "Fork" button on the right-top side of this web page
```

Go to you forked repo webpage and copy the address of it, then clone

How to update your forked repo to follow the original repository
```
git remote add upstream https://github.com/CMSSNU/SKFlatAnalyzer
git fetch upstream
git pull upstream master
```
