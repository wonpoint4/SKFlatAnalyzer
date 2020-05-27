# Plotter
General purpose plotter using ROOT histogram

[1.Basic classes](#basic-classes)  
[2.Quick start using Plotter](#quick-start-using-plotter)  
[3.Draw options](#draw-options)  
[4.Verbosity](#verbosity)  


## Basic Classes
### Style
### Sample
* Define a sample. It can be a file(`Sample::Type::FILE`), a sample, or a collection of samples(`Sample::Type::SUM||Sample::Type::STACK`)
* main variables
  * title: the title of this sample in plots
  * subs: a vector of sub-samples for Collection type or files for Sample type
  * style: a plotting style for this sample
  * sytle_alt: an alternative plotting style for this sample (for systematic uncertainty)
### Plot
* Define a plot. It is passed as argument of `Plotter::GetHist` function
### Plotter
* Base class for plotters
* Draw or save a TCanvas with `vector<Sample> Plotter::entries`
  
## Quick start using Plotter
```c++
root [0] .L Plotter.cc
root [1] Plotter aa
root [2] aa.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/EfficiencyValidation/2016/")
```
This will add ROOT files. You can check added files using `Plotter::PrintFiles` function.
```c++
root [3] aa.PrintFiles(true)
--------------------------
@Key: EfficiencyValidation_WZ_pythia 
  /data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/EfficiencyValidation/2016/EfficiencyValidation_WZ_pythia.root 2020-04-21 02:42:45
--------------------------
@Key: EfficiencyValidation_ZZ_pythia 
  /data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/EfficiencyValidation/2016/EfficiencyValidation_ZZ_pythia.root 2020-04-22 11:07:36
--------------------------
```
Next, setup entries.
```c++
root [4] aa.SetupEntries("EfficiencyValidation_ZZ_pythia EfficiencyValidation_WZ_pythia")
root [5] aa.PrintEntries(true)
--------------------------
+Title:EfficiencyValidation_ZZ_pythia Type:UNDEF  FC:0 LC:1 MC:1 DO:'e hist' 
  /data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/EfficiencyValidation/2016/EfficiencyValidation_ZZ_pythia.root 2020-04-22 11:07:36
--------------------------
+Title:EfficiencyValidation_WZ_pythia Type:UNDEF  FC:0 LC:2 MC:2 DO:'e hist' 
  /data6/Users/hsseo/SKFlatOutput/Run2Legacy_v4/EfficiencyValidation/2016/EfficiencyValidation_WZ_pythia.root 2020-04-21 02:42:45
--------------------------
```
For Draw or save a plot
```c++
root [6] aa.DrawPlot("ee2016/m80to100/dimass_MediumID_Q","norm")
root [7] aa.SavePlot("ee2016/m80to100/dimass_MediumID_Q","norm")
```
## Draw options
* norm
* noleg
* type:INT
* xmin,xmax
* widey, widewidey
  
## Verbosity
You can set verbosity using `Verbosity` variable (QUIET, ERROR, WARNING, INFO, DEBUG, ALL)
```c++
root [7] Verbosity=VERBOSITY::ERROR
```


