LIBRARY = Analyzers
OBJDIR = obj
DEPDIR = $(OBJDIR)/dep
SRCDIR = src
INCDIR = include
#ADDITIONAL_ROOTMAPLIBRARY = -rml libAnalyzerTools.so -rml libGEScaleSyst.so -rml libDataFormats.so -rml libTH4D.so -rml libTree.so -rml libKinematicFitter.so
ADDITIONAL_ROOTMAPLIBRARY = -rml libAnalyzerTools.so -rml libGEScaleSyst.so -rml libDataFormats.so -rml libTH4D.so -rml libTree.so

include $(SKFlat_WD)/makefile.common

INCLUDES += -I$(SKFlat_WD)/DataFormats/include/
INCLUDES += -I$(SKFlat_WD)/AnalyzerTools/include/
INCLUDES += -I$(SKFlat_WD)/external/GEScaleSyst/
INCLUDES += -I$(SKFlat_WD)/external/TH4D/
INCLUDES += -I$(SKFlat_WD)/external/KinematicFitter/include/
INCLUDES += -I/cvmfs/cms.cern.ch/$(SCRAM_ARCH)/cms/$(cmsswrel)/src/
#INCLUDES += -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1-fmblme/include/
