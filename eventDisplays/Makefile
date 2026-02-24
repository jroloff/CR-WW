
CC = g++
CFLAGS = -c -g -Wall `root-config --cflags` -std=c++17
LDFLATS =`root-config --glibs` 
GLIBS = $(shell root-config --glibs)

INCLUDE += `root-config --cflags` -D_HAS_AUTO_PTR_ETC=1


MYOBJS =
RootClass=myClass.h
RootDict=Dict.cpp
RootLinkDef=LinkDef.h
INC = -l

fileConversion : fileConversion.o $(MYOBJS) 
	$(CC) $(LDFLATS) fileConversion.o $(MYOBJS) -o fileConversion


fileConversion.o: fileConversion.cxx $<
	$(CC) $(INCLUDE) $(CFLAGS) fileConversion.cxx


$(RootDict): $(RootClass) $(RootLinkDef)
	rootcint -v4 -f $(RootDict) -c $(INC) $(RootClass) $(RootLinkDef)




clean:
	rm *o


