# Automatically finds all .cxx files in a directory and compiles it with correct flags
CXX = g++
CXXFLAGS = -g -Wall -fPIC
ROOTFLAGS = `root-config --cflags --glibs`
INCLUDE += -I $(PLOTUTILSROOT)/../include
#INCLUDE += -I $(UNFOLDUTILSROOT) -I $(PLOTUTILSROOT)/ -I $(CCQENU_ANA)/ -I$(CCQENU_INCLUDE_PATH) -I$(MCREWEIGHTROOT)
LDLIBS += -L$(PLOTUTILSROOT) -lMAT
#LDLIBS += -L$(UNFOLDUTILSROOT)/$(CMTCONFIG) -lUnfoldUtils -L$(PLOTUTILSROOT)/$(CMTCONFIG) -lplotutils -L$(MCREWEIGHTROOT)/$(CMTCONFIG) -lMCReweight -L$(CCQENU_ANA)/lib -lCCQENU

TARGETS   = $(wildcard *.o)
BINARIES  = $(TARGETS:.o=)
all : $(patsubst %.cxx, %.o, $(wildcard *.cxx))

#--- if making individual targets ---#
%.o : %.cxx $(DEPS)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $(LDLIBS) -c $< #compile
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $@        #link 
.PHONY: clean

clean:
	rm -f $(BINARIES) $(TARGETS)
