# Makefile copied from organism.

# This Makefile should automatically keeps track of dependencies in the
# code. It saves dependencies in *.d files and includes them whenever
# necessary. The dependency part is based on Scott McPeak's Makefile at
# http://www.cs.berkeley.edu/~smcpeak/autodepend/autodepend.html. The
# other parts are based on Simon's Makefile and tons of obscure
# webpages I no longer recall.  
#
#'make' build executable files
# 	'bin/simulator' and 'bin/optimizer' 
#
#'make debug'	Compiles with -g and no optimization 
#
#'make gprof'	Compiles with -pg and no optimization 
#
#'make clean' 	removes all .o, executables and .d 

SHELL = /bin/sh
# Compilers and flags
CXX = g++


# Make flags OS independent (note, not yet for Windows or debug or gprof options)

UNAME_S := $(shell uname -s)	# os
UNAME_A := $(shell uname -p)	# architecture

CXXFLAGS = -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
LDFLAGS = -O3 -DNDEBUG -Wall -pedantic -std=c++11

# Linux (assume default)
ifeq ($(UNAME_S),Linux)
	# CXXFLAGS = -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
	# LDFLAGS = -O3 -DNDEBUG -Wall -pedantic -std=c++11
endif

# Cygwin 10.0
ifeq ($(UNAME_S),CYGWIN_NT-10.0)
	# CXXFLAGS = -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
	# LDFLAGS = -O3 -DNDEBUG -Wall -pedantic -std=c++11
endif

# Mac OS X
ifeq ($(UNAME_S),Darwin)
	# CXXFLAGS = -O3 -DNDEBUG -Wall -pedantic -ansi -std=c++11
	# LDFLAGS = -O3 -DNDEBUG -Wall -pedantic -std=c++11

	CXXFLAGS += -stdlib=libc++ -I/sw/include/ -I/opt/local/include -I/opt/homebrew/include
	LDFLAGS += -L/sw/lib/ -L/opt/local/lib/ -L/opt/homebrew/lib/
endif

# if on an M1 mac, use clang to fix linker bug in gcc
ifeq ($(UNAME_A), arm)
	CXX = clang
endif

# Adding anaconda to library paths (for conda installed boost)
CXXFLAGS += -I${CONDA_PREFIX}/include
LDFLAGS  += -L${CONDA_PREFIX}/lib

# add HomeBrew to library paths (for homebrew installed boost (may need to change))
CXXFLAGS += -I/opt/homebrew/include
LDFLAGS  += -L/opt/homebrew/lib

#Sourcefiles etc.
SIM_SRC = simulator/simulator.cc
SIM_OBJ = $(SIM_SRC:.cc=.o)
CONV_SRC = tools/converter.cc
CONV_OBJ = $(CONV_SRC:.cc=.o)
INITMANIP_SRC = tools/manipulateInit.cc
INITMANIP_OBJ = $(INITMANIP_SRC:.cc=.o)
SRC_DIR = ./ ./ply/
SRCS := $(wildcard $(SRC_DIR:/=/*.cc))
OBJS = $(SRCS:.cc=.o)
DEP = makedep.rule

#Binaries
SIMULATOR = ../bin/simulator
CONVERTER = ../bin/converter
INITMANIP = ../bin/manipulateInit

all: $(SIMULATOR) $(CONVERTER) $(INITMANIP)

$(SIMULATOR): $(OBJS) $(SIM_OBJ)
	$(CXX) $(SIM_OBJ) $(OBJS) $(LDFLAGS) -o $(SIMULATOR) 

$(CONVERTER): $(OBJS) $(CONV_OBJ)
	$(CXX) $(CONV_OBJ) $(OBJS) $(LDFLAGS) -o $(CONVERTER) 

$(INITMANIP): $(OBJS) $(INITMANIP_OBJ)
	$(CXX) $(INITMANIP_OBJ) $(OBJS) $(LDFLAGS) -o $(INITMANIP) 

$(SIM_OBJ) : $(SIM_SRC)

$(CONV_OBJ) : $(CONV_SRC)

$(INITMANIP_OBJ) : $(INITMANIP_SRC)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

# Compile and generate dependency info.  The two first lines fixes a
# bug in gcc(?), and prints dir/foo.o: dir/foo.cc... in the dependency
# file instead of foo.o: dir/foo.cc... as is done by gcc. The other
# two lines use a trick to get make to detect name changes etc.
%.o: %.cc
	$(CXX) -o $*.o -c $(CXXFLAGS) $*.cc
	$(CXX) -MM $(CXXFLAGS) $*.cc > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

debug:
	make "CXXFLAGS = -g -pedantic -Wall -ansi" "LDFLAGS = -g -pedantic -Wall -ansi" 		
gprof:
	make "CXXFLAGS = -pg -pedantic -Wall -ansi" "LDFLAGS = -pg -pedantic -Wall -ansi" 		


.PHONY: clean

clean:
	rm -f $(OBJS)
	rm -f $(SIM_OBJ)
	rm -f $(CONV_OBJ)
	rm -f $(INITMANIP_OBJ)
	rm -f $(SIMULATOR)
	rm -f $(CONVERTER)
	rm -f $(INITMANIP)
	rm -f *.d
	rm -f */*.d
