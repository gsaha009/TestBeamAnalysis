GCC=g++
CXXFLAGS=`root-config --libs --cflags` -g -std=c++17 -fPIC -Wall -I./
SOFLAGS=-shared

SRCDIR=src
BINDIR=bin

SRC=src/BTAnalyzer.cc src/main.cc src/ModuleHistos.cc
OBJ=$(patsubst $(SRCDIR)/%.cc, $(BINDIR)/%.o , $(SRC)  )
PYOBJ=$(filter-out $(BINDIR)/main.o,$(OBJ))

.PHONY: all
all:
	$(info, "--- Full compilation --- ")	
	$(MAKE) libBTAnalysis.so
	$(MAKE) baseline.exe

libBTAnalysis.so: $(PYOBJ) Dict | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(RPATH) $(SOFLAGS) -o $(BINDIR)/$@ $(PYOBJ) $(BINDIR)/dict.o

baseline.exe:  $(BINDIR)/main.o |  $(BINDIR)
	$(GCC) $(CXXFLAGS) $(RPATH) -o $(BINDIR)/$@  $(OBJ) `root-config --libs`

$(OBJ) : $(BINDIR)/%.o : $(SRCDIR)/%.cc | $(BINDIR)
	$(GCC) $(CXXFLAGS) $(RPATH) -c -o $(BINDIR)/$*.o $<

.PHONY: Dict
Dict: $(BINDIR)/dict.o

$(BINDIR)/dict.o: $(SRC) | $(BINDIR)
	genreflex $(SRCDIR)/classes.h -s $(SRCDIR)/classes_def.xml -o $(BINDIR)/dict.cc --fail_on_warnings --rootmap=$(BINDIR)/dict.rootmap --rootmap-lib=libBTAnalysis.so -I./$(SRCDIR)/
	$(GCC) -c -o $(BINDIR)/dict.o $(CXXFLAGS) $(RPATH) -I./$(SRCDIR)/ $(BINDIR)/dict.cc

$(BINDIR):
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	-rm $(OBJ)
	-rm $(BINDIR)/dict*
	-rm $(BINDIR)/libBTAnalysis.so
	-rm $(BINDIR)/baseline.exe
	-rmdir $(BINDIR)
