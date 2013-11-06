# rename is dos program
# dynarray.* stack.* are not used

include Makefile.common

UTILS_DIR = ./

all: brkpnts

brkpnts: $(UTILS_DIR)brkpnts.c 
	$(LINK) -o $(UTILS_DIR)brkpnts -DUNIX $(UTILS_DIR)brkpnts.c -lm

.cpp.o:
	$(CRULES)

clean:
	$(RMCOM) *.o
	$(RMCOM) brkpnts

