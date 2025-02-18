TARGETS := chess-cover chess-cover-locking

default: $(TARGETS)
chess: default

include ../config.mk

CXXFLAGS=$(OPT) $(DBG) $(ARCH) $(CILKFLAG)
CFLAGS=$(OPT) $(DBG) $(ARCH) $(CILKFLAG)

chess-cover: chess-cover.o ktiming.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $@ chess-cover.o ktiming.o

chess-cover-locking: chess-cover-locking.o ktiming.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -pthread -o $@ chess-cover-locking.o ktiming.o

chess-cover-locking.o chess-cover.o ktiming.o: ktiming.h

clean:
	rm -f $(TARGETS) *.o
