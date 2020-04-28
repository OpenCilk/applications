TARGETS := chess-cover chess-cover-locking

default: $(TARGETS)
chess: default

include ../common.mk

NOINLINE=1
ifeq ($(NOINLINE),1)
	CFLAGS += -Wfatal-errors -fcilk-no-inline # -flto
	CXXFLAGS += -Wfatal-errors -fcilk-no-inline # -flto
	LDFLAGS += # -flto
endif

chess-cover: chess-cover.o ktiming.o
chess-cover-locking: chess-cover-locking.o ktiming.o

chess-cover-locking: chess-cover-locking.o $(CILKRTS_SLIB) $(PORR_SLIB)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

%: %.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f $(TARGETS) *.o
