include ../config.mk

AR = ar
RTSDEF = -DCILKVIEW=0 -DTIMES_TO_RUN=10
LIBS = 

COMPUTE_INTENSE = -DCOMPUTE_INTENSE
COMPUTE_PARENTS = -DCOMPUTE_PARENTS
REDUCER_ARRAY = -DUSE_REDUCER_ARRAY
PARALLEL_REDUCE = -DPARALLEL_REDUCE

CILKFLAGS= $(OPT) $(ARCH) $(DBG) $(CILKFLAG)

OTHER= Makefile

BFS_SRC= bfs.cpp
GRAPH_SRC= graph.cpp graph.h reducer_min_uint.h
BAG_SRC= bag.cpp
BAG_INC= bag.h

TARGETS = bfs.out bfs-sp.out bfs-pp.out bfs-p-aor.out \
          bfs-c.out bfs-c-sp.out bfs-c-pp.out bfs-c-p-aor.out

.PHONY: clean cilkm cilkplus

all: $(TARGETS)

bfs.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(BAG_INC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) -o $@

bag.o : bag.cpp bag.h
	$(CXX) $(CILKFLAGS) -c bag.cpp -o $@

bfs-sp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(PARALLEL_REDUCE)=0 $(REDUCER_ARRAY)=1 $< -o $@

bfs-pp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(PARALLEL_REDUCE)=1 $(REDUCER_ARRAY)=1 $< -o $@

bfs-p-aor.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(REDUCER_ARRAY)=0 $< -o $@

bfs-c.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_INTENSE)=1 $< -o $@

bfs-c-sp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(PARALLEL_REDUCE)=0 $(REDUCER_ARRAY)=1 $< -o $@

bfs-c-pp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(PARALLEL_REDUCE)=1 $(REDUCER_ARRAY)=1 $< -o $@

bfs-c-p-aor.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CXX) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(REDUCER_ARRAY)=0 $< -o $@

.o:
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) -c $<

clean :
	rm -f $(TARGETS) *.o
