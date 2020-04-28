cilkm_dir = /scratch/angelee/cilk_dev/cilk-M/src_for_IntelCC

AR = ar
RTSDEF = -DCILKM_RTS=0 -DCILKVIEW=0 -DTIMES_TO_RUN=10
LIBS = -ltbbmalloc -ltbb
# CILKUTIL_INC = /home/angelee/cilkutil/include/cilktools
INC = -I$(cilkm_dir)/cilkplus_tests/hyperobjs \
      -I$(cilkm_dir)/cilkplus_tests/hyperobjs/reducer_lib # $(CILKUTIL_INC)

COMPUTE_INTENSE = -DCOMPUTE_INTENSE
COMPUTE_PARENTS = -DCOMPUTE_PARENTS
REDUCER_ARRAY = -DUSE_REDUCER_ARRAY
PARALLEL_REDUCE = -DPARALLEL_REDUCE

CILK=icpc
CILKFLAGS= -O2 -ggdb
#CILKUTIL_LD = -L/scratch/angelee/cilkutil/lib64
LDFLAGS= -L$(CURDIR)

OTHER= Makefile

BFS_SRC= bfs.cpp
GRAPH_SRC= graph.cpp graph.h reducer_min_uint.h
BAG_SRC= $(cilkm_dir)/cilkplus_tests/hyperobjs/reducer_lib/bag.cpp \
         $(cilkm_dir)/cilkplus_tests/hyperobjs/reducer_lib/bag.h

TARGETS = bfs.out bfs-sp.out bfs-pp.out bfs-p-aor.out \
          bfs-c.out bfs-c-sp.out bfs-c-pp.out bfs-c-p-aor.out

.PHONY: clean cilkm cilkplus

all: $(TARGETS)

bfs.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $< -o $@

bfs-sp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(PARALLEL_REDUCE)=0 $(REDUCER_ARRAY)=1 $< -o $@

bfs-pp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(PARALLEL_REDUCE)=1 $(REDUCER_ARRAY)=1 $< -o $@

bfs-p-aor.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS) $(REDUCER_ARRAY)=0 $< -o $@

bfs-c.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_INTENSE)=1 $< -o $@

bfs-c-sp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(PARALLEL_REDUCE)=0 $(REDUCER_ARRAY)=1 $< -o $@

bfs-c-pp.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(PARALLEL_REDUCE)=1 $(REDUCER_ARRAY)=1 $< -o $@

bfs-c-p-aor.out : $(BFS_SRC) $(GRAPH_SRC) $(BAG_SRC) $(OTHER)
	$(CILK) $(CILKFLAGS) $(INC) $(LDFLAGS) $(LIBS) $(RTSDEF) $(COMPUTE_PARENTS)  $(COMPUTE_INTENSE)=1 $(REDUCER_ARRAY)=0 $< -o $@

.o:
	$(CC) $(CFLAGS) $(LDFLAGS) $(LIBS) -c $<

clean :
	rm -f $(TARGETS)

cilkm :
	$(cilkm_dir)/setup_icc.sh cilkm
	make clean
	make RTSDEF=-DCILKM_RTS=1

cilkplus :
	$(cilkm_dir)/setup_icc.sh cilkplus
	make clean
	make RTSDEF=-DCILKM_RTS=0
