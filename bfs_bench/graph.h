// -*- Mode:C++ -*-
/* graph.h by TB Schardl
   July 13, 2011
*/

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <sys/types.h>
#include <time.h>

#include <cilk/cilk.h>


#if CILKM_RTS
#include <reducer_lib/reducer_ptr.h>
#else 
#include <cilk/reducer_ptr.h>
#endif

#if PARALLEL_REDUCE
#include <reducer_lib/parallel_reducer_array.h>
#else 
#include <reducer_lib/reducer_array.h>
#endif

#include "reducer_lib/bag.h"
#include "reducer_min_uint.h" // user defined special reducer


#ifndef COMPUTE_INTENSE
#define COMPUTE_INTENSE 0
#endif 

#if COMPUTE_INTENSE
#define TIMES_TO_CALL_TRIVIAL 1000
#define WHEN_COMPUTE_INTENSE(ex) ex
extern "C" void trivial(void);
#else
#define WHEN_COMPUTE_INTENSE(ex)
#endif


using namespace cilk;

struct results {
    unsigned long long usec;
    unsigned long long usec_mem;
    unsigned long long usec_mutex;
    unsigned long long usec_mem_mutex;
    unsigned long long usec_prof;
    unsigned long long usec_s;
    unsigned long long usec_mem_s;
    unsigned long long usec_prof_s;
    int nodesAdded;
    int nodesAdded_s;
    int edgesTouched;
    int edgesTouched_s;
    int span;
};

class Graph {
private:
    typedef Bag< uint, scalable_allocator<void> > BagType;
    typedef Bag_reducer< uint, scalable_allocator<void> > BagRdcerType;
    typedef Pennant<uint, BagType::TAlloc> PennantType;
#if USE_REDUCER_ARRAY
    typedef reducer_array<typename reducer_min_uint::Monoid> MinRdcerArray;
#else 
    typedef reducer_ptr<typename reducer_min_uint::Monoid> MinRdcerPtr;
#endif

    uint nNodes; // Number of nodes
    uint nEdges; // Number of edges
    
    uint * nodes;
    uint * edges;
  
#ifdef COMPUTE_PARENTS
#if USE_REDUCER_ARRAY
    inline void pbfs_setup_helper( uint, uint, const uint, uint[], 
                                   MinRdcerArray&, 
                                   reducer_ptr<BagRdcerType::Monoid>& ) const; 
    void pbfs_walk_Bag( BagType *, reducer_ptr<BagRdcerType::Monoid>&, 
                        uint, uint[], MinRdcerArray& ) const;
    void pbfs_walk_Pennant( PennantType*, reducer_ptr<BagRdcerType::Monoid>&, 
                            uint, uint[], MinRdcerArray& ) const;
    void pbfs_proc_Node( const uint[], reducer_ptr<BagRdcerType::Monoid>&, 
                         uint, uint[],  MinRdcerArray&, uint ) const;
#else
    inline void pbfs_setup_helper( uint, uint, const uint, uint[], 
                                   MinRdcerPtr[], 
                                   reducer_ptr<BagRdcerType::Monoid>& ) const; 
    void pbfs_walk_Bag( BagType *, reducer_ptr<BagRdcerType::Monoid >&, 
                        uint, uint[], MinRdcerPtr[] ) const;
    void pbfs_walk_Pennant( PennantType*, reducer_ptr<BagRdcerType::Monoid>&, 
                            uint, uint[], MinRdcerPtr[] ) const;
    void pbfs_proc_Node( const uint[], reducer_ptr<BagRdcerType::Monoid>&, 
                         uint, uint[],  MinRdcerPtr[], uint ) const;
#endif // USE_REDUCER_ARRAY
#else
    void pbfs_walk_Bag( BagType *, BagRdcerType *, uint, uint[] ) const;
    void pbfs_walk_Pennant( PennantType*, BagRdcerType *, uint, uint[] ) const;
    void pbfs_proc_Node(const uint[], BagRdcerType *, uint, uint[], uint) const;
#endif // COMPUTE_PARENTS
    unsigned long long todval(struct timeval*) const;

public:
    // Constructor/Destructor
    Graph(const int nNodes, const int nEdges, std::vector<int>* edges[]);
    Graph(int *ir, int *jc, int m, int n, int nnz);
    ~Graph();

    // Various BFS versions
#ifdef COMPUTE_PARENTS
    void bfs(const uint s, uint distances[], uint parents[]) const;
    void pbfs(const uint s, uint distances[], uint parents[]) const;
#else
    int bfs(const uint s, uint distances[]) const;
    int pbfs(const uint s, uint distances[]) const;
#endif
};

#include "graph.cpp"
#endif

