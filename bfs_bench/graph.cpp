// Graph Representation
#ifndef GRAPH_CPP
#define GRAPH_CPP

#include <assert.h>
#include <vector>
#include <cstdlib>
#include <sys/types.h>

#include "graph.h"

#define GraphDebug 0

#ifndef GRAINSIZE
#define GRAINSIZE BLK_SIZE 
#endif


//static inline uint
//MIN(uint a, uint b) { return b ^ ( (a ^ b) & -(uint)(a < b) ); }
#define MIN(a, b) \
  ({ typeof(a) _a = (a); \
    typeof(b) _b = (b); \
    _b ^ ( (_a ^ _b) & -(typeof(a))(_a < _b) ); })

unsigned long long
Graph::todval(struct timeval *tp) const {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}

// Constructor
Graph::Graph(const int nNodes, const int nEdges, std::vector<int>* edges[]) {

    uint lastEdge = 0;
    this->nNodes = nNodes;
    this->nEdges = nEdges;

    this->nodes = (uint *) scalable_malloc( sizeof(uint) * (nNodes+1) );
    this->edges = (uint *) scalable_malloc( sizeof(uint) * (nEdges*2) );

    for(uint i = 0; i < nNodes; ++i) {
        this->nodes[i] = lastEdge;
        for(uint j = 0; j < edges[i]->size(); ++j, ++lastEdge)
            this->edges[lastEdge] = edges[i]->at(j);
    }
    this->nodes[nNodes] = lastEdge;

    if(GraphDebug) {
        assert(lastEdge == 2*nEdges);
        for(uint i = 0; i < nNodes; i++)
            printf("%u -> %d, ", i, this->nodes[i]);
        printf("\n");

        for(uint i = 0; i < 2*nEdges; i++)
            printf("%d ", this->edges[i]);
        printf("\n");
    }
}

// Constructor
Graph::Graph(int *ir, int *jc, int m, int n, int nnz) {

    this->nNodes = m;
    this->nEdges = nnz;

    this->nodes = (uint *) scalable_malloc( sizeof(uint) * (m+1) );
    this->edges = (uint *) scalable_malloc( sizeof(uint) * nnz );

    int *w = (int *) scalable_malloc( sizeof(int)*m );

    for(uint i = 0; i < m; ++i) {
        w[i] = 0;
    }

    for(uint i = 0; i < jc[n]; ++i)
        w[ir[i]]++;

    int prev;
    int tempnz = 0;
    for(uint i = 0; i < m; ++i) {
        prev = w[i];
        w[i] = tempnz;
        tempnz += prev;
    }
    this->nodes[m] = tempnz;
    for(uint i = 0; i < m; ++i)
        this->nodes[i] = w[i];

    for(uint i = 0; i < n; ++i) {
        for(uint j = jc[i]; j < jc[i+1]; j++)
            this->edges[w[ir[j]]++] = i;
    }

    scalable_free(w);
}

// Destructor
Graph::~Graph() {
    scalable_free(this->nodes);
    scalable_free(this->edges);
}

// Serial BFS
#ifdef COMPUTE_PARENTS
void Graph::bfs(const uint s, uint distances[], uint parents[]) const {

    uint newdist;
    uint head, tail, current;
    uint *queue = (uint *) scalable_malloc( sizeof(uint) * nNodes );

    if(s > nNodes) {
        fprintf(stderr, "Starting node is greater than number of nodes.\n");
        return;
    }

    for(uint i = 0; i < nNodes; ++i) {
        distances[i] = UINT_MAX;
        parents[i] = UINT_MAX;
    }

    current = s;
    distances[s] = 0;
    head = 0;
    tail = 0;

    do {
        newdist = distances[current] + 1;

        for(uint i = nodes[current]; i < nodes[current+1]; i++) {
            if(newdist < distances[edges[i]]) {
                queue[tail++] = edges[i];
                distances[edges[i]] = newdist;
                parents[edges[i]] = current;

            } else if(newdist == distances[edges[i]]) {
                parents[edges[i]] = MIN(current, parents[edges[i]]);
            }
        }
        current = queue[head++];
    } while(head <= tail);

    scalable_free(queue);
}

#else 
int Graph::bfs(const uint s, uint distances[]) const {

    uint *queue = new uint[nNodes];
    uint head, tail, current;
    uint newdist;

    if(s > nNodes)
        return -1;

    for(uint i = 0; i < nNodes; ++i) {
        distances[i] = UINT_MAX;
    }

    current = s;
    distances[s] = 0;
    head = 0;
    tail = 0;

    do {
        newdist = distances[current] + 1;

        for(uint i = nodes[current]; i < nodes[current+1]; i++) {
            if(newdist < distances[edges[i]]) {
                queue[tail++] = edges[i];
                distances[edges[i]] = newdist;
            }
        }
        current = queue[head++];
    } while (head <= tail);

    delete[] queue;

    return 0;
}
#endif

// Helper routine to walk a bag for PBFS
#ifdef COMPUTE_PARENTS
#if USE_REDUCER_ARRAY
inline void
Graph::pbfs_walk_Bag( BagType *b, reducer_ptr<BagRdcerType::Monoid> &next, 
                      uint newdist, uint distances[], 
                      MinRdcerArray &parents ) const {

    PennantType *p = NULL;

    if(b->getFill() > 0) { 
        // Split the bag and recurse
        b->split(&p); 
        cilk_spawn pbfs_walk_Bag(b, next, newdist, distances, parents);
        pbfs_walk_Pennant(p, next, newdist, distances, parents);
    } else {
        pbfs_proc_Node( b->getFilling(), next, newdist, distances, 
                        parents, b->getFillingSize() );
    }
    cilk_sync;
}

// Helper routine to walk a pennant for PBFS
inline void
Graph::pbfs_walk_Pennant( PennantType *p, 
                          reducer_ptr<BagRdcerType::Monoid> &next, 
                          uint newdist, uint distances[], 
                          MinRdcerArray &parents ) const {
    if(p->getLeft() != NULL) {
        cilk_spawn pbfs_walk_Pennant(
                        p->getLeft(), next, newdist, distances, parents);
        if(p->getRight() != NULL)
            cilk_spawn pbfs_walk_Pennant(
                        p->getRight(), next, newdist, distances, parents);
    }
    // Process the current element
    pbfs_proc_Node( p->getElements(), next, newdist, distances, 
                    parents, BLK_SIZE );
    cilk_sync;

    BagType::PAlloc().deallocate(p, 1);
}

// Helper routine to process a node for PBFS
inline void
Graph::pbfs_proc_Node( const uint n[], reducer_ptr<BagRdcerType::Monoid> &next, 
                       uint newdist, uint distances[],  
                       MinRdcerArray &parents, uint fillSize ) const {

    uint* localFilling = BagType::TAlloc().allocate(BLK_SIZE);
    uint localFillingSize = 0;
    BagType *bnext = &(*next);
    // ANGE: has to expose the SPA view instead of the underlying array;
    // otherwise we will miss the logging 
    MinRdcerArray::View &p = parents.view_ref();

    for(uint j = 0; j < fillSize; j++) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        for(uint i = nodes[n[j]]; i < nodes[n[j]+1]; ++i) {

            WHEN_COMPUTE_INTENSE({ 
                for(uint k = 0; k < TIMES_TO_CALL_TRIVIAL; k++) 
                    trivial();
            })

            if(newdist < distances[edges[i]]) {
                distances[edges[i]] = newdist;
                // parents[edges[i]].set(n[j]);  
                p[edges[i]] = n[j];

                localFilling[localFillingSize++] = edges[i];
                if(localFillingSize == BLK_SIZE) {
                    bnext->insert_fblk(localFilling);
                    localFilling = BagType::TAlloc().allocate(BLK_SIZE);
                    localFillingSize = 0;
                }
            } else if (newdist == distances[edges[i]]) {
                //parents[edges[i]].take_min(n[j]);
                if(p[edges[i]] > n[j])
                    p[edges[i]] = n[j];
            }
        }
    }
    bnext->insert_blk(localFilling, localFillingSize);
}

#else // !USE_REDUCER_ARRAY
inline void
Graph::pbfs_walk_Bag( BagType *b, reducer_ptr<BagRdcerType::Monoid> &next, 
                      uint newdist, uint distances[], 
                      MinRdcerPtr parents[] ) const {

    PennantType *p = NULL;

    if(b->getFill() > 0) { 
        // Split the bag and recurse
        b->split(&p); 
        cilk_spawn pbfs_walk_Bag(b, next, newdist, distances, parents);
        pbfs_walk_Pennant(p, next, newdist, distances, parents);
    } else {
        pbfs_proc_Node( b->getFilling(), next, newdist, distances, 
                        parents, b->getFillingSize() );
    }
    cilk_sync;
}

// Helper routine to walk a pennant for PBFS
inline void
Graph::pbfs_walk_Pennant( PennantType *p, 
                          reducer_ptr<BagRdcerType::Monoid> &next, 
                          uint newdist, uint distances[], 
                          MinRdcerPtr parents[] ) const {
    if(p->getLeft() != NULL) {
        cilk_spawn pbfs_walk_Pennant(
                        p->getLeft(), next, newdist, distances, parents);
        if(p->getRight() != NULL)
            cilk_spawn pbfs_walk_Pennant(
                        p->getRight(), next, newdist, distances, parents);
    }
    // Process the current element
    pbfs_proc_Node( p->getElements(), next, newdist, distances, 
                    parents, BLK_SIZE );
    cilk_sync;

    BagType::PAlloc().deallocate(p, 1);
}

// Helper routine to process a node for PBFS
inline void
Graph::pbfs_proc_Node( const uint n[], reducer_ptr<BagRdcerType::Monoid> &next, 
                       uint newdist, uint distances[],  
                       MinRdcerPtr parents[], uint fillSize ) const {

    uint* localFilling = BagType::TAlloc().allocate(BLK_SIZE);
    uint localFillingSize = 0;
    BagType *bnext = &(*next);

    for(uint j = 0; j < fillSize; j++) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        for(uint i = nodes[n[j]]; i < nodes[n[j]+1]; ++i) {

            WHEN_COMPUTE_INTENSE({
                for(uint k = 0; k < TIMES_TO_CALL_TRIVIAL; k++) 
                    trivial(); 
            })

            if(newdist < distances[edges[i]]) {
                distances[edges[i]] = newdist;
	            // parents[edges[i]]->set(n[j]);
	            *(parents[edges[i]]) = n[j];

                localFilling[localFillingSize++] = edges[i];
                if(localFillingSize == BLK_SIZE) {
                    bnext->insert_fblk(localFilling);
                    localFilling = BagType::TAlloc().allocate(BLK_SIZE);
                    localFillingSize = 0;
                }
            } else if (newdist == distances[edges[i]]) {
                // parents[edges[i]]->take_min(n[j]);
                if( *parents[edges[i]] > n[j] ) 
                    *parents[edges[i]] = n[j];
            }
        }
    }

    bnext->insert_blk(localFilling, localFillingSize);
}
#endif // USE_REDUCER_ARRAy
#else // !COMPUTE_PARENTS
inline void
Graph::pbfs_walk_Bag( BagType *b, BagRdcerType *next, 
                      uint newdist, uint distances[] ) const {

    PennantType *p = NULL;

    if(b->getFill() > 0) { 
        // Split the bag and recurse
        b->split(&p); 
        cilk_spawn pbfs_walk_Bag(b, next, newdist, distances);
        pbfs_walk_Pennant(p, next, newdist, distances);
    } else {
        pbfs_proc_Node( b->getFilling(), next, newdist, distances, 
                        b->getFillingSize() );
    }
    cilk_sync;
}

// Helper routine to walk a pennant for PBFS
inline void
Graph::pbfs_walk_Pennant( PennantType *p, BagRdcerType *next, 
                          uint newdist, uint distances[] ) const {

    if(p->getLeft() != NULL) {
        cilk_spawn pbfs_walk_Pennant(p->getLeft(), next, newdist, distances);
        if(p->getRight() != NULL)
            cilk_spawn pbfs_walk_Pennant(p->getRight(), next, newdist, 
                                         distances);
    }
    // Process the current element
    pbfs_proc_Node( p->getElements(), next, newdist, distances, BLK_SIZE );
    cilk_sync;

    BagType::PAlloc().deallocate(p, 1);
}

// Helper routine to process a node for PBFS
inline void
Graph::pbfs_proc_Node( const uint n[], BagRdcerType *next, 
                       uint newdist, uint distances[], uint fillSize ) const {

    uint* localFilling = BagType::TAlloc().allocate(BLK_SIZE);
    uint localFillingSize = 0;
    BagType *bnext = &(next->get_reference());

    for(uint j = 0; j < fillSize; j++) { 
        // Scan the edges of the current node and add untouched
        // neighbors to the opposite bag
        for(uint i = nodes[n[j]]; i < nodes[n[j]+1]; ++i) {

            WHEN_COMPUTE_INTENSE(
                for(uint k = 0; k < TIMES_TO_CALL_TRIVIAL; k++) 
                    trivial(); 
            )

            if(newdist < distances[edges[i]]) {
                // benign race on distances array
                distances[edges[i]] = newdist;

                localFilling[localFillingSize++] = edges[i];
                if(localFillingSize == BLK_SIZE) {
                    bnext->insert_fblk(localFilling);
                    localFilling = BagType::TAlloc().allocate(BLK_SIZE);
                    localFillingSize = 0;
                }
            }
        }
    }
    bnext->insert_blk(localFilling, localFillingSize);
}
#endif // COMPUTE_PARENTS


// PBFS main routine
#ifdef COMPUTE_PARENTS

// ANGE: helper method to parallelize setting up the distances and parents 
// from source node 
#if USE_REDUCER_ARRAY
inline void 
Graph::pbfs_setup_helper(uint begin, uint end, const uint s, 
                         uint distances[], MinRdcerArray &parents,
                         reducer_ptr<BagRdcerType::Monoid> &bagR) const {
tail_recur:
    if( (end - begin) > GRAINSIZE ) {
        uint mid = begin + ((uint) (end - begin)/2);
        cilk_spawn pbfs_setup_helper(begin, mid, s, distances, parents, bagR); 
        begin = mid;
        goto tail_recur;
    } else {
        uint i;
        // ANGE: lift the lookup out of the loop
        BagType *bag = bagR.view_ptr();
        MinRdcerArray::View &p = parents.view_ref();
 
        for(i = begin; i < end; i++) {
            if(edges[i] != s) {
                bag->insert(edges[i]);
                distances[edges[i]] = 1;
                p[edges[i]] = s; 
            }
        }
    }
    cilk_sync;
}

#else // !USE_REDUCER_ARRAY
inline void 
Graph::pbfs_setup_helper(uint begin, uint end, const uint s, 
                         uint distances[], MinRdcerPtr parents[],
                         reducer_ptr<BagRdcerType::Monoid> &bagR) const {
tail_recur:
    if( (end - begin) > GRAINSIZE ) {
        uint mid = begin + ((uint) (end - begin)/2);
        cilk_spawn pbfs_setup_helper(begin, mid, s, distances, parents, bagR); 
        begin = mid;
        goto tail_recur;
    } else {
        uint i;
        // ANGE: lift the lookup out of the loop
        BagType *bag = bagR.view_ptr();
 
        for(i = begin; i < end; i++) {
            if(edges[i] != s) {
                bag->insert(edges[i]);
                distances[edges[i]] = 1;
	            *(parents[edges[i]]) = s;
            }
        }
    }
    cilk_sync;
}
#endif

// Assumes distances array has been initialized
void Graph::pbfs( const uint s, uint distances[], uint parents[]) const {

    BagType bags[2];
    reducer_ptr<BagRdcerType::Monoid> bag_ptr[2];
    bag_ptr[0] = reducer_ptr<BagRdcerType::Monoid>(&(bags[0]));
    bag_ptr[1] = reducer_ptr<BagRdcerType::Monoid>(&(bags[1]));

    bool queuei = 1;
    uint current, newdist;

    if(s > nNodes) {
        fprintf(stderr, "Starting node is greater than number of nodes.\n");
        return;
    }
    
    #pragma cilk grainsize=GRAINSIZE
     cilk_for(uint i = 0; i < nNodes; ++i) {
        distances[i] = UINT_MAX;
        parents[i] = UINT_MAX;
    }
    distances[s] = 0;
    uint begin = nodes[s], end = nodes[s+1];

#if USE_REDUCER_ARRAY
    MinRdcerArray p(nNodes, parents); 
    MinRdcerArray::View *pView = p.view_ptr();
#else 
    MinRdcerPtr *p = (MinRdcerPtr *) 
                        scalable_malloc( sizeof(MinRdcerPtr) * nNodes );
    for(uint i = 0; i < nNodes; ++i) {
        p[i] = MinRdcerPtr( &(parents[i]) );
    }
#endif

    if( (end - begin) <= GRAINSIZE ) {
        for(uint i = begin; i < end; i++) {
            if(edges[i] != s) {
                bags[queuei].insert(edges[i]);
                distances[edges[i]] = 1;
                // ANGE: properly speaking, we should be accessing the View 
                // here if we are using RdcerMinArray, so that the access is
                // actually logged, but the left most view never logs
                // anyway, so it's ok. 
                parents[edges[i]] = s; 
            }
        }
    } else {
        pbfs_setup_helper(begin, end, s, distances, p, bag_ptr[queuei]);
    }
    newdist = 2;

    while( !(bags[queuei].isEmpty()) ) {
        bags[!queuei].clear();
        pbfs_walk_Bag(&(bags[queuei]), bag_ptr[!queuei], newdist, distances, p);
        queuei = !queuei;
        ++newdist;
    }

#if !USE_REDUCER_ARRAY
    scalable_free(p);
#endif
}

#else // !COMPUTE_PARENTS
int Graph::pbfs(const uint s, uint distances[]) const {

    BagRdcerType bags[2];

    bool queuei = 1;
    uint current, newdist;

    if(s > nNodes)
        return -1;

    #pragma cilk grainsize=GRAINSIZE
    cilk_for(uint i = 0; i < nNodes; ++i) {
        distances[i] = UINT_MAX;
    }
    distances[s] = 0;

    #pragma cilk grainsize=GRAINSIZE
    cilk_for(uint i = nodes[s]; i < nodes[s+1]; ++i) {
        if(edges[i] != s) {
            bags[queuei].insert(edges[i]);
            distances[edges[i]] = 1;
        }
    }
    newdist = 2;

    while( !(bags[queuei].isEmpty()) ) {
        bags[!queuei].clear();
        pbfs_walk_Bag(&(bags[queuei].get_reference()), 
                      &(bags[!queuei]), newdist, distances);
        queuei = !queuei;
        ++newdist;
    }

    return 0;
}
#endif // COMPUTE_PARENTS 

#endif // #ifndef GRAPH_CPP
