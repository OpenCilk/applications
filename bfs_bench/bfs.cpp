#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/types.h>
#include <climits>
#include <sys/time.h>
#include <iostream>

#include <cilk/cilk.h>
#include "graph.h"
#include "cilktime.h"

#ifndef TIMES_TO_RUN
#define TIMES_TO_RUN 1
#endif


#if CILKVIEW 
#include <cilkview.h>
#define WHEN_CILKVIEW(ex) ex
#else
#define WHEN_CILKVIEW(ex)
#endif
/* To use cilkview to measure certain section of code:
 * cilkview_data_t start, end;
 * __cilkview_query(start);
 * __cilkview_query(end);
 * __cilkview_report(&start, &end, "BFS", 0);
*/


#if COMPUTE_INTENSE
extern "C" 
void trivial() {
    asm volatile("incq %%rax" : : : "memory");
    return;
}
#endif

// Get Cilk-M to work properly
#if CILKM_RTS
extern "C" __attribute__((visibility("default")))
int __cilkrts_main(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    return __cilkrts_main(argc, argv);
}

#else
extern "C" __attribute__((visibility("default")))
int cilk_main(int argc, char *argv[]);

int main(int argc, char *argv[]) {
    int ret;

    /*
    ret = __cilkrts_set_param("nworkers", "2");
    if(ret) {
        fprintf(stderr, "Setting nworker to 1 failed.\n");
    } */
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
    
/*
    ret = __cilkrts_set_param("nworkers", "2");
    if(ret) {
        fprintf(stderr, "Setting nworker to 2 failed.\n");
    }
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
    
    ret = __cilkrts_set_param("nworkers", "4");
    if(ret) {
        fprintf(stderr, "Setting nworker to 4 failed.\n");
    }
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
    
    ret = __cilkrts_set_param("nworkers", "8");
    if(ret) {
        fprintf(stderr, "Setting nworker to 8 failed.\n");
    }
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
    
    ret = __cilkrts_set_param("nworkers", "12");
    if(ret) {
        fprintf(stderr, "Setting nworker to 12 failed.\n");
    }
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
    
    ret = __cilkrts_set_param("nworkers", "16");
    if(ret) {
        fprintf(stderr, "Setting nworker to 16 failed.\n");
    }
    __cilkrts_init();
    ret = cilk_main(argc, argv);
    __cilkrts_end_cilk(); 
*/    

    return ret;
}
#endif



using namespace std;

static inline unsigned long long rdtsc() {
    unsigned long long r;
    __asm__ __volatile__ ("rdtsc\n"
                          "shl $32,%%rdx\n"
                          "or %%rdx,%%rax\n"
			  "movq %%rax,%0" : "=r"(r) : : "edx", "eax", "rdx", "rax");
    return r;
}

unsigned long long todval (struct timeval *tp) {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}

#define DEBUG false

#define CACHE_NUKE_SIZE 16*1048576

volatile char cache_nuke[CACHE_NUKE_SIZE];

static void
purge_cache() {

    for(int i = 0; i < CACHE_NUKE_SIZE; i++)
        cache_nuke[i]++;
}

static void
init_cache_nuke() {

    for(int i = 0; i < CACHE_NUKE_SIZE; i++)
        cache_nuke[i] = 0;
}

static int
CumulativeSum(int* arr, int size) {

    int prev;
    int tempnz = 0;

    for(int i = 0; i < size; ++i) {
        prev = arr[i];
        arr[i] = tempnz;
        tempnz += prev;
    }
    return tempnz;
}

static inline void
init_distances(uint distances[], int nodes) {

    cilk_for(size_t i = 0; i < nodes; ++i)
        distances[i] = UINT_MAX;
}

static bool
check(const uint distances[], const uint distverf[], int nodes) {

    for(int i = 0; i < nodes; i++) {
        if(distances[i] != distverf[i]) {
            fprintf(stderr,
                    "distances[%d] = %d; distverf[%d] = %d\n",
                    i, distances[i], i, distverf[i]);
            return false;
        }
    }
    return true;
}

extern "C"
int cilk_main(int argc, char** argv) {

    int m, n, nnz;
    string inputname;

    int bfs_sel;

    // Read arguments
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "Usage: [CILK_NWORKERS=<num_workers>] %s "
                "[<bfs_selection>] <inputname>\n\tbfs_selection=1 by default\n",
                argv[0]);
        return -1;
    } else {
        if(argc > 2) {
            bfs_sel = atoi(argv[1]);
            inputname = argv[2];
        } else {
            bfs_sel = 1; // PBFS is default selection
            inputname = argv[1];
        }
    }

    if(DEBUG)
        printf("bfs_sel = %d\n", bfs_sel);

    // Read binary CSB matrix input
    // Code and matrices adapted from oskitest.cpp by Aydin Buluc
    if(DEBUG)
        printf("Reading input file %s\n", inputname.c_str());

    FILE *f = fopen(inputname.c_str(), "r");
    if(!f) {
        fprintf(stderr, "Problem reading binary input file\n");
        return -1;
    }

    fread(&m, sizeof(int), 1, f);
    fread(&n, sizeof(int), 1, f);
    fread(&nnz, sizeof(int), 1, f);

    if(m <= 0 || n <= 0 || nnz <= 0) {
        fprintf(stderr,
                "Problem with matrix size in binary input file\n");
        return -1;
    }

    if(m != n) {
        fprintf(stderr,
                "Input file does not describe a graph\n");
        return -1;
    }

    if (DEBUG)
       fprintf(stderr, "Reading %d-by-%d matrix having %d nonzeros\n",
               m, n, nnz); 
    int *rowindices = (int *) scalable_malloc( sizeof(int) * nnz );
    int *colindices = (int *) scalable_malloc( sizeof(int) * nnz );
    double *vals = (double *) scalable_malloc( sizeof(double) * nnz );

    size_t rows = fread(rowindices, sizeof(int), nnz, f);
    size_t cols = fread(colindices, sizeof(int), nnz, f);
    size_t nums = fread(vals, sizeof(double), nnz, f);
    fclose(f);

    if(rows != nnz || cols != nnz || nums != nnz) {
        fprintf(stderr, "Problem with FREAD. Aborting.\n");
        return -1;
    }
  
    double *num = (double *) scalable_malloc( sizeof(double) * nnz );
    int *ir = (int *) scalable_malloc( sizeof(int) * nnz );
    int *jc = (int *) scalable_malloc( sizeof(int) * (n+1) );
    int *w = (int *) scalable_malloc( sizeof(int) * n );

    for(int k = 0; k < n; ++k)
        w[k] = 0;

    for(int k = 0; k < nnz; ++k)
        w[colindices[k]]++;

    jc[n] = CumulativeSum(w,n);
    for(int k = 0; k < n; ++k)
        jc[k] = w[k];

    int last;
    for(int k = 0; k < nnz; ++k) {
        ir[last = w[colindices[k]]++] = rowindices[k];
        num[last] = vals[k];
    }

    scalable_free(w);
    scalable_free(rowindices);
    scalable_free(colindices);
    scalable_free(vals);

    // Initialize extra data structures
    uint s;
    uint errors = 0;
    uint *distances = (uint *) scalable_malloc( sizeof(uint) * m );
#if COMPUTE_PARENTS
    uint *parents = (uint *) scalable_malloc( sizeof(uint) * m );
#endif

    //unsigned long long t1, t2, t_total = 0;
    //clock_t t1, t2, t_total = 0;
    struct timeval t1, t2;
    double time[TIMES_TO_RUN];

    struct results res;

    if(DEBUG)
        fprintf(stderr, "Making graph\n");
    // Make graph
    Graph graph(ir, jc, m, n, nnz);

    scalable_free(ir);
    scalable_free(jc);
    scalable_free(num);

    if(DEBUG)
        fprintf(stderr, "Initializing cache nuke\n");
    // Initialize cache nuke
    init_cache_nuke();

    s = 0;

    // Purge the cache
    if(DEBUG)
        fprintf(stderr, "Purging cache\n");
    purge_cache();

    if(DEBUG)
        fprintf(stderr, "Computing BFS\n");

    // Execute BFS
    for(int i=0; i < TIMES_TO_RUN; i++) {
        switch(bfs_sel) {

            case 0:
                gettimeofday(&t1,0);
#ifdef COMPUTE_PARENTS
                // graph.bfs(s, distances, parents);
                graph.bfs(s, distances, parents);
#else
                graph.bfs(s, distances);
#endif
                gettimeofday(&(t2),0);
                break;

            case 1:
                WHEN_CILKVIEW( cilkview_data_t start; )
                    WHEN_CILKVIEW( cilkview_data_t end; )

                    gettimeofday(&t1,0);

                WHEN_CILKVIEW( __cilkview_query(start); )

#ifdef COMPUTE_PARENTS
                    graph.pbfs(s, distances, parents);
#else
                graph.pbfs(s, distances);
#endif
                WHEN_CILKVIEW( __cilkview_query(end); )
                    gettimeofday(&t2,0);

                WHEN_CILKVIEW( __cilkview_report(&start, &end, 
                            "BFS", CV_REPORT_WRITE_TO_LOG); )
                    break;

            default:
                gettimeofday(&t1,0);
                gettimeofday(&t2,0);
                break;
        }
        time[i] = (todval(&t2)-todval(&t1)) / 1000000.0;
    }

    // Get verification distances
    if(DEBUG)
        fprintf(stderr, "Getting verification distances\n");

    uint *distverf  = (uint *) scalable_malloc( sizeof(uint) * m );
#ifdef COMPUTE_PARENTS
    uint *parentverf = (uint *) scalable_malloc( sizeof(uint) * m );
    graph.bfs(s, distverf, parentverf);
#else
    graph.bfs(s, distverf);
#endif

    // Verify correctness
    if(!check(distances, distverf, m))
        errors++;

#if COMPUTE_PARENTS
    if(parents && parentverf) {
        for(uint i = 0; i < m; ++i) {
            if(distverf[i] < UINT_MAX && parents[i] != parentverf[i]) {
                fprintf(stderr, "ERROR: parents[%u] = %u, ", i, parents[i]);
                fprintf(stderr, "parentverf[%u] = %u.\n", i, parentverf[i]);
                ++errors;
                break;
            }
        }
        scalable_free(parents);
        scalable_free(parentverf);
    }
#endif
    scalable_free(distverf);
    scalable_free(distances);

    printf( "BFS algorithm %d, %d iter on %2d workers, %u errors.\n",
            bfs_sel, cache_nuke[0], __cilkrts_get_nworkers(), errors );
    print_runtime(time, TIMES_TO_RUN);

    return 0;
}

