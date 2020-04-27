// -*- C++ -*-

/*
 * qsort.cilk
 *
 * An implementation of quicksort using Cilk parallelization.
 *
 * Copyright (c) 2007-2008 Cilk Arts, Inc.  55 Cambridge Street,
 * Burlington, MA 01803.  Patents pending.  All rights reserved. You may
 * freely use the sample code to guide development of your own works,
 * provided that you reproduce this notice in any works you make that
 * use the sample code.  This sample code is provided "AS IS" without
 * warranty of any kind, either express or implied, including but not
 * limited to any implied warranty of non-infringement, merchantability
 * or fitness for a particular purpose.  In no event shall Cilk Arts,
 * Inc. be liable for any direct, indirect, special, or consequential
 * damages, or any other damages whatsoever, for any use of or reliance
 * on this sample code, including, without limitation, any lost
 * opportunity, lost profits, business interruption, loss of programs or
 * data, even if expressly advised of or otherwise aware of the
 * possibility of such damages, whether in an action of contract,
 * negligence, tort, or otherwise.
 *
 * Modified for OpenCilk April 2020.
 */
#include <algorithm>
#include <iostream>
#include <time.h>

typedef long Datum;

// Sort the range between bidirectional iterators begin and end.
// end is one past the final element in the range.
// Use the Quick Sort algorithm, using recursive divide and conquer.
// This function is NOT the same as the Standard C Library qsort() function.
// This implementation is pure C++ code before Cilk++ conversion.
template<typename T>
void sample_qsort(T * begin, T * end)
{
    if (begin != end) {
        // Exclude last element (pivot) from partition
        --end;
        T * middle = std::partition(begin, end,
            [pivot = *end](T a) -> bool { return a < pivot; });
        using std::swap;
        swap(*end, *middle);    // move pivot to middle
        _Cilk_spawn sample_qsort(begin, middle);
        sample_qsort(++middle, ++end); // Exclude pivot and restore end
        _Cilk_sync;
    }
}

static void print_time(std::ostream &out, struct timespec start,
                       struct timespec end)
{
    time_t s = end.tv_sec - start.tv_sec;
    long ns = end.tv_nsec - start.tv_nsec;
    if (ns < 0) {
	ns += 1000000000;
	s -= 1;
    }
    out << (s * 1000 + ((ns + 500000) / 1000000));
}

// A simple test harness 
int qmain(int n)
{
    Datum* a = new Datum[n];
    Datum* b = new Datum[n];

#if 0 /* slow */
    for (int i = 0; i < n; ++i)
        a[i] = i;

    std::random_shuffle(a, a + n);
#else
    for (int i = 0; i < n; ++i) {
	long r = random();
        a[i] = r;
	b[i] = r;
    }
#endif
    std::cout << "Sorting " << n << " integers" << std::endl;

    struct timespec start_stl, end_stl;
    if (clock_gettime(CLOCK_MONOTONIC, &start_stl) < 0) {
	perror("clock_gettime");
	return 1;
    }
    std::sort(&b[0], &b[n]);
    if (clock_gettime(CLOCK_MONOTONIC, &end_stl) < 0) {
	perror("clock_gettime");
	return 1;
    }

    struct timespec start, end;
    if (clock_gettime(CLOCK_MONOTONIC, &start) < 0) {
	perror("clock_gettime");
	return 1;
    }
    sample_qsort(a, a + n);
    if (clock_gettime(CLOCK_MONOTONIC, &end) < 0) {
	perror("clock_gettime");
	return 1;
    }

    std::cout << "Sorted in ";
    print_time(std::cout, start, end);
    std::cout << " ms with Cilk quicksort, ";
    print_time(std::cout, start_stl, end_stl);
    std::cout << " with STL sort\n";

    bool pass = true;
    // Confirm that a is sorted
    for (int i = 0; i < n - 1; ++i) {
        if (a[i] > a[i + 1] || a[i] != b[i]) {
            std::cout << "Sort failed at location i=" << i << " a[i] = "
                      << a[i] << " a[i+1] = " << a[i + 1] << std::endl;
	    pass = false;
	    break;
        }
    }
    if (pass) {
	std::cout << "Sort succeeded." << std::endl;
    }
    delete[] a;
    delete[] b;
    return pass ? 0 : 1;
}

int main(int argc, char* argv[])
{
    int n = 10 * 1000 * 1000;
    if (argc > 1) {
         n = std::atoi(argv[1]);
         if (n <= 0) {
              std::cerr << "Invalid argument" << std::endl;
              std::cerr << "Usage: qsort N" << std::endl;
              std::cerr << "       N = number of elements to sort" << std::endl;
              return 1;
         }
    }

    return qmain(n);
}
