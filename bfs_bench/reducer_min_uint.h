#ifndef REDUCER_MIN_UINT_H
#define REDUCER_MIN_UINT_H

#include <cilk/reducer.h>

#ifdef __cplusplus
#include <cassert>

using namespace cilk;

class reducer_min_uint {

public:
    typedef uint value_type;
    // View of reducer data
    struct Monoid: monoid_base<uint> {
        static void reduce(uint *left, uint *right) {
            if( *left > *right ) *left = *right; 
        }
        void identity(uint *p) const { *p = UINT_MAX; }
    };
private:
    // Hyperobject to serve up views
    reducer<Monoid> imp_;
    // Not copyable
    reducer_min_uint(const reducer_min_uint&);
    reducer_min_uint& operator=(const reducer_min_uint&);

public:
    // Construct a 'reducer_min' object with a value of 'Type()'.
    reducer_min_uint() : imp_() {}
    reducer_min_uint(const uint& init_val) : imp_(init_val) {}
};
#endif // __cplusplus


#endif // defined REDUCER_MIN_UINT_H
