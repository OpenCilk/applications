#include "opencilk_reducer.hpp"

#define NN __attribute__((nonnull))

void *opencilk_reducer_base::allocate(void *reducer, size_t bytes) {
    return malloc(bytes);
}

void opencilk_reducer_base::deallocate(void *reducer, void *view) {
    free(view);
}

namespace {

extern "C" void reduce_fn(void *reducer NN, void *left NN, void *right NN) {
    __cilkrts_hyperobject_base *h =
        static_cast<__cilkrts_hyperobject_base *>(reducer);
    opencilk_reducer_base::cast(h)->reduce(reducer, left, right);
}

extern "C" void identity_fn(void *reducer NN, void *view NN) {
    __cilkrts_hyperobject_base *h =
        static_cast<__cilkrts_hyperobject_base *>(reducer);
    opencilk_reducer_base::cast(h)->identity(reducer, view);
}

extern "C" void destroy_fn(void *reducer NN, void *view NN) {
    __cilkrts_hyperobject_base *h =
        static_cast<__cilkrts_hyperobject_base *>(reducer);
    opencilk_reducer_base::cast(h)->destroy(reducer, view);
}

extern "C" void *allocate_fn(void *reducer NN, size_t bytes) {
    __cilkrts_hyperobject_base *h =
        static_cast<__cilkrts_hyperobject_base *>(reducer);
    return opencilk_reducer_base::cast(h)->allocate(reducer, bytes);
}

extern "C" void deallocate_fn(void *reducer NN, void *view) {
    __cilkrts_hyperobject_base *h =
        static_cast<__cilkrts_hyperobject_base *>(reducer);
    opencilk_reducer_base::cast(h)->deallocate(reducer, view);
}

} // namespace

opencilk_reducer_base::opencilk_reducer_base(size_t view_size,
                                             void *left_view) {
    char *monoid_addr = static_cast<char *>(
        static_cast<void *>(static_cast<__cilkrts_hyperobject_base *>(this)));
    char *view_addr = static_cast<char *>(left_view);

    __c_monoid.reduce_fn = &reduce_fn;
    __c_monoid.identity_fn = &identity_fn;
    __c_monoid.destroy_fn = &destroy_fn;
    __c_monoid.allocate_fn = &allocate_fn;
    __c_monoid.deallocate_fn = &deallocate_fn;
    __id_num = ~0;
    __view_offset = view_addr - monoid_addr;
    __view_size = view_size;
    __cilkrts_hyper_create(static_cast<__cilkrts_hyperobject_base *>(this));
}

opencilk_reducer_base::~opencilk_reducer_base() {
    __cilkrts_hyper_destroy(this);
}
