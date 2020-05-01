#ifndef OPENCILK_REDUCER_H
#include <cilk/hyperobject_base.h>
#include <new>
#include <stdlib.h>

#define NN __attribute__((nonnull))

struct opencilk_reducer_base : protected __cilkrts_hyperobject_base {
    static opencilk_reducer_base *cast(__cilkrts_hyperobject_base *h NN) {
        return static_cast<opencilk_reducer_base *>(h);
    }
    opencilk_reducer_base(size_t view_size, void *left_view NN);
    virtual ~opencilk_reducer_base();
    virtual void identity(void *NN, void *NN) = 0;
    virtual void destroy(void *NN, void *NN) = 0;
    virtual void reduce(void *NN, void *NN, void *NN) = 0;
    static void *allocate(void *reducer NN, size_t bytes) NN;
    static void deallocate(void *reducer NN, void *view NN);
};

template <typename View>
struct opencilk_reducer final : protected opencilk_reducer_base {
  private:
    void identity(void *reducer NN, void *view NN) { new (view) View(); }
    void destroy(void *reducer NN, void *view NN) {
        static_cast<View *>(view)->~View();
    }
    void reduce(void *reducer NN, void *left_v NN, void *right_v NN) {
        View *left = static_cast<View *>(left_v);
        View *right = static_cast<View *>(right_v);
        left->reduce(right);
    }

  public:
    template <typename... Args>
    opencilk_reducer(Args... args)
        : opencilk_reducer_base(sizeof(View), &leftmost), leftmost(args...) {}
    ~opencilk_reducer() {}

    View &view() { return *static_cast<View *>(__cilkrts_hyper_lookup(this)); }

  private:
    View leftmost;
};

#undef NN

#endif /*  OPENCILK_REDUCER_H */
