#ifndef OPENCILK_REDUCER_H
#define OPENCILK_REDUCER_H
#include <new>
#include <stdlib.h>

#define NN __attribute__((nonnull))

// Wrap a class with a default constructor and reduce method.
template <typename View>
struct opencilk_reducer final {
  private:
    static void identity(void *view NN) { new (view) View(); }
    static void reduce(void *left_v NN, void *right_v NN) {
        View *left = static_cast<View *>(left_v);
        View *right = static_cast<View *>(right_v);
        left->reduce(right);
        right->~View();
    }

  public:
    template <typename... Args>
    opencilk_reducer(Args... args)
        : leftmost(args...) {}
    ~opencilk_reducer() {}

    View &view() { return leftmost; }

  private:
    View cilk_reducer(identity, reduce) leftmost;
};

#undef NN

#endif /*  OPENCILK_REDUCER_H */
