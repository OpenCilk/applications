// -*- C++ -*-

#ifndef BAG_H
#define BAG_H

#include <stdlib.h>
#include <assert.h>
#include <cilk/cilk.h>
#include <cilk/reducer.h>
#include <cilk/cilk_api.h>

#include <new>
#include <tbb/scalable_allocator.h>

#ifndef BAG_PROFILE
#define BAG_PROFILE 0
#endif

#if BAG_PROFILE
#define WHEN_BAG_PROFILE(ex) ex
#else
#define WHEN_BAG_PROFILE(ex)
#endif // BAG_PROFILE


const uint BAG_SIZE = 64;
const uint BLK_SIZE = 2048;

template <typename T, typename Alloc> class Bag;
template <typename T, typename Alloc> class Bag_reducer;

template <typename T, typename TAlloc>
class Pennant {
private:
    Pennant<T,TAlloc> *l, *r;
    T* els;
  
public:
    Pennant();
    Pennant(T*);
    ~Pennant();
    
    inline const T* getElements();
    inline Pennant<T,TAlloc>* getLeft();
    inline Pennant<T,TAlloc>* getRight();
    
    inline void clearChildren() {
        l = NULL;
        r = NULL;
    }
  
    inline Pennant<T,TAlloc>* combine(Pennant<T,TAlloc>*);
    inline Pennant<T,TAlloc>* split();
};

template < typename T, typename Alloc = tbb::scalable_allocator<void> >
class Bag {

public:
    typedef typename Alloc::template rebind<T>::other TAlloc;
    typedef typename Alloc::template rebind< Pennant<T,TAlloc> >::other PAlloc;
    typedef typename Alloc::template rebind< Pennant<T,TAlloc>* >::other BAlloc;


    Bag();
    // ANGE: BUGFIX --- make it private
    // Bag(Bag<T,Alloc>*);
    
    ~Bag();
    
    inline void allocBag() {
        this->bag = BAlloc().allocate(BAG_SIZE);
        for(uint i = 0; i < BAG_SIZE; ++i)
            this->bag[i] = NULL;
    }
    
    inline void allocFilling() {
        // ANGE: BUGFIX
        this->filling = TAlloc().allocate(BLK_SIZE); 
    }
    
    inline void insert(T);
    inline void insert_fblk(T* fblk);
    inline void insert_blk(T* blk, uint size);
    void merge(Bag<T,Alloc>*);
    inline bool split(Pennant<T,TAlloc>**);
    int split(Pennant<T,TAlloc>**, int);
    
    inline uint numPennants() const;
    inline uint numElements() const;
    inline uint getFill() const;
    inline bool isEmpty() const;
    inline Pennant<T,TAlloc>* getFirst() const;
    inline T* getFilling() const;
    inline uint getFillingSize() const;
   
    WHEN_BAG_PROFILE(
        inline unsigned int getViewsCreated() { return num_of_views; }
        inline unsigned int getExtraFillingsCreated() { return extra_fillings; }
    )
    
    void clear();
    
private:
    // fill points to the bag entry one position beyond the MSP
    uint fill;
    Pennant<T,TAlloc>* *bag;
    T* filling;

    // ANGE: BUGFIX --- private copy constructor
    Bag(Bag<T,Alloc>*);
  
    // size points to the filling array entry on position
    // beyond last valid element.
    uint size;
    WHEN_BAG_PROFILE( unsigned int num_of_views; )
    WHEN_BAG_PROFILE( unsigned int extra_fillings; )
  
    inline void insert_h();
};

template < typename T, typename Alloc = tbb::scalable_allocator<void> >
class Bag_reducer {

public:
    struct Monoid:cilk::monoid_base< Bag<T,Alloc> > {
        static void reduce (Bag<T,Alloc> *left, Bag<T,Alloc> *right) {
            left->merge(right);
        }
    };
  
private:
    cilk::reducer<Monoid> imp_;
    
public:
    typedef typename Bag<T,Alloc>::TAlloc TAlloc;

    Bag_reducer();
        
    void insert(T);
    void merge(Bag_reducer<T,Alloc>*);
    inline bool split(Pennant<T,TAlloc>**);
    int split(Pennant<T,TAlloc>**, int);
    
    inline Bag<T,Alloc> &get_reference();
    
    inline uint numElements() const;
    inline uint getFill() const;
    inline bool isEmpty() const;
    inline Pennant<T,TAlloc>* getFirst() const;
    inline T* getFilling() const;
    inline uint getFillingSize() const;
    
    inline void clear();
};

template <typename T, typename Alloc>
Bag_reducer<T,Alloc>::Bag_reducer() : imp_() { }

template <typename T, typename Alloc>
void Bag_reducer<T,Alloc>::insert(T el) {
    imp_.view().insert(el);
}

template <typename T, typename Alloc>
void Bag_reducer<T,Alloc>::merge(Bag_reducer<T,Alloc>* that) {
    this->imp_.view().merge(that->imp_.view());
}

template <typename T, typename Alloc>
inline bool Bag_reducer<T,Alloc>::split(Pennant<T,TAlloc>* *p) {
    return imp_.view().split(p);
}

template <typename T, typename Alloc>
int Bag_reducer<T,Alloc>::split(Pennant<T,TAlloc>* *p, int pos) {
    return imp_.view().split(p, pos);
}

template <typename T, typename Alloc>
inline Bag<T,Alloc>& Bag_reducer<T,Alloc>::get_reference() {
    return imp_.view();
}

template <typename T, typename Alloc>
inline uint Bag_reducer<T,Alloc>::numElements() const {
    return imp_.view().numElements();
}

template <typename T, typename Alloc>
inline uint Bag_reducer<T,Alloc>::getFill() const {
    return imp_.view().getFill();
}

template <typename T, typename Alloc>
inline bool Bag_reducer<T,Alloc>::isEmpty() const {
    return imp_.view().isEmpty();
}

template <typename T, typename Alloc>
inline Pennant< T, typename Bag<T,Alloc>::TAlloc >* 
Bag_reducer<T,Alloc>::getFirst() const {
    return imp_.view().getFirst();
}

template <typename T, typename Alloc>
inline T* Bag_reducer<T,Alloc>::getFilling() const {
    return imp_.view().getFilling();
}

template <typename T, typename Alloc>
inline uint Bag_reducer<T,Alloc>::getFillingSize() const {
    return imp_.view().getFillingSize();
}

template <typename T, typename Alloc>
inline void Bag_reducer<T,Alloc>::clear() {
    imp_.view().clear();
}

#include "bag.cpp"

#endif

