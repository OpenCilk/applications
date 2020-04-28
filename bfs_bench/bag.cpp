#ifndef BAG_CPP
#define BAG_CPP

#include <stdlib.h>
#include <string.h>
#include "bag.h"

#define MAX(a, b) \
    ({ typeof(a) _a = (a); \
       typeof(b) _b = (b); \
       _a ^ ( (_a ^ _b) & -(_a < _b) ); })

//////////////////////////////////
/// Pennant method definitions ///
//////////////////////////////////
template <typename T, typename TAlloc> 
Pennant<T, TAlloc>::Pennant()
    : l(NULL), r(NULL), els( TAlloc().allocate(BLK_SIZE) ) {}

// els_array must have size BLK_SIZE
template <typename T, typename TAlloc> 
Pennant<T,TAlloc>::Pennant(T els_array[])
    : l(NULL), r(NULL), els(els_array) {}

template <typename T, typename TAlloc> 
Pennant<T,TAlloc>::~Pennant() {
    TAlloc().deallocate(els, BLK_SIZE);
}

template <typename T, typename TAlloc> 
inline const T* Pennant<T,TAlloc>::getElements() {
    return this->els;
}

template <typename T, typename TAlloc> 
inline Pennant<T,TAlloc>* Pennant<T,TAlloc>::getLeft() {
    return l;
}

template <typename T, typename TAlloc> 
inline Pennant<T,TAlloc>* Pennant<T,TAlloc>::getRight() {
    return r;
}

/*
 * This method assumes that the pennant rooted at <that> has the
 * same number of nodes as the pennant rooted at <this>.
 */
template <typename T, typename TAlloc> 
inline Pennant<T,TAlloc>* Pennant<T,TAlloc>::combine(Pennant<T,TAlloc>* that) {
    that->r = this->l;
    this->l = that;

    return this;
}

/*
 * This method performs the inverse operation of
 * Pennant<T>::combine and places half of the split here and the
 * the other half at the returned pointer
 */
template <typename T, typename TAlloc> 
inline Pennant<T,TAlloc>* Pennant<T,TAlloc>::split() {

    Pennant<T,TAlloc>* that;

    that = this->l;
    this->l = that->r;
    that->r = NULL;

    return that;
}

//////////////////////////////
/// Bag method definitions ///
//////////////////////////////
template <typename T, typename Alloc> 
Bag<T,Alloc>::Bag() : fill(0), size(0) {
    this->bag = BAlloc().allocate(BAG_SIZE);
    this->filling = TAlloc().allocate(BLK_SIZE);
    WHEN_BAG_PROFILE( num_of_views = 1; )
    WHEN_BAG_PROFILE( extra_fillings = 0; )
}

/*
 * Copy Constructor. Performs a shallow bag copy.
 * Useful for verifying Bag implementation correctness, since
 * the split method destroys the Bag.

ANGE: BUGFIX I don't like copy contructor
template <typename T, typename Alloc> 
Bag<T,Alloc>::Bag(Bag<T,Alloc> *that) : fill(that->fill), size(that->size) {

    // ANGE: BUGFIX ? Why shallow copy the filling but not also the bag
    this->bag = BAlloc().allocate(BAG_SIZE);
    for(uint i = 0; i < BAG_SIZE; i++)
        this->bag[i] = that->bag[i];

    this->filling = that->filling;
}
*/

template <typename T, typename Alloc> 
Bag<T,Alloc>::~Bag() {

    if(this->bag != NULL) {
        for(uint i = 0; i < this->fill; i++) {
            if(this->bag[i] != NULL) {
                this->bag[i]->~Pennant<T,TAlloc>();
                PAlloc().deallocate(this->bag[i], 1);
                this->bag[i] = NULL;
            }
        }

        BAlloc().deallocate(this->bag, BAG_SIZE);
    }
    TAlloc().deallocate(this->filling, BLK_SIZE);
}

template <typename T, typename Alloc> 
inline uint Bag<T,Alloc>::numPennants() const {

    uint count = 0;
    uint k = 1;

    for(uint i = 0; i < this->fill; i++) {
        if(this->bag[i] != NULL)
            count += k;
        k = k * 2;
    }

    return count;
}

template <typename T, typename Alloc> 
inline uint Bag<T,Alloc>::numElements() const {

    uint count = this->size;
    count += this->numPennants() * BLK_SIZE;

    return count;
}

// helper routine to perform bag-insert with filled filling array
template <typename T, typename Alloc> 
inline void Bag<T,Alloc>::insert_h() {

    Pennant<T,TAlloc> *c = PAlloc().allocate(1);
    new (c) Pennant<T,TAlloc>(this->filling);
    this->filling = TAlloc().allocate(BLK_SIZE);
    this->size = 0;
    uint i = 0;

    do {
        if(i < this->fill && this->bag[i] != NULL) {
            c = this->bag[i]->combine(c);
            this->bag[i] = NULL;
        } else {
            this->bag[i] = c;
            this->fill = MAX(this->fill, i+1);
            return;
        }
        ++i;
    } while(i < BAG_SIZE);

    this->fill = BAG_SIZE;
}

// helper routine to perform bag-insert with filled filling array
template <typename T, typename Alloc> 
inline void Bag<T,Alloc>::insert_fblk(T fblk[]) {

    Pennant<T,TAlloc> *c = PAlloc().allocate(1);
    new (c) Pennant<T,TAlloc>(fblk);

    uint i = 0;

    do {
        if(i < this->fill && this->bag[i] != NULL) {
            c = this->bag[i]->combine(c);
            this->bag[i] = NULL;
        } else {
            this->bag[i] = c;
            this->fill = MAX(this->fill, i+1);
            return;
        }
        ++i;
    } while(i < BAG_SIZE);

    this->fill = BAG_SIZE;
}

template <typename T, typename Alloc> 
inline void Bag<T,Alloc>::insert_blk(T blk[], uint size) {

    int i;

    // Deal with the partially-filled Pennants
    if(this->size < size) {
        // Copy contents of this->filling into blk
        i = this->size - (BLK_SIZE - size);

        if(i >= 0) {
            // Contents of this->filling fill blk
            memcpy( blk+size, this->filling+i, (BLK_SIZE-size) * sizeof(T) );
            //carry = blk;
            insert_fblk(blk);
            this->size = i;

        } else {
            // Contents of this->filling do not fill blk
            memcpy( blk+size, this->filling, this->size * sizeof(T) );
            TAlloc().deallocate(this->filling, BLK_SIZE);
            WHEN_BAG_PROFILE( this->extra_fillings += 1; )

            this->filling = blk;
            this->size += size;
        }

    } else {
        // Copy contents of blk into this->filling
        T* carry;
        i = size - (BLK_SIZE - this->size);

        if(i >= 0) {
            // Contents of blk fill this->filling
            memcpy( this->filling + this->size, 
                    blk + i, (BLK_SIZE - this->size) * sizeof(T) );

            carry = this->filling;
            this->filling = blk;

            this->size = i;
            insert_fblk(carry);

        } else {
            // Contents of blk do not fill this->filling
            memcpy( this->filling + this->size, blk, size * sizeof(T) );
            this->size += size;
        }
    }
}

template <typename T, typename Alloc> 
inline void Bag<T,Alloc>::insert(T el) {

    this->filling[this->size++] = el;
    if(this->size < BLK_SIZE) {
        return;
    }

    Pennant<T,TAlloc> *c = PAlloc().allocate(1);
    new (c) Pennant<T,TAlloc>(this->filling);
    this->filling = TAlloc().allocate(BLK_SIZE);
    this->size = 0;
    uint i = 0;

    do {
        if(i < this->fill && this->bag[i] != NULL) {
            c = this->bag[i]->combine(c);
            this->bag[i] = NULL;
        } else {
            this->bag[i] = c;
            //this->fill = this->fill > i+1 ? this->fill : i+1;
            this->fill = MAX(this->fill, i+1);
            return;
        }
        ++i;
    } while(i < BAG_SIZE);

    this->fill = BAG_SIZE;
}

template <typename T, typename Alloc> 
void Bag<T,Alloc>::merge(Bag<T,Alloc>* that) {

    Pennant<T,TAlloc> *c = NULL;
    T* carry = NULL;
    char x;
    int i;

    WHEN_BAG_PROFILE( this->num_of_views += that->num_of_views; )
    WHEN_BAG_PROFILE( this->extra_fillings += that->extra_fillings; )

    // Deal with the partially-filled Pennants
    if(this->size < that->size) {
        i = this->size - (BLK_SIZE - that->size);

        if(i >= 0) {
            memcpy( that->filling + that->size,
                    this->filling + i, (BLK_SIZE - that->size) * sizeof(T) );

            carry = that->filling;
            this->size = i;
        } else {
            memcpy( that->filling + that->size,
                    this->filling, this->size * sizeof(T) );
            TAlloc().deallocate(this->filling, BLK_SIZE);
            WHEN_BAG_PROFILE( this->extra_fillings += 1; )

            this->filling = that->filling;
            this->size += that->size;
        }

    } else {
        i = that->size - (BLK_SIZE - this->size);

        if(i >= 0) {
            memcpy( this->filling + this->size,
                    that->filling + i, (BLK_SIZE - this->size) * sizeof(T) );
            carry = this->filling;

            this->filling = that->filling;
            this->size = i;

        } else {
            memcpy( this->filling + this->size,
                    that->filling, that->size * sizeof(T) );
            this->size += that->size;
        }
    }

    that->filling = NULL;
    // Update this->fill (assuming no final carry)
    uint min, max;
    if(this->fill > that->fill) {
        min = that->fill;
        max = this->fill;
    } else {
        min = this->fill;
        max = that->fill;
    }

    if(carry != NULL) {
        c = PAlloc().allocate(1);
        new (c) Pennant<T,TAlloc>(carry);
    }

    // Merge
    for(i = 0; i < min; ++i) {
        x = (this->bag[i] != NULL) << 2 |
            (i < that->fill && that->bag[i] != NULL) << 1 |
            (c != NULL);

        switch(x) {
            case 0x0:
                this->bag[i] = NULL;
                c = NULL;
                that->bag[i] = NULL;
                break;
            case 0x1:
                this->bag[i] = c;
                c = NULL;
                that->bag[i] = NULL;
                break;
            case 0x2:
                this->bag[i] = that->bag[i];
                that->bag[i] = NULL;
                c = NULL;
                break;
            case 0x3:
                c = that->bag[i]->combine(c);
                that->bag[i] = NULL;
                this->bag[i] = NULL;
                break;
            case 0x4:
                /* this->bag[i] = this->bag[i]; */
                c = NULL;
                that->bag[i] = NULL;
                break;
            case 0x5:
                c = this->bag[i]->combine(c);
                this->bag[i] = NULL;
                that->bag[i] = NULL;
                break;
            case 0x6:
                c = this->bag[i]->combine(that->bag[i]);
                that->bag[i] = NULL;
                this->bag[i] = NULL;
                break;
            case 0x7:
                /* this->bag[i] = this->bag[i]; */
                c = that->bag[i]->combine(c);
                that->bag[i] = NULL;
                break;
            default: break;
        }
    }

    that->fill = 0;

    if(this->fill == max) {
        if(c == NULL) return;

        do {
            if(i < max && this->bag[i] != NULL) {
                c = this->bag[i]->combine(c);
                this->bag[i] = NULL;
            } else {
                this->bag[i] = c;
                this->fill = MAX(max, i+1);
                return;
            }
            ++i;
        } while(i < BAG_SIZE);

    } else { // that->fill == max
        int j;

        if(c == NULL) {
            this->fill = max;
            for(j = i; j < this->fill; ++j) {
                this->bag[j] = that->bag[j];
                that->bag[j] = NULL;
            }
            return;
        }

        do {
            if(i < max && that->bag[i] != NULL) {
                c = that->bag[i]->combine(c);
                this->bag[i] = NULL;
                that->bag[i] = NULL;
            } else {
                this->bag[i] = c;
                this->fill = MAX(max, i+1);

                for(j = i+1; j < this->fill; ++j) {
                    this->bag[j] = that->bag[j];
                    that->bag[j] = NULL;
                }
                return;
            }

            ++i;
        } while(i < BAG_SIZE);    
    }

    this->fill = BAG_SIZE;
}

template <typename T, typename Alloc>
inline bool Bag<T,Alloc>::split(Pennant<T,TAlloc>* *p) {

    if(this->fill == 0) {
        *p = NULL;
        return false;
    }

    this->fill--;
    *p = this->bag[this->fill];
    this->bag[this->fill] = NULL;

    for(this->fill; this->fill > 0; this->fill--) {
        if(this->bag[this->fill-1] != NULL)
            break;
    }
    return true;
}

template <typename T, typename Alloc>
int Bag<T,Alloc>::split(Pennant<T,TAlloc>* *p, int pos) {

    if(pos < 0 || this->fill <= pos) {
        *p = NULL;
        return this->fill - 1;
    }

    *p = this->bag[pos];

    for(int i = pos-1; i >= 0; i--) {
        if(this->bag[i] != NULL)
            return i;
    }
    return -1;
}

template <typename T, typename Alloc>
inline uint Bag<T,Alloc>::getFill() const {
    return this->fill;
}

template <typename T, typename Alloc>
inline bool Bag<T,Alloc>::isEmpty() const {
    return this->fill == 0 && this->size == 0;
}

template <typename T, typename Alloc>
inline Pennant<T, typename Bag<T,Alloc>::TAlloc>* 
Bag<T,Alloc>::getFirst() const {
    return this->bag[0];
}

template <typename T, typename Alloc>
inline T* Bag<T,Alloc>::getFilling() const {
    return this->filling;
}

template <typename T, typename Alloc>
inline uint Bag<T,Alloc>::getFillingSize() const {
    return this->size;
}

template <typename T, typename Alloc>
inline void Bag<T,Alloc>::clear() {

    // reset the profiling counter
    WHEN_BAG_PROFILE({
        this->num_of_views = 0; 
        this->extra_fillings = 0; 
    })
    
    // ANGE: BUGFIX I added this if; memory leak otherwise
    if(this->bag) { 
        for(uint i = 0; i < this->fill; i++) {
            if(this->bag[i]) {
                this->bag[i]->~Pennant<T,TAlloc>();
                PAlloc().deallocate(this->bag[i], 1);
                this->bag[i] = NULL;
            }
        }
    }
    this->fill = 0;
    this->size = 0;
}

#endif // BAG_CPP

