/*
Copyright (c) 2013, Intel Corporation
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, 
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list 
  of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this 
  list of conditions and the following disclaimer in the documentation and/or 
  other materials provided with the distribution.
* Neither the name of the Intel Corporation nor the names of its contributors may 
  be used to endorse or promote products derived from this software without 
 specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
OF THE POSSIBILITY OF SUCH DAMAGE.
*/ 

/* Cilk Plus implementation for algorithm described in the paper:

         "Eight Pieces Cannot Cover a Chess Board",
         A. D. Robison, B. J. Hafner, and S. S. Skiena,
         The Computer Journal, 32(6), 1989, p.567-570.
         
   Implementation by Arch D. Robison. */

#define NOMINMAX /* Required on Windows to use std::max */

//! Set to nonzero if parallel implementation is desired.
#ifndef PARALLEL
#define PARALLEL 1
#endif

//! Set to nonzero if bishops of the same color are allowed.
#ifndef BISHOPS_CAN_BE_ON_SAME_COLOR
#define BISHOPS_CAN_BE_ON_SAME_COLOR 1
#endif

#if PARALLEL
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/opadd_reducer.h>
#include <cilk/ostream_reducer.h>
#else
#define cilk_spawn
#endif

#include <iostream>
#include <cassert>
#include <cstdlib>

#include "ktiming.h"


//! Number of squared along an edge of a chess board
const long BoardSize = 8;

//! Number of squares on a chess board
const long N_Square = BoardSize*BoardSize;

//! Enumerate kinds of pieces
enum Piece {
    King,
    Queen,
    Bishop,       //!< Bishop on any square
    WhiteBishop,  //!< Bishop restricted to white squares in bottom-right triangl.
    BlackBishop,  //!< Bishop restricted to black squares.
    Knight,
    Rook
};

//! Array indexed by Piece
static const char Symbol[] = {'K','Q','B','B','B','N','R'};

//! A subset of squares on a chess board.
class SquareSet {
    typedef int64_t maskType;
    //! Bitmask representation of the set
    maskType mask;
    //! Make bitmask for square at position (x,y)
    static maskType makeMask( unsigned x, unsigned y ) {
        assert( x<BoardSize );
        assert( y<BoardSize );
        return maskType(1) << (y*BoardSize+x);
    }
    //! Add square at position (x,y) if it is on the board 
    void addIfOnBoard( size_t x, size_t y ) {
        if( x<BoardSize && y<BoardSize ) 
            add(x,y);
    }
    //! A direction vector.  Each component is drawn from {-1,0,1}.
    struct dirType {
        int dx;
        int dy;
    };
    //! Add all squares along rays starting at (x,y) exclusive, until an occupied square or edge of board is encountered.
    void addRay( unsigned x, unsigned y, const dirType* first, const dirType* last, SquareSet occupied );
public:
    //! Construct empty set
    SquareSet() : mask(0) {}
    //! Construct set of squares strongly attacked by piece p, given a set of occupied squares.
    SquareSet( Piece p, unsigned x, unsigned y, SquareSet occupied=SquareSet() );
    //! True if *this has square with coordinates (x,y)
    bool has( unsigned x, unsigned y ) const {
        return (mask & makeMask(x,y)) !=0;
    }
    //! Add square at (x,y) to *this.
    void add( size_t x, size_t y ) {
        mask |= makeMask(x,y);
    }
    //! Add squares in b to this set.
    void operator|=( const SquareSet& b ) {mask|=b.mask;}
    //! True if this set contains all squares
    bool isAll() const {return mask==~maskType(0);}
    //! Return number of squares in this set
    unsigned long count() const {
        unsigned long c = 0;
        for( maskType m=mask; m; m&=m-1 )
            ++c;
        return c;
    }
};

void SquareSet::addRay( unsigned x, unsigned y, const dirType* first, const dirType* last, SquareSet occupied ) {
    //! Loop over directions
    for( const dirType* d=first; d!=last; ++d ) {
        // Start at (x,y)
        unsigned u = x;
        unsigned v = y;
        do {
            // Move in direction *d
            u += d->dx;
            v += d->dy;
            if( u>=BoardSize || v>=BoardSize ) 
                // Walked off board
                break;
            add(u,v);
        } while( !occupied.has(u,v) );
    }
}

SquareSet::SquareSet( Piece p, unsigned x, unsigned y, SquareSet occupied ) : mask(0) {
    static const dirType dir[8] = {
        {0,1},{1,0},{-1,0},{0,-1},      // Compass points
        {1,1},{-1,1},{-1,-1},{1,-1}     // Diagnonals
    };
    switch( p ) {
        case King:
            for( size_t k=0; k<8; ++k )
                addIfOnBoard(x+dir[k].dx, y+dir[k].dy);
            break;
        case Queen:
            addRay( x, y, dir, dir+8, occupied );
            break;
        case Bishop:
        case WhiteBishop:
        case BlackBishop:
            addRay( x, y, dir+4, dir+8, occupied );
            break;
        case Knight: {
            for( int m=-2; m<=2; m+=4 )
                for( int n=-1; n<=1; n+=2 ) {
                    addIfOnBoard( x+m, y+n );
                    addIfOnBoard( x+n, y+m );
                }
            break;
        }
        case Rook: {
            addRay( x, y, dir, dir+4, occupied );
            break;
        }
    }
    assert(!has(x,y));
}

//! Return z with odd-indexed bits squashed out.
unsigned EvenBits( unsigned z ) {
    unsigned x = 0;
    unsigned mask = 1;
    for( ;z; z>>=2 ) {
        if( z&1 ) x|=mask;
        mask<<=1;
    }
    return x;
}

//! X-Y coordinates of a piece.
/** Each coordinate is in [0,BoardSize). */
struct Position {
    unsigned x;
    unsigned y;
};

//! A node in an axis tree as described in the paper.
struct Node {
    SquareSet weakAttacks;
    unsigned weight;            //!< 1 at a leaf, >1 for internal node
    unsigned char level;        //!< Distance to furthest leaf below this node.
    unsigned char leftmost;     //!< Index of leftmost leaf below this node
    unsigned char rightmost;    //!< Index of rightmost leaf below this node
    union {
        Node* child[2];         //!< Pointers to halves after unfolding this Node.
        Position pos;           //!< Position of leaf on board
    };
};

//! An axis tree as described in the paper.
class Tree {
    //! Pointer to root of tree
    Node* myRoot;
    //! Number of leaves in tree
    unsigned long myLeafCount;
    //! Initialize kth leaf to represent piece of kind p at square (x,y)
    void setLeaf( size_t k, Piece p, size_t x, size_t y ) {
        assert( x<BoardSize );
        assert( y<BoardSize );
        Node& n = myRoot[k+myLeafCount-1];
        n.weight = 1;
        n.pos.x = x;
        n.pos.y = y;
        n.weakAttacks = SquareSet( p, x, y );
        n.level = 0;
        n.leftmost = n.rightmost = k;
    };
public:
    //! Construct folding tree for piece of kind p
    Tree( Piece p );
    //! Return pointer to root 
    Node* root() const {return myRoot;}
};

Tree::Tree( Piece p ) {
    switch( p ) {
        case Queen: 
#if BISHOPS_CAN_BE_ON_SAME_COLOR
            // Queens occupy lower triangle of lower left quadrant (including diagonal)
            myLeafCount = BoardSize*(BoardSize+2)/8;
#else
            // Queens occupy only lower left corner
            myLeafCount = N_Square/4; 
#endif /* BISHOPS_CAN_BE_ON_SAME_COLOR */
            break;
        case WhiteBishop:
            myLeafCount = N_Square/4; 
            break;
        case BlackBishop: 
            myLeafCount = N_Square/2; 
            break;
        default:
            myLeafCount = N_Square;
            break;
    }
    myRoot = new Node[2*myLeafCount-1];
    for( unsigned long k=0; k<myLeafCount; ++k ) {
        int x, y;
        switch( p ) {
            default:
                assert(0);
            case King:
                x = EvenBits(k);
                y = EvenBits(k>>1);
                break;
            case Queen:
#if BISHOPS_CAN_BE_ON_SAME_COLOR
                static const int pos[] = {23,33,2,3,12,13,22,11,1,0};
                assert( sizeof(pos)/sizeof(pos[0])==myLeafCount );
                x = pos[k]%10;
                y = pos[k]/10;
#else
                x = EvenBits(k);
                y = EvenBits(k>>1);
#endif
                break;
            case WhiteBishop: {
                static const int pos[] = {34,25,23,14,45,36,56,47,16,7,27,5,1,12,67,3};
                assert( sizeof(pos)/sizeof(pos[0])==myLeafCount );
                x = pos[k]%10;
                y = pos[k]/10;
                assert( x+y & 1 );
                break;
            }
            case BlackBishop: 
                x = EvenBits(k<<1);
                y = EvenBits(k);
                x += y&1;
                break;
            case Bishop:
            case Knight: {
                // Make topmost split partition by color of square.
                unsigned a = (k<<1)%N_Square;
                x = EvenBits(a);
                y = EvenBits(a>>1);
                x += (y&1) ^ (k >= N_Square/2);
                break;
            }
            case Rook: {
                x = EvenBits(k);
                y = EvenBits(k>>1);
                break;
            }
        }
        setLeaf( k, p, x, y );
    }

    // Verify that pieces are on distinct squares.
    SquareSet tally;
    for( unsigned long j=0; j<myLeafCount; ++j ) {
        Node& n = myRoot[j+myLeafCount-1];
        tally.add( n.pos.x, n.pos.y );
    }
    assert(tally.count()==myLeafCount);

    // Initialize internal nodes
    for( unsigned long k=myLeafCount-1; k-->0; ) {
        for( size_t i=0; i<2; ++i ) 
            myRoot[k].child[i] = &myRoot[2*k+1+i];
        Node& n = myRoot[k];
        assert( n.child[0]->level==n.child[1]->level || p==Queen && myLeafCount==10 && k<=1 );
        n.level = std::max( n.child[0]->level, n.child[1]->level ) + 1;
        n.weakAttacks = n.child[0]->weakAttacks;
        n.weakAttacks |= n.child[1]->weakAttacks;
        n.weight = n.weakAttacks.count();
        n.leftmost = n.child[0]->leftmost;
        n.rightmost = n.child[1]->rightmost;
    }
    assert( myRoot[0].level<=6 );
    assert( myRoot[0].level>=4 );
}

//! Number of pieces that we are trying to cover the board with.
static const long N_Piece = 8;

//! An arrangement of pieces on the board, possibly with multiple pieces on a square. 
class Board {
    //! kth element describes kind of kth piece.
    static Piece thePiece[N_Piece];
    //! kth element describes possible locations of kth piece
    Node* myPositions[N_Piece];    
public:
    //! Construct undefined board.
    Board() {}
    //! Construct arrangement that splits axis i and takes the "whichChild" half of it.
    Board( const Board& p, int i, int whichChild ) {
        assert( whichChild==0||whichChild==1 ); 
        assert( p.myPositions[i]->weight>1 ); // Axis i must be splittable
        *this = p;
        myPositions[i] = p.myPositions[i]->child[whichChild];
    }
    //! Return set of squares weakly attacked by this arrangement.
    SquareSet weakAttacks() const {
        SquareSet a;
        for( long i=0; i<N_Piece; ++i ) 
            a |= myPositions[i]->weakAttacks; 
        return a;
    }
    //! Return set of squares strongly attacked by this arrangement of pieces.
    SquareSet strongAttacks() const;
    //! Return true if board has square with more than one piece
    bool hasSuperposition() const;  
    //! Return true if search should reject this arrangment of pieces.
    bool reject() const;
    //! Choose axis to split.
    /** Return index of axis if splittable; -1 -> no axis is splittable. */
    int chooseAxis() const;
    //! Initialize ith axis
    void initializeAxis( size_t i, Piece p, const Tree& t ) {
        assert( i<N_Piece );
        myPositions[i] = t.root();
        thePiece[i] = p;
    }
    //! Return true if each piece has a unique position.
    bool isLeaf() const {
        for( size_t i=0; i<N_Piece; ++i ) 
            if( myPositions[i]->level > 0 )
                return false;
        return true;
    }
    friend std::ostream& operator<<( std::ostream& o, const Board& b );
};

Piece Board::thePiece[N_Piece];

bool Board::hasSuperposition() const {
    SquareSet o;
    for( size_t i=0; i<N_Piece; ++i ) {
        unsigned x = myPositions[i]->pos.x;
        unsigned y = myPositions[i]->pos.y;
        if( o.has(x,y) )
            return true;
        o.add(x,y);
    }
    return false;
}

SquareSet Board::strongAttacks() const {
    SquareSet o;
    for( size_t i=0; i<N_Piece; ++i )
        o.add(myPositions[i]->pos.x,myPositions[i]->pos.y);
    SquareSet s;
    for( size_t i=0; i<N_Piece; ++i ) 
        s |= SquareSet(thePiece[i],myPositions[i]->pos.x,myPositions[i]->pos.y,o);
    return s;
}

//! Print a board
std::ostream& operator<<( std::ostream& o, const Board& b ) {
    if( b.isLeaf() ) {
        // Print the board arrangement
        for( long y=BoardSize; --y>=0; ) {
            int extra[N_Piece];
            int* e=extra;
            int columns=0;
            for( long x=0; x<BoardSize; ++x ) {
                char filled = '\0';
                for( int i=0; i<N_Piece; ++i ) {
                    assert( b.myPositions[i]->level==0 );
                    Position pos = b.myPositions[i]->pos;
                    if( pos.x==x && pos.y==y ) {
                        if( filled ) {
                            // Already printed another piece for current square, so defer printing piece i.
                            *e++ = i;
                        } else {
                            o << Symbol[Board::thePiece[i]];
                            filled = true;
                        }
                    }
                }
                if( filled ) 
                    ++columns;
                else
                    o << char(((x^y)&1) ? '.' : ' '); 
            }
            // Print deferred pieces
            for( int* f=extra; f!=e; ++f ) {
                int i = *f;
                o << " " << Symbol[Board::thePiece[i]];
                if( columns>1 ) {
                    // Need to distinguish where extra pieces go
                    o << b.myPositions[i]->pos.x;
                } 
            } 
            o << std::endl;
        }
    } else {
        // Some pieces are replicated.  
        // Just print the tree level of each piece.
        o << "[";
        for( size_t j=0; j<N_Piece; ++j ) {
            if( j>0 ) o << " ";
            o << int(b.myPositions[j]->level);
        }
        o << "]";
    }
    return o;
}

bool Board::reject() const {
    // Reject if some square is not weakly attacked.
    if( !weakAttacks().isAll() )
        return true;
    // Clip rooks
    if( myPositions[0]->leftmost > myPositions[7]->rightmost )
        return true;
    // Clip knights
    if( myPositions[1]->leftmost > myPositions[6]->rightmost )
        return true;
#if BISHOPS_CAN_BE_ON_SAME_COLOR
    if( myPositions[2]->leftmost > myPositions[5]->rightmost )
        return true;
#endif /* BISHOPS_CAN_BE_ON_SAME_COLOR */
    return false;
}

int Board::chooseAxis() const {
    int maxWeight = -1;
    int maxI = -1;
    for( int i=0; i<N_Piece; ++i ) {
        if( myPositions[i]->weight>1 ) {
            int w = 0;
            for( int c=0; c<2; ++c )  {
                if( !Board(*this,i,c).reject() ) {
                    assert( myPositions[i]->child[c]->weight>0 );
                    w += myPositions[i]->child[c]->weight;
                } 
            }
            if( w==0 ) {
                // Found axis that splits board into two not-fully covered boards
                return i;
            }
            if( w>maxWeight ) {
                maxWeight = w;
                maxI = i;
            }
        }
    }
    return maxI;
}

#if PARALLEL
cilk::opadd_reducer<int> WeakCount;
cilk::ostream_reducer<char> Output(std::cout);
#else
int WeakCount;
std::ostream& Output(std::cout);
#endif

//! Recursively search whether Board b or sub-boards are fully covered.
void Search( const Board& b ) {
    if( !b.reject() ) {
        // Find an axis to unfold the search space
        int i = b.chooseAxis();
        if( i<0 ) {
            // Found a weak solution
            ++WeakCount;
            if( b.strongAttacks().isAll() ) { 
                // Found a strong solution 
                if( !b.hasSuperposition() ) {
                    // Found solution with no superposition.  Print it.
                    Output << b << std::endl;
                }
            }
        } else {
            // Unfold on axis i and search both halves in parallel
            cilk_spawn Search( Board(b,i,0) );
            Search( Board(b,i,1) );
        }
    } 
    cilk_sync;
}

//! Report statistics.
void Report() {
#if PARALLEL
    int p = __cilkrts_get_nworkers();
    int w = WeakCount;
#else
    int p = 1;
    int w = WeakCount;
#endif
#if BISHOPS_CAN_BE_ON_SAME_COLOR
    int expected_w = 28790;
#else
    int expected_w = 8715;
#endif 
    if( w!=expected_w ) {
        Output << "Error: weak attacks = " << w << " but should be " << expected_w << std::endl;
        std::abort(); 
    }
    Output << "Threads: " << p << std::endl;
    Output << "Weak attacks: " << w << std::endl;
}

//! Solve the puzzle described in the paper.
void DoPuzzle() {
    // Construct folding trees for each kind of piece
    Tree king(King);
    Tree queen(Queen);
#if BISHOPS_CAN_BE_ON_SAME_COLOR
    Tree bishop(Bishop);
#else
    Tree whiteBishop(WhiteBishop);
    Tree blackBishop(BlackBishop);
#endif
    Tree knight(Knight);
    Tree rook(Rook);

    // Construct board representing fully folded search space.
    Board root;
    root.initializeAxis(0,Rook,rook);
    root.initializeAxis(1,Knight,knight);
#if BISHOPS_CAN_BE_ON_SAME_COLOR
    root.initializeAxis(2,Bishop,bishop);
    root.initializeAxis(5,Bishop,bishop);
#else
    root.initializeAxis(2,BlackBishop,blackBishop);
    root.initializeAxis(5,WhiteBishop,whiteBishop);
#endif
    root.initializeAxis(3,Queen,queen);
    root.initializeAxis(4,King,king);
    root.initializeAxis(6,Knight,knight);
    root.initializeAxis(7,Rook,rook);

    // Do search
    Search(root);
}

int main() {
    // Measure time to solve puzzle from a cold start.
    clockmark_t begin, end;
    begin = ktiming_getmark();
    DoPuzzle();
    end = ktiming_getmark();
    Report();
    double elapsed_time = ktiming_diff_sec(&begin, &end);
    Output << "Benchmark time in second: " << elapsed_time << std::endl;
    return 0;
}
