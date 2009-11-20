#ifndef DUNE_SPGRID_DECOMPOSITION_HH
#define DUNE_SPGRID_DECOMPOSITION_HH

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPDecomposition
  // ---------------

  template< int dim >
  class SPDecomposition
  {
    typedef SPDecomposition< dim > This;

    struct Node;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

  public:
    SPDecomposition ( const unsigned int size, const MultiIndex &width )
    : node_( new Node( size, width ) )
    {}

  private:
    Node *root_;
  };


  template< int dim >
  struct SPDecomposition< dim >::Node
  {

    Node ( const unsigned int size, const MultiIndex &width,
           const MultiIndex &origin = MultiIndex::zero() );

    ~Node ();

  private:
    unsigned int size_;
    MultiIndex width_, origin_;
    Node *left_, *right_;
  };


  template< int dim >
  SPDecomposition< dim >::Node
    ::Node ( const unsigned int size, const MultiIndex &width, const MultiIndex &origin )
  : size_( size ),
    width_( width ),
    origin_( origin ),
    left_( 0 ),
    right_( 0 )
  {
    if( size_ > 1 )
    {
      const int dir = argmax( width_ );

      MultiIndex leftWidth = width_;
      leftWidth[ dir ] = ((size/2) * width_[ dir ]) / size;
      MultiIndex rightWidth = width_;
      rightWidth[ dir ] -= leftWidth[ dir ];

      MultiIndex leftOrigin = origin;
      MultiIndex rightOrigin = origin;
      rightOrigin[ dir ] += leftWidth[ dir ];

      left_ = new Node( size/2, leftWidth, leftOrigin );
      right_ = new Node( size - size/2, rightWidth, rightOrigin );
    }
  }


  template< int dim >
  SPDecomposition< dim >::Node::~Node ()
  {
    delete left_;
    delete right_;
  }

}

#endif // #ifndef DUNE_SPGRID_DECOMPOSITION_HH
