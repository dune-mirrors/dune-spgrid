#ifndef DUNE_SPGRID_DECOMPOSITION_HH
#define DUNE_SPGRID_DECOMPOSITION_HH

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPPartition
  // -----------

  template< int dim >
  class SPPartition
  {
    typedef SPPartition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    SPPartition ( const MultiIndex &origin, const MultiIndex &width )
    : origin_( origin ),
      width_( width )
    {}

    const MultiIndex &origin () const
    {
      return origin_;
    }

    const MultiIndex &width () const
    {
      return width_;
    }

  private:
    MultiIndex origin_, width_;
  };



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
    typedef SPPartition< dimension > Partition;

  public:
    SPDecomposition ( const unsigned int size, const MultiIndex &width )
    : node_( new Node( size, Partition( MultiIndex::zero(), width ) ) )
    {}

  private:
    Node *root_;
  };


  template< int dim >
  struct SPDecomposition< dim >::Node
  {

    Node ( const Partition &partition, const unsigned int size )

    ~Node ();

  private:
    Partition partition_;
    unsigned int size_;
    Node *left_, *right_;
  };


  template< int dim >
  SPDecomposition< dim >::Node::Node ( const Partition &partition, const unsigned int size )
  : partition_( partition ),
    size_( size ),
    left_( 0 ),
    right_( 0 )
  {
    if( size_ > 1 )
    {
      const MultiIndex &width = partition_.width();
      const int dir = argmax( width );

      MultiIndex leftWidth = width;
      leftWidth[ dir ] = ((size/2) * width[ dir ]) / size;
      MultiIndex rightWidth = width;
      rightWidth[ dir ] -= leftWidth[ dir ];

      MultiIndex leftOrigin = partition_.origin();
      MultiIndex rightOrigin = partition_.origin();
      rightOrigin[ dir ] += leftWidth[ dir ];

      left_ = new Node( size/2, Partition( leftOrigin, leftWidth ) );
      right_ = new Node( size - size/2, Partition( rightOrigin, rightWidth ) );
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
