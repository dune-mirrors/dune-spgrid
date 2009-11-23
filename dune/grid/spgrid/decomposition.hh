#ifndef DUNE_SPGRID_DECOMPOSITION_HH
#define DUNE_SPGRID_DECOMPOSITION_HH

#include <dune/common/iostream.hh>
#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

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

    template< class ctype, SPRefinementStrategy strategy >
    SPPartition ( const This &other, const SPRefinement< ctype, dimension, strategy > &refinement )
    : origin_( other.origin_ )
    {
      for( int i = 0; i < dimension; ++i )
        width_[ i ] = refinement.factor( i ) * other.width_[ i ];
    }

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


  template< int dim >
  inline std::ostream &
  operator<< ( std::ostream &out, const SPPartition< dim > &partition )
  {
    return out << "[ " << partition.origin() << ", " << (partition.origin() + partition.width()) << " [";
  }



  // SPDecomposition
  // ---------------

  template< int dim >
  class SPDecomposition
  {
    typedef SPDecomposition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;
    typedef SPPartition< dimension > Partition;

  private:
    struct Node
    {
      Node ( const Partition &partition, const unsigned int size );
      ~Node ();

      const Partition &interiorPartition ( const unsigned int rank ) const;
      const Partition &allPartition () const;

      unsigned int size () const;

    private:
      Partition partition_;
      unsigned int size_;
      Node *left_, *right_;
    };

  public:
    SPDecomposition ( const MultiIndex &width, const unsigned int size,
                      const unsigned int periodic = 0 )
    : root_( Partition( MultiIndex::zero(), width ), size ),
      periodic_( periodic )
    {}

    Partition partition ( const unsigned int rank, const int overlap ) const;

    bool periodic ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return ((periodic_ & (1 << i)) != 0);
    }

    unsigned int size () const
    {
      return root_.size();
    }

  private:
    Node root_;
    unsigned int periodic_;
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

      left_ = new Node( Partition( leftOrigin, leftWidth ), size/2 );
      right_ = new Node( Partition( rightOrigin, rightWidth ), size - size/2 );
    }
  }


  template< int dim >
  SPDecomposition< dim >::Node::~Node ()
  {
    delete left_;
    delete right_;
  }


  template< int dim >
  const typename SPDecomposition< dim >::Partition &
  SPDecomposition< dim >::Node::interiorPartition ( const unsigned int rank ) const
  {
    assert( rank < size_ );
    if( size_ > 1 )
    {
      assert( (left_ != 0) && (right_ != 0) );
      if( rank < size_/2 )
        return left_->interiorPartition( rank );
      else
        return right_->interiorPartition( rank - size_/2 );
    }
    else
      return partition_;
  }


  template< int dim >
  inline const typename SPDecomposition< dim >::Partition &
  SPDecomposition< dim >::Node::allPartition () const
  {
    return partition_;
  }


  template< int dim >
  inline unsigned int SPDecomposition< dim >::Node::size () const
  {
    return size_;
  }


  template< int dim >
  inline typename SPDecomposition< dim >::Partition
  SPDecomposition< dim >::partition ( const unsigned int rank, const int overlap ) const
  {
    const Partition &interior = root_.interiorPartition( rank );
    const Partition &allPartition = root_.allPartition();

    const MultiIndex &allOrigin = allPartition.origin();
    const MultiIndex &allWidth = allPartition.width();

    MultiIndex origin = interior.origin();
    MultiIndex width = interior.width();
    for( int i = 0; i < dimension; ++i )
    {
      origin[ i ] -= overlap;
      width[ i ] += 2*overlap;

      if( !periodic( i ) )
      {
        const int lMargin = origin[ i ] - allOrigin[ i ];
        if( lMargin < 0 )
        {
          origin[ i ] -= lMargin;
          width[ i ] += lMargin;
        }

        const int rMargin = (allOrigin[ i ] + allWidth[ i ]) - (origin[ i ] + width[ i ]);
        if( rMargin < 0 )
          width[ i ] += rMargin;
      }
    }
    return Partition( origin, width );
  }

}

#endif // #ifndef DUNE_SPGRID_DECOMPOSITION_HH
