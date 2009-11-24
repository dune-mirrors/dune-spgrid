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

    SPPartition ( const MultiIndex &width )
    : begin_( MultiIndex::zero() ),
      end_( width )
    {}

    template< SPRefinementStrategy strategy >
    SPPartition ( const This &other, const SPRefinement< dimension, strategy > &refinement );

  private:
    SPPartition ( const MultiIndex &begin, const MultiIndex &end )
    : begin_( begin ),
      end_( end )
    {}

  public:
    const MultiIndex &begin () const
    {
      return begin_;
    }

    const MultiIndex &end () const
    {
      return end_;
    }

    bool empty () const;

    This grow ( const int amount, const unsigned int dir = ((1 << dimension)-1) ) const;
    This intersect ( const This &other ) const;

    std::pair< This, This > split ( const int dir, const int leftWeight, const int rightWeight ) const;

    int volume () const;
    MultiIndex width () const;

  private:
    MultiIndex begin_, end_;
  };



  // Implementation of SPPartition
  // -----------------------------

  template< int dim >
  template< SPRefinementStrategy strategy >
  inline SPPartition< dim >::SPPartition ( const This &other, const SPRefinement< dimension, strategy > &refinement )
  {
    for( int i = 0; i < dimension; ++i )
    {
      const int factor = refinement.factor( i );
      begin_[ i ] = factor * other.begin_[ i ];
      end_[ i ] = factor * other.end_[ i ];
    }
  }


  template< int dim >
  inline bool SPPartition< dim >::empty () const
  {
    bool empty = false;
    for( int i = 0; i < dimension; ++i )
      empty |= (end[ i ] < begin[ i ]);
    return empty;
  }


  template< int dim >
  inline SPPartition< dim >
  SPPartition< dim >::grow ( const int amount, const unsigned int dir ) const
  {
    MultiIndex b, e;
    for( int i = 0; i < dimension; ++i )
    {
      const int s = ((dir >> i) & 1) * amount;
      b[ i ] = begin()[ i ] - s;
      e[ i ] = end()[ i ] + s;
    }
    return This( b, e );
  }


  template< int dim >
  inline SPPartition< dim >
  SPPartition< dim >::intersect ( const This &other ) const
  {
    MultiIndex b, e;
    for( int i = 0; i < dimension; ++i )
    {
      b[ i ] = std::max( begin()[ i ], other.begin()[ i ] );
      e[ i ] = std::min( end()[ i ], other.end()[ i ] );
    }
    return This( b, e );
  }


  template< int dim >
  inline std::pair< SPPartition< dim >, SPPartition< dim > >
  SPPartition< dim >::split ( const int dir, const int leftWeight, const int rightWeight ) const
  {
    const MultiIndex &lbegin = begin();
    const MultiIndex &rend = end();

    assert( (dir >= 0) && (dir < dimension) );
    const int width = (rend[ dir ] - lbegin[ dir ]);
    const int leftWidth = (leftWeight * width) / (leftWeight + rightWeight);

    MultiIndex lend = rend;
    MultiIndex rbegin = lbegin;
    rbegin[ dir ] = lend[ dir ] = lbegin[ dir ] + leftWidth;

    return std::make_pair( This( lbegin, lend ), This( rbegin, rend ) );
  }


  template< int dim >
  inline int SPPartition< dim >::volume () const
  {
    const MultiIndex &w = width();
    int volume = 1;
    for( int i = 0; i < dimension; ++i )
      volume *= w[ i ];
    return volume;
  }


  template< int dim >
  inline typename SPPartition< dim >::MultiIndex SPPartition< dim >::width () const
  {
    MultiIndex width;
    for( int i = 0; i < dimension; ++i )
      width[ i ] = std::max( end()[ i ] - begin()[ i ], 0 );
    return width;
  }



  // Auxilliary functions for SPPartition
  // ------------------------------------

  template< int dim >
  inline std::ostream &
  operator<< ( std::ostream &out, const SPPartition< dim > &partition )
  {
    return out << "[ " << partition.begin() << ", " << partition.end() << " [";
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

      Partition partition ( const int overlap = 0 ) const;
      Partition partition ( const unsigned int rank, const int overlap = 0 ) const;
      void partitions ( std::vector< Partition > &partitions, const int overlap = 0 ) const;

      unsigned int size () const;

    private:
      Partition partition_;
      unsigned int size_;
      Node *left_, *right_;
    };

  public:
    SPDecomposition ( const MultiIndex &width, const unsigned int size,
                      const unsigned int periodic = 0 )
    : root_( Partition( width ), size ),
      periodic_( periodic )
    {}

    Partition partition ( const unsigned int rank, const int overlap = 0 ) const;
    void partitions ( std::vector< Partition > &partitions, const int overlap = 0 ) const;

    unsigned int size () const
    {
      return root_.size();
    }

  private:
    Node root_;
    unsigned int periodic_;
  };



  // Implementation of SPDecomposition::Node
  // ---------------------------------------

  template< int dim >
  inline SPDecomposition< dim >::Node::Node ( const Partition &partition, const unsigned int size )
  : partition_( partition ),
    size_( size ),
    left_( 0 ),
    right_( 0 )
  {
    if( size_ > 1 )
    {
      const int leftWeight = size_/2;
      const int rightWeight = size_ - leftWeight;

      const MultiIndex &width = partition_.width();
      const std::pair< Partition, Partition > split
        = partition_.split( argmax( width ), leftWeight, rightWeight );
      left_ = new Node( split.first, leftWeight );
      right_ = new Node( split.second, rightWeight );
    }
  }


  template< int dim >
  inline SPDecomposition< dim >::Node::~Node ()
  {
    delete left_;
    delete right_;
  }


  template< int dim >
  inline typename SPDecomposition< dim >::Partition
  SPDecomposition< dim >::Node::partition ( const int overlap ) const
  {
    return partition_.grow( overlap );
  }


  template< int dim >
  inline typename SPDecomposition< dim >::Partition
  SPDecomposition< dim >::Node::partition ( const unsigned int rank, const int overlap ) const
  {
    assert( rank < size_ );
    if( size_ > 1 )
    {
      assert( (left_ != 0) && (right_ != 0) );
      if( rank < size_/2 )
        return left_->partition( rank, overlap );
      else
        return right_->partition( rank - size_/2, overlap );
    }
    else
      return partition( overlap );
  }


  template< int dim >
  inline void
  SPDecomposition< dim >::Node::partitions ( std::vector< Partition > &partitions, const int overlap ) const
  {
    if( size_ > 1 )
    {
      assert( (left_ != 0) && (right_ != 0) );
      left_->partitions( partitions, overlap );
      right_->partitions( partitions, overlap );
    }
    else
      partitions.push_back( partition( overlap ) );
  }


  template< int dim >
  inline unsigned int SPDecomposition< dim >::Node::size () const
  {
    return size_;
  }



  // Implementation of SPDecomposition
  // ---------------------------------

  template< int dim >
  inline typename SPDecomposition< dim >::Partition
  SPDecomposition< dim >::partition ( const unsigned int rank, const int overlap ) const
  {
    const Partition partition = root_.partition( rank, overlap );
    const Partition allPartition = root_.partition().grow( overlap, periodic_ );
    return partition.intersect( allPartition );
  }


  template< int dim >
  inline void
  SPDecomposition< dim >::partitions ( std::vector< Partition > &partitions, const int overlap ) const
  {
    typedef typename std::vector< Partition >::iterator Iterator;

    partitions.reserve( root_.size() );
    root_.partitions( partitions, overlap );

    const Partition allPartition = root_.partition().grow( overlap, periodic_ );
    const Iterator end = partitions.end();
    for( Iterator it = partitions.begin(); it != end; ++it )
      *it = it->intersect( allPartition );
  }

}

#endif // #ifndef DUNE_SPGRID_DECOMPOSITION_HH
