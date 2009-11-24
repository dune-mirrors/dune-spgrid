#ifndef DUNE_SPGRID_DECOMPOSITION_HH
#define DUNE_SPGRID_DECOMPOSITION_HH

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/partition.hh>

namespace Dune
{

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
