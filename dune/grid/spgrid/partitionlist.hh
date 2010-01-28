#ifndef DUNE_SPGRID_PARTITIONLIST_HH
#define DUNE_SPGRID_PARTITIONLIST_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/exceptions.hh>

#include <dune/grid/spgrid/partition.hh>

namespace Dune
{

  // SPPartitionList
  // ---------------

  template< int dim >
  class SPPartitionList
  {
    typedef SPPartitionList< dim > This;

    class Node;

  public:
    typedef SPPartition< dim > Partition;

    typedef typename Partition::MultiIndex MultiIndex;
    typedef typename Partition::Mesh Mesh;
    
    struct Iterator;

    SPPartitionList ()
    : head_( 0 )
    {}

  private:
    SPPartitionList ( const Mesh &mesh )
    : head_( new Node( Partition( mesh ) ) )
    {}

  public:
    SPPartitionList ( const This &other )
    : head_( other.head_ != 0 ? new Node( *other.head_ ) : 0 )
    {}

    ~SPPartitionList ()
    {
      delete head_;
    }

    This &operator= ( const This &other )
    {
      delete head_;
      head_ = (other.head_ != 0 ? new Node( *other.head_ ) : 0);
      return *this;
    }

    This &operator += ( const Partition &partition );

    Iterator begin () const
    {
      return Iterator( head_ );
    }

    Iterator end () const
    {
      return Iterator( 0 );
    }

  private:
    Node *head_;
  };



  // SPPartitionList::Node
  // ---------------------

  template< int dim >
  struct SPPartitionList< dim >::Node
  : public SmallObject
  {
    explicit Node ( const Partition &partition )
    : partition_( partition ),
      next_( 0 )
    {}

    Node ( const Node &other )
    : partition_( other.partition_ ),
      next_( other.next_ != 0 ? new Node( *other.next_ ) : 0 )
    {}

    ~Node ()
    {
      delete next_;
    }

    void append ( Node *other )
    {
      if( next_ != 0 )
        next_->append( other );
      else
        next_ = other;
    }

    const Partition &partition () const
    {
      return partition_;
    }

    const Node *next () const
    {
      return next_;
    }

  private:
    Partition partition_;
    Node *next_;
  };



  // SPPartitionList::Iterator
  // -------------------------

  template< int dim >
  struct SPPartitionList< dim >::Iterator
  {
    explicit Iterator ( const Node *node )
    : node_( node )
    {}

    Iterator &operator++ ()
    {
      assert( node_ != 0 );
      node_ = node_->next();
      return *this;
    }

    operator bool () const
    {
      return (node_ != 0);
    }

    bool operator== ( const Iterator &other ) const
    {
      return (node_ == other.node_);
    }

    bool operator!= ( const Iterator &other ) const
    {
      return (node_ != other.node_);
    }

    const Partition &operator* () const
    {
      assert( node_ != 0 );
      return node_->partition();
    }

    const Partition *operator-> () const
    {
      assert( node_ != 0 );
      return &(node_->partition());
    }

  private:
    const Node *node_;
  };



  // SPPartitionPool
  // ---------------

  template< int dim >
  class SPPartitionPool
  {
    typedef SPPartitionPool< dim > This;

  public:
    typedef SPPartitionList< dim > PartitionList;

    typedef typename PartitionList::Partition Partition;
    typedef typename PartitionList::MultiIndex MultiIndex;
    typedef typename PartitionList::Mesh Mesh;
    
    SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                      const MultiIndex &overlap, unsigned int periodic = 0 );

    template< PartitionIteratorType pitype >
    const PartitionList &get () const;

  private:
    static Partition
    removeBorder ( const Mesh &localMesh, const Mesh &globalMesh );

    PartitionList interior_;
    PartitionList interiorBorder_;
    PartitionList overlap_;
    PartitionList overlapFront_;
    PartitionList all_;
    PartitionList ghost_;
  };



  // Implementation of SPPartitionList
  // ---------------------------------

  template< int dim >
  inline typename SPPartitionList< dim >::This &
  SPPartitionList< dim >::operator+= ( const Partition &partition )
  {
    if( head_ != 0 )
      head_->append( new Node( partition ) );
    else
      head_ = new Node( partition );
    return *this;
  }



  // Implementation of SPPartitionPool
  // ---------------------------------

  template< int dim >
  inline SPPartitionPool< dim >
    ::SPPartitionPool ( const Mesh &localMesh, const Mesh &globalMesh,
                        const MultiIndex &overlap, unsigned int periodic )
  {
    interior_ += removeBorder( localMesh, globalMesh );
    interiorBorder_ += Partition( localMesh );

    Mesh overlapMesh = localMesh.grow( overlap );
    overlap_ += removeBorder( globalMesh.intersect( overlapMesh ), globalMesh );
    overlapFront_ += Partition( globalMesh.intersect( overlapMesh ) );
    
    Mesh allMesh = localMesh.grow( overlap ).grow( 1 );
    all_ += Partition( globalMesh.intersect( overlapMesh ) );
  }


  template< int dim >
  template< PartitionIteratorType pitype >
  inline const typename SPPartitionPool< dim >::PartitionList &
  SPPartitionPool< dim >::get () const
  {
    switch( pitype )
    {
    case Interior_Partition:
      return interior_;

    case InteriorBorder_Partition:
      return interiorBorder_;

    case Overlap_Partition:
      return overlap_;

    case OverlapFront_Partition:
      return overlapFront_;

    case All_Partition:
      return all_;

    case Ghost_Partition:
      return ghost_;

    default:
      DUNE_THROW( GridError, "No such PartitionIteratorType." );
    }
  }


  template< int dim >
  inline typename SPPartitionPool< dim >::Partition
  SPPartitionPool< dim >::removeBorder ( const Mesh &localMesh, const Mesh &globalMesh )
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh.begin();
    const MultiIndex &gend = globalMesh.end();

    MultiIndex begin, end;
    for( int i = 0; i < dim; ++i )
    {
      begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
      end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
    }
    return Partition( begin, end );
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
