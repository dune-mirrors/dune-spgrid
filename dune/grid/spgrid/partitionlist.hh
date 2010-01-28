#ifndef DUNE_SPGRID_PARTITIONLIST_HH
#define DUNE_SPGRID_PARTITIONLIST_HH

#include <dune/common/smallobject.hh>

#include <dune/grid/common/gridenums.hh>

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

    SPPartitionList ( const Mesh &mesh )
    : head_( new Node( Partition( mesh ) ) )
    {}

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
      head_ = (other.head_ != 0 ? new Node( other.head_ ) : 0);
    }

    Iterator begin () const
    {
      return Iterator( head_ );
    }

    Iterator end () const
    {
      return Iterator( 0 );
    }

    template< PartitionIteratorType pitype >
    static This
    create ( const Mesh &localMesh, const Mesh &globalMesh, const int overlap = 0 );

  private:
    void append ( const Partition &partition )
    {
      if( head_ != 0 )
        head_->append( new Node( partition ) );
      else
        head_ = new Node( partition );
    }

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

    void append ( const This &other )
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



  // Implementation of SPPartitionList
  // ---------------------------------

  template< int dim >
  template< PartitionIteratorType pitype >
  inline typename SPPartitionList< dim >::This SPPartitionList< dim >
    ::create ( const Mesh &localMesh, const Mesh &globalMesh, const int overlap )
  {
    const MultiIndex &lbegin = localMesh.begin();
    const MultiIndex &lend = localMesh.end();
    const MultiIndex &gbegin = globalMesh.begin();
    const MultiIndex &gend = globalMesh.end();

    This list;

    if( pitype == Interior_Partition )
    {
      MultiIndex begin, end;
      for( int i = 0; i < dim; ++i )
      {
        begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
        end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
      }
      list.append( Partition( begin, end ) );
    }
    else if( pitype == InteriorBorder_Partition )
      list.append( Partition( 2*lbegin, 2*lend ) );
    else if( pitype == Overlap_Partition )
    {
      MultiIndex begin, end;
      for( int i = 0; i < dim; ++i )
      {
        begin[ i ] = 2*lbegin[ i ] + int( lbegin[ i ] != gbegin[ i ] );
        end[ i ] = 2*lend[ i ] - int( lend[ i ] != gend[ i ] );
      }
      list.append( Partition( begin, end ) );
    }
    else
      list.append( Partition( 2*lbegin, 2*lend ) );

    return list;
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
