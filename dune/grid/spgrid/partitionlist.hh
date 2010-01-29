#ifndef DUNE_SPGRID_PARTITIONLIST_HH
#define DUNE_SPGRID_PARTITIONLIST_HH

#include <dune/common/smallobject.hh>

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

#if 0
  private:
    SPPartitionList ( const Mesh &mesh )
    : head_( new Node( Partition( mesh ) ) )
    {}

  public:
#endif

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

    This &operator+= ( const Partition &partition );

    Iterator begin () const
    {
      return Iterator( head_ );
    }

    Iterator end () const
    {
      return Iterator( 0 );
    }

    bool contains ( const MultiIndex &id ) const
    {
      for( const Node *it = head_; it != 0; it = it->next() )
      {
        if( it->partition().contains( id ) )
          return true;
      }
      return false;
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

}

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
