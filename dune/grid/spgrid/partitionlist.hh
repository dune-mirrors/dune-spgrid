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
    typedef typename Partition::Mesh Mesh;
    
    struct Iterator;

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
      head_ = new Node( other.head_ );
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
    static This create ( const Mesh &localMesh, const Mesh &globalMesh );

  private:
    Node *head_;
  };


  template< int dim >
  class SPPartitionList< dim >::Node
  : public SmallObject
  {
    friend struct SPPartitionList< dim >::Iterator;

  public:
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

  private:
    Partition partition_;
    Node *next_;
  };


  template< int dim >
  struct SPPartitionList< dim >::Iterator
  {
    explicit Iterator ( const Node *node )
    : node_( node )
    {}

    Iterator &operator++ ()
    {
      assert( !!(*this) );
      node_ = node_->next_;
      return *this;
    }

    bool operator== ( const Iterator &other ) const
    {
      return (node_ == other.node_);
    }

    bool operator!= ( const Iterator &other ) const
    {
      return (node_ != other.node_);
    }

    bool operator! () const
    {
      return (node_ == 0 );
    }

    const Partition &operator* () const
    {
      assert( !!(*this) );
      return node_->partition_;
    }

    const Partition *operator-> () const
    {
      assert( !!(*this) );
      return &(node_->partition_);
    }

  private:
    const Node *node_;
  };



  // Implementation of SPPartitionList
  // ---------------------------------

  template< int dim >
  template< PartitionIteratorType pitype >
  inline typename SPPartitionList< dim >::This
  SPPartitionList< dim >::create ( const Mesh &localMesh, const Mesh &globalMesh )
  {
    return This( localMesh );
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
