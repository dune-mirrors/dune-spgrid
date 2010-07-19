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

  protected:
    class Node;

  public:
    typedef SPPartition< dim > Partition;

    typedef typename Partition::MultiIndex MultiIndex;
    typedef typename Partition::Mesh Mesh;
    
    struct Iterator;

    SPPartitionList ()
    : head_( 0 )
    {}

    SPPartitionList ( const This &other )
    : head_( other.head_ != 0 ? new Node( *other.head_ ) : 0 )
    {}

    ~SPPartitionList ()
    {
      delete head_;
    }

    const This &operator= ( const This &other )
    {
      delete head_;
      head_ = (other.head_ != 0 ? new Node( *other.head_ ) : 0);
      return *this;
    }

    const This &operator+= ( const Partition &partition );

    Iterator begin () const;
    Iterator end () const;

    bool contains ( const MultiIndex &id, const unsigned int number ) const;
    const Partition *findPartition ( const MultiIndex &id ) const;
    int volume () const;

    bool empty () const;
    unsigned int size () const;

  protected:
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
  inline const typename SPPartitionList< dim >::This &
  SPPartitionList< dim >::operator+= ( const Partition &partition )
  {
    if( head_ != 0 )
      head_->append( new Node( partition ) );
    else
      head_ = new Node( partition );
    return *this;
  }


  template< int dim >
  inline typename SPPartitionList< dim >::Iterator
  SPPartitionList< dim >::begin () const
  {
    return Iterator( head_ );
  }


  template< int dim >
  inline typename SPPartitionList< dim >::Iterator
  SPPartitionList< dim >::end () const
  {
    return Iterator( 0 );
  }


  template< int dim >
  inline bool
  SPPartitionList< dim >
    ::contains ( const MultiIndex &id, const unsigned int number ) const
  {
    const Partition *partition = findPartition( id );
    assert( !partition || (partition->number() == number) );
    return bool( partition );
  }


  template< int dim >
  inline const typename SPPartitionList< dim >::Partition *
  SPPartitionList< dim >::findPartition ( const MultiIndex &id ) const
  {
    for( const Node *it = head_; it != 0; it = it->next() )
    {
      if( it->partition().contains( id ) )
        return &(it->partition());
    }
    return 0;
  }


  template< int dim >
  inline int SPPartitionList< dim >::volume () const
  {
    int volume = 0;
    for( const Node *it = head_; it != 0; it = it->next() )
      volume += it->partition().volume();
    return volume;
  }


  template< int dim >
  inline bool SPPartitionList< dim >::empty () const
  {
    return (head_ == 0);
  }


  template< int dim >
  inline unsigned int SPPartitionList< dim >::size () const
  {
    unsigned int size = 0;
    for( const Node *it = head_; it != 0; it = it->next() )
      ++size;
    return size;
  }



  // Auxilliary Functions for SPPartitionList
  // ----------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out,
               const SPPartitionList< dim > &partition )
  {
    typedef typename SPPartitionList< dim >::Iterator Iterator;
    std::string separator = "";
    for( Iterator it = partition.begin(); it; ++it )
    {
      out << separator << *it;
      separator = "; ";
    }
    return out;
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITIONLIST_HH
