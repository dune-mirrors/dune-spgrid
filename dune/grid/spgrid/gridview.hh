#ifndef DUNE_SPGRID_GRIDVIEW_HH
#define DUNE_SPGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridview.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/capabilities.hh>
#include <dune/grid/spgrid/communication.hh>
#include <dune/grid/spgrid/indexset.hh>
#include <dune/grid/spgrid/intersection.hh>
#include <dune/grid/spgrid/intersectioniterator.hh>
#include <dune/grid/spgrid/iterator.hh>
#include <dune/grid/spgrid/superentityiterator.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ViewTraits >
  class SPGridView;



  // SPLevelGridViewTraits
  // ---------------------

  template< class G, PartitionIteratorType pitype >
  struct SPLevelGridViewTraits
  {
    typedef SPGridView< SPLevelGridViewTraits< G, pitype > > GridViewImp;

    typedef typename remove_const< G >::type Grid;

    typedef SPIndexSet< const Grid > IndexSet;
    typedef Dune::Intersection< const Grid, SPIntersection > Intersection;
    typedef Dune::IntersectionIterator< const Grid, SPIntersectionIterator, SPIntersection >
      IntersectionIterator;

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    static const bool conforming = Capabilities::isLevelwiseConforming< Grid >::v;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< codim >::Entity Entity;
      typedef typename Grid::Traits::template Codim< codim >::EntityPointer EntityPointer;

      typedef typename Grid::Traits::template Codim< codim >::Geometry Geometry;
      typedef typename Grid::Traits::template Codim< codim >::LocalGeometry LocalGeometry;

      template< PartitionIteratorType pit >
      struct Partition
      {
        typedef Dune::LevelIterator< codim, pit, const Grid, SPIterator > Iterator;
        typedef SPIterator< codim, pit, const Grid > IteratorImpl;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
      typedef typename Partition< pitype >::IteratorImpl IteratorImpl;

      static const bool hasSuperEntityIterator = true;
      typedef Dune::SuperEntityIterator< const Grid, SPSuperEntityIterator > SuperEntityIterator;
    };
  };



  // SPLeafGridViewTraits
  // --------------------

  template< class G, PartitionIteratorType pitype >
  struct SPLeafGridViewTraits
  {
    typedef SPGridView< SPLeafGridViewTraits< G, pitype > > GridViewImp;

    typedef typename remove_const< G >::type Grid;

    typedef SPIndexSet< const Grid > IndexSet;
    typedef Dune::Intersection< const Grid, SPIntersection > Intersection;
    typedef Dune::IntersectionIterator< const Grid, SPIntersectionIterator, SPIntersection >
      IntersectionIterator;

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    static const bool conforming = Capabilities::isLeafwiseConforming< Grid >::v;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< codim >::Entity Entity;
      typedef typename Grid::Traits::template Codim< codim >::EntityPointer EntityPointer;

      typedef typename Grid::Traits::template Codim< codim >::Geometry Geometry;
      typedef typename Grid::Traits::template Codim< codim >::LocalGeometry LocalGeometry;

      template< PartitionIteratorType pit >
      struct Partition
      {
        typedef Dune::LeafIterator< codim, pit, const Grid, SPIterator > Iterator;
        typedef SPIterator< codim, pit, const Grid > IteratorImpl;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
      typedef typename Partition< pitype >::IteratorImpl IteratorImpl;

      static const bool hasSuperEntityIterator = true;
      typedef Dune::SuperEntityIterator< const Grid, SPSuperEntityIterator > SuperEntityIterator;
    };
  };



  // SPGridView
  // ----------

  template< class ViewTraits >
  class SPGridView
  {
    typedef SPGridView< ViewTraits > This;

    template< class, int, SPRefinementStrategy, class > friend class SPGrid;
    template< class > friend class SPGridView;

  public:
    typedef typename ViewTraits::Grid Grid;
    typedef typename ViewTraits::IndexSet IndexSet;
    typedef typename ViewTraits::IntersectionIterator IntersectionIterator;
    typedef typename ViewTraits::CollectiveCommunication CollectiveCommunication;

    typedef SPGridLevel< const Grid > GridLevel;

    template< int codim >
    struct Codim
    : public ViewTraits::template Codim< codim >
    {};

  private:
    typedef std::pair< IndexSet, unsigned int > IndexSetPair;

    typedef SPIntersectionIterator< const Grid > IntersectionIteratorImpl;

    SPGridView ();

    explicit SPGridView ( const GridLevel &gridLevel );

  public:
    template< class VT >
    SPGridView ( const SPGridView< VT > &other );

    SPGridView ( const This &other );

    ~SPGridView ();

    This &operator= ( const This &other );

    const Grid &grid () const;

    const IndexSet &indexSet () const;

    int size ( int codim ) const;
    int size ( const GeometryType &type ) const;

    int overlapSize ( const int codim ) const;
    int ghostSize ( const int codim ) const;

    template< int codim >
    typename Codim< codim >::Iterator
    begin ( const unsigned int sweepDir = 0 ) const;

    template< int codim >
    typename Codim< codim >::Iterator
    end ( const unsigned int sweepDir = 0 ) const;

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator
    begin ( const unsigned int sweepDir = 0 ) const;

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator
    end ( const unsigned int sweepDir = 0 ) const;

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const;
    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const;

    template< int codim >
    typename Codim< codim >::SuperEntityIterator
    superEntityBegin ( const typename Codim< codim >::Entity &entity ) const;

    template< int codim >
    typename Codim< codim >::SuperEntityIterator
    superEntityEnd ( const typename Codim< codim >::Entity &entity ) const;

    const CollectiveCommunication &comm () const;

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType iftype, CommunicationDirection dir ) const;

    const GridLevel &gridLevel () const;

    void update ( const GridLevel &gridLevel );

  private:
    IndexSetPair *indexSet_;
  };



  // Implementation of SPGridView
  // ----------------------------

  template< class ViewTraits >
  inline SPGridView< ViewTraits >::SPGridView ()
  : indexSet_( new IndexSetPair )
  {
    indexSet_->second = 1;
  }


  template< class ViewTraits >
  inline SPGridView< ViewTraits >::SPGridView ( const GridLevel &gridLevel )
  : indexSet_( new IndexSetPair )
  {
    indexSet_->first.update( gridLevel );
    indexSet_->second = 1;
  }


  template< class ViewTraits >
  template< class VT >
  inline SPGridView< ViewTraits >::SPGridView ( const SPGridView< VT > &other )
  : indexSet_( other.indexSet_ )
  {
    ++indexSet_->second;
  }


  template< class ViewTraits >
  inline SPGridView< ViewTraits >::SPGridView ( const This &other )
  : indexSet_( other.indexSet_ )
  {
    ++indexSet_->second;
  }


  template< class ViewTraits >
  inline SPGridView< ViewTraits >::~SPGridView ()
  {
    if( --indexSet_->second == 0 )
      delete indexSet_;
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::This &
  SPGridView< ViewTraits >::operator= ( const This &other )
  {
    ++other.indexSet_->second;
    if( --indexSet_->second == 0 )
      delete indexSet_;
    indexSet_ = other.indexSet_;
    return *this;
  }


  template< class ViewTraits >
  inline const typename SPGridView< ViewTraits >::Grid &
  SPGridView< ViewTraits >::grid () const
  {
    return gridLevel().grid();
  }


  template< class ViewTraits >
  inline const typename SPGridView< ViewTraits >::IndexSet &
  SPGridView< ViewTraits >::indexSet () const
  {
    return indexSet_->first;
  }


  template< class ViewTraits >
  inline int SPGridView< ViewTraits >::size ( int codim ) const
  {
    return indexSet().size( codim );
  }


  template< class ViewTraits >
  inline int SPGridView< ViewTraits >::size ( const GeometryType &type ) const
  {
    return indexSet().size( type );
  }


  template< class ViewTraits >
  inline int SPGridView< ViewTraits >::overlapSize ( const int codim ) const
  {
    if( codim != 0 )
      DUNE_THROW( NotImplemented, "overlapSize not implemented for codim > 0." );
    int volume = gridLevel().template partition< OverlapFront_Partition >().volume();
    volume -= gridLevel().template partition< InteriorBorder_Partition >().volume();
    return volume;
  }


  template< class ViewTraits >
  inline int SPGridView< ViewTraits >::ghostSize ( const int codim ) const
  {
    if( codim != 0 )
      DUNE_THROW( NotImplemented, "ghostSize not implemented for codim > 0." );
    return gridLevel().template partition< Ghost_Partition >().volume();
  }


  template< class ViewTraits >
  template< int codim >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::Iterator
  SPGridView< ViewTraits >::begin ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::IteratorImpl IteratorImpl;
    typename IteratorImpl::Begin begin;
    return IteratorImpl( gridLevel(), begin, sweepDir );
  }


  template< class ViewTraits >
  template< int codim >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::Iterator
  SPGridView< ViewTraits >::end ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::IteratorImpl IteratorImpl;
    typename IteratorImpl::End end;
    return IteratorImpl( gridLevel(), end, sweepDir );
  }


  template< class ViewTraits >
  template< int codim, PartitionIteratorType pitype >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::template Partition< pitype >::Iterator
  SPGridView< ViewTraits >::begin ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
    typename IteratorImpl::Begin begin;
    return IteratorImpl( gridLevel(), begin, sweepDir );
  }


  template< class ViewTraits >
  template< int codim, PartitionIteratorType pitype >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::template Partition< pitype >::Iterator
  SPGridView< ViewTraits >::end ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
    typename IteratorImpl::End end;
    return IteratorImpl( gridLevel(), end, sweepDir );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::IntersectionIterator
  SPGridView< ViewTraits >::ibegin ( const typename Codim< 0 >::Entity &entity ) const
  {
    return IntersectionIteratorImpl( entity, 0 );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::IntersectionIterator
  SPGridView< ViewTraits >::iend ( const typename Codim< 0 >::Entity &entity ) const
  {
    return IntersectionIteratorImpl( entity, GridLevel::ReferenceCube::numFaces );
  }


  template< class ViewTraits >
  template< int codim >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::SuperEntityIterator
  SPGridView< ViewTraits >::superEntityBegin ( const typename Codim< codim >::Entity &entity ) const
  {
    typedef SPSuperEntityIterator< const Grid > Impl;
    return Impl( Grid::getRealImplementation( entity ), typename Impl::Begin() );
  }


  template< class ViewTraits >
  template< int codim >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::SuperEntityIterator
  SPGridView< ViewTraits >::superEntityEnd ( const typename Codim< codim >::Entity &entity ) const
  {
    typedef SPSuperEntityIterator< const Grid > Impl;
    return Impl( Grid::getRealImplementation( entity ), typename Impl::End() );
  }


  template< class ViewTraits >
  inline const typename SPGridView< ViewTraits >::CollectiveCommunication &
  SPGridView< ViewTraits >::comm () const
  {
    return grid().comm();
  }


  template< class ViewTraits >
  template< class DataHandle, class Data >
  inline void SPGridView< ViewTraits >
    ::communicate ( CommDataHandleIF< DataHandle, Data > &data,
                    InterfaceType iftype, CommunicationDirection dir ) const
  {
    typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
    typedef typename GridLevel::CommInterface Interface;

    SPCommunication< Grid, DataHandleIF > communication( gridLevel(), data );
    const Interface &interface = gridLevel().commInterface( iftype );

    const typename Interface::Iterator &end = interface.end();
    for( typename Interface::Iterator it = interface.begin(); it != end; ++it )
      communication.gather( it->rank(), it->sendList( dir ) );
    for( typename Interface::Iterator it = interface.begin(); it != end; ++it )
      communication.scatter( it->rank(), it->receiveList( dir ) );
  }


  template< class ViewTraits >
  inline const typename SPGridView< ViewTraits >::GridLevel &
  SPGridView< ViewTraits >::gridLevel () const
  {
    return indexSet().gridLevel();
  }


  template< class ViewTraits >
  inline void SPGridView< ViewTraits >::update ( const GridLevel &gridLevel )
  {
    indexSet_->first.update( gridLevel );
  }

}

#endif // #ifndef DUNE_SPGRID_GRIDVIEW_HH
