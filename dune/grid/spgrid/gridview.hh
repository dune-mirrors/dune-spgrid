#ifndef DUNE_SPGRID_GRIDVIEW_HH
#define DUNE_SPGRID_GRIDVIEW_HH

#include <memory>
#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/boundarysegmentiterator.hh>
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



  // SPGridViewTraits
  // ----------------

  template< class G >
  struct SPGridViewTraits
  {
    typedef SPGridView< SPGridViewTraits< G > > GridViewImp;

    typedef typename std::remove_const< G >::type Grid;

    typedef SPIndexSet< const Grid > IndexSet;
    typedef Dune::Intersection< const Grid, SPIntersection< const Grid > > Intersection;
    typedef Dune::IntersectionIterator< const Grid, SPIntersectionIterator< const Grid >, SPIntersection< const Grid > > IntersectionIterator;

    static const bool hasBoundarySegmentIterator = true;
    typedef Dune::IntersectionIterator< const Grid, SPBoundarySegmentIterator< const Grid >, SPIntersection< const Grid > > BoundarySegmentIterator;

    typedef typename Grid::Communication Communication;
    typedef Communication CollectiveCommunication;

    static const bool conforming = true;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::Traits::template Codim< codim >::Entity Entity;

      typedef typename Grid::Traits::template Codim< codim >::Geometry Geometry;
      typedef typename Grid::Traits::template Codim< codim >::LocalGeometry LocalGeometry;

      template< PartitionIteratorType pit >
      struct Partition
      {
        typedef SPPartitionIterator< codim, const Grid > IteratorImpl;
        typedef Dune::EntityIterator< codim, const Grid, IteratorImpl > Iterator;
      };

      typedef typename Partition< All_Partition >::Iterator Iterator;
      typedef typename Partition< All_Partition >::IteratorImpl IteratorImpl;

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

    template< class, int, template< int > class, class > friend class SPGrid;
    template< class > friend class SPGridView;

  public:
    typedef typename ViewTraits::Grid Grid;
    typedef typename ViewTraits::IndexSet IndexSet;
    typedef typename ViewTraits::IntersectionIterator IntersectionIterator;
    typedef typename ViewTraits::BoundarySegmentIterator BoundarySegmentIterator;
    typedef typename ViewTraits::Communication Communication;
    typedef Communication CollectiveCommunication;

    typedef SPGridLevel< Grid > GridLevel;

    template< int codim >
    struct Codim
      : public ViewTraits::template Codim< codim >
    {};

  private:
    typedef std::pair< IndexSet, unsigned int > IndexSetPair;

    typedef SPIntersectionIterator< const Grid > IntersectionIteratorImpl;

    SPGridView () : indexSet_( new IndexSet ) {}

    explicit SPGridView ( const GridLevel &gridLevel ) : indexSet_( new IndexSet( gridLevel ) ) {}

  public:
    const Grid &grid () const;

    const IndexSet &indexSet () const;

    bool isConforming() const { return bool(ViewTraits::conforming); }

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

    template< class Entity >
    typename Codim< Entity::codimension >::SuperEntityIterator
    superEntityBegin ( const Entity &entity ) const;

    template< class Entity >
    typename Codim< Entity::codimension >::SuperEntityIterator
    superEntityEnd ( const Entity &entity ) const;

    BoundarySegmentIterator boundarySegmentBegin ( int face = 0 ) const;
    BoundarySegmentIterator boundarySegmentEnd ( int face = GridLevel::numFaces-1 ) const;

    const Communication &comm () const { return grid().comm(); }

    template< class DataHandle, class Data >
    SPCommunication< Grid, CommDataHandleIF< DataHandle, Data > >
    communicate ( CommDataHandleIF< DataHandle, Data > &data, InterfaceType iftype, CommunicationDirection dir ) const
    {
      return SPCommunication< Grid, CommDataHandleIF< DataHandle, Data > >( gridLevel(), data, iftype, dir );
    }

    const GridLevel &gridLevel () const { return indexSet().gridLevel(); }

    void update ( const GridLevel &gridLevel ) { assert( indexSet_ ); indexSet_->update( gridLevel ); }

  private:
    std::shared_ptr< IndexSet > indexSet_;
  };



  // Implementation of SPGridView
  // ----------------------------

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
    return *indexSet_;
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
    return IteratorImpl( gridLevel(), gridLevel().template partition< All_Partition >(), begin, sweepDir );
  }


  template< class ViewTraits >
  template< int codim >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::Iterator
  SPGridView< ViewTraits >::end ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::IteratorImpl IteratorImpl;
    typename IteratorImpl::End end;
    return IteratorImpl( gridLevel(), gridLevel().template partition< All_Partition >(), end, sweepDir );
  }


  template< class ViewTraits >
  template< int codim, PartitionIteratorType pitype >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::template Partition< pitype >::Iterator
  SPGridView< ViewTraits >::begin ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
    typename IteratorImpl::Begin begin;
    return IteratorImpl( gridLevel(), gridLevel().template partition< pitype >(), begin, sweepDir );
  }


  template< class ViewTraits >
  template< int codim, PartitionIteratorType pitype >
  inline typename SPGridView< ViewTraits >::template Codim< codim >::template Partition< pitype >::Iterator
  SPGridView< ViewTraits >::end ( const unsigned int sweepDir ) const
  {
    typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
    typename IteratorImpl::End end;
    return IteratorImpl( gridLevel(), gridLevel().template partition< pitype >(), end, sweepDir );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::IntersectionIterator
  SPGridView< ViewTraits >::ibegin ( const typename Codim< 0 >::Entity &entity ) const
  {
    return IntersectionIteratorImpl( entity.impl().entityInfo(), 0 );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::IntersectionIterator
  SPGridView< ViewTraits >::iend ( const typename Codim< 0 >::Entity &entity ) const
  {
    return IntersectionIteratorImpl( entity.impl().entityInfo(), GridLevel::ReferenceCube::numFaces );
  }


  template< class ViewTraits >
  template< class Entity >
  inline typename SPGridView< ViewTraits >::template Codim< Entity::codimension >::SuperEntityIterator
  SPGridView< ViewTraits >::superEntityBegin ( const Entity &entity ) const
  {
    typedef SPSuperEntityIterator< const Grid > Impl;
    return Impl( entity.impl().entityInfo(), typename Impl::Begin() );
  }


  template< class ViewTraits >
  template< class Entity >
  inline typename SPGridView< ViewTraits >::template Codim< Entity::codimension >::SuperEntityIterator
  SPGridView< ViewTraits >::superEntityEnd ( const Entity &entity ) const
  {
    typedef SPSuperEntityIterator< const Grid > Impl;
    return Impl( entity.impl().entityInfo(), typename Impl::End() );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::BoundarySegmentIterator
  SPGridView< ViewTraits >::boundarySegmentBegin ( int face ) const
  {
    typedef SPBoundarySegmentIterator< const Grid > Impl;
    return Impl( gridLevel(), face, typename Impl::Begin() );
  }


  template< class ViewTraits >
  inline typename SPGridView< ViewTraits >::BoundarySegmentIterator
  SPGridView< ViewTraits >::boundarySegmentEnd ( int face ) const
  {
    typedef SPBoundarySegmentIterator< const Grid > Impl;
    return Impl( gridLevel(), face, typename Impl::End() );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GRIDVIEW_HH
