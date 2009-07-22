#ifndef DUNE_SPGRID_GRIDVIEW_HH
#define DUNE_SPGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridview.hh>

#include <dune/grid/spgrid/capabilities.hh>
#include <dune/grid/spgrid/indexset.hh>
#include <dune/grid/spgrid/intersection.hh>
#include <dune/grid/spgrid/intersectioniterator.hh>
#include <dune/grid/spgrid/iterator.hh>

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
    };
  };



  // SPLeqafGridViewTraits
  // ---------------------

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
    };
  };



  // SPGridView
  // ----------

  template< class ViewTraits >
  class SPGridView
  {
    typedef SPGridView< ViewTraits > This;

    template< class, int > friend class SPGrid;
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

    SPGridView ()
    : indexSet_( new IndexSetPair )
    {
      indexSet_->second = 1;
    }

    explicit SPGridView ( const GridLevel &gridLevel )
    : indexSet_( new IndexSetPair )
    {
      indexSet_->first.update( gridLevel );
      indexSet_->second = 1;
    }

  public:
    template< class VT >
    SPGridView ( const SPGridView< VT > &other )
    : indexSet_( other.indexSet_ )
    {
      ++indexSet_->second;
    }

    SPGridView ( const This &other )
    : indexSet_( other.indexSet_ )
    {
      ++indexSet_->second;
    }

    ~SPGridView ()
    {
      if( --indexSet_->second == 0 )
        delete indexSet_;
    }

    This &operator= ( const This &other )
    {
      ++other.indexSet_->second;
      if( --indexSet_->second == 0 )
        delete indexSet_;
      indexSet_ = other.indexSet_;
      return *this;
    }

    const Grid &grid () const
    {
      return gridLevel().grid();
    }

    const IndexSet &indexSet () const
    {
      return indexSet_->first;
    }

    int size ( int codim ) const
    {
      return indexSet().size( codim );
    }

    int size ( const GeometryType &type ) const
    {
      return indexSet().size( type );
    }

    template< int codim >
    typename Codim< codim >::Iterator
    begin ( const unsigned int sweepDir = 0 ) const
    {
      typedef typename Codim< codim >::IteratorImpl IteratorImpl;
      typename IteratorImpl::Begin begin;
      return IteratorImpl( gridLevel(), begin, sweepDir );
    }

    template< int codim >
    typename Codim< codim >::Iterator
    end ( const unsigned int sweepDir = 0 ) const
    {
      typedef typename Codim< codim >::IteratorImpl IteratorImpl;
      typename IteratorImpl::End end;
      return IteratorImpl( gridLevel(), end, sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator
    begin ( const unsigned int sweepDir = 0 ) const
    {
      typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
      typename IteratorImpl::Begin begin;
      return IteratorImpl( gridLevel(), begin, sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator
    end ( const unsigned int sweepDir = 0 ) const
    {
      typedef typename Codim< codim >::template Partition< pitype >::IteratorImpl IteratorImpl;
      typename IteratorImpl::End end;
      return IteratorImpl( gridLevel(), end, sweepDir );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      return IntersectionIteratorImpl( entity, 0 );
    }

    IntersectionIterator iend ( const typename Codim< 0 >::Entity &entity ) const
    {
      return IntersectionIteratorImpl( entity, GridLevel::Cube::numFaces );
    }

    const CollectiveCommunication &comm () const
    {
      return grid().comm();
    }

    int overlapSize ( const int codim ) const
    {
      return 0;
    }

    int ghostSize ( const int codim ) const
    {
      return 0;
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType interface, CommunicationDirection dir ) const
    {}

    const GridLevel &gridLevel () const
    {
      return indexSet().gridLevel();
    }

    void update ( const GridLevel &gridLevel )
    {
      indexSet_->first.update( gridLevel );
    }

  private:
    IndexSetPair *indexSet_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDVIEW_HH
