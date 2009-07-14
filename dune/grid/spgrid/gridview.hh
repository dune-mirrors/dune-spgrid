#ifndef DUNE_SPGRID_GRIDVIEW_HH
#define DUNE_SPGRID_GRIDVIEW_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridview.hh>

#include <dune/grid/spgrid/indexset.hh>
#include <dune/grid/spgrid/intersection.hh>
#include <dune/grid/spgrid/intersectioniterator.hh>

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

    typedef SPIndexSet< Grid > IndexSet;
    typedef Dune::Intersection< Grid, SPIntersection > Intersection;
    typedef Dune::IntersectionIterator< Grid, SPIntersectionIterator, SPIntersection >
      IntersectionIterator;

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    static const bool conforming = Capabilibies::isLevelwiseConforming< Grid >::v;

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
        typedef Dune::LevelIterator< codim, pit, Grid, SPIterator > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };
  };



  // SPLeqafGridViewTraits
  // ---------------------

  template< class G, PartitionIteratorType pitype >
  struct SPLeafGridViewTraits
  {
    typedef SPGridView< SPLeafGridViewTraits< G, pitype > > GridViewImp;

    typedef typename remove_const< G >::type Grid;

    typedef SPIndexSet< Grid > IndexSet;
    typedef Dune::Intersection< Grid, SPIntersection > Intersection;
    typedef Dune::IntersectionIterator< Grid, SPIntersectionIterator, SPIntersection >
      IntersectionIterator;

    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    static const bool conforming = Capabilibies::isLeafwiseConforming< Grid >::v;

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
        typedef Dune::LeafIterator< codim, pit, Grid, SPIterator > Iterator;
      };

      typedef typename Partition< pitype >::Iterator Iterator;
    };
  };



  // SPGridView
  // ----------

  template< class ViewTraits >
  class SPLevelGridView
  {
  public:
    const Grid &grid () const
    {
      return *grid_;
    }

    const IndexSet &indexSet () const
    {
      return *indexSet_;
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
    typename Codim< codim >::Iterator begin () const
    {
      // ...
    }

    template< int codim >
    typename Codim< codim >::Iterator end () const
    {
      // ...
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator begin () const
    {
      // ...
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::Iterator end () const
    {
      // ...
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      return IntersectionIteratorImpl( *this, 0 );
    }

    IntersectionIterator ibegin ( const typename Codim< 0 >::Entity &entity ) const
    {
      return IntersectionIteratorImpl( *this, GridLevel::Cube::numFaces );
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

  private:
    const Grid *grid_;
    const IndexSet *indexSet_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDVIEW_HH
