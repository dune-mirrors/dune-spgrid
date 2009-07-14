#ifndef DUNE_SPGRID_GRID_HH
#define DUNE_SPGRID_GRID_HH

#include <dune/grid/common/grid.hh>

#include <dune/grid/spgrid/intersectioniterator.hh>
#include <dune/grid/spgrid/iterator.hh>
#include <dune/grid/spgrid/idset.hh>
#include <dune/grid/spgrid/indexset.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ct, int dim >
  class SPGrid;



  // SPGridFamily
  // ------------

  template< class ct, int dim >
  struct SPGridFamily
  {
    struct Traits
    {
      typedef SPGrid< ct, dim > Grid;

      typedef ct ctype;

      static const int dimension = dim;
      static const int dimensionworld = dimworld;

      typedef Dune::CollectiveCommunication< Grid > CollectiveCommunication;

      typedef Dune::Intersection< const Grid, SPIntersection >
        LevelIntersection;
      typedef LevelIntersection LeafIntersection;

      typedef Dune::IntersectionIterator
        < const Grid, SPIntersectionIterator, SPIntersection >
        LevelIntersectionIterator;
      typedef LevelIntersectionIterator LeafIntersectionIterator;

      typedef Dune::HierarchicIterator< const Grid, SPHierarchicIterator >
        HierarchicIterator;

      typedef SPIndexSet< const Grid > LevelIndexSet;
      typedef LevelIndexSet LeafIndexSet;

      typedef SPGlobalIdSet< const Grid > GlobalIdSet;
      typedef SPLocalIdSet< const Grid > LocalIdSet;
      typedef typename GlobalIdSet::IdType GlobalIdType;
      typedef typename LocalIdSet::IdType LocalIdType;

      template< int codim >
      struct Codim
      {
        typedef Dune::Entity< codim, dimension, const Grid, SPEntity > Entity;
        typedef Dune::EntityPointer< const Grid, SPEntityPointer > EntityPointer;

        typedef Dune::Geometry
          < dimension - codim, dimensionworld, const Grid, SPGeometry >
          Geometry;
        typedef Dune::Geometry
          < dimension - codim, dimension, const Grid, SPGeometry >
          LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::LevelIterator< codim, pitype, const Grid, SPIterator >
            LevelIterator;
          typedef Dune::LeafIterator< codim, pitype, const Grid, SPIterator >
            LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune::GridView< SPGridViewTraits< const Grid, pitype > >
          LevelGridView;
        typedef LevelGridView LeafGridView;
      };
    };
  };



  template< class ct, int dim >
  class SPGrid
  : public GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim > >
  {
    typedef SPGrid< ct, dim > This;
    typedef GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim > > Base;

    typedef SPGridFamily< ct, dim > GridFamily;

  public:
    typedef typename GridFamiily::Traits Traits;

    typedef typename Traits::ctype ctype;

    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;

    typedef typename Traits::GlobalIdSet GlobalIdSet;
    typedef typename Traits::LocalIdSet LocalIdSet;

    typedef typename Traits::template Partition< pitype >::LevelGridView LevelGridView;
    typedef typename Traits::template Partition< pitype >::LeafGridView LeafGridView;

    const std::string &name () const
    {
      return name_;
    }

    int maxLevel () const
    {
      return levelViews_.size()-1;
    }

    int size ( const int level, const int codim ) const
    {
      return levelView( level ).size( codim );
    }

    int size ( const int codim ) const
    {
      return leafView().size( codim );
    }

    int size ( const int level, const GeometryType &type ) const
    {
      return levelView( level ).size( type );
    }

    int size ( const GeometryType &type ) const
    {
      return leafView().size( type );
    }

    template< PartitionIteratorType pitype >
    typename Traits::template Partition< pitype >::LevelGridView
    levelView ( const int level ) const
    {
      assert( (level >= 0) && (level <= maxLevel()) );
      return levelViews_[ level ];
    }

    template< PartitionIteratorType pitype >
    typename Traits::template Partition< pitype >::LeafGridView
    leafView ( const int level ) const
    {
      return levelView< pitype >( maxLevel() );
    }

    LevelGridView levelView ( const int level ) const
    {
      assert( (level >= 0) && (level <= maxLevel()) );
      return levelViews_[ level ];
    }

    LeafGridView leafView () const
    {
      return levelView( maxLevel() );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lbegin ( const int level ) const
    {
      return levelView< pitype >( level ).template begin< codim >();
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lend ( const int level ) const
    {
      return levelView< pitype >( level ).template end< codim >();
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lbegin ( const int level ) const
    {
      return levelView( level ).template begin< codim >();
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lend ( const int level ) const
    {
      return levelView( level ).template end< codim >();
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    lbafbegin () const
    {
      return leafView< pitype >().template begin< codim >();
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafend () const
    {
      return leafView< pitype >().template end< codim >();
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafbegin () const
    {
      return leafView().template begin< codim >();
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafend () const
    {
      return leafView().template end< codim >();
    }

    const GlobalIdSet &globalIdSet () const
    {
      return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
      return localIdSet_;
    }

    const LevelIndexSet &levelIndexSet ( const int level ) const
    {
      return levelView( level ).indexSet();
    }

    const LeafIndexSet &leafIndexSet () const
    {
      return leafView().indexSet();
    }

    void globalRefine ( int refCount )
    {
      const int maxLevel = maxLevel();
      for( int i = 0; i < refCount; ++i )
        levelViews_.push_back( LevelGridView( levelViews_[ maxLevel+i ], numDirections-1 ) );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType interface, CommunicationDirection dir,
                       int level ) const
    {
      levelView( level ).communicate( data, interface, dir );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType interface, CommunicationDirection dir ) const
    {
      leafView().communicate( data, interface, dir );
    }

    const CollectiveCommunication &comm () const
    {
      return comm_;
    }

  private:
    std::string name_;
    std::vector< LevelGridView > levelViews_;
    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;
    CollectiveCommunication comm_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRID_HH
