#ifndef DUNE_SPGRID_GRID_HH
#define DUNE_SPGRID_GRID_HH

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

}

#endif // #ifndef DUNE_SPGRID_GRID_HH
