#ifndef DUNE_SPGRID_HIERARCHICSEARCH_HH
#define DUNE_SPGRID_HIERARCHICSEARCH_HH

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/entity.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  class SPIndexSet;

  template< class Grid, class IndexSet >
  class HierarchicSearch;



  // SPBasicHierarchicSearch
  // -----------------------

  template< class Grid >
  class SPBasicHierarchicSearch
  {
  public:
    typedef typename Grid::ctype ctype;
    typedef typename Grid::template Codim< 0 >::Entity Entity;

    static const int dimension = Grid::dimension;

    typedef FieldVector< ctype, dimension > GlobalVector;

    SPBasicHierarchicSearch ( const Grid &grid )
    : grid_( grid )
    {}

    Entity findEntity ( const GlobalVector &global, int level ) const
    {
      typedef SPEntity< 0, dimension, const Grid > EntityImpl;
      typedef typename Grid::GridLevel GridLevel;
      typedef typename GridLevel::PartitionList PartitionList;

      assert( grid_.domain().contains( global ) );
      const GridLevel &gridLevel = grid_.gridLevel( level );
      const PartitionList &partitionList = gridLevel.template partition< All_Partition >();

      const GlobalVector y = global - grid_.domain().cube().origin();
      GlobalVector z;
      SPDirectionIterator< dimension, 0 > dirIt;
      gridLevel.template geometryCache< 0 >( *dirIt ).jacobianInverseTransposed().mv( y, z );

      typename GridLevel::MultiIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*int( z[ i ] ) + 1;

      const typename PartitionList::Partition *partition = partitionList.findPartition( id );
      if( partition )
        return Entity( EntityImpl( gridLevel, id, partition->number() ) );
      else
      {
        for( int i = 0; i < dimension; ++i )
        {
          // check upper bound 
          if( id[ i ] - 1 == 2*gridLevel.localMesh().bound( 1 )[ i ] ) 
            id[ i ] = 2*int( z[ i ] ) - 1;
        }
        const typename PartitionList::Partition *leftPartition = partitionList.findPartition( id );

        if( leftPartition ) 
          return Entity( EntityImpl( gridLevel, id, leftPartition->number() ) );
        else  
          DUNE_THROW( GridError, "Coordinate " << global << " is outside the grid." );
      }
    }
    
  protected:
    const Grid &grid_;
  };



  // SPHierarchicSearch
  // ------------------

  template< class Grid, class IndexSet >
  class SPHierarchicSearch
  : protected SPBasicHierarchicSearch< Grid >
  {
    typedef SPBasicHierarchicSearch< Grid > Base;

  public:
    typedef typename Base::ctype ctype;
    typedef typename Base::Entity Entity;

    using Base::dimension;

    typedef typename Base::GlobalVector GlobalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid ),
      indexSet_( indexSet )
    {}

    Entity findEntity ( const GlobalVector &global ) const
    {
      const Entity e = Base::findEntity( global, 0 );
      return (indexSet_.contains( e ) ? e : hFindEntity( e, global ));
    }
  
  private:
    typedef GlobalVector LocalVector;
    typedef typename Entity::HierarchicIterator HierarchicIterator;

    Entity hFindEntity ( const Entity &e, const GlobalVector &global ) const
    {
      // To Do: This method should use the Cartesian structure, too
      const HierarchicIterator end = e.hend( e.level()+1 );
      for( HierarchicIterator it = e.hbegin( e.level()+1 ); it != end; ++it )
      {
        const Entity &child = *it;
        LocalVector local = child.geometry().local( global );
        if( ReferenceElements< ctype, dimension >::cube().checkInside( local ) )
          return (indexSet_.contains( child ) ? child : hFindEntity( child, global ));
      }
      DUNE_THROW( Exception, "Unexpected internal Error" );
    }

  protected:
    using Base::grid_;

  private:
    const IndexSet &indexSet_;
  };



  // SPHierarchicSearch for SPIndexSet
  // ---------------------------------

  template< class Grid >
  class SPHierarchicSearch< Grid, SPIndexSet< Grid > >
  : protected SPBasicHierarchicSearch< Grid >
  {
    typedef SPBasicHierarchicSearch< Grid > Base;
    typedef SPIndexSet< Grid > IndexSet;

  public:
    typedef typename Base::Entity Entity;

    typedef typename Base::GlobalVector GlobalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
      : Base( grid ),
        indexSet_( indexSet )
    {}

    Entity findEntity ( const GlobalVector &global ) const
    {
      return Base::findEntity( global, indexSet_.gridLevel().level() );
    }

  private:
    const IndexSet &indexSet_;
  };



  // SPHierarchicSearch for IndexSet< SPIndexSet >
  // ---------------------------------------------

  template< class Grid >
  class SPHierarchicSearch< Grid, IndexSet< Grid, SPIndexSet< Grid >, typename SPIndexSet< Grid >::IndexType > >
  : public SPHierarchicSearch< Grid, SPIndexSet< Grid > >
  {
    typedef SPHierarchicSearch< Grid, SPIndexSet< Grid > > Base;
    typedef Dune::IndexSet< Grid, SPIndexSet< Grid >, typename SPIndexSet< Grid >::IndexType > IndexSet;

  public:
    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid, static_cast< const SPIndexSet< Grid > & >( indexSet ) )
    {}
  };



  // HierarchicSearch for SPGrid
  // ---------------------------

  template< class ct, int dim, template< int > class Ref, class Comm, class IndexSet >
  class HierarchicSearch< SPGrid< ct, dim, Ref, Comm >, IndexSet >
    : public SPHierarchicSearch< SPGrid< ct, dim, Ref, Comm >, IndexSet >
  {
    typedef SPHierarchicSearch< SPGrid< ct, dim, Ref, Comm >, IndexSet > Base;
    typedef SPGrid< ct, dim, Ref, Comm > Grid;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::GlobalVector GlobalVector;

    HierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid, indexSet )
    {}

    Entity findEntity ( const GlobalVector &global ) const
    {
      return Base::findEntity( global );
    }

    template< PartitionIteratorType pitype >
    Entity findEntity ( const GlobalVector &global ) const
    {
      const Entity entity = Base::findEntity( global );
      if( !contains< pitype >( entity.partitionType() ) )
        DUNE_THROW( GridError, "Coordinate " << global << " does not belong to partition " << pitype );
      return entity;
    }

  private:
    template< PartitionIteratorType pitype >
    static bool contains ( const PartitionType partitionType )
    {
      switch( pitype )
      {
      case Interior_Partition:
        return ( partitionType == InteriorEntity );

      case InteriorBorder_Partition:
        return ( partitionType == InteriorEntity || partitionType == BorderEntity );

      case Overlap_Partition:
        return ( partitionType != FrontEntity && partitionType != GhostEntity );

      case OverlapFront_Partition:
        return ( partitionType != GhostEntity );

      case All_Partition:
        return true;

      case Ghost_Partition:
        return ( partitionType == GhostEntity );

      default:
        DUNE_THROW( InvalidStateException, "wrong PartitionIteratorType: " << pitype << "." );
      }
    }
  };

} // namespace Dune

#endif // #infdef DUNE_SPGRID_HIERARCHICSEARCH_HH
