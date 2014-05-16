#ifndef DUNE_SPGRID_HIERARCHICSEARCH_HH
#define DUNE_SPGRID_HIERARCHICSEARCH_HH

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/indexidset.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, class Grid >
  class SPEntityPointer;

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid;

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
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    static const int dimension = Grid::dimension;

    typedef FieldVector< ctype, dimension > GlobalVector;

    SPBasicHierarchicSearch ( const Grid &grid )
    : grid_( grid )
    {}

    EntityPointer findEntity ( const GlobalVector &global, int level ) const
    {
      typedef SPEntityPointer< 0, const Grid > EntityPointerImpl;
      typedef typename Grid::GridLevel GridLevel;
      typedef typename GridLevel::PartitionList PartitionList;

      assert( grid_.domain().contains( global ) );
      const GridLevel &gridLevel = grid_.gridLevel( level );
      const PartitionList &partitionList = gridLevel.template partition< All_Partition >();

      const GlobalVector y = global - grid_.domain().cube().origin();
      GlobalVector z;
      gridLevel.template geometryCache< 0 >( (1 << dimension) - 1 ).jacobianInverseTransposed().mv( y, z );

      typename GridLevel::MultiIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*int( z[ i ] ) + 1;

      const typename PartitionList::Partition *partition = partitionList.findPartition( id );
      if( partition )
        return EntityPointerImpl( gridLevel, id, partition->number() );
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
          return EntityPointerImpl( gridLevel, id, leftPartition->number() );
        else  
          DUNE_THROW( GridError, "Coordinate " << global << " is outside the grid." );
      }
    }
    
  private:
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
    typedef typename Base::EntityPointer EntityPointer;

    using Base::dimension;

    typedef typename Base::GlobalVector GlobalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid ),
      indexSet_( indexSet )
    {}

    EntityPointer findEntity ( const GlobalVector &global ) const
    {
      EntityPointer ep = Base::findEntity( global, 0 );
      return (indexSet_.contains( *ep ) ? ep : hFindEntity( *ep, global ));
    }
  
  private:
    typedef GlobalVector LocalVector;
    typedef typename Entity::HierarchicIterator HierarchicIterator;

    EntityPointer hFindEntity ( const Entity &e, const GlobalVector &global ) const
    {
      // To Do: This method should use the Cartesian structure, too
      const HierarchicIterator end = e.hend( e.level()+1 );
      for( HierarchicIterator it = e.hbegin( e.level()+1 ); it != end; ++it )
      {
        const Entity &child = *it;
        LocalVector local = child.geometry().local( global );
        if( ReferenceElements< ctype, dimension >::cube().checkInside( local ) )
          return (indexSet_.contains( *it ) ? EntityPointer( it ) : hFindEntity( *it, global ));
      }
      DUNE_THROW( Exception, "Unexpected internal Error" );
    }

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
    typedef typename Base::EntityPointer EntityPointer;

    typedef typename Base::GlobalVector GlobalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid ),
      indexSet_( indexSet )
    {}

    EntityPointer findEntity ( const GlobalVector &global ) const
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

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  class HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
  : public SPHierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
  {
    typedef SPHierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet > Base;
    typedef SPGrid< ct, dim, strategy, Comm > Grid;

  public:
    typedef typename Base::EntityPointer EntityPointer;
    typedef typename Base::GlobalVector GlobalVector;

    HierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
    : Base( grid, indexSet )
    {}

    EntityPointer findEntity ( const GlobalVector &global ) const
    {
      return Base::findEntity( global );
    }

    template< PartitionIteratorType pitype >
    EntityPointer findEntity ( const GlobalVector &global ) const
    {
      const EntityPointer entityPointer = Base::findEntity( global );
      if( !contains< pitype >( entityPointer->partitionType() ) )
        DUNE_THROW( GridError, "Coordinate " << global << " does not belong to partition " << pitype );
      return entityPointer;
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
