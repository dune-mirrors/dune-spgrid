#ifndef DUNE_SPGRID_SPHIERARCHICSEARCH_HH
#define DUNE_SPGRID_SPHIERARCHICSEARCH_HH

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid;

  template< class Grid, class IndexSet >
  class HierarchicSearch;



  // HierarchicSearch for SPGrid
  // ---------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  class HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
  {
    typedef SPGrid< ct, dim, strategy, Comm > Grid;

  public:
    typedef FieldVector< ct, dim > GlobalVector;
    typedef FieldVector< ct, dim > LocalVector;
    
    typedef typename Grid::template Codim< 0 >::Entity Entity;
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    HierarchicSearch ( const Grid &grid, const IndexSet &indexSet );

    EntityPointer findEntity ( const GlobalVector &global ) const;
    
  private:
    EntityPointer hFindEntity ( const Entity &e, const GlobalVector &global ) const;

    const Grid &grid_;
    const IndexSet &indexSet_;
  };



  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  inline HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
    ::HierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : grid_( grid ),
    indexSet_( indexSet )
  {}


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  inline typename HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >::EntityPointer
  HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
    ::findEntity ( const GlobalVector &global ) const
  {
    typedef typename Grid::template Codim< 0 >::LevelIterator LevelIterator;
    const LevelIterator end = grid_.template lend< 0 >( 0 );
    for( LevelIterator it = grid_.template lbegin< 0 >( 0 ); it != end; ++it )
    {
      const Entity &e = *it;
      const typename Entity::Geometry &geo = e.geometry();

      LocalVector local = geo.local( global );
      if( !GenericReferenceElements< double, dim >::general( geo.type() ).checkInside( local ) )
        continue;

      if( indexSet_.contains( *it ) )
        return EntityPointer( it );
      else
        return hFindEntity( *it, global );
    }
    DUNE_THROW( GridError, "Coordinate " << global << " is outside the grid." );
  }



  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  inline typename HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >::EntityPointer
  HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
    ::hFindEntity ( const Entity &e, const GlobalVector &global ) const
  {
    typedef typename Grid::template Codim< 0 >::HierarchicIterator HierarchicIterator;
    const HierarchicIterator end = e.hend( e.level()+1 );
    for( HierarchicIterator it = e.hbegin( e.level()+1 ); it != end; ++it )
    {
      const Entity &child = *it;
      LocalVector local = child.geometry().local( global );
      if( GenericReferenceElements< double, dim >::general( child.type() ).checkInside( local ) )
      {
        if( indexSet_.contains( *it ) )
          return EntityPointer( it );
        else
          return hFindEntity( *it, global );
      }
    }
    DUNE_THROW( Exception, "Unexpected internal Error" );
  }

}

#endif // #infdef DUNE_SPGRID_SPHIERARCHICSEARCH_HH
