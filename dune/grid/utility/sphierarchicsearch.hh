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



  // SPHierarchicSearch
  // ------------------

  template< class Grid, class IndexSet >
  struct SPHierarchicSearch
  {
    typedef typename Grid::ctype ctype;

    static const int dimension = Grid::dimension;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef GlobalVector LocalVector;
    
    typedef typename Grid::template Codim< 0 >::Entity Entity;
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet );

    EntityPointer findEntity ( const GlobalVector &global ) const;
    
  private:
    EntityPointer hFindEntity ( const Entity &e, const GlobalVector &global ) const;

    const Grid &grid_;
    const IndexSet &indexSet_;
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
    HierarchicSearch ( const Grid &grid, const IndexSet &indexSet );
  };



  // Implementation of SPHierarchicSearch
  // ------------------------------------

  template< class Grid, class IndexSet >
  inline SPHierarchicSearch< Grid, IndexSet >
    ::SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : grid_( grid ),
    indexSet_( indexSet )
  {}


  template< class Grid, class IndexSet >
  inline typename SPHierarchicSearch< Grid, IndexSet >::EntityPointer
  SPHierarchicSearch< Grid, IndexSet >::findEntity ( const GlobalVector &global ) const
  {
    EntityPointer ep = grid_.findEntity( global, 0 );
    if( indexSet_.contains( *ep ) )
      return ep;
    else
      return hFindEntity( *ep, global );
  }



  template< class Grid, class IndexSet >
  inline typename SPHierarchicSearch< Grid, IndexSet >::EntityPointer
  SPHierarchicSearch< Grid, IndexSet >
    ::hFindEntity ( const Entity &e, const GlobalVector &global ) const
  {
    typedef typename Grid::template Codim< 0 >::HierarchicIterator HierarchicIterator;

    const HierarchicIterator end = e.hend( e.level()+1 );
    for( HierarchicIterator it = e.hbegin( e.level()+1 ); it != end; ++it )
    {
      const Entity &child = *it;
      LocalVector local = child.geometry().local( global );
      if( GenericReferenceElements< ctype, dimension >::cube().checkInside( local ) )
      {
        if( indexSet_.contains( *it ) )
          return EntityPointer( it );
        else
          return hFindEntity( *it, global );
      }
    }
    DUNE_THROW( Exception, "Unexpected internal Error" );
  }



  // Implementation of HierarchicSearch for SPGrid
  // ---------------------------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  inline HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
    ::HierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : Base( grid, indexSet )
  {}

}

#endif // #infdef DUNE_SPGRID_SPHIERARCHICSEARCH_HH
