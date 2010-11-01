#ifndef DUNE_SPGRID_HIERARCHICSEARCH_HH
#define DUNE_SPGRID_HIERARCHICSEARCH_HH

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/indexidset.hh>

#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid;

  template< class Grid >
  class SPIndexSet;

  template< class Grid, class IndexSet >
  class HierarchicSearch;



  // SPHierarchicSearch
  // ------------------

  template< class Grid, class IndexSet >
  class SPHierarchicSearch
  {
  public:
    typedef typename Grid::ctype ctype;
    typedef typename Grid::template Codim< 0 >::Entity Entity;
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    static const int dimension = Grid::dimension;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef GlobalVector LocalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet );

    EntityPointer findEntity ( const GlobalVector &global ) const;
    
  private:
    EntityPointer hFindEntity ( const Entity &e, const GlobalVector &global ) const;

    const Grid &grid_;
    const IndexSet &indexSet_;
  };



  // SPHierarchicSearch for SPIndexSet
  // ---------------------------------

  template< class Grid >
  class SPHierarchicSearch< Grid, SPIndexSet< Grid > >
  {
    typedef SPIndexSet< Grid > IndexSet;

  public:
    typedef typename Grid::ctype ctype;
    typedef typename Grid::template Codim< 0 >::Entity Entity;
    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;

    static const int dimension = Grid::dimension;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef GlobalVector LocalVector;

    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet );

    EntityPointer findEntity ( const GlobalVector &global ) const;
    
  private:
    const Grid &grid_;
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
    SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet );
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
    typedef typename Grid::HierarchicIterator HierarchicIterator;

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



  // Implementation of SPHierarchicSearch for SPIndexSet
  // ---------------------------------------------------

  template< class Grid >
  inline SPHierarchicSearch< Grid, SPIndexSet< Grid > >
    ::SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : grid_( grid ),
    indexSet_( indexSet )
  {}


  template< class Grid >
  inline typename SPHierarchicSearch< Grid, SPIndexSet< Grid > >::EntityPointer
  SPHierarchicSearch< Grid, SPIndexSet< Grid > >::findEntity ( const GlobalVector &global ) const
  {
    return grid_.findEntity( global, indexSet_.gridLevel().level() );
  }



  // Implementation of SPHierarchicSearch for IndexSet< SPIndexSet >
  // ---------------------------------------------------------------

  template< class Grid >
  inline SPHierarchicSearch< Grid, IndexSet< Grid, SPIndexSet< Grid >, typename SPIndexSet< Grid >::IndexType > >
    ::SPHierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : Base( grid, static_cast< const SPIndexSet< Grid > & >( indexSet ) )
  {}



  // Implementation of HierarchicSearch for SPGrid
  // ---------------------------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm, class IndexSet >
  inline HierarchicSearch< SPGrid< ct, dim, strategy, Comm >, IndexSet >
    ::HierarchicSearch ( const Grid &grid, const IndexSet &indexSet )
  : Base( grid, indexSet )
  {}

}

#endif // #infdef DUNE_SPGRID_HIERARCHICSEARCH_HH
