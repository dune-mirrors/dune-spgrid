#ifndef DUNE_GRID_TEST_CHECKTREE_HH
#define DUNE_GRID_TEST_CHECKTREE_HH

#include <functional>
#include <iterator>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/rangegenerators.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, class, class >
  class EntityTree;

  template< class, class >
  class IntersectionTree;



  // makeEntityTree
  // --------------

  template< class Grid, class Entity, class Predicate >
  inline static EntityTree< Entity::codimension, Grid, Predicate >
  makeEntityTree ( const Grid &grid, const Entity &entity, const Predicate &predicate )
  {
    return EntityTree< Entity::codimension, Grid, Predicate >( grid, entity, predicate );
  }



  // checkEntityTree
  // ---------------

  template< class Grid, class Entity >
  void checkEntityTree ( const Grid &grid, const Entity &entity )
  {
    std::function< bool( const Entity & ) > allLeaf = [] ( const Entity &e ) { return true; };
    auto tree = makeEntityTree( grid, entity, allLeaf );

    auto noOp = [] ( const Entity &e ) {};
    if( testForwardIterator( tree.begin(), tree.end(), noOp ) != 0 )
      DUNE_THROW( Dune::Exception, "Tree iterator does not fulfill the forward iterator concept." );

    if( std::distance( tree.begin(), tree.end() ) != 1 )
      DUNE_THROW( Dune::Exception, "Tree iterator with all leaf predicate should iterate over exactly 1 element." );
  }

  template< int codim, class GridView >
  void checkEntityTree ( const GridView &gridView, Codim< codim > = Codim< codim >() )
  {
    for( const auto &entity : entities( gridView, Codim< codim >() ) )
      checkEntityTree( gridView.grid(), entity );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_TEST_CHECKTREE_HH
