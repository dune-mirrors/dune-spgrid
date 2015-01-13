#ifndef DUNE_SPGRID_PERSISTENTCONTAINER_HH
#define DUNE_SPGRID_PERSISTENTCONTAINER_HH

#include <vector>

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainervector.hh>

#include <dune/grid/spgrid/declaration.hh>

namespace Dune
{

  // PersistentContainer for SPGrid
  // ------------------------------

  template< class ct, int dim, template< int > class Ref, class Comm, class T >
  class PersistentContainer< SPGrid< ct, dim, Ref, Comm >, T >
    : public PersistentContainerVector< SPGrid< ct, dim, Ref, Comm >, typename SPGrid< ct, dim, Ref, Comm >::HierarchicIndexSet, std::vector< T > >
  {
    typedef PersistentContainerVector< SPGrid< ct, dim, Ref, Comm >, typename SPGrid< ct, dim, Ref, Comm >::HierarchicIndexSet, std::vector< T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid.hierarchicIndexSet(), codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_PERSISTENTCONTAINER_HH
