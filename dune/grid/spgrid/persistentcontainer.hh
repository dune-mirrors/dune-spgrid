#ifndef DUNE_SPGRID_PERSISTENTCONTAINER_HH
#define DUNE_SPGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>

#include <dune/grid/spgrid/declaration.hh>

namespace Dune
{

  // PersistentContainer for SPGrid
  // ------------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm,
            class Data, class Allocator >
  class PersistentContainer< SPGrid< ct, dim, strategy, Comm >, Data, Allocator >
  : public PersistentContainerVector< SPGrid< ct, dim, strategy, Comm >, 
                                      typename SPGrid< ct, dim, strategy, Comm >::HierarchicIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef SPGrid< ct, dim, strategy, Comm > Grid;
    typedef PersistentContainerVector< Grid, typename Grid::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
    : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1, allocator )
    {
    }
  };

} // end namespace Dune

#endif // end DUNE_SPGRID_PERSISTENTCONTAINER_HH
