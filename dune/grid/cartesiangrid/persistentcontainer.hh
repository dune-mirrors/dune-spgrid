#ifndef DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH
#define DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainerwrapper.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class >
  class CartesianGrid;



  // PersistentContainer for CartesianGrid
  // -------------------------------------

  template< class HostGrid, class T, class Allocator >
  class PersistentContainer< CartesianGrid< HostGrid >, T, Allocator >
  : public PersistentContainerWrapper< CartesianGrid< HostGrid >, T, Allocator >
  {
    typedef PersistentContainerWrapper< CartesianGrid< HostGrid >, T, Allocator > Base;

  public:
    typedef typename Base::Grid Grid;

    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
    : Base( grid, codim, allocator )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH
