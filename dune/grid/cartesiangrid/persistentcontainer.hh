#ifndef DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH
#define DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH

#include <dune/grid/freiburg/persistentcontainerwrapper.hh>

namespace Dune
{
  // PersistentContainer for CartesianGrid
  // -------------------------------------

  template< class HostGrid, class DataImp, class Allocator >
  class PersistentContainer< CartesianGrid< HostGrid >, DataImp, Allocator >
    : public PersistentContainerWrapper< CartesianGrid< HostGrid >, DataImp, Allocator >
  {
    typedef CartesianGrid< HostGrid > Grid ;
    typedef PersistentContainerWrapper< Grid,  DataImp, Allocator > Base ;
  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
      : Base( grid, codim, allocator )
    {
    }
  };
} // end namespace Dune

#endif // end DUNE_CARTESIANGRID_PERSISTENTCONTAINER_HH
