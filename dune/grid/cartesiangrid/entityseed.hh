#ifndef DUNE_CARTESIANGRID_ENTITYSEED_HH
#define DUNE_CARTESIANGRID_ENTITYSEED_HH

#include <dune/common/typetraits.hh>

namespace Dune
{

  // CartesianGridEntitySeed
  // -----------------------

  template< int codim, class Grd >
  class CartesianGridEntitySeed
  {
    typedef typename remove_const< Grd >::type::Traits Traits;

  public:
    static const int codimension = codim;
    static const int dimension = Traits::dimension;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = Traits::dimensionworld;

    typedef typename Traits::Grid Grid;
    typedef typename Traits::template Codim< codim >::Entity Entity;

    typedef typename Traits::HostGrid HostGrid;
    typedef typename HostGrid::template Codim< codim >::EntitySeed HostEntitySeed;

    explicit CartesianGridEntitySeed ( const HostEntitySeed &hostEntitySeed )
    : hostEntitySeed_( hostEntitySeed )
    {}

    const HostEntitySeed &hostEntitySeed () const { return hostEntitySeed_; }

  private:
    HostEntitySeed hostEntitySeed_;
  };

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_ENTITYSEED_HH
