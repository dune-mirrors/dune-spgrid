#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <vector>

namespace Dune
{

  template< class ct, int dim >
  class SPGridLevel
  {
    typedef SPGridLevel< ct, dim > This;

  public:
    static const int dimension = dim;

    typedef unsigned int MultiIndex[ dimension ];

    unsigned int level () const
    {
      return level_;
    }

    unsigned int
    n ( const unsigned int codim, const unsigned int direction, int i ) const
    {
      return n_[ i ] + multiDirection_[ codim ][ direction ][ i ];
    }

  private:
    unsigned int level_;
    MultiIndex n_;
    std::vector< MultiIndex > multiDirection_[ dimension+1 ];
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH
