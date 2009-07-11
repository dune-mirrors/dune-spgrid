#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

namespace Dune
{

  template< class ct, int dim, int codim >
  class SPGridLevel
  {

    unsigned int level () const
    {
      return level_;
    }

  private:
    unsigned int level_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH
