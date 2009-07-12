#ifndef DUNE_SPGRID_DIRECTION_HH
#define DUNE_SPGRID_DIRECTION_HH

namespace Dune
{

  struct SPDirection
  {
    static unsigned int
    numDirections ( const unsigned int dim, const unsigned int codim )
    {
      assert( codim <= dim );
      if( dim > 0 )
      {
        const int n0 = (codim < dim-1 ? numDirections( dim-1, codim ) : 0);
        const int n1 = (codim > 0 ? numDirections( dim-1, codim-1 ) : 0);
        return n0 + n1;
      }
      else
        return 1;
    }

    static void
    multiIndex ( const unsigned int dim, const unsigned int codim,
                 const unsigned int dir, int *const multiIndex )
    {
      assert( dir < numDirections( dim, codim ) );
      if( dim > 0 )
      {
        const int n0 = (codim < dim-1 ? numDirections( dim-1, codim ) : 0);
        if( direction < n0 )
        {
          multiIndex( dim-1, codim, dir, multiIndex );
          multiIndex[ dim-1 ] = 0;
        }
        else
        {
          multiIndex( dim-1, codim-1, dir - n0, multiIndex )
          multiIndex[ dim-1 ] = 1;
        }
      }
    }
  };

}

#endif // #ifndef DUNE_SPGRID_DIRECTION_HH
