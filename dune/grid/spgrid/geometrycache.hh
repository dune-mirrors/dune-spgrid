#ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
#define DUNE_SPGRID_GEOMETRYCACHE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  // SPGeometryCache
  // ---------------

  template< class ct, int dim, int codim >
  class SPGeometryCache
  {
    typedef SPGeometryCache< ct, dim, codim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldMatrix< ctype, dimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, dimension > JacobianTransposed;

    SPGeometryCache ( const GlobalVector &h, const unsigned int dir )
    : volume_( 1 ),
      jacobianTransposed_( 0 ),
      jacobianInverseTransposed_( 0 )
    {
      int k = 0;
      for( int j = 0; j < dimension; ++j )
      {
        if( ((dir >> j) & 1) == 0 )
          continue;
        volume_ *= h[ j ];
        jacobianTransposed_[ j ][ k ] = h[ j ];
        jacobianInverseTransposed_[ k ][ j ] = ctype( 1 ) / h[ j ];
        ++k;
      }
    }

    const ctype &volume () const
    {
      return volume_;
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      return jacobianTransposed_;
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      return jacobianInverseTransposed_;
    }

  private:
    ctype volume_;
    JacobianTransposed jacobianTransposed_;
    Jacobian jacobianInverseTransposed_;
  };

}

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
