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

    SPGeometryCache ( const GlobalVector &h, const unsigned int dir );

    const ctype &volume () const;
    const JacobianTransposed &jacobianTransposed () const;
    const Jacobian &jacobianInverseTransposed () const;

    GlobalVector global ( const GlobalVector &origin, const LocalVector &local ) const;
    LocalVector local ( const GlobalVector &origin, const GlobalVector &global ) const;

  private:
    ctype volume_;
    JacobianTransposed jacobianTransposed_;
    Jacobian jacobianInverseTransposed_;
  };



  // Implementation of SPGeometryCache
  // ---------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >
    ::SPGeometryCache ( const GlobalVector &h, const unsigned int dir )
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
      jacobianTransposed_[ k ][ j ] = h[ j ];
      jacobianInverseTransposed_[ j ][ k ] = ctype( 1 ) / h[ j ];
      ++k;
    }
  }


  template< class ct, int dim, int codim >
  inline const typename SPGeometryCache< ct, dim, codim >::ctype &
  SPGeometryCache< ct, dim, codim >::volume () const
  {
    return volume_;
  }


  template< class ct, int dim, int codim >
  inline const typename SPGeometryCache< ct, dim, codim >::JacobianTransposed &
  SPGeometryCache< ct, dim, codim >::jacobianTransposed () const
  {
    return jacobianTransposed_;
  }


  template< class ct, int dim, int codim >
  inline const typename SPGeometryCache< ct, dim, codim >::Jacobian &
  SPGeometryCache< ct, dim, codim >::jacobianInverseTransposed () const
  {
    return jacobianInverseTransposed_;
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::GlobalVector
  SPGeometryCache< ct, dim, codim >
    ::global ( const GlobalVector &origin, const LocalVector &local ) const
  {
    // this can be optimized
    GlobalVector y( origin );
    jacobianTransposed().umtv( local, y );
    return y;
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::LocalVector
  SPGeometryCache< ct, dim, codim >
    ::local ( const GlobalVector &origin, const GlobalVector &global ) const
  {
    // this can be optimized
    GlobalVector y = global - origin;
    LocalVector x;
    jacobianInverseTransposed().mtv( y, x );
    return x;
  }

}

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
