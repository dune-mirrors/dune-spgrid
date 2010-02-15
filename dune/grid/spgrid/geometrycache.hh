#ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
#define DUNE_SPGRID_GEOMETRYCACHE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  // SPGeometryCacheBase
  // -------------------

  template< int dim, int codim >
  class SPGeometryCacheBase
  {
    typedef SPGeometryCacheBase< dim, codim > This;

  protected:
    explicit SPGeometryCacheBase ( const unsigned int dir );

    int direction ( const int k ) const;

  private:
    int direction_[ dim - codim ];
  };



  template< int dim >
  class SPGeometryCacheBase< dim, 0 >
  {
    typedef SPGeometryCacheBase< dim, 0 > This;

  protected:
    explicit SPGeometryCacheBase ( const unsigned int dir );

    int direction ( const int k ) const;
  };



  template< int dim >
  class SPGeometryCacheBase< dim, dim >
  {
    typedef SPGeometryCacheBase< dim, dim > This;

  protected:
    explicit SPGeometryCacheBase ( const unsigned int dir );

    int direction ( const int k ) const;
  };



  // SPGeometryCache
  // ---------------

  template< class ct, int dim, int codim >
  class SPGeometryCache
  : public SPGeometryCacheBase< dim, codim >
  {
    typedef SPGeometryCache< ct, dim, codim > This;
    typedef SPGeometryCacheBase< dim, codim > Base;

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

  protected:
    using Base::direction;

  private:
    LocalVector h_;
    ctype volume_;
    JacobianTransposed jacobianTransposed_;
    Jacobian jacobianInverseTransposed_;
  };



  // Implementation of SPGeometryCacheBase
  // -------------------------------------

  template< int dim, int codim >
  inline SPGeometryCacheBase< dim, codim >
    ::SPGeometryCacheBase ( const unsigned int dir )
  {
    int k = 0;
    for( int j = 0; j < dim; ++j )
    {
      if( ((dir >> j) & 1) != 0 )
        direction_[ k++ ] = j;
    }
    assert( k == dim - codim );
  }


  template< int dim >
  inline SPGeometryCacheBase< dim, 0 >
    ::SPGeometryCacheBase ( const unsigned int dir )
  {}


  template< int dim >
  inline SPGeometryCacheBase< dim, dim >
    ::SPGeometryCacheBase ( const unsigned int dir )
  {}


  template< int dim, int codim >
  inline int
  SPGeometryCacheBase< dim, codim >::direction ( const int k ) const
  {
    assert( (k >= 0) && (k < dim - codim) );
    return direction_[ k ];
  }


  template< int dim >
  inline int
  SPGeometryCacheBase< dim, 0 >::direction ( const int k ) const
  {
    assert( (k >= 0) && (k < dim) );
    return k;
  }


  template< int dim >
  inline int
  SPGeometryCacheBase< dim, dim >::direction ( const int k ) const
  {
    assert( false );
    return k;
  }



  // Implementation of SPGeometryCache
  // ---------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >
    ::SPGeometryCache ( const GlobalVector &h, const unsigned int dir )
  : Base( dir ),
    volume_( 1 ),
    jacobianTransposed_( 0 ),
    jacobianInverseTransposed_( 0 )
  {
    for( int k = 0; k < mydimension; ++k )
    {
      h_[ k ] = h[ direction( k ) ];
      volume_ *= h_[ k ];
      jacobianTransposed_[ k ][ direction( k ) ] = h_[ k ];
      jacobianInverseTransposed_[ direction( k ) ][ k ] = ctype( 1 ) / h_[ k ];
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
    GlobalVector y( origin );
    for( int k = 0; k < mydimension; ++k )
      y[ direction( k ) ] += h_[ k ] * local[ k ];
    return y;
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::LocalVector
  SPGeometryCache< ct, dim, codim >
    ::local ( const GlobalVector &origin, const GlobalVector &global ) const
  {
    LocalVector x;
    for( int k = 0; k < mydimension; ++k )
      x[ k ] = (global[ direction( k ) ] - origin[ direction( k ) ]) / h_[ k ];
    return x;
  }

}

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
