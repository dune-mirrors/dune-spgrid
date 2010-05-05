#ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
#define DUNE_SPGRID_GEOMETRYCACHE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  // SPGeometryPattern
  // -----------------

  template< int dim, int codim >
  struct SPGeometryPattern
  {
    explicit SPGeometryPattern ( const unsigned int dir );

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;

  private:
    int nonzero_[ dim - codim ];
    int zero_[ codim ];
  };

  template< int dim >
  struct SPGeometryPattern< dim, 0 >
  {
    explicit SPGeometryPattern ( const unsigned int dir );

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };

  template< int dim >
  struct SPGeometryPattern< dim, dim >
  {
    explicit SPGeometryPattern ( const unsigned int dir );

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };



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
    typedef FieldMatrix< ctype, dimension, mydimension > JacobianMatrix;
    typedef FieldMatrix< ctype, mydimension, dimension > JacobianMatrixTransposed;

    struct JacobianTransposed;
    struct JacobianInverseTransposed;

    SPGeometryCache ( const GlobalVector &h, const unsigned int dir );

    const ctype &volume () const;
    JacobianTransposed jacobianTransposed () const;
    JacobianInverseTransposed jacobianInverseTransposed () const;

    GlobalVector global ( const GlobalVector &origin, const LocalVector &local ) const;
    LocalVector local ( const GlobalVector &origin, const GlobalVector &global ) const;

  private:
    SPGeometryPattern< dimension, codimension > pattern_;
    LocalVector h_;
    LocalVector hInv_;
    ctype volume_;
    JacobianMatrixTransposed jacobianTransposed_;
    JacobianMatrix jacobianInverseTransposed_;
  };



  // SPGeometryCache::JacobianTransposed
  // -----------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianTransposed
  {
    typedef SPGeometryCache< ct, dim, codim > GeometryCache;

    JacobianTransposed ( const GeometryCache &geometryCache );

    operator const JacobianMatrixTransposed & () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;

  private:
    const GeometryCache &geometryCache_;
  };



  // SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
  {
    typedef SPGeometryCache< ct, dim, codim > GeometryCache;

    JacobianInverseTransposed ( const GeometryCache &geometryCache );

    operator const JacobianMatrix & () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;

  private:
    const GeometryCache &geometryCache_;
  };



  // Implementation of SPGeometryPattern
  // -----------------------------------

  template< int dim, int codim >
  inline SPGeometryPattern< dim, codim >::SPGeometryPattern ( const unsigned int dir )
  {
    int k = 0;
    for( int j = 0; j < dim; ++j )
    {
      if( ((dir >> j) & 1) != 0 )
        nonzero_[ k++ ] = j;
      else
        zero_[ j - k ] = j;
    }
    assert( k == dim - codim );
  }


  template< int dim >
  inline SPGeometryPattern< dim, 0 >::SPGeometryPattern ( const unsigned int dir )
  {}


  template< int dim >
  inline SPGeometryPattern< dim, dim >::SPGeometryPattern ( const unsigned int dir )
  {}


  template< int dim, int codim >
  inline int SPGeometryPattern< dim, codim >::nonzero ( const int k ) const
  {
    assert( (k >= 0) && (k < dim - codim) );
    return nonzero_[ k ];
  }


  template< int dim >
  inline int SPGeometryPattern< dim, 0 >::nonzero ( const int k ) const
  {
    assert( (k >= 0) && (k < dim) );
    return k;
  }


  template< int dim >
  inline int SPGeometryPattern< dim, dim >::nonzero ( const int k ) const
  {
    assert( false );
    return k;
  }


  template< int dim, int codim >
  inline int SPGeometryPattern< dim, codim >::zero ( const int k ) const
  {
    assert( (k >= 0) && (k < codim) );
    return zero_[ k ];
  }


  template< int dim >
  inline int SPGeometryPattern< dim, 0 >::zero ( const int k ) const
  {
    assert( false );
    return k;
  }


  template< int dim >
  inline int SPGeometryPattern< dim, dim >::zero ( const int k ) const
  {
    assert( (k >= 0) && (k < dim) );
    return k;
  }



  // Implementation of SPGeometryCache
  // ---------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >
    ::SPGeometryCache ( const GlobalVector &h, const unsigned int dir )
  : pattern_( dir ),
    volume_( 1 ),
    jacobianTransposed_( ctype( 0 ) ),
    jacobianInverseTransposed_( ctype( 0 ) )
  {
    for( int k = 0; k < mydimension; ++k )
    {
      h_[ k ] = h[ pattern_.nonzero( k ) ];
      hInv_[ k ] = ctype( 1 ) / h_[ k ];
      volume_ *= h_[ k ];
      jacobianTransposed_[ k ][ pattern_.nonzero( k ) ] = h_[ k ];
      jacobianInverseTransposed_[ pattern_.nonzero( k ) ][ k ] = hInv_[ k ];
    }
  }


  template< class ct, int dim, int codim >
  inline const typename SPGeometryCache< ct, dim, codim >::ctype &
  SPGeometryCache< ct, dim, codim >::volume () const
  {
    return volume_;
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::JacobianTransposed
  SPGeometryCache< ct, dim, codim >::jacobianTransposed () const
  {
    return JacobianTransposed( *this );
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
  SPGeometryCache< ct, dim, codim >::jacobianInverseTransposed () const
  {
    return JacobianInverseTransposed( *this );
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::GlobalVector
  SPGeometryCache< ct, dim, codim >
    ::global ( const GlobalVector &origin, const LocalVector &local ) const
  {
    GlobalVector y( origin );
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] += h_[ k ] * local[ k ];
    return y;
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::LocalVector
  SPGeometryCache< ct, dim, codim >
    ::local ( const GlobalVector &origin, const GlobalVector &global ) const
  {
    LocalVector x;
    for( int k = 0; k < mydimension; ++k )
    {
      const int j = pattern_.nonzero( k );
      x[ k ] = (global[ j ] - origin[ j ]) / h_[ j ];
    }
    return x;
  }



  // Implementation of SPGeometryCache::JacobianTransposed
  // -----------------------------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed
    ::JacobianTransposed ( const GeometryCache &geometryCache )
  : geometryCache_( geometryCache )
  {}


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed
    ::operator const JacobianMatrixTransposed & () const
  {
    return geometryCache_.jacobianTransposed_;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = geometryCache_.h_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mtv ( const X &x, Y &y ) const
  {
    if( mydimension < dimension )
    {
      for( int k = 0; k < dimension; ++k )
        y[ k ] = 0;
    }
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] = geometryCache_.h_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += geometryCache_.h_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] += geometryCache_.h_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= geometryCache_.h_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] -= geometryCache_.h_[ k ] * x[ k ];
  }



  // Implementation of SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
    ::JacobianInverseTransposed ( const GeometryCache &geometryCache )
  : geometryCache_( geometryCache )
  {}


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
    ::operator const JacobianMatrix & () const
  {
    return geometryCache_.jacobianInverseTransposed_;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mv ( const X &x, Y &y ) const
  {
    if( mydimension < dimension )
    {
      for( int k = 0; k < dimension; ++k )
        y[ k ] = 0;
    }
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] = geometryCache_.hInv_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = geometryCache_.hInv_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] += geometryCache_.hInv_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += geometryCache_.hInv_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ geometryCache_.pattern_.nonzero( k ) ] -= geometryCache_.hInv_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= geometryCache_.hInv_[ k ] * x[ geometryCache_.pattern_.nonzero( k ) ];
  }

}

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
