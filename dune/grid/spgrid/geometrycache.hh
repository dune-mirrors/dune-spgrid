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

  protected:
    using Base::direction;

  private:
    LocalVector h_;
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
    geometryCache_.jacobianTransposed_.mv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianTransposed_.mtv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianTransposed_.umv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianTransposed_.umtv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianTransposed_.mmv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianTransposed_.mmtv( x, y );
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
    geometryCache_.jacobianInverseTransposed_.mv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianInverseTransposed_.mtv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianInverseTransposed_.umv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianInverseTransposed_.umtv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianInverseTransposed_.mmv( x, y );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmtv ( const X &x, Y &y ) const
  {
    geometryCache_.jacobianInverseTransposed_.mmtv( x, y );
  }

}

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
