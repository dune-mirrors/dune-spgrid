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

  template<>
  struct SPGeometryPattern< 0, 0 >
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

    struct JacobianTransposed;
    struct JacobianInverseTransposed;

    SPGeometryCache ( const GlobalVector &h, const unsigned int dir );

    const ctype &volume () const;
    const JacobianTransposed &jacobianTransposed () const;
    const JacobianInverseTransposed &jacobianInverseTransposed () const;

  private:
    JacobianTransposed jacobianTransposed_;
    JacobianInverseTransposed jacobianInverseTransposed_;
    ctype volume_;
  };



  // SPGeometryCache::JacobianTransposed
  // -----------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianTransposed
  {
    static const int rows = mydimension;
    static const int cols = dimension;

    typedef Dune::FieldMatrix< ctype, rows, cols > FieldMatrix;

    JacobianTransposed ( const GlobalVector &h, const unsigned int dir );

    operator const FieldMatrix & () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;

    ctype det () const;

  private:
    friend class SPGeometryCache< ct, dim, codim >;

    JacobianTransposed ( const JacobianTransposed &other );

    // prohibit assignment
    const JacobianTransposed &operator= ( const JacobianTransposed & );

    SPGeometryPattern< dimension, codimension > pattern_;
    LocalVector h_;
    FieldMatrix matrix_;
  };



  // SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
  {
    static const int rows = dimension;
    static const int cols = mydimension;

    typedef Dune::FieldMatrix< ctype, rows, cols > FieldMatrix;

    JacobianInverseTransposed ( const GlobalVector &h, const unsigned int dir );

    operator const FieldMatrix & () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;

    ctype det () const;

  private:
    friend class SPGeometryCache< ct, dim, codim >;

    JacobianInverseTransposed ( const JacobianInverseTransposed &other );

    // prohibit assignment
    const JacobianInverseTransposed &operator= ( const JacobianInverseTransposed & );

    SPGeometryPattern< dimension, codimension > pattern_;
    LocalVector hInv_;
    FieldMatrix matrix_;
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


  inline SPGeometryPattern< 0, 0 >::SPGeometryPattern ( const unsigned int dir )
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


  inline int SPGeometryPattern< 0, 0 >::nonzero ( const int k ) const
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


  inline int SPGeometryPattern< 0, 0 >::zero ( const int k ) const
  {
    assert( false );
    return k;
  }



  // Implementation of SPGeometryCache
  // ---------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >
    ::SPGeometryCache ( const GlobalVector &h, const unsigned int dir )
  : jacobianTransposed_( h, dir ),
    jacobianInverseTransposed_( h, dir ),
    volume_( jacobianTransposed_.det() )
  {}


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
  inline const typename SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed &
  SPGeometryCache< ct, dim, codim >::jacobianInverseTransposed () const
  {
    return jacobianInverseTransposed_;
  }



  // Implementation of SPGeometryCache::JacobianTransposed
  // -----------------------------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed
    ::JacobianTransposed ( const GlobalVector &h, const unsigned int dir )
  : pattern_( dir ),
    matrix_( ctype( 0 ) )
  {
    for( int k = 0; k < mydimension; ++k )
    {
      const int j = pattern_.nonzero( k );
      h_[ k ] = h[ j ];
      matrix_[ k ][ j ] = h_[ k ];
    }
  }


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed
    ::JacobianTransposed ( const JacobianTransposed &other )
  : pattern_( other.pattern_ ),
    h_( other.h_ ),
    matrix_( other.matrix_ )
  {}


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed
    ::operator const FieldMatrix & () const
  {
    return matrix_;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = h_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] = h_[ k ] * x[ k ];
    for( int k = 0; k < codimension; ++k )
      y[ pattern_.zero( k ) ] = ctype( 0 );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += h_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] += h_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= h_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] -= h_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::ctype
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::det () const
  {
    ctype det( 1 );
    for( int k = 0; k < mydimension; ++k )
      det *= h_[ k ];
    return det;
  }



  // Implementation of SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
    ::JacobianInverseTransposed ( const GlobalVector &h, const unsigned int dir )
  : pattern_( dir ),
    matrix_( ctype( 0 ) )
  {
    for( int k = 0; k < mydimension; ++k )
    {
      const int j = pattern_.nonzero( k );
      hInv_[ k ] = ctype( 1 ) / h[ j ];
      matrix_[ j ][ k ] = hInv_[ k ];
    }
  }


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
    ::JacobianInverseTransposed ( const JacobianInverseTransposed &other )
  : pattern_( other.pattern_ ),
    hInv_( other.hInv_ ),
    matrix_( other.matrix_ )
  {}


  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
    ::operator const FieldMatrix & () const
  {
    return matrix_;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] = hInv_[ k ] * x[ k ];
    for( int k = 0; k < codimension; ++k )
      y[ pattern_.zero( k ) ] = ctype( 0 );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = hInv_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] += hInv_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += hInv_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern_.nonzero( k ) ] -= hInv_[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= hInv_[ k ] * x[ pattern_.nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::ctype
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::det () const
  {
    ctype det( 1 );
    for( int k = 0; k < mydimension; ++k )
      det *= hInv_[ k ];
    return det;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
