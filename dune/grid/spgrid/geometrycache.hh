#ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
#define DUNE_SPGRID_GEOMETRYCACHE_HH

#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/normal.hh>

namespace Dune
{

  // SPGeometryPattern
  // -----------------

  template< int dim, int codim >
  struct SPGeometryPattern
  {
    typedef SPDirection< dim > Direction;

    SPGeometryPattern ();
    explicit SPGeometryPattern ( Direction dir );

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;

  private:
    int nonzero_[ dim - codim ];
    int zero_[ codim ];
  };

  template< int dim >
  struct SPGeometryPattern< dim, 0 >
  {
    typedef SPDirection< dim > Direction;

    SPGeometryPattern () = default;
    explicit SPGeometryPattern ( Direction /*dir*/ ) {}

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };

  template< int dim >
  struct SPGeometryPattern< dim, dim >
  {
    typedef SPDirection< dim > Direction;

    SPGeometryPattern () = default;
    explicit SPGeometryPattern ( Direction /*dir*/ ) {}

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };

  template<>
  struct SPGeometryPattern< 0, 0 >
  {
    typedef SPDirection< 0 > Direction;

    SPGeometryPattern () = default;
    explicit SPGeometryPattern ( Direction /*dir*/ ) {}

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };



  // SPJacobianTransposed
  // --------------------

  template< class ct, int dim, int mydim >
  class SPJacobianTransposed
    : private SPGeometryPattern< dim, dim - mydim >
  {
    typedef SPGeometryPattern< dim, dim - mydim > Pattern;

  public:
    typedef ct field_type;
    typedef ct value_type;

    typedef std::size_t size_type;

    static const int rows = mydim;
    static const int cols = dim;

    typedef FieldVector< field_type, dim > GlobalVector;
    typedef FieldVector< field_type, mydim > LocalVector;

    typedef Dune::FieldMatrix< field_type, rows, cols > FieldMatrix;
    typedef typename Dune::FieldTraits< field_type >::real_type real_type;

    SPJacobianTransposed () = default;

    SPJacobianTransposed ( const GlobalVector &h, SPDirection< dim > dir )
      : Pattern( std::move( dir ) )
    {
      for( int k = 0; k < rows; ++k )
        h_[ k ] = h[ pattern().nonzero( k ) ];
    }

    operator FieldMatrix () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;
    template< class X, class Y > void mhv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;
    template< class X, class Y > void umhv ( const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmhv ( const X &x, Y &y ) const;

    template< class X, class Y > void usmv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y > void usmtv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y > void usmhv ( const field_type &alpha, const X &x, Y &y ) const;

    field_type det () const;

    field_type determinant () const
    {
      if( rows == cols )
        return det();
      else
        DUNE_THROW( FMatrixError, "There is no determinant for a " << rows << "x" << cols << " matrix" );
    }

    real_type frobenius_norm () const { return h().two_norm(); }
    real_type frobenius_norm2 () const { return h().two_norm2(); }
    real_type infinity_norm () const { return h().infinity_norm(); }
    real_type infinity_norm_real () const { return h().infinity_norm_real(); }

  private:
    const Pattern &pattern () const { return static_cast< const Pattern & >( *this ); }
    const LocalVector &h () const { return h_; }

    LocalVector h_;
  };



  // FieldTraits for SPJacobianTransposed
  // ------------------------------------

  template< class ct, int dim, int mydim >
  struct FieldTraits< SPJacobianTransposed< ct, dim, mydim > >
  {
    typedef typename FieldTraits< ct >::field_type field_type;
    typedef typename FieldTraits< ct >::real_type real_type;
  };



  // SPJacobianInverseTransposed
  // ---------------------------

  template< class ct, int dim, int mydim >
  class SPJacobianInverseTransposed
    : private SPGeometryPattern< dim, dim - mydim >
  {
    typedef SPGeometryPattern< dim, dim - mydim > Pattern;

  public:
    typedef ct field_type;
    typedef ct value_type;

    typedef std::size_t size_type;

    static const int rows = dim;
    static const int cols = mydim;

    typedef FieldVector< field_type, dim > GlobalVector;
    typedef FieldVector< field_type, mydim > LocalVector;

    typedef Dune::FieldMatrix< field_type, rows, cols > FieldMatrix;
    typedef typename Dune::FieldTraits< field_type >::real_type real_type;

    SPJacobianInverseTransposed () = default;

    SPJacobianInverseTransposed ( const GlobalVector &h, SPDirection< dim > dir )
      : Pattern( dir )
    {
      for( int k = 0; k < cols; ++k )
        hInv_[ k ] = field_type( 1 ) / h[ pattern().nonzero( k ) ];
    }

    operator FieldMatrix () const;

    template< class X, class Y > void mv ( const X &x, Y &y ) const;
    template< class X, class Y > void mtv ( const X &x, Y &y ) const;
    template< class X, class Y > void mhv ( const X &x, Y &y ) const;

    template< class X, class Y > void umv ( const X &x, Y &y ) const;
    template< class X, class Y > void umtv ( const X &x, Y &y ) const;
    template< class X, class Y > void umhv ( const X &x, Y &y ) const;

    template< class X, class Y > void usmv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y > void usmtv ( const field_type &alpha, const X &x, Y &y ) const;
    template< class X, class Y > void usmhv ( const field_type &alpha, const X &x, Y &y ) const;

    template< class X, class Y > void mmv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmtv ( const X &x, Y &y ) const;
    template< class X, class Y > void mmhv ( const X &x, Y &y ) const;

    field_type det () const;

    field_type determinant () const
    {
      if( rows == cols )
        return det();
      else
        DUNE_THROW( FMatrixError, "There is no determinant for a " << rows << "x" << cols << " matrix" );
    }

    real_type frobenius_norm () const { return hInv().two_norm(); }
    real_type frobenius_norm2 () const { return hInv().two_norm2(); }
    real_type infinity_norm () const { return hInv().infinity_norm(); }
    real_type infinity_norm_real () const { return hInv().infinity_norm_real(); }

  private:
    const Pattern &pattern () const { return static_cast< const Pattern & >( *this ); }
    const LocalVector &hInv () const { return hInv_; }

    LocalVector hInv_;
  };



  // FieldTraits for SPJacobianInverseTransposed
  // -------------------------------------------

  template< class ct, int dim, int mydim >
  struct FieldTraits< SPJacobianInverseTransposed< ct, dim, mydim > >
  {
    typedef typename FieldTraits< ct >::field_type field_type;
    typedef typename FieldTraits< ct >::real_type real_type;
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

    typedef SPDirection< dimension > Direction;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef FieldVector< ctype, mydimension > LocalVector;

    typedef SPJacobianTransposed< ctype, dimension, mydimension > JacobianTransposed;
    typedef SPJacobianInverseTransposed< ctype, dimension, mydimension > JacobianInverseTransposed;

    SPGeometryCache ( const GlobalVector &h, Direction dir )
      : jacobianTransposed_( h, dir ), jacobianInverseTransposed_( h, dir ), volume_( jacobianTransposed_.det() )
    {}

    const ctype &volume () const { return volume_; }
    const JacobianTransposed &jacobianTransposed () const { return jacobianTransposed_; }
    const JacobianInverseTransposed &jacobianInverseTransposed () const { return jacobianInverseTransposed_; }

    JacobianTransposed jacobianTransposed_;
    JacobianInverseTransposed jacobianInverseTransposed_;
    ctype volume_;
  };



  // Implementation of SPGeometryPattern
  // -----------------------------------

  template< int dim, int codim >
  inline SPGeometryPattern< dim, codim >::SPGeometryPattern ()
  {
    const int mydim = dim - codim;
    for( int k = 0; k < mydim; ++k )
      nonzero_[ k ] = k;
    for( int k = 0; k < codim; ++k )
      zero_[ k ] = mydim + k;
  }

  template< int dim, int codim >
  inline SPGeometryPattern< dim, codim >::SPGeometryPattern ( Direction dir )
  {
    int k = 0;
    for( int j = 0; j < dim; ++j )
    {
      if( dir[ j ] != 0 )
        nonzero_[ k++ ] = j;
      else
        zero_[ j - k ] = j;
    }
    assert( k == dim - codim );
  }


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



  // Implementation of SPJacobianTransposed
  // --------------------------------------

  template< class ct, int dim, int mydim >
  inline SPJacobianTransposed< ct, dim, mydim >::operator FieldMatrix () const
  {
    FieldMatrix matrix( field_type( 0 ) );
    for( int k = 0; k < rows; ++k )
      matrix[ k ][ pattern().nonzero( k ) ] = h()[ k ];
    return matrix;
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ k ] = h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] = h()[ k ] * x[ k ];
    for( int k = 0; k < cols - rows; ++k )
      y[ pattern().zero( k ) ] = field_type( 0 );
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] = conjugateComplex( h()[ k ] ) * x[ k ];
    for( int k = 0; k < cols - rows; ++k )
      y[ pattern().zero( k ) ] = field_type( 0 );
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ k ] += h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] += h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::umhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] += conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::usmv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ k ] += alpha * h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::usmtv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] += alpha * h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::usmhv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] += alpha * conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ k ] -= h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] -= h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianTransposed< ct, dim, mydim >::mmhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < rows; ++k )
      y[ pattern().nonzero( k ) ] -= conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int mydim >
  inline typename SPJacobianTransposed< ct, dim, mydim >::field_type SPJacobianTransposed< ct, dim, mydim >::det () const
  {
    field_type det( 1 );
    for( int k = 0; k < rows; ++k )
      det *= h()[ k ];
    return det;
  }



  // Implementation of SPJacobianInverseTransposed
  // ---------------------------------------------

  template< class ct, int dim, int mydim >
  inline SPJacobianInverseTransposed< ct, dim, mydim >::operator FieldMatrix () const
  {
    FieldMatrix matrix( field_type( 0 ) );
    for( int k = 0; k < cols; ++k )
      matrix[ pattern().nonzero( k ) ][ k ] = hInv()[ k ];
    return matrix;
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ pattern().nonzero( k ) ] = hInv()[ k ] * x[ k ];
    for( int k = 0; k < rows - cols; ++k )
      y[ pattern().zero( k ) ] = field_type( 0 );
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] = hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] = conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ pattern().nonzero( k ) ] += hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] += hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::umhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] += conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::usmv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ pattern().nonzero( k ) ] += alpha * hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::usmtv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] += alpha * hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::usmhv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] += alpha * conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ pattern().nonzero( k ) ] -= hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] -= hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  template< class X, class Y >
  inline void SPJacobianInverseTransposed< ct, dim, mydim >::mmhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < cols; ++k )
      y[ k ] -= conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int mydim >
  inline typename SPJacobianInverseTransposed< ct, dim, mydim >::field_type SPJacobianInverseTransposed< ct, dim, mydim >::det () const
  {
    field_type det( 1 );
    for( int k = 0; k < cols; ++k )
      det *= hInv()[ k ];
    return det;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
