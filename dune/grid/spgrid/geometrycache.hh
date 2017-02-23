#ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
#define DUNE_SPGRID_GEOMETRYCACHE_HH

#include <dune/common/exception.hh>
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

    explicit SPGeometryPattern ( Direction dir ) {}

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };

  template< int dim >
  struct SPGeometryPattern< dim, dim >
  {
    typedef SPDirection< dim > Direction;

    explicit SPGeometryPattern ( Direction dir ) {}

    int nonzero ( const int k ) const;
    int zero ( const int k ) const;
  };

  template<>
  struct SPGeometryPattern< 0, 0 >
  {
    typedef SPDirection< 0 > Direction;

    explicit SPGeometryPattern ( Direction dir ) {}

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

    typedef SPDirection< dimension > Direction;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef FieldVector< ctype, mydimension > LocalVector;

  private:
    typedef SPGeometryPattern< dimension, codimension > Pattern;

    struct MatrixStorage
      : public Pattern
    {
      explicit MatrixStorage ( Direction dir ) : Pattern( dir ) {}

      LocalVector h;
    };

  public:
    struct JacobianTransposed;
    struct JacobianInverseTransposed;

    SPGeometryCache ( const GlobalVector &h, Direction dir )
      : jacobianTransposed_( h, dir ),
        jacobianInverseTransposed_( h, dir ),
        volume_( jacobianTransposed_.det() )
    {}

    const ctype &volume () const;
    const JacobianTransposed &jacobianTransposed () const;
    const JacobianInverseTransposed &jacobianInverseTransposed () const;

    JacobianTransposed jacobianTransposed_;
    JacobianInverseTransposed jacobianInverseTransposed_;
    ctype volume_;
  };



  // FieldTraits for SPGeometryCache::JacobianTransposed
  // ---------------------------------------------------

  template< class ct, int dim, int codim >
  struct FieldTraits< SPGeometryCache< ct, dim, codim >::JacobianTransposed >
  {
    typedef typename FieldTraits< ct >::field_type field_type;
    typedef typename FieldTraits< ct >::real_type real_type;
  };



  // SPGeometryCache::JacobianTransposed
  // -----------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianTransposed
  {
    typedef ct field_type;
    typedef ct value_type;

    typedef std::size_t size_type;

    static const int rows = mydimension;
    static const int cols = dimension;

    typedef Dune::FieldMatrix< field_type, rows, cols > FieldMatrix;
    typedef typename Dune::FieldTraits< field_type >::real_type real_type;

    JacobianTransposed ( const GlobalVector &h, Direction dir )
      : storage_( dir )
    {
      for( int k = 0; k < mydimension; ++k )
        storage_.h[ k ] = h[ pattern().nonzero( k ) ];
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
    const SPGeometryPattern< dimension, codimension > &pattern () const { return storage_; }
    const LocalVector h () const { return storage_.h; }

    MatrixStorage storage_;
  };



  // FieldTraits for SPGeometryCache::JacobianInverseTransposed
  // ----------------------------------------------------------

  template< class ct, int dim, int codim >
  struct FieldTraits< SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed >
  {
    typedef typename FieldTraits< ct >::field_type field_type;
    typedef typename FieldTraits< ct >::real_type real_type;
  };



  // SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------

  template< class ct, int dim, int codim >
  struct SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed
  {
    typedef ct field_type;
    typedef ct value_type;

    typedef std::size_t size_type;

    static const int rows = dimension;
    static const int cols = mydimension;

    typedef Dune::FieldMatrix< field_type, rows, cols > FieldMatrix;
    typedef typename Dune::FieldTraits< field_type >::real_type real_type;

    JacobianInverseTransposed ( const GlobalVector &h, Direction dir )
      : storage_( dir )
    {
      for( int k = 0; k < mydimension; ++k )
        storage_.h[ k ] = ctype( 1 ) / h[ pattern().nonzero( k ) ];
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
    const SPGeometryPattern< dimension, codimension > &pattern () const { return storage_; }
    const LocalVector hInv () const { return storage_.h; }

    MatrixStorage storage_;
  };



  // Implementation of SPGeometryPattern
  // -----------------------------------

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



  // Implementation of SPGeometryCache
  // ---------------------------------

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
  inline SPGeometryCache< ct, dim, codim >::JacobianTransposed::operator FieldMatrix () const
  {
    FieldMatrix matrix( ctype( 0 ) );
    for( int k = 0; k < mydimension; ++k )
      matrix[ k ][ pattern().nonzero( k ) ] = h()[ k ];
    return matrix;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] = h()[ k ] * x[ k ];
    for( int k = 0; k < codimension; ++k )
      y[ pattern().zero( k ) ] = ctype( 0 );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] = conjugateComplex( h()[ k ] ) * x[ k ];
    for( int k = 0; k < codimension; ++k )
      y[ pattern().zero( k ) ] = ctype( 0 );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::umhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::usmv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += alpha * h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::usmtv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += alpha * h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::usmhv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += alpha * conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= h()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] -= h()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::mmhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] -= conjugateComplex( h()[ k ] ) * x[ k ];
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::ctype
  SPGeometryCache< ct, dim, codim >::JacobianTransposed::det () const
  {
    ctype det( 1 );
    for( int k = 0; k < mydimension; ++k )
      det *= h()[ k ];
    return det;
  }



  // Implementation of SPGeometryCache::JacobianInverseTransposed
  // ------------------------------------------------------------

  template< class ct, int dim, int codim >
  inline SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::operator FieldMatrix () const
  {
    FieldMatrix matrix( ctype( 0 ) );
    for( int k = 0; k < mydimension; ++k )
      matrix[ pattern().nonzero( k ) ][ k ] = hInv()[ k ];
    return matrix;
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] = hInv()[ k ] * x[ k ];
    for( int k = 0; k < codimension; ++k )
      y[ pattern().zero( k ) ] = ctype( 0 );
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] = conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::umhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::usmv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] += alpha * hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::usmtv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += alpha * hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::usmhv ( const field_type &alpha, const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] += alpha * conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ pattern().nonzero( k ) ] -= hInv()[ k ] * x[ k ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmtv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= hInv()[ k ] * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  template< class X, class Y >
  inline void
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::mmhv ( const X &x, Y &y ) const
  {
    for( int k = 0; k < mydimension; ++k )
      y[ k ] -= conjugateComplex( hInv()[ k ] ) * x[ pattern().nonzero( k ) ];
  }


  template< class ct, int dim, int codim >
  inline typename SPGeometryCache< ct, dim, codim >::ctype
  SPGeometryCache< ct, dim, codim >::JacobianInverseTransposed::det () const
  {
    ctype det( 1 );
    for( int k = 0; k < mydimension; ++k )
      det *= hInv()[ k ];
    return det;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRYCACHE_HH
