#ifndef DUNE_SPGRID_NORMAL_HH
#define DUNE_SPGRID_NORMAL_HH

#include <dune/common/fvector.hh>

namespace Dune
{

  // SPNormalVector
  // --------------

  template< class ct, int dim >
  class SPNormalVector
  {
    typedef SPNormalVector< ct, dim > This;

  public:
    static const int dimension = dim;

    typedef ct field_type;

    typedef Dune::FieldVector< ct, dim > FieldVector;

    SPNormalVector ( const int i, const field_type &p );

    operator FieldVector () const;

    This &operator*= ( const field_type &s );
    This &operator/= ( const field_type &s );

    bool operator== ( const This &other ) const;
    bool operator!= ( const This &other ) const;

    field_type operator* ( const This &other ) const;
    field_type operator* ( const FieldVector &other ) const;

    field_type one_norm () const;
    field_type two_norm () const;
    field_type two_norm2 () const;
    field_type infinity_norm () const;

  private:
    int i_;
    field_type p_;
  };



  // Auxilliary Functions for SPNormalVector
  // ---------------------------------------

  template< class ct, int dim >
  inline ct operator* ( const FieldVector< ct, dim > &a, const SPNormalVector< ct, dim > &b )
  {
    return b*a;
  }


  template< class ct, int dim >
  inline SPNormalVector< ct, dim > operator* ( const ct &a, const SPNormalVector< ct, dim > &b )
  {
    SPNormalVector< ct, dim > c( b );
    c *= a;
    return c;
  }


  template< class ct, int dim >
  inline SPNormalVector< ct, dim > operator* ( const SPNormalVector< ct, dim > &a, const ct &b )
  {
    return b*a;
  }


  template< class ct, int dim >
  inline bool operator== ( const FieldVector< ct, dim > &a, const SPNormalVector< ct, dim > &b )
  {
    return (a == static_cast< FieldVector< ct, dim > >( b ));
  }


  template< class ct, int dim >
  inline bool operator!= ( const FieldVector< ct, dim > &a, const SPNormalVector< ct, dim > &b )
  {
    return (a != static_cast< FieldVector< ct, dim > >( b ));
  }


  template< class ct, int dim >
  inline bool operator== ( const SPNormalVector< ct, dim > &a, const FieldVector< ct, dim > &b )
  {
    return (static_cast< FieldVector< ct, dim > >( a ) == b);
  }


  template< class ct, int dim >
  inline bool operator!= ( const SPNormalVector< ct, dim > &a, const FieldVector< ct, dim > &b )
  {
    return (static_cast< FieldVector< ct, dim > >( a ) != b);
  }



  // Implementation of SPNormalVector
  // --------------------------------

  template< class ct, int dim >
  inline SPNormalVector< ct, dim >
    ::SPNormalVector ( const int i, const field_type &p )
  : i_( i ), p_( p )
  {}


  template< class ct, int dim >
  inline SPNormalVector< ct, dim >::operator FieldVector () const
  {
    FieldVector v( field_type( 0 ) );
    v[ i_ ] = p_;
    return v;
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::This &
  SPNormalVector< ct, dim >::operator*= ( const field_type &s )
  {
    p_ *= s;
    return *this;
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::This &
  SPNormalVector< ct, dim >::operator/= ( const field_type &s )
  {
    p_ /= s;
    return *this;
  }


  template< class ct, int dim >
  inline bool SPNormalVector< ct, dim >::operator== ( const This &other ) const
  {
    return (i_ == other.i_) && (p_ == other.p_);
  }


  template< class ct, int dim >
  inline bool SPNormalVector< ct, dim >::operator!= ( const This &other ) const
  {
    return (i_ != other.i_) || (p_ != other.p_);
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::operator* ( const This &other ) const
  {
    return (i_ == other.i_ ? p_ * other.p_ : field_type( 0 ));
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::operator* ( const FieldVector &other ) const
  {
    return p_ * other[ i_ ];
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::one_norm () const
  {
    return std::abs( p_ );
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::two_norm () const
  {
    return std::abs( p_ );
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::two_norm2 () const
  {
    return p_ * p_;
  }


  template< class ct, int dim >
  inline typename SPNormalVector< ct, dim >::field_type
  SPNormalVector< ct, dim >::infinity_norm () const
  {
    return std::abs( p_ );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_NORMAL_HH
