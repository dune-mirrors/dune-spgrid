#ifndef DUNE_GRID_SPGRID_DIRECTION_HH
#define DUNE_GRID_SPGRID_DIRECTION_HH

#include <cassert>

#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPDirection
  // -----------

  template< int dim >
  class SPDirection
  {
    typedef SPDirection< dim > This;

    static_assert( (dim >= 0) && (dim < 8*sizeof( unsigned int )), "Invalid dimension." );

    template< int, int >
    friend class SPEntityDirection;
    template< int, int >
    friend class SPDirectionIterator;

  public:
    static const int dimension = dim;

    explicit SPDirection ( const SPMultiIndex< dimension > &id );

  private:
    explicit SPDirection ( unsigned int bits ) : bits_( bits ) { assert( (bits_ & ~((1u << dimension) - 1u)) == 0 ); }

  public:
    bool operator== ( const This &other ) const { return (bits_ == other.bits_); }
    bool operator!= ( const This &other ) const { return (bits_ != other.bits_); }

    unsigned int operator[] ( int i ) const { return ((bits_ >> i) & 1u); }

    int mydimension () const { return bitCount( bits_, dim ); }

    unsigned int bits () const { return bits_; }

  private:
    unsigned int bits_;
  };



  // SPBasicEntityDirection
  // ----------------------

  template< int dim, int mydim >
  class SPBasicEntityDirection
  {
    typedef SPBasicEntityDirection< dim, mydim > This;

    static_assert( (mydim >= 0) && (mydim <= dim), "Invalid entity dimension." );

  public:
    static const int dimension = dim;

    typedef SPDirection< dimension > Direction;

    int mydimension () const { return mydim; }
  };



  // SPEntityDirection
  // -----------------

  template< int dim, int mydim >
  class SPEntityDirection
    : public SPBasicEntityDirection< dim, mydim >
  {
    typedef SPEntityDirection< dim, mydim > This;
    typedef SPBasicEntityDirection< dim, mydim > Base;

  public:
    using Base::dimension;

    typedef typename Base::Direction Direction;

    SPEntityDirection () : direction_( 0u ) {}
    explicit SPEntityDirection ( const SPMultiIndex< dimension > &id ) : direction_( id ) {}

    operator Direction () const { return direction_; }

  private:
    Direction direction_;
  };


  template< int dim >
  class SPEntityDirection< dim, 0 >
    : public SPBasicEntityDirection< dim, 0 >
  {
    typedef SPEntityDirection< dim, 0 > This;
    typedef SPBasicEntityDirection< dim, 0 > Base;

  public:
    using Base::dimension;

    typedef typename Base::Direction Direction;

    SPEntityDirection () = default;
    explicit SPEntityDirection ( const SPMultiIndex< dimension > &id ) { assert( Direction( id ) == static_cast< Direction >( This() ) ); }

    operator Direction () const { return Direction( 0u ); }
  };


  template< int dim >
  class SPEntityDirection< dim, dim >
    : public SPBasicEntityDirection< dim, dim >
  {
    typedef SPEntityDirection< dim, dim > This;
    typedef SPBasicEntityDirection< dim, dim > Base;

  public:
    using Base::dimension;

    typedef typename Base::Direction Direction;

    SPEntityDirection () = default;
    explicit SPEntityDirection ( const SPMultiIndex< dimension > &id ) { assert( Direction( id ) == static_cast< Direction >( This() ) ); }

    operator Direction () const { return Direction( (1u << dim) - 1u ); }
  };


  template< >
  class SPEntityDirection< 0, 0 >
    : public SPBasicEntityDirection< 0, 0 >
  {
    typedef SPEntityDirection< 0, 0 > This;
    typedef SPBasicEntityDirection< 0, 0 > Base;

  public:
    using Base::dimension;

    typedef Base::Direction Direction;

    SPEntityDirection () = default;
    explicit SPEntityDirection ( const SPMultiIndex< dimension > &id ) { assert( Direction( id ) == static_cast< Direction >( This() ) ); }

    operator Direction () const { return Direction( 0u ); }
  };



  // SPDirectionIterator
  // -------------------

  template< int dim, int codim >
  class SPDirectionIterator
  {
    typedef SPDirectionIterator< dim, codim > This;

  public:
    typedef SPDirection< dim > Direction;

    SPDirectionIterator () : direction_( (1u << (dim - codim)) - 1u ) {}
    explicit SPDirectionIterator ( const Direction &direction ) : direction_( direction ) {}

    SPDirectionIterator ( const This & ) = default;
    SPDirectionIterator ( This && ) = default;

    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

    operator bool () const { return (direction_.bits_ != (1u << dim) - (1u << codim) + 1u); }

    const Direction &operator* () const { return direction_; }
    const Direction *operator-> () const { return &direction_; }

    const This &operator++ ();

  private:
    Direction direction_;
  };



  // Implementation of SPDirection
  // -----------------------------

  template< int dim >
  const int SPDirection< dim >::dimension;


  template< int dim >
  inline SPDirection< dim >
  ::SPDirection ( const SPMultiIndex< dimension > &id )
    : bits_( 0 )
  {
    for( int i = 0; i < dimension; ++i )
      bits_ |= ((unsigned int)( id[ i ] & 1 ) << i);
  }



  // Implementation of SPBasicEntityDirection
  // ----------------------------------------

  template< int dim, int mydim >
  const int SPBasicEntityDirection< dim, mydim >::dimension;



  // Implementation of SPDirectionIterator
  // -------------------------------------

  template< int dim, int codim >
  inline const typename SPDirectionIterator< dim, codim >::This &
  SPDirectionIterator< dim, codim >::operator++ ()
  {
    assert( *this );
    do {
      ++direction_.bits_;
    } while( *this && (direction_.mydimension() != (dim - codim)) );
    return *this;
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_DIRECTION_HH
