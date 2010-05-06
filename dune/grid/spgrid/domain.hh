#ifndef DUNE_SPGRID_DOMAIN_HH
#define DUNE_SPGRID_DOMAIN_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPDomain
  // --------

  template< class ct, int dim >
  class SPDomain
  {
    typedef SPDomain< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    SPDomain ( const GlobalVector &a, const GlobalVector &b,
               const unsigned int periodic = 0 );

    const GlobalVector &origin () const;
    const GlobalVector &width () const;

    bool contains ( const GlobalVector &x ) const;

    bool periodic ( const int i ) const;
    unsigned int periodic () const;

    static This unitCube ();

  private:
    GlobalVector origin_, width_;
    unsigned int periodic_;
  };



  // Implementation of SPDomain
  // --------------------------

  template< class ct, int dim >
  inline SPDomain< ct, dim >
    ::SPDomain ( const GlobalVector &a, const GlobalVector &b,
                 const unsigned int periodic )
  : periodic_( periodic )
  {
    for( int i = 0; i < dimension; ++i )
    {
      origin_[ i ] = std::min( a[ i ], b[ i ] );
      width_[ i ] = std::max( a[ i ], b[ i ] ) - origin_[ i ];
    }
  }


  template< class ct, int dim >
  inline const typename SPDomain< ct, dim >::GlobalVector &
  SPDomain< ct, dim >::origin () const
  {
    return origin_;
  }


  template< class ct, int dim >
  inline const typename SPDomain< ct, dim >::GlobalVector &
  SPDomain< ct, dim >::width () const
  {
    return width_;
  }


  template< class ct, int dim >
  inline bool SPDomain< ct, dim >::contains ( const GlobalVector &x ) const
  {
    bool contains = true;
    for( int i = 0; i < dimension; ++i )
    {
      const ctype y = x[ i ] - origin()[ i ];
      contains &= ((y >= 0) && (y <= width()[ i ]));
    }
    return contains;
  }


  template< class ct, int dim >
  inline bool SPDomain< ct, dim >::periodic ( const int i ) const
  {
    assert( (i >= 0) && (i < dimension) );
    return ((periodic_ & (1 << i)) != 0);
  }


  template< class ct, int dim >
  inline unsigned int SPDomain< ct, dim >::periodic () const
  {
    return periodic_;
  }


  template< class ct, int dim >
  inline typename SPDomain< ct, dim >::This
  SPDomain< ct, dim >::unitCube ()
  {
    GlobalVector a, b;
    for( int i = 0; i < dimension; ++i )
    {
      a = ctype( 0 );
      b = ctype( 1 );
    }
    return This( a, b );
  }

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH
