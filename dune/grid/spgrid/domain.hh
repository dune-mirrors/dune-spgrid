#ifndef DUNE_SPGRID_DOMAIN_HH
#define DUNE_SPGRID_DOMAIN_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPDomain
  // --------

  /** \class SPDomain
   *  \brief description of the computational domain
   *
   *  \tparam  ct   coordinate type (e.g., double)
   *  \tparam  dim  dimension of the domain
   */
  template< class ct, int dim >
  class SPDomain
  {
    typedef SPDomain< ct, dim > This;

  public:
    /** \brief coordinate type */
    typedef ct ctype;

    /** \brief dimension of the domain */
    static const int dimension = dim;

    /** \brief type of global vectors, i.e., vectors within the domain */
    typedef FieldVector< ctype, dimension > GlobalVector;
    /** \brief type of multiindices */
    typedef SPMultiIndex< dimension > MultiIndex;

    /** \brief constructor
     *
     *  \param[in]  a         one corner of the domain
     *  \param[in]  b         the opposite corner of the domain
     *  \param[in]  periodic  bit field specifying which directions should be
     *                        periodic (defaults to 0)
     *
     *  \note The only restriction on the given corners is that they are
     *        opposite to each other.
     *        It is not guaranteed, that one of the corners will be returned
     *        by the method origin.
     */
    SPDomain ( const GlobalVector &a, const GlobalVector &b,
               const unsigned int periodic = 0 );

    /** \brief obtain lower left corner
     *
     *  \returns a reference to the origin of the domain
     */
    const GlobalVector &origin () const;

    /** \brief obtain width of the domain
     * 
     *  \returns a reference to the width of the domain
     */
    const GlobalVector &width () const;

    /** \brief determine whether the domain contains a point x
     *
     *  \param[in]  x  point to consider
     *
     *  \returns true, if x is contained in the domain
     */
    bool contains ( const GlobalVector &x ) const;

    /** \brief determine whether a direction is periodic
     *
     *  \param[in]  i  direction (0 <= i < dimension)
     *
     *  \returns true, if direction i is periodic
     */
    bool periodic ( const int i ) const;

    /** \brief obtain the periodicity bit field
     *
     *  \returns the bitfield specifying which directions are periodic
     */
    unsigned int periodic () const;

    /** \brief obtain a domain modelling the unit cube
     *
     *  \returns a domain modelling \f$[0,1]^{dim}\f$
     */
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
