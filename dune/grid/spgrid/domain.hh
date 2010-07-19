#ifndef DUNE_SPGRID_DOMAIN_HH
#define DUNE_SPGRID_DOMAIN_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/topology.hh>
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

    typedef SPTopology< dimension > Topology;

    struct Cube;

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

    /** \todo please doc me */
    const Cube &cube () const { return cube_; }

    /** \brief determine whether the domain contains a point x
     *
     *  \param[in]  x  point to consider
     *
     *  \returns true, if x is contained in the domain
     */
    bool contains ( const GlobalVector &x ) const { return cube().contains( x ); }

    /** \brief determine whether a direction is periodic
     *
     *  \param[in]  i  direction (0 <= i < dimension)
     *
     *  \returns true, if direction i is periodic
     */
    bool periodic ( const int i ) const { return topology_.periodic( i ); }

    /** \brief obtain the periodicity bit field
     *
     *  \returns the bitfield specifying which directions are periodic
     */
    unsigned int periodic () const { return topology_.periodic(); }

    /** \brief obtain a domain modelling the unit cube
     *
     *  \returns a domain modelling \f$[0,1]^{dim}\f$
     */
    static This unitCube ();

  private:
    Cube cube_;
    Topology topology_;
  };



  // SPDomain::Cube
  // --------------

  template< class ct, int dim >
  class SPDomain< ct, dim >::Cube
  {
    friend class SPDomain< ct, dim >;

    typedef SPDomain< ct, dim > Domain;

  public:
    /** \brief dimension of the domain */
    static const int dimension = Domain::dimension;

    /** \brief type of global vectors, i.e., vectors within the domain */
    typedef typename Domain::GlobalVector GlobalVector;

    /** \brief constructor
     *
     *  \param[in]  a         one corner of the cube
     *  \param[in]  b         the opposite corner of the cube
     *
     *  \note The only restriction on the given corners is that they are
     *        opposite to each other.
     *        It is not guaranteed, that one of the corners will be returned
     *        by the method origin.
     */
    Cube ( const GlobalVector &a, const GlobalVector &b );

    /** \brief obtain lower left corner
     *
     *  \returns a reference to the origin of the cube
     */
    const GlobalVector &origin () const;

    /** \brief obtain width
     * 
     *  \returns a reference to the width of the cube
     */
    const GlobalVector &width () const;

    /** \brief determine whether the cube contains a point x
     *
     *  \param[in]  x  point to consider
     *
     *  \returns true, if x is contained in the cube
     */
    bool contains ( const GlobalVector &x ) const;

  private:
    GlobalVector origin_, width_;
  };



  // Implementation of SPDomain
  // --------------------------

  template< class ct, int dim >
  inline SPDomain< ct, dim >
    ::SPDomain ( const GlobalVector &a, const GlobalVector &b,
                 const unsigned int periodic )
  : cube_( a, b ),
    topology_( periodic )
  {}


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



  // Implementation of SPDomain::Cube
  // --------------------------------

  template< class ct, int dim >
  inline SPDomain< ct, dim >::Cube
    ::Cube ( const GlobalVector &a, const GlobalVector &b )
  {
    for( int i = 0; i < dimension; ++i )
    {
      origin_[ i ] = std::min( a[ i ], b[ i ] );
      width_[ i ] = std::max( a[ i ], b[ i ] ) - origin_[ i ];
    }
  }


  template< class ct, int dim >
  inline const typename SPDomain< ct, dim >::Cube::GlobalVector &
  SPDomain< ct, dim >::Cube::origin () const
  {
    return origin_;
  }


  template< class ct, int dim >
  inline const typename SPDomain< ct, dim >::Cube::GlobalVector &
  SPDomain< ct, dim >::Cube::width () const
  {
    return width_;
  }


  template< class ct, int dim >
  inline bool SPDomain< ct, dim >::Cube::contains ( const GlobalVector &x ) const
  {
    bool contains = true;
    for( int i = 0; i < dimension; ++i )
    {
      const ctype y = x[ i ] - origin()[ i ];
      contains &= ((y >= 0) && (y <= width()[ i ]));
    }
    return contains;
  }

}

#endif // #ifndef DUNE_SPGRID_DOMAIN_HH
