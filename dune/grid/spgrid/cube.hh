#ifndef DUNE_SPGRID_CUBE_HH
#define DUNE_SPGRID_CUBE_HH

#include <dune/common/fvector.hh>

#include <dune/common/iostream.hh>

namespace Dune
{

  // SPCube
  // ------

  template< class ct, int dim >
  class SPCube
  {
    typedef SPCube< ct, dim > This;

  public:
    /** \brief coordinate type */
    typedef ct ctype;

    /** \brief dimension of the domain */
    static const int dimension = dim;

    /** \brief type of global vectors, i.e., vectors within the domain */
    typedef FieldVector< ctype, dimension > GlobalVector;

    /** \brief default constructor */
    SPCube ();

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
    SPCube ( const GlobalVector &a, const GlobalVector &b );

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

    /** \brief obtain a domain modelling the unit cube
     *
     *  \returns a domain modelling \f$[0,1]^{dim}\f$
     */
    static This unitCube ();

  private:
    GlobalVector origin_, width_;
  };



  // Implementation of SPCube
  // ------------------------

  template< class ct, int dim >
  inline SPCube< ct, dim >::SPCube ()
  {
    for( int i = 0; i < dimension; ++i )
      origin_[ i ] = width_[ i ] = 0;
  }


  template< class ct, int dim >
  inline SPCube< ct, dim >
    ::SPCube ( const GlobalVector &a, const GlobalVector &b )
  {
    for( int i = 0; i < dimension; ++i )
    {
      origin_[ i ] = std::min( a[ i ], b[ i ] );
      width_[ i ] = std::max( a[ i ], b[ i ] ) - origin_[ i ];
    }
  }


  template< class ct, int dim >
  inline const typename SPCube< ct, dim >::GlobalVector &
  SPCube< ct, dim >::origin () const
  {
    return origin_;
  }


  template< class ct, int dim >
  inline const typename SPCube< ct, dim >::GlobalVector &
  SPCube< ct, dim >::width () const
  {
    return width_;
  }


  template< class ct, int dim >
  inline bool SPCube< ct, dim >::contains ( const GlobalVector &x ) const
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
  inline typename SPCube< ct, dim >::This
  SPCube< ct, dim >::unitCube ()
  {
    GlobalVector a, b;
    for( int i = 0; i < dimension; ++i )
    {
      a = ctype( 0 );
      b = ctype( 1 );
    }
    return This( a, b );
  }



  // Auxilliary Functions for SPCube
  // -------------------------------

  template< class char_type, class traits, class ct, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out,
               const SPCube< ct, dim > &cube )
  {
    typedef SPCube< ct, dim > Cube;
    typename Cube::GlobalVector a = cube.origin();
    typename Cube::GlobalVector b = a + cube.width();
    for( int i = 0; i < Cube::dimension; ++i )
      out << (i > 0 ? "x[" : "[") << a[ i ] << "," << b[ i ] << "]";
    return out;
  }


  template< class char_type, class traits, class ct, int dim >
  inline std::basic_istream< char_type, traits > &
  operator>> ( std::basic_istream< char_type, traits > &in,
               SPCube< ct, dim > &cube )
  {
    typedef SPCube< ct, dim > Cube;
    typename Cube::GlobalVector a;
    typename Cube::GlobalVector b;
    for( int i = 0; i < Cube::dimension; ++i )
    {
      if( i > 0 )
        in >> match( 'x' );
      in >> match( '[' ) >> a[ i ] >> match( ',' ) >> b[ i ] >> match( ']' );
    }
    if( !in.fail() )
      cube = Cube( a, b );
    return in;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_CUBE_HH
