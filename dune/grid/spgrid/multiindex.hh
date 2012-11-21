#ifndef DUNE_SPGRID_MULTIINDEX_HH
#define DUNE_SPGRID_MULTIINDEX_HH

#include <limits>

#include <dune/common/array.hh>

#include <dune/common/iostream.hh>

#include <dune/grid/spgrid/misc.hh>

namespace Dune
{

  // SPMultiIndex
  // ------------

  /** \class SPMultiIndex
   *  \brief multiindex
   *
   *  A structured grid is most easily addressed using a multiindex, which is
   *  realized through this class.
   *
   *  \todo describe increment, codimension and direction - i.e., numbering concept
   *
   *  \tparam  dim  dimension of the multiindex
   */
  template< int dim >
  class SPMultiIndex
  {
    typedef SPMultiIndex< dim > This;

  public:
    /** \brief dimension of the multiindex */
    static const int dimension = dim;

    /** \brief default constructor
     *
     *  \note The default constructor does not initialze the multiindex.
     */
    SPMultiIndex ()
    {}

    /** \brief constructor from int array
     *
     *  \note This constructor defines an implicit conversion.
     *
     *  \param[in]  index  int array to copy
     */
    SPMultiIndex ( const int (&index)[ dimension ] )
    {
      *this = index;
    }

    /** \brief constructor from int array
     *
     *  \note This constructor defines an implicit conversion.
     *
     *  \param[in]  index  int array to copy
     */
    SPMultiIndex ( const Dune::array< int, dimension > &index )
    {
      *this = index;
    }

    /** \brief copy constructor */
    SPMultiIndex ( const This &other )
    {
      *this = other;
    }

    /** \brief assignment operator */
    const This &operator= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = other.index_[ i ];
      return *this;
    }

    /** \brief assignment operator from int array */
    const This &operator= ( const int (&index)[ dimension ] )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = index[ i ];
      return *this;
    }

    const This &operator= ( const Dune::array< int, dimension > &index )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = index[ i ];
      return *this;
    }

    /** \brief add another multiindex to this one (vector operation) */
    const This &operator+= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] += other.index_[ i ];
      return *this;
    }

    /** \brief subtract another multiindex from this one (vector operation) */
    const This &operator-= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] -= other.index_[ i ];
      return *this;
    }

    /** \brief scale this multiindex (vector operation) */
    This &operator*= ( const int a )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] *= a;
      return *this;
    }

    /** \brief scale this multiindex (vector operation) */
    This &operator/= ( const int a )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] /= a;
      return *this;
    }

    /** \brief access i-th component */
    const int &operator[] ( const int i ) const
    {
      return index_[ i ];
    }

    /** \brief access i-th component */
    int &operator[] ( const int i )
    {
      return index_[ i ];
    }

    /** \brief compare two multiindices for equality */
    bool operator== ( const This &other ) const
    {
      bool equals = true;
      for( int i = 0; i < dimension; ++i )
        equals &= (index_[ i ] == other.index_[ i ]);
      return equals;
    }
    
    /** \brief compare two multiindices for inequality */
    bool operator!= ( const This &other ) const
    {
      bool equals = false;
      for( int i = 0; i < dimension; ++i )
        equals |= (index_[ i ] != other.index_[ i ]);
      return equals;
    }

    /** \brief add multiple of a multiindex to this one (vector operation) */
    void axpy( const int a, const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] += a*other.index_[ i ];
    }

    /** \brief initialize to zero */
    void clear ()
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = 0;
    }

    /** \todo please doc me */
    void increment ( const This &bound, const int k = 1 )
    {
      for( int i = 0; i < dimension; ++i )
      {
        index_[ i ] += k;
        if( index_[ i ] < bound[ i ] )
          return;
        index_[ i ] = 0;
      }
    }

    /** \todo please doc me */
    int codimension () const
    {
      int codim = dimension;
      for( int i = 0; i < dimension; ++i )
        codim -= (index_[ i ] & 1);
      return codim;
    }

    /** \todo please doc me */
    unsigned int direction () const
    {
      unsigned int dir = 0;
      for( int i = 0; i < dimension; ++i )
        dir |= (index_[ i ] & 1) << i;
      return dir;
    }

    /** \brief obtain the zero multiindex */
    static This zero ()
    {
      This zero;
      zero.clear();
      return zero;
    }

  private:
    int index_[ dimension ];
  };



  // Auxilliary Functions for SPMultiIndex
  // -------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPMultiIndex< dim > &multiIndex )
  {
    out << "( " << multiIndex[ 0 ];
    for( int i = 1; i < dim; ++i )
      out << ", " << multiIndex[ i ];
    return out << " )";
  }


  template< class char_type, class traits, int dim >
  inline std::basic_istream< char_type, traits > &
  operator>> ( std::basic_istream< char_type, traits > &in, SPMultiIndex< dim > &multiIndex )
  {
    SPMultiIndex< dim > m;
    in >> match( '(' ) >> m[ 0 ];
    for( int i = 1; i < dim; ++i )
      in >> match( ',' ) >> m[ i ];
    in >> match( ')' );
    if( !in.fail() )
      multiIndex = m;
    return in;
  }


  template< int dim >
  inline SPMultiIndex< dim >
  operator+ ( const SPMultiIndex< dim > &a, const SPMultiIndex< dim > &b )
  {
    SPMultiIndex< dim > c = a;
    c += b;
    return c;
  }


  template< int dim >
  inline SPMultiIndex< dim >
  operator- ( const SPMultiIndex< dim > &a, const SPMultiIndex< dim > &b )
  {
    SPMultiIndex< dim > c = a;
    c -= b;
    return c;
  }


  template< int dim >
  inline SPMultiIndex< dim >
  operator* ( const SPMultiIndex< dim > &a, const int &b )
  {
    SPMultiIndex< dim > c = a;
    c *= b;
    return c;
  }

  
  template< int dim >
  inline SPMultiIndex< dim >
  operator* ( const int &a, const SPMultiIndex< dim > &b )
  {
    return (b*a);
  }


  template< int dim >
  inline SPMultiIndex< dim >
  operator/ ( const SPMultiIndex< dim > &a, const int &b )
  {
    SPMultiIndex< dim > c = a;
    c /= b;
    return c;
  }


  template< int dim >
  inline int argmax( const SPMultiIndex< dim > &multiIndex )
  {
    int m = 0;
    for( int i = 1; i < dim; ++i )
      m = argmax( multiIndex, i, m );
    return m;
  }


  template< int dim >
  inline int argmin( const SPMultiIndex< dim > &multiIndex )
  {
    int m = 0;
    for( int i = 1; i < dim; ++i )
      m = argmin( multiIndex, i, m );
    return m;
  }

} // namespace Dune


namespace std
{

  // Auxilliary functions for SPMultiIndex
  // -------------------------------------

  template< int dim >
  inline Dune::SPMultiIndex< dim >
  min ( const Dune::SPMultiIndex< dim > &a, const Dune::SPMultiIndex< dim > &b )
  {
    Dune::SPMultiIndex< dim > c;
    for( int i = 0; i < dim; ++i )
      c[ i ] = min( a[ i ], b[ i ] );
    return c;
  }


  template< int dim >
  inline Dune::SPMultiIndex< dim >
  max ( const Dune::SPMultiIndex< dim > &a, const Dune::SPMultiIndex< dim > &b )
  {
    Dune::SPMultiIndex< dim > c;
    for( int i = 0; i < dim; ++i )
      c[ i ] = max( a[ i ], b[ i ] );
    return c;
  }



  // numeric_limits for SPMultiIndex
  // -------------------------------

  template< int dim >
  struct numeric_limits< Dune::SPMultiIndex< dim > >
  {
    typedef Dune::SPMultiIndex< dim > MultiIndex;

    static MultiIndex min ()
    {
      MultiIndex multiIndex;
      for( int i = 0; i < dim; ++i )
        multiIndex[ i ] = numeric_limits< int >::min();
      return multiIndex;
    }

    static MultiIndex max ()
    {
      MultiIndex multiIndex;
      for( int i = 0; i < dim; ++i )
        multiIndex[ i ] = numeric_limits< int >::max();
      return multiIndex;
    }
  };

} // namespace std

#endif // #ifndef DUNE_SPGRID_MULTIINDEX_HH
