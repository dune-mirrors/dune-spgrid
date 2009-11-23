#ifndef DUNE_SPGRID_MULTIINDEX_HH
#define DUNE_SPGRID_MULTIINDEX_HH

#include <dune/common/iostream.hh>
#include <dune/grid/spgrid/misc.hh>

namespace Dune
{

  // SPMultiIndex
  // ------------

  template< int dim >
  class SPMultiIndex
  {
    typedef SPMultiIndex< dim > This;

  public:
    static const int dimension = dim;

    SPMultiIndex ()
    {}

    SPMultiIndex ( const int (&index)[ dimension ] )
    {
      *this = index;
    }

    SPMultiIndex ( const This &other )
    {
      *this = other;
    }

    This &operator= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = other.index_[ i ];
      return *this;
    }

    This &operator= ( const int (&index)[ dimension ] )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = index[ i ];
      return *this;
    }

    This &operator+= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] += other.index_[ i ];
      return *this;
    }

    This &operator-= ( const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] -= other.index_[ i ];
      return *this;
    }

    This &operator*= ( const int a )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] *= a;
      return *this;
    }

    const int &operator[] ( const int i ) const
    {
      return index_[ i ];
    }

    int &operator[] ( const int i )
    {
      return index_[ i ];
    }

    bool operator== ( const This &other ) const
    {
      bool equals = true;
      for( int i = 0; i < dimension; ++i )
        equals &= (index_[ i ] == other.index_[ i ]);
      return equals;
    }
    
    bool operator!= ( const This &other ) const
    {
      bool equals = false;
      for( int i = 0; i < dimension; ++i )
        equals |= (index_[ i ] != other.index_[ i ]);
      return equals;
    }

    void axpy( const int a, const This &other )
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] += a*other.index_[ i ];
    }

    void clear ()
    {
      for( int i = 0; i < dimension; ++i )
        index_[ i ] = 0;
    }

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

    int codimension () const
    {
      int codim = dimension;
      for( int i = 0; i < dimension; ++i )
        codim -= (index_[ i ] & 1);
      return codim;
    }

    unsigned int direction () const
    {
      unsigned int dir = 0;
      for( int i = 0; i < dimension; ++i )
        dir |= (index_[ i ] & 1) << i;
      return dir;
    }

    static This zero ()
    {
      This zero;
      zero.clear();
      return zero;
    }

  private:
    int index_[ dimension ];
  };


  template< class char_type, class Traits, int dim >
  inline std::basic_ostream< char_type, Traits > &
  operator<< ( std::basic_ostream< char_type, Traits > &out, const SPMultiIndex< dim > &multiIndex )
  {
    out << "( " << multiIndex[ 0 ];
    for( int i = 1; i < dim; ++i )
      out << ", " << multiIndex[ i ];
    return out << " )";
  }


  template< class char_type, class Traits, int dim >
  inline std::basic_istream< char_type, Traits > &
  operator>> ( std::basic_istream< char_type, Traits > &in, SPMultiIndex< dim > &multiIndex )
  {
    SPMultiIndex< dim > m;
    in >> match( '(' ) >> m[ 0 ];
    for( int i = 1; i < dim; ++i )
      in >> match( ',' ) >> m[ i ];
    in >> match( ')' );
    if( in.good() )
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

}

#endif // #ifndef DUNE_SPGRID_MULTIINDEX_HH
