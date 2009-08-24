#ifndef DUNE_SPGRID_MISC_HH
#define DUNE_SPGRID_MISC_HH

namespace Dune
{

  inline unsigned int bitCount ( unsigned int i )
  {
    unsigned int bitCount = 0;
    for( ; i != 0; i /= 2 )
      bitCount += (i & 1);
    return bitCount;
  }



  template< class Array, class Index >
  inline Index argmax( const Array &array, Index i, Index j )
  {
    return (array[ i ] > array[ j ] ? i : j);
  }

  template< class Array, class Index >
  inline Index argmin( const Array &array, Index i, Index j )
  {
    return (array[ i ] < array[ j ] ? i : j);
  }

}

#endif // #ifndef DUNE_SPGRID_MISC_HH
