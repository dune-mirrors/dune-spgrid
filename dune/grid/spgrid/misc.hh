#ifndef DUNE_SPGRID_MISC_HH
#define DUNE_SPGRID_MISC_HH

#include <utility>
#include <vector>

/** \file
 *  \author Martin Nolte
 *  \brief  miscellaneous helper functions
 */

namespace Dune
{

  /**
   * \brief count the number of set bits within an integer
   *
   * \param[in]  i  integer whose set bits shall be counted
   *
   * \returns number of set bits in the integer i
   */
  inline unsigned int bitCount ( unsigned int i )
  {
    unsigned int bitCount = 0;
    for( ; i != 0; i /= 2 )
      bitCount += (i & 1);
    return bitCount;
  }



  /**
   * \brief count the number of set bits within an integer
   *
   * \param[in]  i        integer whose set bits shall be counted
   * \param[in]  numBits  maximal number of used bits (starting from bit 0) in i
   *
   *  \returns number of set bits in the integer i
   */
  inline unsigned int bitCount ( unsigned int i, unsigned int numBits )
  {
    unsigned int bitCount = 0;
    for( unsigned int k = 0; k < numBits; ++k )
    {
      bitCount += (i & 1);
      i /= 2;
    }
    assert( i == 0 );
    return bitCount;
  }


  /**
   * \brief copy a vector, performing an operation on each element
   *
   * \param[in]  in  vector of input data
   * \param[in]  op  operation to perform on each element
   */

  template< class T, class Op >
  inline std::vector< decltype( std::declval< Op >()( std::declval< T >() ) ) >
  transform ( const std::vector< T > &in, Op op )
  {
    const std::size_t size = in.size();
    std::vector< T > out;
    out.reserve( size );
    for( const T &v : in )
      out.push_back( op( v ) );
    return std::move( out );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_MISC_HH
