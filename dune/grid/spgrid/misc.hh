#ifndef DUNE_SPGRID_MISC_HH
#define DUNE_SPGRID_MISC_HH

/** \file
 *  \author Martin Nolte
 *  \brief  miscellaneous helper functions
 */

namespace Dune
{

  /** \brief count the number of set bits within an integer
   *
   *  \param[in]  i  integer whose set bits shall be counted
   *
   *  \returns number of set bits in the integer i
   */
  inline unsigned int bitCount ( unsigned int i )
  {
    unsigned int bitCount = 0;
    for( ; i != 0; i /= 2 )
      bitCount += (i & 1);
    return bitCount;
  }



  /** \brief obtain index with greater value in array
   *
   *  \param[in]  array  array of comparable values
   *  \param[in]  i      first index
   *  \param[in]  j      second index
   *
   *  \note The array values must implement the operator <.
   *
   *  \returns i if array[ i ] > array[ j ], j otherwise
   */
  template< class Array, class Index >
  inline Index argmax( const Array &array, Index i, Index j )
  {
    return (array[ j ] < array[ i ] ? i : j);
  }

  /** \brief obtain index with lesser value in array
   *
   *  \param[in]  array  array of comparable values
   *  \param[in]  i      first index
   *  \param[in]  j      second index
   *
   *  \note The array values must implement the operator <.
   *
   *  \returns i if array[ i ] < array[ j ], j otherwise
   */
  template< class Array, class Index >
  inline Index argmin( const Array &array, Index i, Index j )
  {
    return (array[ i ] < array[ j ] ? i : j);
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_MISC_HH
