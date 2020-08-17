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
    return out;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_MISC_HH
