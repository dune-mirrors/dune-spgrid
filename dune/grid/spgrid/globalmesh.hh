#ifndef DUNE_SPGRID_GLOBALMESH_HH
#define DUNE_SPGRID_GLOBALMESH_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  template< class ct, int dim >
  class GlobalMesh
  {
    typedef GlobalMesh< ct, dim > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    explicit GlobalMesh ( const GlobalVector &h, const GlobalVector &origin = GlobalVector( 0 ) )
    : h_( h ),
      origin_( origin )
    {}

    GlobalVector coordinate ( const MultiIndex &id ) const
    {
      GlobalVector coordinate( origin_ );
      for( int i = 0; i < dimension; ++i )
        coordinate[ i ] += (id[ i ] / 2) * h_[ i ];
      return coordinate;
    }

    template< SPRefinementStrategy strategy >
    This refine ( const SPRefinement< dimension, strategy > &refinement ) const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = h_[ i ] / ctype( refinement.factor( i ) );
      return This( h, origin_ );
    }

  private:
    GlobalVector h_;
    GlobalVector origin_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GLOBALMESH_HH
