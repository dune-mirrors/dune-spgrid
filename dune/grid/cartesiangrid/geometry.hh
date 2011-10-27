#ifndef DUNE_CARTESIANGRID_GEOMETRY_HH
#define DUNE_CARTESIANGRID_GEOMETRY_HH

#include <dune/grid/common/geometry.hh>
#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // CartesianGridGeometry
  // ---------------------

  template< int mydim, int cdim, class Grid >
  class CartesianGridGeometry
  : public SPBasicGeometry< mydim, cdim, Grid, CartesianGridGeometry< mydim, cdim, Grid > >
  {
    typedef CartesianGridGeometry< mydim, cdim, Grid > This;
    typedef SPBasicGeometry< mydim, cdim, Grid, This > Base;

  protected:
    typedef typename Base::Traits Traits;

  public:
    typedef typename Base::ctype ctype;

    static const int mydimension = Base::mydimension;
    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;

    typedef typename Base::ReferenceCube ReferenceCube;
    typedef typename Base::GeometryCache GeometryCache;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;

  public:
    template< class GeometricGridLevel, class OriginVector >
    CartesianGridGeometry ( const GeometricGridLevel& gridLevel,
                            const unsigned int dir,
                            const OriginVector &origin )
    : refCube_( &gridLevel.template referenceCube< codimension >() ),
      geometryCache_( &gridLevel.template geometryCache< codimension >( dir ) ),
      origin_()
    {
      // copy by hand, since OriginVector could differ from GlobalVector
      for( int i=0; i < GlobalVector::dimension; ++i )
        origin_[ i ] = origin[ i ]; 
    }

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const ReferenceCube &referenceCube () const
    {
      assert( refCube_ );
      return *refCube_ ;
    }

    const GlobalVector& origin () const
    {
      return origin_;
    }

    const GeometryCache &geometryCache () const
    {
      assert( geometryCache_ );
      return *geometryCache_;
    }

  private:
    const ReferenceCube* refCube_;
    const GeometryCache *geometryCache_; 
    GlobalVector origin_;
  };

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_GEOMETRY_HH
