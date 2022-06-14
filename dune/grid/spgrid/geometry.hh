#ifndef DUNE_SPGRID_GEOMETRY_HH
#define DUNE_SPGRID_GEOMETRY_HH

#include <type_traits>

#include <dune/geometry/type.hh>

#include <dune/common/typetraits.hh>
#include <dune/common/transpose.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/spgrid/referencecube.hh>
#include <dune/grid/spgrid/entityinfo.hh>

namespace Dune
{

  // SPBasicGeometry
  // ---------------

  template< int mydim, int cdim, class Grid, class Impl >
  class SPBasicGeometry
  {
    typedef SPBasicGeometry< mydim, cdim, Grid, Impl > This;

  protected:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::ReferenceCube::ctype ctype;

    static const int mydimension = mydim;
    static const int coorddimension = cdim;
    static const int dimension = Traits::ReferenceCube::dimension;
    static const int codimension = dimension - mydimension;

    typedef SPReferenceCube< ctype, mydimension > ReferenceCube;
    typedef SPGeometryCache< ctype, dimension, codimension > GeometryCache;

    static const int numCorners = (1 << mydimension);

    typedef typename GeometryCache::GlobalVector GlobalVector;
    typedef typename GeometryCache::LocalVector LocalVector;

    // just to make the Dune interface happy
    typedef GlobalVector GlobalCoordinate;
    typedef LocalVector LocalCoordinate;

    typedef typename GeometryCache::JacobianTransposed JacobianTransposed;
    typedef typename GeometryCache::JacobianInverseTransposed JacobianInverseTransposed;

    using Jacobian = std::decay_t<decltype(transposedView(std::declval<const JacobianTransposed&>()))>;
    using JacobianInverse = std::decay_t<decltype(transposedView(std::declval<const JacobianInverseTransposed&>()))>;

  protected:
    SPBasicGeometry ()
    {}

  public:
    GeometryType type () const { return GeometryTypes::cube( mydimension ); }

    int corners () const { return numCorners; }
    GlobalVector corner ( const int i ) const { return global( ReferenceCube::corner( i ) ); }
    GlobalVector center () const { return global( ReferenceCube::center() ); }

    bool affine () const { return true; }

    GlobalVector global ( const LocalVector &local ) const;
    LocalVector local ( const GlobalVector &global ) const;

    ctype volume () const { return asImpl().geometryCache().volume(); }
    ctype integrationElement ( const LocalVector &local ) const { return volume(); }

    const JacobianTransposed &jacobianTransposed ( const LocalVector &local ) const;
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalVector &local ) const;

    auto jacobian ( const LocalVector &local ) const;
    auto jacobianInverse ( const LocalVector &local ) const;

  protected:
    const Impl &asImpl () const { return static_cast< const Impl & >( *this ); }
  };



  // SPGeometry
  // ----------

  template< int mydim, int cdim, class Grid >
  class SPGeometry
  : public SPBasicGeometry< mydim, cdim, Grid, SPGeometry< mydim, cdim, Grid > >
  {
    typedef SPGeometry< mydim, cdim, Grid > This;
    typedef SPBasicGeometry< mydim, cdim, Grid, This > Base;

  protected:
    typedef typename Base::Traits Traits;

  public:
    typedef typename Base::ctype ctype;

    static const int mydimension = Base::mydimension;
    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;

    typedef __SPGrid::EntityInfo< Grid, codimension > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;

    typedef typename Base::ReferenceCube ReferenceCube;
    typedef typename Base::GeometryCache GeometryCache;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    explicit SPGeometry ( const GridLevel &gridLevel )
    : entityInfo_( gridLevel )
    {}

    explicit SPGeometry ( const EntityInfo &entityInfo )
    : entityInfo_( entityInfo ),
      origin_( computeOrigin() )
    {}

    SPGeometry ( const GridLevel &gridLevel, const MultiIndex &id )
    : entityInfo_( gridLevel, id ),
      origin_( computeOrigin() )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    using Base::jacobian;
    using Base::jacobianInverse;

    GlobalVector origin () const { return origin_; }
    const GeometryCache &geometryCache () const { return entityInfo().geometryCache(); }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

  private:
    GlobalVector computeOrigin () const
    {
      const GlobalVector &h = gridLevel().h();
      GlobalVector origin = gridLevel().domain().cube().origin();
      for( int i = 0; i < dimension; ++i )
        origin[ i ] += (entityInfo().id()[ i ] / 2) * h[ i ];
      return origin;
    }

    EntityInfo entityInfo_;
    GlobalVector origin_;
  };



  // SPLocalGeometry
  // ---------------

  template< int mydim, int cdim, class Grid >
  class SPLocalGeometry
  : public SPBasicGeometry< mydim, cdim, Grid, SPLocalGeometry< mydim, cdim, Grid > >
  {
    typedef SPLocalGeometry< mydim, cdim, Grid > This;
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
    SPLocalGeometry ( const GeometryCache &geometryCache,
                      const GlobalVector &origin )
    : geometryCache_( geometryCache ),
      origin_( origin )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    using Base::jacobian;
    using Base::jacobianInverse;

    GlobalVector origin () const { return origin_; }
    const GeometryCache &geometryCache () const { return geometryCache_; }

  private:
    GeometryCache geometryCache_;
    GlobalVector origin_;
  };



  // Implementation of SPBasicGeometry
  // ---------------------------------

  template< int mydim, int cdim, class Grid, class Impl >
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::GlobalVector
  SPBasicGeometry< mydim, cdim, Grid, Impl >::global ( const LocalVector &local ) const
  {
    GlobalVector y( asImpl().origin() );
    asImpl().geometryCache().jacobianTransposed().umtv( local, y );
    return y;
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::LocalVector
  SPBasicGeometry< mydim, cdim, Grid, Impl >::local ( const GlobalVector &global ) const
  {
    LocalVector x;
    GlobalVector y = global - asImpl().origin();
    asImpl().geometryCache().jacobianInverseTransposed().mtv( y, x );
    return x;
  }

  
  template< int mydim, int cdim, class Grid, class Impl >
  inline const typename SPBasicGeometry< mydim, cdim, Grid, Impl >::JacobianTransposed &
  SPBasicGeometry< mydim, cdim, Grid, Impl >::jacobianTransposed ( const LocalVector &local ) const
  {
    return asImpl().geometryCache().jacobianTransposed();
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline const typename SPBasicGeometry< mydim, cdim, Grid, Impl >::JacobianInverseTransposed &
  SPBasicGeometry< mydim, cdim, Grid, Impl >::jacobianInverseTransposed ( const LocalVector &local ) const
  {
    return asImpl().geometryCache().jacobianInverseTransposed();
  }

  template< int mydim, int cdim, class Grid, class Impl >
  inline auto SPBasicGeometry< mydim, cdim, Grid, Impl >::jacobian ( const LocalVector &local ) const
  {
    // Handing out a transposedView is OK, because jacobianTransposed
    // returns a reference to a cached value.
    return transposedView(jacobianTransposed(local));
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline auto SPBasicGeometry< mydim, cdim, Grid, Impl >::jacobianInverse ( const LocalVector &local ) const
  {
    // Handing out a transposedView is OK, because jacobianInverseTransposed
    // returns a reference to a cached value.
    return transposedView(jacobianInverseTransposed(local));
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRY_HH
