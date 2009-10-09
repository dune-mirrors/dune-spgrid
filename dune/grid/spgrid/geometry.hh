#ifndef DUNE_SPGRID_GEOMETRY_HH
#define DUNE_SPGRID_GEOMETRY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/grid/spgrid/cube.hh>
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
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::Cube::ctype ctype;

    static const int mydimension = mydim;
    static const int coorddimension = cdim;
    static const int dimension = Traits::Cube::dimension;
    static const int codimension = dimension - mydimension;

    typedef SPCube< ctype, mydimension > Cube;
    typedef SPGeometryCache< ctype, dimension, codimension > GeometryCache;

    static const int numCorners = (1 << mydimension);

    typedef FieldVector< ctype, coorddimension > GlobalVector;
    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldMatrix< ctype, coorddimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, coorddimension > JacobianTransposed;

  protected:
    SPBasicGeometry ()
    {}

  public:
    GeometryType type () const
    {
      return GeometryType( GeometryType::cube, mydimension );
    }

    int corners () const
    {
      return numCorners;
    }

    const GlobalVector &operator[] ( const int i ) const
    {
      DUNE_THROW( NotImplemented, "SPGrid does not implement Geometry::operator[], use Geometry::corner instead." );
    }

    GlobalVector corner ( const int i ) const
    {
      const Cube &cube = asImpl().cube();
      return global( cube.corner( i ) );
    }

    bool affine () const
    {
      return true;
    }

    GlobalVector global ( const LocalVector &local ) const
    {
      GlobalVector y = asImpl().origin();
      // this can be optimized
      asImpl().geometryCache().jacobianTransposed().umtv( local, y );
      return y;
    }

    LocalVector local ( const GlobalVector &global ) const
    {
      GlobalVector y = global - asImpl().origin();
      LocalVector x;
      // this can be optimized
      asImpl().geometryCache().jacobianInverseTransposed().mtv( y, x );
      return x;
    }

    ctype volume () const
    {
      return asImpl().geometryCache().volume();
    }

    ctype integrationElement ( const LocalVector &local ) const
    {
      return volume();
    }

    const JacobianTransposed &
    jacobianTransposed ( const LocalVector &local ) const
    {
      return asImpl().geometryCache().jacobianTransposed();
    }

    const Jacobian &
    jacobianInverseTransposed ( const LocalVector &local ) const
    {
      return asImpl().geometryCache().jacobianInverseTransposed();
    }

  protected:
    const Impl &asImpl () const
    {
      return static_cast< const Impl & >( *this );
    }
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

    typedef SPEntityInfo< Grid, codimension > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;

    typedef typename Base::Cube Cube;
    typedef typename Base::GeometryCache GeometryCache;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;
    typedef typename Base::Jacobian Jacobian;
    typedef typename Base::JacobianTransposed JacobianTransposed;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    explicit SPGeometry ( const GridLevel &gridLevel )
    : entityInfo_( gridLevel )
    {}

    explicit SPGeometry ( const EntityInfo &entityInfo )
    : entityInfo_( entityInfo )
    {}

    SPGeometry ( const GridLevel &gridLevel, const MultiIndex &id )
    : entityInfo_( gridLevel, id )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const Cube &cube () const
    {
      return gridLevel().template cube< codimension >();
    }

    GlobalVector origin () const
    {
      return entityInfo().origin();
    }

    const GeometryCache &geometryCache () const
    {
      return entityInfo().geometryCache();
    }

    const GridLevel &gridLevel () const
    {
      return entityInfo().gridLevel();
    }

    const EntityInfo &entityInfo () const
    {
      return entityInfo_;
    }

    EntityInfo &entityInfo ()
    {
      return entityInfo_;
    }

  private:
    EntityInfo entityInfo_;
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

    typedef typename Base::Cube Cube;
    typedef typename Base::GeometryCache GeometryCache;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;
    typedef typename Base::Jacobian Jacobian;
    typedef typename Base::JacobianTransposed JacobianTransposed;

  public:
    SPLocalGeometry ( const Cube &cube, const GeometryCache &geometryCache,
                      const GlobalVector &origin )
    : cube_( &cube ),
      geometryCache_( geometryCache ),
      origin_( origin )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const Cube &cube () const
    {
      return *cube_;
    }

    GlobalVector origin () const
    {
      return origin_;
    }

    const GeometryCache &geometryCache () const
    {
      return geometryCache_;
    }

  private:
    const Cube *cube_;
    GeometryCache geometryCache_;
    GlobalVector origin_;
  };

}

#endif // #ifndef DUNE_SPGRID_GEOMETRY_HH
