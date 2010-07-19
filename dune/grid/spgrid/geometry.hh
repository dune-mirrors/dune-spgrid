#ifndef DUNE_SPGRID_GEOMETRY_HH
#define DUNE_SPGRID_GEOMETRY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/genericreferenceelements.hh>

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
    typedef typename remove_const< Grid >::type::Traits Traits;

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

    typedef typename GeometryCache::JacobianTransposed JacobianTransposed;
    typedef typename GeometryCache::JacobianInverseTransposed JacobianInverseTransposed;

  protected:
    SPBasicGeometry ();

  public:
    GeometryType type () const;

    int corners () const;
    const GlobalVector &operator[] ( const int i ) const;
    GlobalVector corner ( const int i ) const;
    GlobalVector center () const;

    bool affine () const;

    GlobalVector global ( const LocalVector &local ) const;
    LocalVector local ( const GlobalVector &global ) const;

    ctype volume () const;
    ctype integrationElement ( const LocalVector &local ) const;

    const JacobianTransposed &jacobianTransposed ( const LocalVector &local ) const;
    const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalVector &local ) const;

  protected:
    const Impl &asImpl () const;
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
    : entityInfo_( entityInfo )
    {}

    SPGeometry ( const GridLevel &gridLevel, const MultiIndex &id )
    : entityInfo_( gridLevel, id )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const ReferenceCube &referenceCube () const
    {
      return gridLevel().template referenceCube< codimension >();
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

    typedef typename Base::ReferenceCube ReferenceCube;
    typedef typename Base::GeometryCache GeometryCache;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;

  public:
    SPLocalGeometry ( const ReferenceCube &refCube,
                      const GeometryCache &geometryCache,
                      const GlobalVector &origin )
    : refCube_( &refCube ),
      geometryCache_( geometryCache ),
      origin_( origin )
    {}

    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const ReferenceCube &referenceCube () const
    {
      return *refCube_;
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
    const ReferenceCube *refCube_;
    GeometryCache geometryCache_;
    GlobalVector origin_;
  };



  // Implementation of SPBasicGeometry
  // ---------------------------------

  template< int mydim, int cdim, class Grid, class Impl >
  inline SPBasicGeometry< mydim, cdim, Grid, Impl >::SPBasicGeometry ()
  {}


  template< int mydim, int cdim, class Grid, class Impl >
  inline GeometryType SPBasicGeometry< mydim, cdim, Grid, Impl >::type () const
  {
    return GeometryType( GeometryType::cube, mydimension );
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline int SPBasicGeometry< mydim, cdim, Grid, Impl >::corners () const
  {
    return numCorners;
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline const typename SPBasicGeometry< mydim, cdim, Grid, Impl >::GlobalVector &
  SPBasicGeometry< mydim, cdim, Grid, Impl >::operator[] ( const int i ) const
  {
    DUNE_THROW( NotImplemented, "SPGrid does not implement Geometry::operator[], use Geometry::corner instead." );
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::GlobalVector
  SPBasicGeometry< mydim, cdim, Grid, Impl >::corner ( const int i ) const
  {
    const ReferenceCube &refCube = asImpl().referenceCube();
    return global( refCube.corner( i ) );
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::GlobalVector
  SPBasicGeometry< mydim, cdim, Grid, Impl >::center () const
  {
    const ReferenceCube &refCube = asImpl().referenceCube();
    return global( refCube.center() );
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline bool SPBasicGeometry< mydim, cdim, Grid, Impl >::affine () const
  {
    return true;
  }


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
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::ctype
  SPBasicGeometry< mydim, cdim, Grid, Impl >::volume () const
  {
    return asImpl().geometryCache().volume();
  }


  template< int mydim, int cdim, class Grid, class Impl >
  inline typename SPBasicGeometry< mydim, cdim, Grid, Impl >::ctype
  SPBasicGeometry< mydim, cdim, Grid, Impl >::integrationElement ( const LocalVector &local ) const
  {
    return volume();
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
  inline const Impl &SPBasicGeometry< mydim, cdim, Grid, Impl >::asImpl () const
  {
    return static_cast< const Impl & >( *this );
  }

}

#endif // #ifndef DUNE_SPGRID_GEOMETRY_HH
