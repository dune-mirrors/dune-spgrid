#ifndef DUNE_SPGRID_GEOMETRY_HH
#define DUNE_SPGRID_GEOMETRY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/grid/spgrid/cube.hh>
#include <dune/grid/spgrid/entityinfo.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;




  // SPBasicGeometry
  // ---------------

  template< int mydim, int cdim, class Grid, class Impl >
  class SPBasicGeometry
  {
    typedef SPBasicGeometry< mydim, cdim, Grid, Impl > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::ctype ctype;

    static const int mydimension = mydim;
    static const int coorddimension = cdim;
    static const int dimension = Traits::dimension;
    static const int codimension = dimension - mydimension;

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
      return global( asImpl().cube().corner( i ) );
    }

    GlobalVector global ( const LocalVector &local ) const
    {
      GlobalVector y = asImpl().origin();
      // this can be optimized
      asImpl().jacobianTransposed().umtv( local, y );
      return y;
    }

    LocalVector local ( const GlobalVector &global ) const
    {
      GlobalVector y = global - asImpl().origin();
      LocalVector x;
      // this can be optimized
      asImpl().jacobianInverseTransposed().mtv( y, x );
      return x;
    }

    ctype integrationElement ( const LocalVector &local ) const
    {
      return asImpl().volume();
    }

    const JacobianTransposed &
    jacobianTransposed ( const LocalVector &local ) const
    {
      return asImpl().jacobianTransposed();
    }

    const Jacobian &
    jacobianInverseTransposed ( const LocalVector &local ) const
    {
      return asImpl().jacobianInverseTransposed();
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

    typedef typename Base::Traits Traits;

    friend class SPEntity< Base::codimension, Base::dimension, Grid >;

  public:
    typedef typename Base::ctype ctype;

    static const int mydimension = Base::mydimension;
    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;

    typedef SPEntityInfo< Grid, codimension > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;

    typedef typename GridLevel::Cube Cube;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;
    typedef typename Base::Jacobian Jacobian;
    typedef typename Base::JacobianTransposed JacobianTransposed;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
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
      return gridLevel().cube();
    }

    GlobalVector origin () const
    {
      return entityInfo_.origin();
    }

    ctype volume () const
    {
      return entityInfo_.volume();
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      return entityInfo_.jacobianTransposed();
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      return entityInfo_.jacobianInverseTransposed();
    }

    const GridLevel &gridLevel () const
    {
      return entityInfo_.gridLevel();
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

    typedef typename Base::Traits Traits;

  public:
    typedef typename Base::ctype ctype;

    static const int mydimension = Base::mydimension;
    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;

    typedef SPCube< ctype, dimension > Cube;

    typedef typename Base::GlobalVector GlobalVector;
    typedef typename Base::LocalVector LocalVector;
    typedef typename Base::Jacobian Jacobian;
    typedef typename Base::JacobianTransposed JacobianTransposed;

  public:
    using Base::jacobianTransposed;
    using Base::jacobianInverseTransposed;

    const Cube &cube () const
    {
      // ...
    }

    GlobalVector origin () const
    {
      // ...
    }

    ctype volume () const
    {
      // ...
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      // ...
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      // ...
    }
  };

}

#endif // #ifndef DUNE_SPGRID_GEOMETRY_HH
