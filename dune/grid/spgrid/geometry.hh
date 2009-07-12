#ifndef DUNE_SPGRID_GEOMETRY_HH
#define DUNE_SPGRID_GEOMETRY_HH

#include <dune/common/typetraits.hh>
#include <dune/common/geometrytype.hh>

#include <dune/grid/spgrid/entityinfo.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;



  // SPGeometry
  // ----------

  template< int mydim, int cdim, class Grid >
  class SPGeometry
  {
    typedef SPGeometry< mydim, cdim, Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    friend class SPEntity< Traits::dimension - mydim, cdim, Grid >;

  public:
    typedef typename Traits::ctype ctype;

    static const int mydimension = mydim;
    static const int coorddimension = cdim;
    static const int dimension = Traits::dimension;
    static const int codimension = dimension - mydimension;

    static const int numCorners = (1 << mydimension);

    typedef SPEntityInfo< ctype, dimension, codimension > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;

    typedef typename EntityInfo::GlobalVector GlobalVector;
    typedef typename EntityInfo::LocalVector LocalVector;
    typedef typename EntityInfo::Jacobian Jacobian;
    typedef typename EntityInfo::JacobianTransposed JacobianTransposed;

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
      // this can be optimized (refCorner[ i ] is 0 or 1)
      // return global( SPCubeCorners< ctype, mydimension >::corner( i ) );
    }

    GlobalVector global ( const LocalVector &local ) const
    {
      GlobalVector y = entityInfo_.origin();
      // this can be optimized
      jacobianTransposed().umtv( local, y );
      return y;
    }

    LocalVector local ( const GlobalVector &global ) const
    {
      GlobalVector y = global - entityInfo_.origin();
      LocalVector x;
      // this can be optimized
      jacobianInverseTransposed().mtv( y, x );
      return x;
    }

    ctype integrationElement ( const LocalVector &local ) const
    {
      return entityInfo_.volume();
    }

    ctype volume () const
    {
      return entityInfo_.volume();
    }

    const JacobianTransposed &
    jacobianTransposed ( const LocalVector &local ) const
    {
      return entityInfo_.jacobianTransposed();
    }

    const Jacobian &
    jacobianInverseTransposed ( const LocalVector &local ) const
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

}

#endif // #ifndef DUNE_SPGRID_GEOMETRY_HH
