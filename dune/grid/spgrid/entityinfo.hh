#ifndef DUNE_SPGRID_ENTITYINFO_HH
#define DUNE_SPGRID_ENTITYINFO_HH

#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  template< class ct, int dim, int codim >
  class SPEntityInfo
  {
    typedef SPEntityInfo< ct, dim, codim > This;

  public:
    typedef SPGridLevel< ct, dim > GridLevel;

    typedef typename GridLevel::ctype ctype;
    
    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    typedef typename GridLevel::template GeometryCache< codim > GeometryCache;
    typedef typename GeometryCache::LocalVector LocalVector;
    typedef typename GeometryCache::Jacobian Jacobian;
    typedef typename GeometryCache::JacobianTransposed JacobianTransposed;

    SPEntityInfo ( const GridLevel &gridLevel )
    : gridLevel_( &gridLevel )
    {}

    SPEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
    : gridLevel_( &gridLevel ),
      id_( id ),
      direction_( id.direction() )
    {}

    unsigned int direction () const
    {
      return direction_;
    }

    const MultiIndex &id () const
    {
      return id_;
    }

    bool equals ( const This &other ) const
    {
      return (gridLevel_ == other.gridLevel_) && (id_ == other.id_);
    }

    GlobalVector origin () const
    {
      const GlobalVector &h = gridLevel().h();
      GlobalVector origin = gridLevel().domain().origin();
      for( int i = 0; i < dimension; ++i )
        origin[ i ] += (id_ / 2) * h[ i ];
      return origin;
    }

    const ctype &volume () const
    {
      return geometryCache().volume( direction_ );
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      return geometryCache().jacobianTransposed( direction_ );
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      return geometryCache().jacobianInverseTransposed( direction_ );
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

    const GeometryCache &geometryCache () const
    {
      return gridLevel().template geometryCache< codim >( direction_ );
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
    unsigned int direction_;
  };

}

#endif
