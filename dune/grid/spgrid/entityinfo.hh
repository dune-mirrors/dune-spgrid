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
      direction_( 0 )
    {
      for( int j = 0; j < dimension; ++j )
        direction_ |= ((id[ j ] & 1) << j);
    }

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
      bool equals = (gridLevel_ == other.gridLevel_);
      for( int i = 0; i < dimension; ++i )
        equals &= (id_[ i ] == other.id_[ i ]);
    }

    This father () const
    {
      This father( gridLevel().father() );
      father.direction_ = direction_;
      for( int i = 0; i < dimension; ++i )
        father.id_[ i ] = (id_[ i ] / 2) | 1;
      return father;
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
      return gridLevel().template geometryCache< codim >();
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
    unsigned int direction_;
  };

}

#endif
