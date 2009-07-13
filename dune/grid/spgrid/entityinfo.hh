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

    unsigned int direction () const
    {
      return direction_;
    }

    const MultiIndex &multiDirection () const
    {
      return gridLevel().multiDirection( codimension, direction_ );
    }

    const MultiIndex &multiIndex () const
    {
      return multiIndex_;
    }

    bool equals ( const This &other ) const
    {
      bool equals = (gridLevel_ == other.gridLevel_) & (direction_ == other.direction_);
      for( int i = 0; i < dimension; ++i )
        equals &= (multiIndex_[ i ] == other.multiIndex_[ i ]);
    }

    This father () const
    {
      This father( gridLevel().father() );
      father.direction_ = direction_;
      for( int i = 0; i < dimension; ++i )
      {
        assert( multiIndex_[ i ] % 1 == 0 );
        father.multiIndex_[ i ] = multiIndex_[ i ] / 2;
      }
      return father;
    }

    GlobalVector origin () const
    {
      const GlobalVector &h = gridLevel().h();
      GlobalVector origin = gridLevel().domain().origin();
      for( int i = 0; i < dimension; ++i )
        origin[ i ] += multiIndex_[ i ] * h[ i ];
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
    unsigned int direction_;
    MultiIndex multiIndex_;
  };

}

#endif
