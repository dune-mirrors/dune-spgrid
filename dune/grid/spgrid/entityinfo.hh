#ifndef DUNE_SPGRID_ENTITYINFO_HH
#define DUNE_SPGRID_ENTITYINFO_HH

#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  // SPBasicEntityInfo
  // -----------------

  template< class ct, int dim, int codim >
  class SPBasicEntityInfo
  {
    typedef SPBasicEntityInfo< ct, dim, codim > This;

  public:
    typedef SPGridLevel< ct, dim > GridLevel;

    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo ( const GridLevel &gridLevel )
    : gridLevel_( &gridLevel )
    {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
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

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
    unsigned int direction_;
  };



  // SPBasicEntityInfo (for codim =  0)
  // ----------------------------------

  template< class ct, int dim >
  class SPBasicEntityInfo< ct, dim, 0 >
  {
    typedef SPBasicEntityInfo< ct, dim, 0 > This;

  public:
    typedef SPGridLevel< ct, dim > GridLevel;

    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo ( const GridLevel &gridLevel )
    : gridLevel_( &gridLevel )
    {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
    : gridLevel_( &gridLevel ),
      id_( id )
    {}

    unsigned int direction () const
    {
      return (1 << dimension) - 1;
    }

    const MultiIndex &id () const
    {
      return id_;
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
  };



  // SPBasicEntityInfo (for codim =  dim)
  // ------------------------------------

  template< class ct, int dim >
  class SPBasicEntityInfo< ct, dim, dim >
  {
    typedef SPBasicEntityInfo< ct, dim, dim > This;

  public:
    typedef SPGridLevel< ct, dim > GridLevel;

    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo ( const GridLevel &gridLevel )
    : gridLevel_( &gridLevel )
    {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
    : gridLevel_( &gridLevel ),
      id_( id )
    {}

    unsigned int direction () const
    {
      return 0;
    }

    const MultiIndex &id () const
    {
      return id_;
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
  };



  // SPEntityInfo
  // ------------

  template< class ct, int dim, int codim >
  class SPEntityInfo
  : public SPBasicEntityInfo< ct, dim, codim >
  {
    typedef SPEntityInfo< ct, dim, codim > This;
    typedef SPBasicEntityInfo< ct, dim, codim > Base;

  public:
    typedef typename Base::GridLevel GridLevel;

    typedef typename GridLevel::ctype ctype;
    
    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    typedef typename GridLevel::template GeometryCache< codim > GeometryCache;
    typedef typename GeometryCache::LocalVector LocalVector;
    typedef typename GeometryCache::Jacobian Jacobian;
    typedef typename GeometryCache::JacobianTransposed JacobianTransposed;

    SPEntityInfo ( const GridLevel &gridLevel )
    : Base( gridLevel )
    {}

    SPEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
    : Base( gridLevel, id )
    {}

    using Base::direction;
    using Base::id;
    using Base::gridLevel;

    bool equals ( const This &other ) const
    {
      return (&gridLevel() == &other.gridLevel()) && (id() == other.id());
    }

    GlobalVector origin () const
    {
      const GlobalVector &h = gridLevel().h();
      GlobalVector origin = gridLevel().domain().origin();
      for( int i = 0; i < dimension; ++i )
        origin[ i ] += (id()[ i ] / 2) * h[ i ];
      return origin;
    }

    const ctype &volume () const
    {
      return geometryCache().volume( direction() );
    }

    const JacobianTransposed &jacobianTransposed () const
    {
      return geometryCache().jacobianTransposed( direction() );
    }

    const Jacobian &jacobianInverseTransposed () const
    {
      return geometryCache().jacobianInverseTransposed( direction() );
    }

    const GeometryCache &geometryCache () const
    {
      return gridLevel().template geometryCache< codim >( direction() );
    }
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITYINFO_HH
