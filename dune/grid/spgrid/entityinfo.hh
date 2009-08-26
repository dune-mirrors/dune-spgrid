#ifndef DUNE_SPGRID_ENTITYINFO_HH
#define DUNE_SPGRID_ENTITYINFO_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  // SPBasicEntityInfo
  // -----------------

  template< class Grid, int dim, int codim >
  class SPBasicEntityInfo
  {
    typedef SPBasicEntityInfo< Grid, dim, codim > This;

  public:
    typedef SPGridLevel< Grid > GridLevel;

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

    MultiIndex &id ()
    {
      return id_;
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

    void update ()
    {
      direction_ = id_.direction();
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
    unsigned int direction_;
  };



  // SPBasicEntityInfo (for codim =  0)
  // ----------------------------------

  template< class Grid, int dim >
  class SPBasicEntityInfo< Grid, dim, 0 >
  {
    typedef SPBasicEntityInfo< Grid, dim, 0 > This;

  public:
    typedef SPGridLevel< Grid > GridLevel;

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

    MultiIndex &id ()
    {
      return id_;
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

    void update ()
    {
      assert( id_.direction() == direction() );
    }

    void down ()
    {
      assert( !gridLevel().isLeaf() );
      gridLevel_ = &gridLevel().childLevel();
      const unsigned int refDir = gridLevel().refinementDirection();
      for( int i = 0; i < dimension; ++i )
        id_[ i ] = (((refDir >> i) & 1) != 0 ? (id_[ i ] << 1) - 1 : id_[ i ]);
    }

    void up ()
    {
      assert( gridLevel().level() > 0 );
      const unsigned int refDir = gridLevel().refinementDirection();
      for( int i = 0; i < dimension; ++i )
        id_[ i ] = (((refDir >> i) & 1) != 0 ? (id_[ i ] >> 1) | 1 : id_[ i ]);
      gridLevel_ = &gridLevel().fatherLevel();
    }

    bool nextChild ()
    {
      const unsigned int refDir = gridLevel_->refinementDirection();
      for( int i = 0; i < dimension; ++i )
      {
        if( ((refDir >> i) & 1) == 0 )
          continue;
        id_[ i ] ^= 2;
        if( (id_[ i ] & 2) != 0 )
          return true;
      }
      return false;
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
  };



  // SPBasicEntityInfo (for codim = dim)
  // -----------------------------------

  template< class Grid, int dim >
  class SPBasicEntityInfo< Grid, dim, dim >
  {
    typedef SPBasicEntityInfo< Grid, dim, dim > This;

  public:
    typedef SPGridLevel< Grid > GridLevel;

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

    MultiIndex &id ()
    {
      return id_;
    }

    const GridLevel &gridLevel () const
    {
      return *gridLevel_;
    }

    void update ()
    {
      assert( id_.direction() == direction() );
    }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
  };



  // SPEntityInfo
  // ------------

  template< class Grid, int codim >
  class SPEntityInfo
  : public SPBasicEntityInfo< Grid, SPGridLevel< Grid >::dimension, codim >
  {
    typedef SPEntityInfo< Grid, codim > This;
    typedef SPBasicEntityInfo< Grid, SPGridLevel< Grid >::dimension, codim > Base;

  public:
    typedef typename Base::GridLevel GridLevel;

    typedef typename GridLevel::Traits Traits;
    typedef typename GridLevel::ctype ctype;

    static const int dimension = GridLevel::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;
    
    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    typedef typename GridLevel::template Codim< codim >::GeometryCache GeometryCache;
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

    PartitionType partitionType () const
    {
      return InteriorEntity;
    }

    GlobalVector origin () const
    {
      const GlobalVector &h = gridLevel().h();
      GlobalVector origin = gridLevel().domain().origin();
      for( int i = 0; i < dimension; ++i )
        origin[ i ] += (id()[ i ] / 2) * h[ i ];
      return origin;
    }

    const GeometryCache &geometryCache () const
    {
      return gridLevel().template geometryCache< codim >( direction() );
    }
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITYINFO_HH
