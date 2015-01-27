#ifndef DUNE_SPGRID_ENTITYINFO_HH
#define DUNE_SPGRID_ENTITYINFO_HH

#include <algorithm>
#include <limits>
#include <type_traits>

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
    typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo () : gridLevel_( nullptr ) {}

    SPBasicEntityInfo ( const GridLevel &gridLevel ) : gridLevel_( &gridLevel ) {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
      : gridLevel_( &gridLevel ),
        id_( id ),
        direction_( id.direction() )
    {}

    unsigned int direction () const { return direction_; }

    const MultiIndex &id () const { return id_; }
    MultiIndex &id () { return id_; }

    const GridLevel &gridLevel () const { assert( gridLevel_ ); return *gridLevel_; }

    void update () { direction_ = id_.direction(); }

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
    typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

    static const int dimension = GridLevel::dimension;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo () : gridLevel_( nullptr ) {}
    SPBasicEntityInfo ( const GridLevel &gridLevel ) : gridLevel_( &gridLevel ) {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
      : gridLevel_( &gridLevel ),
        id_( id )
    {}

    unsigned int direction () const { return (1 << dimension) - 1; }

    const MultiIndex &id () const { return id_; }
    MultiIndex &id () { return id_; }

    const GridLevel &gridLevel () const { assert( gridLevel_ ); return *gridLevel_; }

    void update () { assert( (id_.direction() == direction()) ); }

    void down ()
    {
      const Grid &grid = gridLevel().grid();
      const int level = gridLevel().level();
      gridLevel_ = &grid.gridLevel( level+1 );
      gridLevel().refinement().firstChild( id_ );
    }

    void up ()
    {
      const Grid &grid = gridLevel().grid();
      const int level = gridLevel().level();
      gridLevel().refinement().father( id_ );
      gridLevel_ = &grid.gridLevel( level-1 );
    }

    bool nextChild ()
    {
      return gridLevel().refinement().nextChild( id_ );
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
    typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

    typedef typename GridLevel::MultiIndex MultiIndex;
    typedef typename GridLevel::GlobalVector GlobalVector;

    SPBasicEntityInfo () : gridLevel_( nullptr ) {}
    SPBasicEntityInfo ( const GridLevel &gridLevel ) : gridLevel_( &gridLevel ) {}

    SPBasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
      : gridLevel_( &gridLevel ),
        id_( id )
    {}

    unsigned int direction () const { return 0; }

    const MultiIndex &id () const { return id_; }
    MultiIndex &id () { return id_; }

    const GridLevel &gridLevel () const { assert( gridLevel_ ); return *gridLevel_; }

    void update () { assert( (id_.direction() == direction()) ); }

  private:
    const GridLevel *gridLevel_;
    MultiIndex id_;
  };



  // SPEntityInfo
  // ------------

  template< class Grid, int codim >
  class SPEntityInfo
    : public SPBasicEntityInfo< Grid, std::remove_const< Grid >::type::dimension, codim >
  {
    typedef SPEntityInfo< Grid, codim > This;
    typedef SPBasicEntityInfo< Grid, std::remove_const< Grid >::type::dimension, codim > Base;

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

    SPEntityInfo () : Base(), partitionNumber_( std::numeric_limits< unsigned int >::max() ) {}

    SPEntityInfo ( const GridLevel &gridLevel )
      : Base( gridLevel ),
        partitionNumber_( std::numeric_limits< unsigned int >::max() )
    {}

    SPEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
      : Base( gridLevel, id ),
        partitionNumber_( partitionNumber )
    {}

    using Base::direction;
    using Base::id;
    using Base::gridLevel;

    bool equals ( const This &other ) const
    {
      return (&gridLevel() == &other.gridLevel()) && (id() == other.id());
    }

    unsigned int partitionNumber () const { return partitionNumber_; }

    PartitionType partitionType () const
    {
      return gridLevel().template partitionType< codim >( id(), partitionNumber() );
    }

    const GeometryCache &geometryCache () const
    {
      return gridLevel().template geometryCache< codim >( direction() );
    }

    void update ()
    {
      assert( std::find( id().begin(), id().end(), std::numeric_limits< int >::max() ) == id().end() );
      assert( gridLevel().template partition< All_Partition >().contains( id(), partitionNumber() ) );
      Base::update();
    }

    void update ( unsigned int partitionNumber )
    {
      partitionNumber_ = partitionNumber;
      update();
    }

  private:
    unsigned int partitionNumber_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYINFO_HH
