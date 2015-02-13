#ifndef DUNE_SPGRID_ENTITYINFO_HH
#define DUNE_SPGRID_ENTITYINFO_HH

#include <algorithm>
#include <limits>
#include <type_traits>

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/gridlevel.hh>

namespace Dune
{

  namespace __SPGrid
  {

    // BasicEntityInfo
    // ---------------

    template< class Grid, int dim, int codim >
    class BasicEntityInfo
    {
      typedef BasicEntityInfo< Grid, dim, codim > This;

      typedef SPEntityDirection< dim, dim - codim > EntityDirection;

    public:
      typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

      static const int dimension = GridLevel::dimension;

      typedef typename GridLevel::MultiIndex MultiIndex;
      typedef typename GridLevel::GlobalVector GlobalVector;

      typedef typename EntityDirection::Direction Direction;

      BasicEntityInfo () : gridLevel_( nullptr ) {}

      BasicEntityInfo ( const GridLevel &gridLevel ) : gridLevel_( &gridLevel ) {}

      BasicEntityInfo ( const GridLevel &gridLevel, const MultiIndex &id )
        : gridLevel_( &gridLevel ), id_( id ), direction_( id )
      {}

      Direction direction () const { return direction_; }

      const MultiIndex &id () const { return id_; }
      MultiIndex &id () { return id_; }

      const GridLevel &gridLevel () const { assert( gridLevel_ ); return *gridLevel_; }

      void update () { direction_ = EntityDirection( id() ); }

      // hierarchic traversal

      bool hasFather () const
      {
        return ((codim == 0) || gridLevel().refinement().hasFather( id() ));
      }

      void up ()
      {
        const Grid &grid = gridLevel().grid();
        const int level = gridLevel().level();
        gridLevel().refinement().father( id() );
        gridLevel_ = &grid.gridLevel( level-1 );
      }

      void down ()
      {
        const Grid &grid = gridLevel().grid();
        const int level = gridLevel().level();
        gridLevel_ = &grid.gridLevel( level+1 );
        gridLevel().refinement().firstChild( id() );
      }

      bool nextChild ()
      {
        return gridLevel().refinement().nextChild( id() );
      }

    private:
      const GridLevel *gridLevel_;
      MultiIndex id_;
      EntityDirection direction_;
    };



    // EntityInfo
    // ----------

    template< class Grid, int codim >
    class EntityInfo
      : public BasicEntityInfo< Grid, std::remove_const< Grid >::type::dimension, codim >
    {
      typedef EntityInfo< Grid, codim > This;
      typedef BasicEntityInfo< Grid, std::remove_const< Grid >::type::dimension, codim > Base;

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

      EntityInfo () : Base(), partitionNumber_( std::numeric_limits< unsigned int >::max() ) {}

      EntityInfo ( const GridLevel &gridLevel )
        : Base( gridLevel ),
          partitionNumber_( std::numeric_limits< unsigned int >::max() )
      {}

      EntityInfo ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
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

  } // namespace __SPGrid

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYINFO_HH
