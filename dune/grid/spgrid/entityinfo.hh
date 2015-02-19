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

    // EntityInfo
    // ----------

    template< class Grid, int codim >
    class EntityInfo
    {
      typedef EntityInfo< Grid, codim > This;

    public:
      typedef SPGridLevel< typename std::remove_const< Grid >::type > GridLevel;

      static const int dimension = GridLevel::dimension;
      static const int codimension = codim;
      static const int mydimension = dimension - codimension;

    private:
      typedef SPEntityDirection< dimension, mydimension > EntityDirection;

    public:
      typedef typename EntityDirection::Direction Direction;

      typedef typename GridLevel::MultiIndex MultiIndex;
      typedef typename GridLevel::GlobalVector GlobalVector;

      typedef typename GridLevel::Traits Traits;

      typedef typename GridLevel::template Codim< codimension >::GeometryCache GeometryCache;

      EntityInfo ()
        : gridLevel_( nullptr ),
          partitionNumber_( std::numeric_limits< unsigned int >::max() )
      {}

      EntityInfo ( const GridLevel &gridLevel )
        : gridLevel_( &gridLevel ),
          partitionNumber_( std::numeric_limits< unsigned int >::max() )
      {}

      EntityInfo ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
        : gridLevel_( &gridLevel ),
          id_( id ),
          direction_( id ),
          partitionNumber_( partitionNumber )
      {}

      // access member data

      const GridLevel &gridLevel () const { assert( gridLevel_ ); return *gridLevel_; }

      const MultiIndex &id () const { return id_; }
      MultiIndex &id () { return id_; }

      Direction direction () const { return direction_; }

      unsigned int partitionNumber () const { return partitionNumber_; }

      // convenience to implement entity

      bool equals ( const This &other ) const
      {
        return (gridLevel_ == other.gridLevel_) && (id() == other.id());
      }

      PartitionType partitionType () const
      {
        return gridLevel().template partitionType< codimension >( id(), partitionNumber() );
      }

      const GeometryCache &geometryCache () const
      {
        return gridLevel().template geometryCache< codimension >( direction() );
      }

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

      // manipulation methods

      void update ()
      {
        assert( std::find( id().begin(), id().end(), std::numeric_limits< int >::max() ) == id().end() );
        assert( gridLevel().template partition< All_Partition >().contains( id(), partitionNumber() ) );
        direction_ = EntityDirection( id() );
      }

      void update ( unsigned int partitionNumber )
      {
        partitionNumber_ = partitionNumber;
        update();
      }

    private:
      const GridLevel *gridLevel_;
      MultiIndex id_;
      EntityDirection direction_;
      unsigned int partitionNumber_;
    };

  } // namespace __SPGrid

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYINFO_HH
