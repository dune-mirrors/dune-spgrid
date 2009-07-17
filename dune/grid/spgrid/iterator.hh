#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPIterator
  // ----------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class SPIterator
  : public SPEntityPointer< codim, Grid >
  {
    typedef SPIterator< codim, pitype, Grid > This;
    typedef SPEntityPointer< Grid > Base;

  public:
    typedef typename Base::GridLevel GridLevel;

    static const unsigned int numDirections = GridLevel::numDirections;

    struct Begin {};
    struct End {};

    using Base::gridLevel;

    SPIterator ( const GridLevel &gridLevel, const Begin &begin,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      unsigned int dir = 0;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );
      assert( dir < numDirections );

      MulitIndex &id = entityInfo.id();
      for( int i = 0; i < dimension-1; ++i )
        id[ i ] = begin( i, dir );
      entityInfo.update();
    }

    SPIterator ( const GridLevel &gridLevel, const End &end,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      unsigned int dir = numDirections-1;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); --dir );
      assert( dir < numDirections );

      MulitIndex &id = entityInfo.id();
      for( int i = 0; i < dimension-1; ++i )
        id[ i ] = begin( i, dir );
      id[ dimension-1 ] = end( dimension-1, dir );
      entityInfo.update();
    }

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      MultiIndex &id = entityInfo.id();
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int sweep = (sweepDirection_ >> i) & 1;
        id[ i ] += (2 - 4*sweep);
        if( id[ i ] != end( i, entityInfo.direction() ) )
          return entityInfo.update();
        id[ i ] = begin( i, entityInfo.direction() );
      }

      unsigned int dir = entityInfo.direction();
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );
      if( dir < numDirections )
      {
        for( int i = 0; i < dimension; ++i )
          id[ i ] = begin( i, dir );
      }

      entityInfo.update();
    }

  private:
    const int begin ( const int i, const unsigned int dir ) const
    {
      const MultiIndex &cells = gridLevel().cells();
      const unsigned int sweep = (sweepDirection_ >> i) & 1;
      const unsigned int d = (dir >> i) & 1;
      return d + sweep*2*(cells[ i ]-d);
    }

    const int end ( const int i, const unsigned int dir ) const
    {
      const MultiIndex &cells = gridLevel().cells();
      const unsigned int sweep = (sweepDirection_ >> i) & 1;
      const unsigned int d = (dir >> i) & 1;
      return (d-2) + (1-sweep)*2*(cells[ i ]-d);
    }

    const unsigned int sweepDirection_;
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH
