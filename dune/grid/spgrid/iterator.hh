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
    struct Begin {};
    struct End {};

    SPIterator ( const GridLevel &gridLevel, const Begin &begin,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      unsigned int dir = 0;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );
      assert( dir < numDirections );

      MulitIndex id;
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int sweep = (sweepDirection_ >> i) & 1;
        const unsigned int d = (dir >> i) & 1;
        id[ i ] = sweep*(2*(cells[ i ]-d)+1) + d;
      }
      entityInfo.setId( id );
    }

    SPIterator ( const GridLevel &gridLevel, const End &end,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      // ...
    }

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      MultiIndex &id = entityInfo.id();
      unsigned int dir = entityInfo.direction();

      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int sweep = (sweepDirection_ >> i) & 1;
        id[ i ] += (2 - 4*sweep);
        if( (id[ i ] >= 0) && (id[ i ] <= 2*cells[ i ]) )
          return;

        const unsigned int d = (dir >> i) & 1;
        id[ i ] = sweep*(2*(cells[ i ]-d)+1) + d;
      }

      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );

      MulitIndex id;
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int sweep = (sweepDirection_ >> i) & 1;
        const unsigned int d = (dir >> i) & 1;
        id[ i ] = sweep*(2*(cells[ i ]-d)+1) + d;
      }
      entityInfo.setId( id );
    }

  private:
    const unsigned int sweepDirection_;
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH
