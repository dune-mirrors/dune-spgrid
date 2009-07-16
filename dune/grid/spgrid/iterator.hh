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

    SPIterator ( const GridLevel &gridLevel, const Begin &begin )
    : Base( gridLevel )
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      unsigned int dir = 0;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );
      assert( dir < numDirections );

      MulitIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = (dir >> i) & 1;
      entityInfo.setId( id );
    }

    SPIterator ( const GridLevel &gridLevel, const End &end )
    : Base( gridLevel )
    {
      // ...
    }

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      // ...

      unsigned int dir = entityInfo.direction();
      for( ; (dir < numDirections) && (bitCount( dir ) != mydimension); ++dir );

      // ...
    }
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH
