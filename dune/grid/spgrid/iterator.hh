#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

#include <dune/grid/spgrid/misc.hh>
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
    typedef SPEntityPointer< codim, Grid > Base;

  public:
    typedef typename Base::EntityInfo EntityInfo;
    typedef typename Base::GridLevel GridLevel;

    static const int dimension = Base::dimension;
    static const int codimension = Base::codimension;
    static const int mydimension = Base::mydimension;

    static const unsigned int numDirections = GridLevel::numDirections;

    struct Begin {};
    struct End {};

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    using Base::gridLevel;

    SPIterator ( const GridLevel &gridLevel, const Begin &b,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      assert( sweepDir < numDirections );
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      unsigned int dir = 0;
      const unsigned int mydim = mydimension;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydim); ++dir );
      assert( dir < numDirections );

      MultiIndex &id = entityInfo.id();
      for( int i = 0; i < dimension; ++i )
        id[ i ] = begin( i, dir );
      entityInfo.update();
      assert( entityInfo.direction() == dir );
    }

    SPIterator ( const GridLevel &gridLevel, const End &e,
                 const unsigned int sweepDir = 0 )
    : Base( gridLevel ),
      sweepDirection_( sweepDir )
    {
      assert( sweepDir < numDirections );
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      unsigned int dir = numDirections-1;
      const unsigned int mydim = mydimension;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydim); --dir );
      assert( dir < numDirections );

      MultiIndex &id = entityInfo.id();
      for( int i = 0; i < dimension-1; ++i )
        id[ i ] = begin( i, dir );
      id[ dimension-1 ] = end( dimension-1, dir );
      entityInfo.update();
    }

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();

      MultiIndex &id = entityInfo.id();
      assert( id[ dimension-1 ] != end( dimension-1, entityInfo.direction() ) );
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int sweep = (sweepDirection_ >> i) & 1;
        id[ i ] += (2 - 4*sweep);
        if( id[ i ] != end( i, entityInfo.direction() ) )
          return entityInfo.update();
        id[ i ] = begin( i, entityInfo.direction() );
      }

      unsigned int dir = entityInfo.direction()+1;
      const unsigned int mydim = mydimension;
      for( ; (dir < numDirections) && (bitCount( dir ) != mydim); ++dir );
      if( dir < numDirections )
      {
        for( int i = 0; i < dimension; ++i )
          id[ i ] = begin( i, dir );
      }
      else
        id[ dimension-1 ] = end( dimension-1, entityInfo.direction() );

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
      return (d-2) + (1-sweep)*2*(cells[ i ]-(d-2));
    }

  protected:
    using Base::entity_;

  private:
    unsigned int sweepDirection_;
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH
