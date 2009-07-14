#ifndef DUNE_SPGRID_HITERATOR_HH
#define DUNE_SPGRID_HITERATOR_HH

#include <dune/grid/common/hierarchiciterator.hh>

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  template< class Grid >
  class SPHierarchicIterator
  : public SPEntityPointer< Grid, SPHierarchicIterator< Grid > >
  {
    typedef SPHierarchicIterator< Grid > This;
    typedef SPEntityPointer< Grid > Base;

  public:
    SPHierarchicIterator ( const Entity &entity, int maxLevel )
    : Base( entity ),
      maxLevel_( maxLevel )
    {
      increment();
    }

    using Base::level;
    using Base::dereference();

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      MultiIndex = entityInfo.id();
      if( level() < maxLevel_ )
      {
        if( entity_.isLeaf() )
          maxLevel_ = level();

      }

    }

  protected:
    using Base::entity_;

  private:
    int maxLevel_;
  };

}

#endif // #ifndef DUNE_SPGRID_HITERATOR_HH
