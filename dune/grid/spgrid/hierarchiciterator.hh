#ifndef DUNE_SPGRID_HITERATOR_HH
#define DUNE_SPGRID_HITERATOR_HH

#include <dune/grid/common/hierarchiciterator.hh>

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  template< class Grid >
  class SPHierarchicIterator
  : public SPEntityPointer< 0, Grid >
  {
    typedef SPHierarchicIterator< Grid > This;
    typedef SPEntityPointer< 0, Grid > Base;

    friend class SPEntity< 0, Base::dimension, Grid >;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::EntityInfo EntityInfo;

  protected:
    typedef typename Base::EntityImpl EntityImpl;

  private:
    SPHierarchicIterator ( const EntityImpl &entityImpl, int maxLevel )
    : Base( entityImpl ),
      minLevel_( entityImpl.level() ),
      maxLevel_( maxLevel )
    {
      increment();
    }

  public:
    using Base::level;
    using Base::dereference;

    void increment ()
    {
      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      if( entity_.isLeaf() || (level() >= maxLevel_) )
      {
        while( (level() > minLevel_) && !entityInfo.nextChild() )
          entityInfo.up();
      }
      else
        entityInfo.down();
    }

  protected:
    using Base::entity_;

  private:
    int minLevel_, maxLevel_;
  };

}

#endif // #ifndef DUNE_SPGRID_HITERATOR_HH
