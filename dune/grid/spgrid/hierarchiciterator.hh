#ifndef DUNE_SPGRID_HITERATOR_HH
#define DUNE_SPGRID_HITERATOR_HH

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPHierarchicIterator
  // --------------------

  template< class Grid, int codim >
  class SPHierarchicIterator
    : public SPEntityPointer< codim, Grid >
  {
    typedef SPHierarchicIterator< Grid, codim > This;
    typedef SPEntityPointer< codim, Grid > Base;

  public:
    typedef typename Base::Entity Entity;
    typedef typename Base::EntityInfo EntityInfo;

  protected:
    typedef typename Base::EntityImpl EntityImpl;

  public:
    SPHierarchicIterator () = default;

    SPHierarchicIterator ( const EntityInfo &entityInfo, int maxLevel )
    : Base( entityInfo ),
      minLevel_( entityInfo.gridLevel().level() ),
      maxLevel_( std::min( maxLevel, entityInfo.gridLevel().grid().maxLevel() ) )
    {
      increment();
    }

  public:
    using Base::entityInfo;
    using Base::level;

    void increment ()
    {
      if( level() >= maxLevel_ )
      {
        while( (level() > minLevel_) && !entityInfo().nextChild() )
          entityInfo().up();
      }
      else
        entityInfo().down();
      entityInfo().update();
    }

  private:
    int minLevel_, maxLevel_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_HITERATOR_HH
