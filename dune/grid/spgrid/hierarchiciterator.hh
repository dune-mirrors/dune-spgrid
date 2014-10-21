#ifndef DUNE_SPGRID_HITERATOR_HH
#define DUNE_SPGRID_HITERATOR_HH

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPHierarchicIterator
  // --------------------

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
    SPHierarchicIterator () = default;

    SPHierarchicIterator ( const EntityImpl &entityImpl, int maxLevel )
    : Base( entityImpl ),
      minLevel_( entityImpl.level() ),
      maxLevel_( std::min( maxLevel, entityImpl.grid().maxLevel() ) )
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
