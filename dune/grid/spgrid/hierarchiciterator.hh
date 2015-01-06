#ifndef DUNE_SPGRID_HITERATOR_HH
#define DUNE_SPGRID_HITERATOR_HH

#include <type_traits>

#include <dune/grid/common/entityiterator.hh>

#include <dune/grid/spgrid/entity.hh>

namespace Dune
{

  // SPHierarchicIterator
  // --------------------

  template< class Grid, int codim >
  class SPHierarchicIterator
  {
    typedef SPHierarchicIterator< Grid, codim > This;

  public:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

    static const int dimension = Traits::ReferenceCube::dimension;
    static const int codimension = 0;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;

  private:
    typedef SPEntity< codimension, dimension, Grid > EntityImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    SPHierarchicIterator () = default;

    SPHierarchicIterator ( const EntityInfo &entityInfo, int maxLevel )
      : entityInfo_( entityInfo ),
        minLevel_( entityInfo.gridLevel().level() ),
        maxLevel_( std::min( maxLevel, entityInfo.gridLevel().grid().maxLevel() ) )
    {
      increment();
    }

    Entity dereference () const { return EntityImpl( entityInfo() ); }

    bool equals ( const This &other ) const { return entityInfo().equals( other.entityInfo() ); }

    void increment ()
    {
      if( gridLevel().level() >= maxLevel_ )
      {
        while( (gridLevel().level() > minLevel_) && !entityInfo().nextChild() )
          entityInfo().up();
      }
      else
        entityInfo().down();
      entityInfo().update();
    }

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

  private:
    EntityInfo entityInfo_;
    int minLevel_, maxLevel_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_HITERATOR_HH
