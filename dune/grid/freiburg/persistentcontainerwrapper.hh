#ifndef DUNE_PERSISTENTCONTAINERWRAPPER_HH
#define DUNE_PERSISTENTCONTAINERWRAPPER_HH

#include <dune/grid/utility/persistentcontainer.hh>

namespace Dune
{

  // PersistentContainerWrapper
  // --------------------------

  template< class G, class T, class Allocator >
  class PersistentContainerWrapper 
  {
    typedef PersistentContainerWrapper< G, T, Allocator > This;

    typedef typename G::HostGrid HostGrid;
    typedef PersistentContainer< HostGrid, T, Allocator > PersistentContainerHostGrid;

  public:
    typedef G Grid;
    typedef T Data;

    typedef typename Grid::template Codim< 0 >::Entity Element;

    typedef typename PersistentContainerHostGrid::Iterator Iterator;
    typedef typename PersistentContainerHostGrid::ConstIterator ConstIterator;

    PersistentContainerWrapper ( const Grid &grid, const int codim, const Allocator &allocator = Allocator() )
    : hostContainer_( grid.hostGrid(), codim, allocator )
    {}

    template< class Entity >
    Data &operator[] ( const Entity &entity )
    {
      return hostContainer_[ hostEntity( entity ) ];
    }

    template< class Entity >
    const Data &operator[] ( const Entity &entity ) const
    {
      return hostContainer_[ hostEntity( entity ) ];
    }

    Data &operator() ( const Element &element, const int subEntity )
    {
      return hostContainer_( hostEntity( element ), subEntity );
    }

    const Data &operator() ( const Element &element, const int subEntity ) const
    {
      return hostContainer_( hostEntity( element ), subEntity );
    }

    Iterator begin () { return hostContainer_.begin(); }
    ConstIterator begin () const { return hostContainer_.begin(); }

    Iterator end () { return hostContainer_.end(); }
    ConstIterator end () const { return hostContainer_.end(); }

    size_t size () const { return hostContainer_.size(); }

    void clear () { hostContainer_.clear(); }
    void reserve () { hostContainer_.reserve(); }
    void update () { hostContainer_.update(); }

  private:
    template< class Entity >
    const typename HostGrid::template Codim< Entity::codimension >::Entity &
    hostEntity ( const Entity &entity ) const
    {
      return Grid::getRealImplementation( entity ).hostEntity();
    }

    PersistentContainerHostGrid hostContainer_ ;
  };

} // namespace Dune

#endif // #ifndef DUNE_PERSISTENTCONTAINERWRAPPER_HH
