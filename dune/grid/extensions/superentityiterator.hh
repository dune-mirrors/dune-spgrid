#ifndef DUNE_GRID_SUPERENTITYITERATOR_HH
#define DUNE_GRID_SUPERENTITYITERATOR_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  template< class Object >
  class SuperEntityIteratorExtension;



  // SuperEntityIterator
  // -------------------

  template< class Grid, template< class > class SuperEntityIteratorImp >
  class SuperEntityIterator
  : public EntityPointer< Grid, SuperEntityIteratorImp< Grid > >
  {
    typedef SuperEntityIterator< Grid, SuperEntityIteratorImp > This;
    typedef EntityPointer< Grid, SuperEntityIteratorImp< Grid > > Base;

    typedef SuperEntityIteratorImp< Grid > Implementation;

  public:
    typedef typename Grid::template Codim< 0 >::Entity Entity;

    SuperEntityIterator ( const Implementation &implementation );

    const This &operator++ ();

    int index () const;

  protected:
    using Base::realIterator;
  };



  // SuperEntityIteratorExtension for GridView
  // -----------------------------------------

  template< class ViewTraits >
  class SuperEntityIteratorExtension< Dune::GridView< ViewTraits > >
  {
    typedef SuperEntityIteratorExtension< Dune::GridView< ViewTraits > > This;

    typedef Dune::GridView< ViewTraits > GridView;

    typedef typename GridView::Grid Grid;

  private:
    typedef typename ViewTraits::GridViewImp GridViewImp;

  public:
    template< int codim >
    struct Codim
    {
      static const bool hasSuperEntityIterator
        = ViewTraits::template Codim< codim >::hasSuperEntityIterator;

      typedef typename ViewTraits::template Codim< codim >::SuperEntityIterator
        SuperEntityIterator;
    };

    SuperEntityIteratorExtension ( const GridView &gridView )
    : gridView_( Grid::getRealImplementation( gridView ) )
    {}

    template< class Entity >
    typename Codim< Entity::codimension >::SuperEntityIterator
    superEntityBegin ( const Entity &entity ) const
    {
      return gridView_.superEntityBegin< Entity::codimension >( entity );
    }

    template< class Entity >
    typename Codim< Entity::codimension >::SuperEntityIterator
    superEntityEnd ( const Entity &entity ) const
    {
      return gridView_.superEntityEnd< Entity::codimension >( entity );
    }

  private:
    const GridViewImp &gridView_;
  };


  // Implementation of SuperEntityIterator
  // -------------------------------------

  template< class Grid, template< class > class SuperEntityIteratorImp >
  inline SuperEntityIterator< Grid, SuperEntityIteratorImp >
    ::SuperEntityIterator ( const Implementation &implementation )
  : Base( implementation )
  {}


  template< class Grid, template< class > class SuperEntityIteratorImp >
  inline const typename SuperEntityIterator< Grid, SuperEntityIteratorImp >::This &
  SuperEntityIterator< Grid, SuperEntityIteratorImp >::operator++ ()
  {
    realIterator.increment();
    return *this;
  }


  template< class Grid, template< class > class SuperEntityIteratorImp >
  inline int SuperEntityIterator< Grid, SuperEntityIteratorImp >::index () const
  {
    return realIterator.index();
  }



  // Extensions
  // ----------

  namespace Extensions
  {

    template< class Grid, int codim >
    struct SuperEntityIterator
    {
      static const bool v = false;
    };

    template< class Grid, int codim >
    struct SuperEntityIterator< const Grid, codim >
    {
      static const bool v = SuperEntityIterator< Grid, codim >::v;
    };

  }

}

#endif // #ifndef DUNE_GRID_EXTENSIONS_SUPERENTITYITERATOR_HH
