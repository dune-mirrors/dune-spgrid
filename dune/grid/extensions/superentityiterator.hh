#ifndef DUNE_GRID_SUPERENTITYITERATOR_HH
#define DUNE_GRID_SUPERENTITYITERATOR_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

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
