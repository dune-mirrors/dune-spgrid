#ifndef DUNE_GRID_SUPERENTITYITERATOR_HH
#define DUNE_GRID_SUPERENTITYITERATOR_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  template< class Grid, template< class > class SuperEntityIteratorImp >
  class SuperEntityIterator
  : public EntityPointer< Grid, SuperEntityIteratorImp< Grid > >
  {
    typedef SuperEntityIterator< Grid, SuperEntityIteratorImp > This;
    typedef EntityPointer< Grid, SuperEntityIteratorImp< Grid > > Base;

    typedef SuperEntityIteratorImp< Grid > Implementation;

  public:
    typedef typename Grid::template Codim< 0 >::Entity Entity;

    This &operator++ ()
    {
      realIterator.increment();
      return *this;
    }
    
    SuperEntityIterator ( const Implementation &implementation )
    : Base( implementation )
    {}

  protected:
    using Base::realIterator;
  };

}

#endif // #ifndef DUNE_GRID_SUPERENTITYITERATOR_HH
