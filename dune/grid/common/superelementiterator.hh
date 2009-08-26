#ifndef DUNE_GRID_SUPERELEMENTITERATOR_HH
#define DUNE_GRID_SUPERELEMENTITERATOR_HH

#include <dune/grid/common/entitypointer.hh>

namespace Dune
{

  template< class Grid, template< class > class SuperElementIteratorImp >
  class SuperElementIterator
  : public EntityPointer< Grid, SuperElementIteratorImp< Grid > >
  {
    typedef SuperElementIterator< Grid, SuperElementIteratorImp > This;
    typedef EntityPointer< Grid, SuperElementIteratorImp< Grid > > Base;

    typedef SuperElementIteratorImp< Grid > Implementation;

  public:
    typedef typename Grid::template Codim< 0 >::Entity Entity;

    This &operator++ ()
    {
      realIterator.increment();
      return *this;
    }
    
    SuperElementIterator ( const Implementation &implementation )
    : Base( implementation )
    {}

  protected:
    using Base::realIterator;
  };

}

#endif // #ifndef DUNE_GRID_SUPERELEMENTITERATOR_HH
