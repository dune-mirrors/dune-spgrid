#ifndef DUNE_SPGRID_ITERATOR_HH
#define DUNE_SPGRID_ITERATOR_HH

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPIterator
  // ----------

  template< int codim, PartitionIteratorType pitype, class Grid >
  class SPIterator
  : public SPEntityPointer< codim, Grid >
  {
    typedef SPIterator< codim, pitype, Grid > This;
    typedef SPEntityPointer< Grid > Base;

  public:
    struct Begin {};
    struct End {};

    SPIterator ( const GridLevel &gridLevel, const Begin &begin )
    : Base( gridLevel )
    {
      // ...
    }

    SPIterator ( const GridLevel &gridLevel, const End &end )
    : Base( gridLevel )
    {
      // ...
    }

    void increment ()
    {
      // ...
    }
  };

}

#endif // #ifndef DUNE_SPGRID_ITERATOR_HH
