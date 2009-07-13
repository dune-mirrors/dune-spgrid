#ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
#define DUNE_SPGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/spgrid/intersectioniterator.hh>

namespace Dune
{

  template< class Grid >
  class SPIntersectionIterator
  {
    typedef SPIntersectionIterator< Grid > This;

  public:
    typedef Dune::Intersection< Grid, SPIntersection > Intersection;

    const Intersection &dereference () const
    {
      return intersection_;
    }

    bool equals ( const This &other ) const
    {
      return intersection_.equals( other.intersection_ );
    }

    void increment ()
    {
      // ...
    }

  private:
    Intersection intersection_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
