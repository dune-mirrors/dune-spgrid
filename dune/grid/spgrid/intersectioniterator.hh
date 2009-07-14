#ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
#define DUNE_SPGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/spgrid/intersection.hh>

#include <dune/grid/common/intersectioniterator.hh>

namespace Dune
{

  template< class Grid >
  class SPIntersectionIterator
  {
    typedef SPIntersectionIterator< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef Dune::Intersection< Grid, SPIntersection > Intersection;
    
    typedef typename Traits::Entity Entity;

  private:
    typedef SPIntersection< Grid > IntersectionImpl;

  public:
    SPIntersectionIterator ( const Entity &entity, const unsigned int face )
    : intersection_( IntersectionImpl( entity, face ) )
    {}

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
      const int face = intersection_.indexInInside();
      assert( face < numFaces );
      Grid::getRealImplementation( intersection_ ).setFace( face+1 );
    }

  private:
    Intersection intersection_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
