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

    typedef typename Intersection::Entity Entity;
    
  private:
    typedef SPEntity< 0, Traits::ReferenceCube::dimension, Grid > EntityImpl;
    typedef SPIntersection< Grid > IntersectionImpl;

  public:
    typedef typename IntersectionImpl::GridLevel GridLevel;

  public:
    SPIntersectionIterator ( const Entity &entity, const int face )
    : intersection_( IntersectionImpl( Grid::getRealImplementation( entity ), face ) )
    {}

    SPIntersectionIterator ( const EntityImpl &entityImpl, const int face )
    : intersection_( IntersectionImpl( entityImpl, face ) )
    {}

    SPIntersectionIterator ( const This &other )
    : intersection_( Grid::getRealImplementation( other.intersection_ ) )
    {}

    const This &operator= ( const This &other )
    {
      Grid::getRealImplementation( intersection_ ) = Grid::getRealImplementation( other.intersection_ );
      return *this;
    }

    const Intersection &dereference () const
    {
      return intersection_;
    }

    bool equals ( const This &other ) const
    {
      return Grid::getRealImplementation( intersection_ ).equals( Grid::getRealImplementation( other.intersection_ ) );
    }

    void increment ()
    {
      const int face = intersection_.indexInInside();
      assert( face < GridLevel::ReferenceCube::numFaces );
      Grid::getRealImplementation( intersection_ ).setFace( face+1 );
    }

  private:
    Intersection intersection_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
