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

    typedef SPIntersection< Grid > IntersectionImpl;

  public:
    typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

    typedef typename Intersection::Entity Entity;
    
    typedef typename IntersectionImpl::EntityInfo EntityInfo;
    typedef typename IntersectionImpl::GridLevel GridLevel;

    SPIntersectionIterator ( const EntityInfo &entityInfo, const int face )
    : intersection_( IntersectionImpl( entityInfo, face ) )
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

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
