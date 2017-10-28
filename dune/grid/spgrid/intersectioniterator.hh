#ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
#define DUNE_SPGRID_INTERSECTIONITERATOR_HH

#include <type_traits>

#include <dune/grid/spgrid/intersection.hh>

#include <dune/grid/common/intersectioniterator.hh>

namespace Dune
{

  // SPIntersectionIterator
  // ----------------------

  template< class Grid >
  class SPIntersectionIterator
  {
    typedef SPIntersectionIterator< Grid > This;

    typedef typename std::remove_const< Grid >::type::Traits Traits;

    typedef SPIntersection< Grid > IntersectionImpl;

  public:
    typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

    typedef typename Intersection::Entity Entity;
    
    typedef typename IntersectionImpl::ElementInfo ElementInfo;
    typedef typename IntersectionImpl::GridLevel GridLevel;

    SPIntersectionIterator () = default;

    SPIntersectionIterator ( const ElementInfo &insideInfo, int face )
      : insideInfo_( insideInfo ), face_( face )
    {}

    Intersection dereference () const { return IntersectionImpl( insideInfo(), face_ ); }

    bool equals ( const This &other ) const
    {
      return (face_ == other.face_) && insideInfo().equals( other.insideInfo() );
    }

    void increment () { assert( face_ < GridLevel::ReferenceCube::numFaces ); ++face_; }

    const ElementInfo &insideInfo () const { return insideInfo_; }

  private:
    ElementInfo insideInfo_;
    int face_ = 0;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_INTERSECTIONITERATOR_HH
