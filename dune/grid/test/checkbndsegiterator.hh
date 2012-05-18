#ifndef DUNE_SPGRID_CHECKBNDSEGITERATOR_HH
#define DUNE_SPGRID_CHECKBNDSEGITERATOR_HH

//- C++ includes
#include <cassert>

//- dune-common includes
#include <dune/common/exceptions.hh>

//- dune-grid includes
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridview.hh>


namespace Dune
{

  template< class VT >
  void checkBoundarySegmentIterator ( const Dune::GridView< VT > &gridView )
  {
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    typedef Dune::GridView< VT > GridView;
    typedef typename GridView::Implementation::BoundarySegmentIterator BoundarySegmentIterator;
    typedef typename BoundarySegmentIterator::Intersection Intersection;
    typedef typename Intersection::EntityPointer EntityPointer;
    typedef typename EntityPointer::Entity Entity;

    const BoundarySegmentIterator bend = gridView.impl().boundarySegmentEnd();
    for( BoundarySegmentIterator bit = gridView.impl().boundarySegmentBegin(); bit != bend; ++bit )
    {
      // get boundary segment
      const Intersection &boundarySegment = *bit;

      // assert boundary segment is on boundary
      if( !boundarySegment.boundary() )
        DUNE_THROW( GridError, "Boundary segment not on boundary" );

      // get index in inside
      const int indexInInside = boundarySegment.indexInInside();

      // get inside entity
      EntityPointer insidePtr = boundarySegment.inside();
      const Entity &inside = *insidePtr;

      // find boundary segment in inside intersections 
      bool foundIntersection = false;
      typedef typename GridView::IntersectionIterator IntersectionIterator;
      const IntersectionIterator iend = gridView.iend( inside );
      for( IntersectionIterator iit = gridView.ibegin( inside ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;

        // check for equality should be symmetrical
        if( intersection.impl().equals( boundarySegment.impl() ) 
              && boundarySegment.impl().equals( intersection.impl() ) )
        {
          // assert boundary segment is not found twice
          if( foundIntersection )
            DUNE_THROW( GridError, "Boundary segment found twice in inside entity" );

          // assert that index in inside conincides
          if( intersection.indexInInside() != indexInInside )
            DUNE_THROW( GridError, "Index in inside of boundary segment and intersection in inside entity do not conincide" );
          else
            foundIntersection = true;
        }
      }
      if( !foundIntersection )
        DUNE_THROW( GridError, "Boundary segment not found in inside entity" );
    }
#else // #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    std::cerr << ">>> Skipping check for BoundarySegmentIterator"
              << "because experimental grid extensions are disabled."
              << std::endl;
#endif // #else // #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
  }

}

#endif // #ifndef DUNE_SPGRID_CHECKBNDSEGITERATOR_HH
