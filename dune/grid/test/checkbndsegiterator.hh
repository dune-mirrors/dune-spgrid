#ifndef DUNE_SPGRID_CHECKBNDSEGITERATOR_HH
#define DUNE_SPGRID_CHECKBNDSEGITERATOR_HH

//- C++ includes
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

//- dune-common includes
#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>

//- dune-grid includes
#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridview.hh>


namespace Dune
{

  template< class VT >
  void checkBoundarySegmentIterator ( const Dune::GridView< VT > &gridView )
  {
    // get grid view types
    typedef Dune::GridView< VT > GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename IntersectionIterator::Intersection Intersection;
    typedef typename Intersection::EntityPointer EntityPointer;
    typedef typename EntityPointer::Entity Entity;

    // boundary segment iterator to be tested
    typedef typename GridView::Implementation::BoundarySegmentIterator BoundarySegmentIterator;

    // get index set
    typedef typename GridView::IndexSet IndexSet;
    const IndexSet &indexSet = gridView.indexSet();

    // storage for all intersections in grid view
    static const int numFaces = 2*GridView::dimension;
    typedef Dune::array< int, numFaces > Data;
    std::vector< Data > boundaryFace( indexSet.size( 0 ) );
    for( typename std::vector< Data >::iterator vit = boundaryFace.begin(); vit != boundaryFace.end(); ++vit )
      vit->fill( 0 );

    // visit all boundary intersections
    int bndIntersections = 0;
    typedef typename GridView::template Codim< 0 >::Iterator Iterator;
    const Iterator end = gridView.template end< 0 >();
    for( Iterator it = gridView.template begin< 0 >(); it != end; ++it )
    {
      // get entity
      const Entity &entity = *it;
      if( !entity.hasBoundaryIntersections() )
        continue;

      // visit boundary intersections
      const typename IndexSet::IndexType index = indexSet.index( entity );
      const IntersectionIterator iend = gridView.iend( entity );
      for( IntersectionIterator iit = gridView.ibegin( entity ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;
        if( intersection.boundary() )
        {
          boundaryFace[ index ][ intersection.indexInInside() ] = 1;
          bndIntersections++;
        }
      }
    }
    std::cerr << "Found " << bndIntersections << " boundary intersections" << std::endl;

    // now iterate over all boundary segments
    int bndSegments = 0;
    const BoundarySegmentIterator bend = gridView.impl().boundarySegmentEnd();
    for( BoundarySegmentIterator bit = gridView.impl().boundarySegmentBegin(); bit != bend; ++bit )
    {
      // dereference boundary segment iterator
      const Intersection &boundarySegment = *bit;
      ++bndSegments;

      // assert boundary segment is on boundary
      if( !boundarySegment.boundary() )
        DUNE_THROW( GridError, "Boundary segment not on boundary" );

      // get index in inside
      const int indexInInside = boundarySegment.indexInInside();

      // get inside entity
      EntityPointer insidePtr = boundarySegment.inside();
      const Entity &inside = *insidePtr;

      // find boundary segment in inside intersections 
      const IntersectionIterator iend = gridView.iend( inside );
      for( IntersectionIterator iit = gridView.ibegin( inside ); iit != iend; ++iit )
      {
        const Intersection &intersection = *iit;

        // equality check should be symmetrical
        if( intersection.impl().equals( boundarySegment.impl() ) != boundarySegment.impl().equals( intersection.impl() ) )
          DUNE_THROW( GridError, "Comparison of intersection and boundary segment failed" );

        // equality check should be symmetrical
        if( intersection.impl().equals( boundarySegment.impl() ) )
        {
          // assert that index in inside conincides
          if( intersection.indexInInside() != indexInInside )
            DUNE_THROW( GridError, "Index in inside of boundary segment and intersection in inside entity do not conincide" );

          // assert boundary segment is not found twice
          if( (boundaryFace[ indexSet.index( inside ) ][ indexInInside ] -= 1) < 0 )
            DUNE_THROW( GridError, "Boundary segment found twice in inside entity" );
        }
      }
    }
    std::cerr << "Found " << bndSegments << " boundary segments" << std::endl;

    // check we found all boundary intersections
    for( typename std::vector< Data >::iterator vit = boundaryFace.begin(); vit != boundaryFace.end(); ++vit )
    {
      for( size_t face = 0; face < vit->size(); ++face )
        if( (*vit)[ face ] != 0 )
          DUNE_THROW( GridError, "Not all boundary intersections were found" );
    }
  }

}

#endif // #ifndef DUNE_SPGRID_CHECKBNDSEGITERATOR_HH
