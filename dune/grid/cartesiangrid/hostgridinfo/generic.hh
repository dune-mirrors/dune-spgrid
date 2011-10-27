#ifndef DUNE_CARTESIANGRID_HOSTGRIDINFO_GENERIC_HH
#define DUNE_CARTESIANGRID_HOSTGRIDINFO_GENERIC_HH

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Dune
{

  template< class HostGrid >
  struct CartesianGridHostGridInfo
  {
    //////////////////////////////////////////////////
    //  direction methods 
    //////////////////////////////////////////////////
    
    //! default direction for given dimension 
    static unsigned int defaultDirection ( const int mydimension )
    {
      return (1 << mydimension)-1;
    }

    //! direction for face i 
    static unsigned int direction ( const int i, const int dimension )
    {
      assert( (i >= 0) && (i < 2*dimension) );
      return (1 << (i / 2)) ^ ((1 << dimension) - 1);
    }

    //! direction for given host entity
    template< class HostEntity >
    static unsigned int direction ( const HostEntity &hostEntity )
    {
      typedef typename HostEntity::Geometry HostGeometry;
      unsigned int direction = 0;
      const HostGeometry &geo = hostEntity.geometry();
      const typename HostGeometry::GlobalCoordinate origin = geo.corner( 0 );
      for( int d = 0; d <HostGeometry::mydimension; ++d )
      {
        const typename HostGeometry::GlobalCoordinate point = geo.corner( 1 << d );
        for( int i = 0; i < HostGeometry::coorddimension; ++i )
        {
          if( std::abs( point[ i ] - origin[ i ] ) > 1e-8 )
            direction |= (1 << i);
        }
      }
      return direction;
    }

    //////////////////////////////////////////////////////
    //  origin methods 
    //////////////////////////////////////////////////////

    template< int codim >
    struct Codim
    {
      // maybe we have to change this 
      typedef typename HostGrid::template Codim< codim >::Geometry::GlobalCoordinate OriginReturnType;
    };

    //! default origin 
    static typename Codim< 0 >::OriginReturnType
    defaultOrigin ()
    {
      return typename Codim< 0 >::OriginReturnType( 0 );
    }

    //! origin for given entity or intersection
    template< class HostItem >
    static typename Codim< HostItem::codimension >::OriginReturnType
    origin ( const HostItem &hostItem )
    {
      return hostItem.geometry().corner( 0 );
    }

    //! origin for given intersection
    template< class HostIntersection >
    static typename Codim< HostIntersection::codimension >::OriginReturnType
    originIntersection ( const HostIntersection &hostIntersection )
    {
      typedef typename HostIntersection::Geometry Geometry;
      const Geometry &geo = hostIntersection.geometry();
      assert( geo.corners() == (1 << Geometry::mydimension) );
      const typename Geometry::GlobalCoordinate a = geo.corner( 0 );
      const typename Geometry::GlobalCoordinate b = geo.corner( (1 << Geometry::mydimension)-1 );
      typename Geometry::GlobalCoordinate origin;
      for( int i = 0; i < Geometry::coorddimension; ++i )
        origin[ i ] = std::min( a[ i ], b[ i ] );
      return origin;
    }

    //////////////////////////////////////////////////////////
    //  child index methods 
    //////////////////////////////////////////////////////////

    //! return child index for entity 
    template< class HostEntity >
    static int childIndex ( const HostEntity &hostEntity )
    {
      int childIndex = 0;
      typedef typename HostGrid::template Codim< HostEntity::codimension >::LocalGeometry LocalGeometry;
      typedef  typename LocalGeometry::GlobalCoordinate GlobalCoordinate;
      const GlobalCoordinate centerInFather = hostEntity.geometryInFather().center();
      for( int i = 0; i < GlobalCoordinate::dimension; ++i )
        childIndex |= (centerInFather[ i ] < 0.5 ? 0 : (1 << i));
      return childIndex;
    }

    //////////////////////////////////////////////////////////
    //  level methods 
    //////////////////////////////////////////////////////////

    //! return inside level for intersection 
    template< class HostIntersection >
    static int insideLevel ( const HostIntersection &hostIntersection )
    {
      return hostIntersection.inside().level();
    }

    //! return outside level for intersection 
    template< class HostIntersection >
    static int outsideLevel ( const HostIntersection &hostIntersection,
                              const int insideLevel )
    {
      return hostIntersection.neighbor() ? hostIntersection.outside().level() : insideLevel;
    }

    template< class HostIntersection >
    static int childIndexInInside ( const HostIntersection &hostIntersection,
                                    const int insideLevel,
                                    const int outsideLevel )
    {
      return getIntersectionChildLevel( hostIntersection.geometryInInside(),
                                        hostIntersection.indexInInside(),
                                        insideLevel,
                                        outsideLevel );
    }

    template< class HostIntersection >
    static int childIndexInOutside ( const HostIntersection &hostIntersection,
                                     const int insideLevel,
                                     const int outsideLevel ) 
    {
      return getIntersectionChildLevel( hostIntersection.geometryInOutside(),
                                        hostIntersection.indexInOutside(),
                                        outsideLevel, 
                                        insideLevel );
    }

  protected:  
    template< class LocalGeometry >
    static int
    getIntersectionChildLevel ( const LocalGeometry &localGeo, const int index,
                                const int myLevel, const int otherLevel )
    {
      typedef typename LocalGeometry::GlobalCoordinate GlobalCoordinate;
      // make sure non-confoming level difference is at most 1 
      assert( std::abs( myLevel - otherLevel ) <= 1 );
      if( myLevel < otherLevel )
      {
        int childIndex = 0;
        const GlobalCoordinate center = localGeo.center();
        for( int i = 0; i < LocalGeometry::mydimension; ++i )
        {
          const int j = (i < index/2 ? i : i+1 );
          childIndex |= (center[ j ] < 0.5 ? 0 : (1 << i));
        }
        return childIndex;
      }
      else
        return -1;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_HOSTGRIDINFO_HH
