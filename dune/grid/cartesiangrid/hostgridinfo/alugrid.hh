#ifndef DUNE_CARTESIANGRID_HOSTGRIDINFO_ALUGRID_HH
#define DUNE_CARTESIANGRID_HOSTGRIDINFO_ALUGRID_HH

#if ENABLE_ALUGRID
#include <dune/common/fvector.hh>

#include <dune/grid/alugrid/3d/grid.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< class HostGrid >
  struct CartesianGridHostGridInfo;



  // ALU3dCubeOrigin
  // ---------------

  template< int dimension, int codimension >
  struct ALU3dCubeOrigin
  {
    typedef ALUCubeGrid< dimension, dimension > HostGrid;

    typedef typename HostGrid::template Codim< codimension >::Geometry::GlobalCoordinate OriginReturnType;

    template< class Item >
    static OriginReturnType get ( const Item &item )
    {
      return item.geometry().corner( 0 );
    }
  };

  // for vertices 
  template< int dimension >
  struct ALU3dCubeOrigin< dimension, dimension >
  {
    typedef ALUCubeGrid< dimension, dimension > HostGrid;

    typedef double coord_t[ dimension ];
    typedef const coord_t &OriginReturnType;

    template< class Item >
    static OriginReturnType get ( const Item &item )
    {
      return HostGrid::getRealImplementation( item ).getItem().Point();
    }
  };

  // for elements 
  template< int dimension >
  struct ALU3dCubeOrigin< dimension, 0 >
  {
    typedef ALUCubeGrid< dimension, dimension > HostGrid;

    typedef double coord_t[ dimension ];
    typedef const coord_t &OriginReturnType;

    template< class Item >
    static OriginReturnType get ( const Item &item )
    {
      return HostGrid::getRealImplementation( item ).getItem().myvertex( 0 )->Point();
    }
  };

  // specialization for intersections 
  template< int dimension >
  struct ALU3dCubeOrigin< dimension, -1 >
  {
    typedef ALUCubeGrid< dimension, dimension > HostGrid;

    typedef double coord_t[ dimension ];
    typedef const coord_t &OriginReturnType;

    template< class Intersection >
    static OriginReturnType get( const Intersection &intersection )
    {
#ifndef NDEBUG 
      typedef FaceTopologyMapping< hexa > CubeFaceMapping;
      const int duneTwist = HostGrid::getRealImplementation( intersection ).twistInInside();
      const int twistedIndex = CubeFaceMapping::twistedDuneIndex( 0, duneTwist );
      assert( twistedIndex == 0 );
#endif
      return HostGrid::getRealImplementation( intersection ).it().getItem().myvertex( 0 )->Point();
    }
  };



  // CartesianGridHostGridInfo for ALUCubeGrid< 3, 3 >
  // -------------------------------------------------

  template<> 
  struct CartesianGridHostGridInfo< ALUCubeGrid< 3, 3 > >
  {
    static const int dimension = 3;

    typedef ALUCubeGrid< 3, 3 > HostGrid; 

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
    static unsigned int direction ( const HostEntity& hostEntity )
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

    //! default origin
    static FieldVector< double, dimension > defaultOrigin ()
    {
      return FieldVector< double, dimension >( 0 );
    }

    //! origin for given entity or intersection 
    template< class HostItem >
    static typename ALU3dCubeOrigin< HostItem::dimension, HostItem::codimension >::OriginReturnType
    origin ( const HostItem &hostItem )
    {
      return ALU3dCubeOrigin< HostItem::dimension, HostItem::codimension >::get( hostItem );
    }

    //! origin for given entity or intersection 
    template< class HostItem >
    static typename ALU3dCubeOrigin< HostItem::dimension, -HostItem::codimension >::OriginReturnType
    originIntersection ( const HostItem &hostItem )
    {
      return ALU3dCubeOrigin< HostItem::dimension, -HostItem::codimension >::get( hostItem );
    }

    //////////////////////////////////////////////////////////
    //  child index methods 
    //////////////////////////////////////////////////////////

    //! return child index for entity 
    template< class HostEntity >
    static int childIndex ( const HostEntity &hostEntity )
    {
      // apply the same change as for the vertices of the hexa 
      typedef ElementTopologyMapping< hexa > ElemTopo;
      return ElemTopo::alu2duneVertex( HostGrid::getRealImplementation( hostEntity ).getItem().nChild() );
    }

    //////////////////////////////////////////////////////////
    //  level methods 
    //////////////////////////////////////////////////////////

    //! return inside level for intersection 
    template <class HostIntersection>
    static int insideLevel ( const HostIntersection &hostIntersection )
    {
      assert( HostGrid::getRealIntersection( hostIntersection ).level() == hostIntersection.inside().level() );
      return HostGrid::getRealIntersection( hostIntersection ).level();
    }

    //! return outside level for intersection 
    template< class HostIntersection >
    static int outsideLevel ( const HostIntersection &hostIntersection, 
                              const int insideLevel )
    {
      // outsideLevel might be less than insideLevel if there is no level neighbor
      const int outsideLevel = HostGrid::getRealIntersection( hostIntersection ).it().outsideLevel();
      assert( !hostIntersection.neighbor() || (outsideLevel == hostIntersection.outside().level()) );
      return outsideLevel;
    }

    template< class HostIntersection >
    static int childIndexInInside ( const HostIntersection &hostIntersection,
                                    const int insideLevel,
                                    const int outsideLevel ) 
    {
      return getIntersectionChildLevel( HostGrid::getRealIntersection( hostIntersection ).twistInInside(),
                                        HostGrid::getRealIntersection( hostIntersection ).it().getItem().nChild(),
                                        insideLevel, outsideLevel );
    }

    template< class HostIntersection >
    static int childIndexInOutside ( const HostIntersection &hostIntersection,
                                     const int insideLevel,
                                     const int outsideLevel ) 
    {
      return getIntersectionChildLevel( HostGrid::getRealIntersection( hostIntersection ).twistInOutside(),
                                        HostGrid::getRealIntersection( hostIntersection ).it().getItem().nChild(),
                                        outsideLevel, insideLevel );
    }

  protected:  
    static int
    getIntersectionChildLevel( const int duneTwist, const int child,
                               const int myLevel, const int otherLevel )
    {
      // make sure non-confoming level difference is at most 1 
      assert( std::abs( myLevel - otherLevel ) <= 1 );
      if( myLevel < otherLevel )
      {
        // swap children 2 and 3 
        static const int map[4] = { 0, 1, 3, 2 }; 
        return map[ child ];
      }
      else
        return -1;
    }
  };

} // namespace Dune 

#endif // #if ENABLE_ALUGRID

#endif // #ifndef DUNE_CARTESIANGRID_HOSTGRIDINFO_HH
