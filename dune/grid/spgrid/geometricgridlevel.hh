#ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
#define DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH

#include <cassert>

#include <dune/common/forloop.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/geometrycache.hh>
#include <dune/grid/spgrid/referencecube.hh>

namespace Dune
{

  // SPGeometricGridLevel
  // --------------------

  template< class ct, int dim >
  class SPGeometricGridLevel
  {
    typedef SPGeometricGridLevel< ct, dim > This;

  public:
    typedef SPReferenceCubeContainer< ct, dim > ReferenceCubeContainer;

    typedef typename ReferenceCubeContainer::ReferenceCube ReferenceCube;

    typedef typename ReferenceCube::ctype ctype;
    static const int dimension = ReferenceCube::dimension;

    typedef typename ReferenceCube::GlobalVector GlobalVector;

    static const unsigned int numDirections = ReferenceCube::numCorners;

    typedef SPDirection< dimension > Direction;

    template< int codim >
    struct Codim
    {
      typedef typename ReferenceCubeContainer::template Codim< codim >::ReferenceCube ReferenceCube;
      typedef SPGeometryCache< ctype, dimension, codim > GeometryCache;
    };

  private:
    template< int codim >
    struct BuildGeometryCache;
    template< int codim >
    struct DestroyGeometryCache;

  public:
    SPGeometricGridLevel ( const ReferenceCubeContainer &refCubes, const GlobalVector &h );
    SPGeometricGridLevel ( const This &other );

    ~SPGeometricGridLevel ();

    const ReferenceCube &referenceCube () const { return refCubes_.get(); }

    template< int codim >
    const typename Codim< codim >::ReferenceCube &
    referenceCube () const { return refCubes_.template get< codim >(); }

    const GlobalVector &h () const { return h_; }

    template< int codim >
    const typename Codim< codim >::GeometryCache &geometryCache ( Direction dir ) const;

    ctype faceVolume ( const int i ) const;
    const GlobalVector &volumeNormal ( const int i ) const;

  private:
    void buildGeometry ();

    const ReferenceCubeContainer &refCubes_;

    GlobalVector h_;
    void *geometryCache_[ numDirections ];
    ctype faceVolume_[ ReferenceCube::numFaces ];
    GlobalVector normal_[ ReferenceCube::numFaces ];
  };



  // SPGeometricGridLevel::BuildGeometryCache
  // ----------------------------------------

  template< class ct, int dim >
  template< int codim >
  struct SPGeometricGridLevel< ct, dim >::BuildGeometryCache
  {
    static void apply ( const GlobalVector &h, void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( SPDirectionIterator< dimension, codim > dirIt; dirIt; ++dirIt )
        geometryCache[ (*dirIt).bits() ] = new GeometryCache( h, *dirIt );
    }
  };



  // SPGeometricGridLevel::DestroyGeometryCache
  // ------------------------------------------

  template< class ct, int dim >
  template< int codim >
  struct SPGeometricGridLevel< ct, dim >::DestroyGeometryCache
  {
    static void apply ( void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( SPDirectionIterator< dimension, codim > dirIt; dirIt; ++dirIt )
      {
        delete static_cast< GeometryCache * >( geometryCache[ (*dirIt).bits() ] );
        geometryCache[ (*dirIt).bits() ] = nullptr;
      }
    }
  };



  // Implementation of SPGeometricGridLevel
  // --------------------------------------

  template< class ct, int dim >
  inline SPGeometricGridLevel< ct, dim >
    ::SPGeometricGridLevel ( const ReferenceCubeContainer &refCubes, const GlobalVector &h )
  : refCubes_( refCubes ),
    h_( h )
  {
    buildGeometry();
  }


  template< class ct, int dim >
  inline SPGeometricGridLevel< ct, dim >::SPGeometricGridLevel ( const This &other )
  : refCubes_( other.refCubes_ ),
    h_( other.h_ )
  {
    buildGeometry();
  }


  template< class ct, int dim >
  inline SPGeometricGridLevel< ct, dim >::~SPGeometricGridLevel ()
  {
    ForLoop< DestroyGeometryCache, 0, dimension >::apply( geometryCache_ );
  }


  template< class ct, int dim >
  template< int codim >
  inline const typename SPGeometricGridLevel< ct, dim >::template Codim< codim >::GeometryCache &
  SPGeometricGridLevel< ct, dim >::geometryCache ( Direction dir ) const
  {
    typedef typename Codim< codim >::GeometryCache GeometryCache;
    assert( dir.mydimension() == dimension - codim );
    return *static_cast< const GeometryCache * >( geometryCache_[ dir.bits() ] );
  }


  template< class ct, int dim >
  inline typename SPGeometricGridLevel< ct, dim >::ctype
  SPGeometricGridLevel< ct, dim >::faceVolume ( const int i ) const
  {
    assert( (i >= 0) && (i < ReferenceCube::numFaces) );
    return faceVolume_[ i ];
  }


  template< class ct, int dim >
  inline const typename SPGeometricGridLevel< ct, dim >::GlobalVector &
  SPGeometricGridLevel< ct, dim >::volumeNormal ( const int i ) const
  {
    assert( (i >= 0) && (i < ReferenceCube::numFaces) );
    return normal_[ i ];
  }


  template< class ct, int dim >
  inline void SPGeometricGridLevel< ct, dim >::buildGeometry ()
  {
    ForLoop< BuildGeometryCache, 0, dimension >::apply( h_, geometryCache_ );
    
    SPDirectionIterator< dimension, 0 > dirIt;
    const ctype volume = geometryCache< 0 >( *dirIt ).volume();
    for( int face = 0; face < ReferenceCube::numFaces; ++face )
    {
      normal_[ face ] = referenceCube().normal( face );
      faceVolume_[ face ] = std::abs( volume / (normal_[ face ] * h_) );
      normal_[ face ] *= faceVolume_[ face ];
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
