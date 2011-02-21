#ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
#define DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH

#include <cassert>

#include <dune/common/forloop.hh>

#include <dune/grid/spgrid/referencecube.hh>
#include <dune/grid/spgrid/geometrycache.hh>

namespace Dune
{

  // SPGeometricGridLevel
  // --------------------

  template< class Grid >
  class SPGeometricGridLevel
  {
    typedef SPGeometricGridLevel< Grid > This;

  public:
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ReferenceCube ReferenceCube;

    typedef typename ReferenceCube::ctype ctype;
    static const int dimension = ReferenceCube::dimension;

    typedef typename ReferenceCube::GlobalVector GlobalVector;

    static const unsigned int numDirections = ReferenceCube::numCorners;

    template< int codim >
    struct Codim
    {
      typedef SPReferenceCube< ctype, dimension-codim > ReferenceCube;
      typedef SPGeometryCache< ctype, dimension, codim > GeometryCache;
    };

  private:
    template< int codim >
    struct BuildGeometryCache;
    template< int codim >
    struct DestroyGeometryCache;

  public:
    SPGeometricGridLevel ( const Grid &grid, const GlobalVector &h );
    SPGeometricGridLevel ( const This &other );

    ~SPGeometricGridLevel ();

    const Grid &grid () const { return *grid_; }

    const ReferenceCube &referenceCube () const { return grid().referenceCube(); }

    template< int codim >
    const typename Codim< codim >::ReferenceCube &
    referenceCube () const { return grid().template referenceCube< codim >(); }

    const GlobalVector &h () const { return h_; }

    template< int codim >
    const typename Codim< codim >::GeometryCache &
    geometryCache ( const unsigned int dir ) const;

    ctype faceVolume ( const int i ) const;
    const GlobalVector &volumeNormal ( const int i ) const;

  private:
    void buildGeometry ();

    const Grid *grid_;

    GlobalVector h_;
    void *geometryCache_[ numDirections ];
    ctype faceVolume_[ ReferenceCube::numFaces ];
    GlobalVector normal_[ ReferenceCube::numFaces ];
  };



  // SPGeometricGridLevel::BuildGeometryCache
  // ----------------------------------------

  template< class Grid >
  template< int codim >
  struct SPGeometricGridLevel< Grid >::BuildGeometryCache
  {
    static void
    apply ( const GlobalVector &h, void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        const int mydim = bitCount( dir );
        if( mydim == dimension - codim )
          geometryCache[ dir ] = new GeometryCache( h, dir );
      }
    }
  };



  // SPGeometricGridLevel::DestroyGeometryCache
  // ------------------------------------------

  template< class Grid >
  template< int codim >
  struct SPGeometricGridLevel< Grid >::DestroyGeometryCache
  {
    static void
    apply ( void *(&geometryCache)[ 1 << dimension ] )
    {
      typedef typename Codim< codim >::GeometryCache GeometryCache;
      for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
      {
        const int mydim = bitCount( dir );
        if( mydim == dimension - codim )
        {
          delete (GeometryCache *)geometryCache[ dir ];
          geometryCache[ dir ] = 0;
        }
      }
    }
  };



  // Implementation of SPGeometricGridLevel
  // --------------------------------------

  template< class Grid >
  inline SPGeometricGridLevel< Grid >
    ::SPGeometricGridLevel ( const Grid &grid, const GlobalVector &h )
  : grid_( &grid ),
    h_( h )
  {
    buildGeometry();
  }


  template< class Grid >
  inline SPGeometricGridLevel< Grid >::SPGeometricGridLevel ( const This &other )
  : grid_( other.grid_ ),
    h_( other.h_ )
  {
    buildGeometry();
  }


  template< class Grid >
  inline SPGeometricGridLevel< Grid >::~SPGeometricGridLevel ()
  {
    ForLoop< DestroyGeometryCache, 0, dimension >::apply( geometryCache_ );
  }


  template< class Grid >
  template< int codim >
  inline const typename SPGeometricGridLevel< Grid >::template Codim< codim >::GeometryCache &
  SPGeometricGridLevel< Grid >::geometryCache ( const unsigned int dir ) const
  {
    typedef typename Codim< codim >::GeometryCache GeometryCache;
    assert( bitCount( dir ) == dimension - codim );
    return *((const GeometryCache *)geometryCache_[ dir ]);
  }


  template< class Grid >
  inline typename SPGeometricGridLevel< Grid >::ctype
  SPGeometricGridLevel< Grid >::faceVolume ( const int i ) const
  {
    assert( (i >= 0) && (i < ReferenceCube::numFaces) );
    return faceVolume_[ i ];
  }


  template< class Grid >
  inline const typename SPGeometricGridLevel< Grid >::GlobalVector &
  SPGeometricGridLevel< Grid >::volumeNormal ( const int i ) const
  {
    assert( (i >= 0) && (i < ReferenceCube::numFaces) );
    return normal_[ i ];
  }


  template< class Grid >
  inline void SPGeometricGridLevel< Grid >::buildGeometry ()
  {
    ForLoop< BuildGeometryCache, 0, dimension >::apply( h_, geometryCache_ );
    
    const ctype volume = geometryCache< 0 >( numDirections-1 ).volume();
    for( int face = 0; face < ReferenceCube::numFaces; ++face )
    {
      normal_[ face ] = referenceCube().normal( face );
      faceVolume_[ face ] = std::abs( volume / (normal_[ face ] * h_) );
      normal_[ face ] *= faceVolume_[ face ];
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
