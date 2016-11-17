#ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
#define DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH

#include <array>
#include <cassert>

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/geometrycache.hh>
#include <dune/grid/spgrid/normal.hh>
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

    SPGeometricGridLevel ( const ReferenceCubeContainer &refCubes, const GlobalVector &h );
    SPGeometricGridLevel ( const This &other );

    ~SPGeometricGridLevel ();

    const ReferenceCube &referenceCube () const { return refCubes_.get(); }

    template< int codim >
    const typename Codim< codim >::ReferenceCube &
    referenceCube () const { return refCubes_.template get< codim >(); }

    const GlobalVector &h () const { return h_; }

    template< int codim >
    const typename Codim< codim >::GeometryCache &
    geometryCache ( Direction dir, Dune::Codim< codim > = Dune::Codim< codim >() ) const
    {
      assert( dir.mydimension() == dimension - codim );
      return *static_cast< const typename Codim< codim >::GeometryCache * >( geometryCache_[ dir.bits() ] );
    }

    const typename Codim< 0 >::GeometryCache &geometryCache ( Dune::Codim< 0 > ) const
    {
      SPDirectionIterator< dimension, 0 > dirIt;
      return geometryCache< 0 >( *dirIt );
    }

    ctype faceVolume ( int i ) const { assert( (i >= 0) && (i < ReferenceCube::numFaces) ); return faceVolume_[ i ]; }

  private:
    void buildGeometry ();

    const ReferenceCubeContainer &refCubes_;

    GlobalVector h_;
    std::array< void *, numDirections > geometryCache_;
    std::array< ctype, ReferenceCube::numFaces > faceVolume_;
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
    Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this ] ( auto codim ) {
        for( SPDirectionIterator< dimension, codim > dirIt; dirIt; ++dirIt )
        {
          delete static_cast< typename Codim< codim >::GeometryCache * >( geometryCache_[ (*dirIt).bits() ] );
          geometryCache_[ (*dirIt).bits() ] = nullptr;
        }
      } );
  }


  template< class ct, int dim >
  inline void SPGeometricGridLevel< ct, dim >::buildGeometry ()
  {
    Hybrid::forEach( std::make_integer_sequence< int, dimension+1 >(), [ this ] ( auto codim ) {
        typedef typename Codim< codim >::GeometryCache GeometryCache;
        for( SPDirectionIterator< dimension, codim > dirIt; dirIt; ++dirIt )
          geometryCache_[ (*dirIt).bits() ] = new GeometryCache( h_, *dirIt );
      } );

    const ctype volume = geometryCache( Dune::Codim< 0 >() ).volume();
    for( int face = 0; face < ReferenceCube::numFaces; ++face )
    {
      SPNormalId< dimension > normalId( face );
      SPNormalVector< ctype, dimension > normal( normalId );
      faceVolume_[ face ] = std::abs( volume / (normal * h_) );
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GEOMETRICGRIDLEVEL_HH
