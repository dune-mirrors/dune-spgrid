#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <cassert>
#include <vector>

#include <dune/grid/genericgeometry/misc.hh>
#include <dune/grid/genericgeometry/codimtable.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/domain.hh>

namespace Dune
{

  template< class ct, int dim >
  class SPGridLevel
  {
    typedef SPGridLevel< ct, dim > This;

  public:
    typedef SPDomain< ct, dim > Domain;
    typedef SPGridLevel< ct, dim > GridLevel;

    typedef typename Domain::ctype ctype;
    static const int dimension = Domain::dimension;

    typedef typename Domain::GlobalVector GlobalVector;
    typedef unsigned int MultiIndex[ dimension ];

  public:
    template< int codim >
    struct GeometryCache;
    
  private:
    typedef GenericGeometry::CodimTable< GeometryCache, dimension > GeometryCacheTable;

  public:
    SPGridLevel ( const Domain &domain, const MultiIndex &n )
    : father_( 0 ),
      domain_( &domain ),
      level_( 0 )
    {
      const GlobalVector &width  = domain.width();
      for( int i = 0; i < dimension; ++i )
      {
        n_[ i ] = n[ i ];
        h_[ i ] = width[ i ] / (ctype)n[ i ];
      }

      for( int codim = 0; codim <= dimension; ++codim )
      {
        const unsigned int numDirections = SPDirection::numDirections( dimension, codim );
        multiDirection_[ codim ].resize( numDirections );
        for( unsigned int dir = 0; dir < numDirections; ++dir )
          SPDirection::multiIndex( dimension, codim, dir, multiDirection_[ codim ][ dir ] );
      }
      GenericGeometry::ForLoop< GeometryCache, 0, dimension >
        ::apply( h_, multiDirection_, geometryCache_ );
    }

    SPGridLevel ( const GridLevel &father )
    : father_( &father ),
      domain_( father.domain_ ),
      level_( father.level_+1 )
    {
      for( int i = 0; i < dimension; ++i )
      {
        n_[ i ] = 2*father.n_[ i ];
        h_[ i ] = 0.5*father.h_[ i ];
      }

      for( int codim = 0; codim <= dimension; ++codim )
        multiDirection_[ codim ] = father.multiDirection_[ codim ];
      GenericGeometry::ForLoop< GeometryCache, 0, dimension >
        ::apply( h_, multiDirection_, geometryCache_ );
    }

    const Domain &domain () const
    {
      return *domain_;
    }

    const GlobalVector &h () const
    {
      return h_;
    }

    unsigned int level () const
    {
      return level_;
    }

    const MultiIndex &
    MultiDirection ( const unsigned int codim, const unsigned int dir ) const
    {
      return multiDirection_[ codim ][ dir ];
    }

    const MultiIndex &n () const
    {
      return n_;
    }

    template< int codim >
    const GeometryCache< codim > geometryCache () const
    {
      Int2Type< codim > codimVariable;
      return geometryCache_[ codimVariable ];
    }

  private:
    const GridLevel *father_;
    const Domain *domain_;
    unsigned int level_;
    MultiIndex n_;
    GlobalVector h_;
    std::vector< MultiIndex > multiDirection_[ dimension+1 ];
    GeometryCacheTable geometryCache_;
  };



  template< class ct, int dim >
  template< int codim >
  class SPGridLevel< ct, dim >::GeometryCache
  {
    typedef GeometryCache< codim > This;

    friend class SPGridLevel< ct, dim >;

  public:
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldMatrix< ctype, dimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, dimension > JacobianTransposed;

  private:
    static void
    apply ( const GlobalVector &h, const std::vector< MultiIndex > &multiDirection,
            GeometryCacheTable &cache )
    {
      Int2Type< codim > codimVariable;
      cache[ codimVariable ].initialize( h, multiDirection );
    }

    void initialize ( const GlobalVector &h, const std::vector< MultiIndex > &multiDirection )
    {
      const unsigned int numDirections = multiDirection.size();

      volume_.resize( numDirections );
      jacobianTransposed_.resize( numDirections );
      jacobianInverseTransposed_.resize( numDirections );

      for( unsigned int dir = 0; dir < numDirections; ++dir )
      {
        volume_[ dir ] = ctype( 1 );
        jacobianTransposed_[ dir ] = ctype( 0 );
        jacobianInverseTransposed_[ dir ] = ctype( 0 );

        int j = 0;
        for( int i = 0; i < dimension; ++i )
        {
          if( multiDirection[ dir ][ i ] != 0 )
            continue;
          volume_ *= h[ i ];
          jacobianTransposed_[ dir ][ i ][ j ] = h[ i ];
          jacobianInverseTransposed_[ dir ][ j ][ i ] = ctype( 1 ) / h[ i ];
          ++j;
        }
      }
    }

  public:
    const ctype &volume ( const unsigned int dir ) const
    {
      assert( dir < volume_.size() );
      return volume_[ dir ];
    }

    const JacobianTransposed &jacobianTransposed ( const unsigned int dir ) const
    {
      assert( dir < jacobianTransposed_.size() );
      return jacobianTransposed_[ dir ];
    }

    const Jacobian &jacobianInverseTransposed ( const unsigned int dir ) const
    {
      assert( dir < jacobianInverseTransposed_.size() );
      return jacobianInverseTransposed_[ dir ];
    }

  private:
    std::vector< ctype > volume_;
    std::vector< JacobianTransposed > jacobianTransposed_;
    std::vector< Jacobian > jacobianInverseTransposed_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH
