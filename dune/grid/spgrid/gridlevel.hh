#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <vector>

#include <dune/grid/genericgeometry/codimtable.hh>

namespace Dune
{

  template< class ct, int dim >
  class SPGridLevel
  {
    typedef SPGridLevel< ct, dim > This;

  public:
    typedef SPDomain< ct, dim > Domain;
    typedef SPGridLevel< dt, dim > GridLevel;

    typedef typename Domain::ctype ctype;
    static const int dimension = Domain::dimension;

    typedef typename Domain::GlobalVector GlobalVector;
    typedef unsigned int MultiIndex[ dimension ];

  private:
    template< int codim >
    struct GeometryCache;

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

    unsigned int
    n ( const unsigned int codim, const unsigned int dir, int i ) const
    {
      return n_[ i ] + multiDirection_[ codim ][ dir ][ i ];
    }

  private:
    const GridLevel *father_;
    const Domain *domain_;
    unsigned int level_;
    MultiIndex n_;
    GlobalVector h_;
    std::vector< MultiIndex > multiDirection_[ dimension+1 ];
  };



  template< class ct, int dim >
  template< int codim >
  class SPGridLevel< ct, dim >::GeometryCache
  {
    typedef GeometryCache< codim > This;

  public:
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef FieldVector< ctype, mydimension > LocalVector;
    typedef FieldMatrix< ctype, dimension, mydimension > Jacobian;
    typedef FieldMatrix< ctype, mydimension, dimension > JacobianTransposed;

    GeometryCache ( const GlobalVector &h, const std::vector< MultiIndex > &multiDirection )
    : volume_( multiDirection.size() ),
      jacobianTransposed( multiDirection.size() ),
      jacobianInverseTransposed( multiDirection.size() )
    {
      // ...
    }

  private:
    std::vector< ctype > volume_;
    std::vector< JacobianTransposed > jacobianTransposed_;
    std::vector< Jacobian > jacobianInverseTransposed_;
  };

}

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH
