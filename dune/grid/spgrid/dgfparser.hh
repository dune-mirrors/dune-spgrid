#ifndef DUNE_SPGRID_DGFPARSER_HH
#define DUNE_SPGRID_DGFPARSER_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/spgrid.hh>

namespace Dune
{

  namespace dgf
  {

    // SPGridParameterBlock
    // --------------------

    template< int dim >
    struct SPGridParameterBlock
    : public GridParameterBlock
    {
      typedef int MultiIndex[ dim ];

      explicit SPGridParameterBlock( std::istream &in );

      const MultiIndex &overlap () const
      {
        return overlap_;
      }

    private:
      MultiIndex overlap_;
    };



    // Implementation of SPGridParameterBlock
    // --------------------------------------

    template< int dim >
    inline SPGridParameterBlock< dim >::SPGridParameterBlock( std::istream &in )
    : GridParameterBlock( in )
    {
      for( int i = 0; i < dim; ++i )
        overlap_[ i ] = 0;

      if( findtoken( "overlap" ) )
      {
        int i = 0;
        for( int x; getnextentry( x ); ++i )
        {
          if( x < 0 )
            DUNE_THROW( DGFException, "Negative overlap specified." );
          if( i < dim )
            overlap_[ i ] = x;
        }

        if( i < 1 )
          dwarn << "GridParameterBlock: Found keyword 'overlap' without valid value, defaulting to no overlap." << std::endl;
        else if( i == 1 )
        {
          for( int i = 1; i < dim; ++i )
            overlap_[ i ] = overlap_[ 0 ];
        }
        else if( i != dim )
          DUNE_THROW( DGFException, "Invalid argument for parameter 'overlap' specified." );
      }
      else
        dwarn << "GridParameterBlock: Parameter 'overlap' not specified, defaulting to no overlap." << std::endl;
    }

  } // namespace dgf



  // DGFGridFactory< SPGrid >
  // ------------------------

  template< class ct, int dim, template< int > class Ref, class Comm >
  class DGFGridFactory< SPGrid< ct, dim, Ref, Comm > >
  {
  public:
    typedef SPGrid< ct, dim, Ref, Comm > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::Communication Communication;

    static const int dimension = Grid::dimension;
    typedef typename Grid::template Codim< 0 >::Entity Element;
    typedef typename Grid::template Codim< dimension >::Entity Vertex;

  private:
    typedef FieldVector< double, dimension > Point;

    typedef dgf::BoundaryDomBlock BoundaryDomainBlock;

  public:
    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

    ~DGFGridFactory ()
    {
      delete boundaryDomainBlock_;
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      std::vector< Point > corners;
      getCorners( intersection.geometry(), corners );
      const dgf::DomainData *data = boundaryDomainBlock_->contains( corners );
      if( data )
        return data->id();
      else
        return intersection.indexInInside();
    }

    bool haveBoundaryParameters () const
    {
      return boundaryDomainBlock_->hasParameter();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
      std::vector< Point > corners;
      getCorners( intersection.geometry(), corners );
      const dgf::DomainData *data = boundaryDomainBlock_->contains( corners );
      if( data )
        return data->parameter();
      else
        return DGFBoundaryParameter::defaultValue();
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      DUNE_THROW( InvalidStateException,
                  "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
    }

  private:
    void generate ( std::istream &input, const Communication &comm );

    template< class Geometry >
    static void getCorners ( const Geometry &geometry, std::vector< Point > &corners )
    {
      corners.resize( geometry.corners() );
      for( int i = 0; i < geometry.corners(); ++i )
      {
        const typename Geometry::GlobalCoordinate corner = geometry.corner( i );
        for( int j = 0; j < dimension; ++j )
          corners[ i ][ j ] = corner[ j ];
      }
    }

    Grid *grid_;
    BoundaryDomainBlock *boundaryDomainBlock_;
  };



  // Implementation of DGFGridFactory< SPGrid >
  // ------------------------------------------

  template< class ct, int dim, template< int > class Ref, class Comm >
  inline DGFGridFactory< SPGrid< ct, dim, Ref, Comm > >
    ::DGFGridFactory ( std::istream &input, MPICommunicatorType comm )
  {
    generate( input, SPCommunicationTraits< Comm >::comm( comm ) );
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline DGFGridFactory< SPGrid< ct, dim, Ref, Comm > >
    ::DGFGridFactory ( const std::string &filename, MPICommunicatorType comm )
  {
    std::ifstream input( filename.c_str() );
    if( !input )
      DUNE_THROW( DGFException, "Unable to open file: " << filename << "." );
    generate( input, SPCommunicationTraits< Comm >::comm( comm ) );
    input.close();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void
  DGFGridFactory< SPGrid< ct, dim, Ref, Comm > >
    ::generate ( std::istream &input, const Communication &comm )
  {
    dgf::IntervalBlock intervalBlock( input );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "DGF stream must contain an interval block to be used with SPGrid< ct, " << dim << " >." );
    if( intervalBlock.numIntervals() != 1 )
      DUNE_THROW( DGFException, "SPGrid< ct, " << dim << " > can only handle 1 interval block." );

    if( intervalBlock.dimw() != dim )
      DUNE_THROW( DGFException, "SPGrid< ct, " << dim << " cannot handle an interval of dimension " << intervalBlock.dimw() << "." );
    const dgf::IntervalBlock::Interval &interval = intervalBlock.get( 0 );

    FieldVector< ct, dim > a, b;
    int cells[ dim ];
    for( int i = 0; i < dim; ++i )
    {
      a[ i ] = interval.p[ 0 ][ i ];
      b[ i ] = interval.p[ 1 ][ i ];
      cells[ i ] = interval.n[ i ];
    }

    typedef dgf::PeriodicFaceTransformationBlock::AffineTransformation Transformation;
    dgf::PeriodicFaceTransformationBlock trafoBlock( input, dim );

    unsigned int periodic = 0;

    const int numTrafos = trafoBlock.numTransformations();
    for( int k = 0; k < numTrafos; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      bool identity = true;
      for( int i = 0; i < dim; ++i )
        for( int j = 0; j < dim; ++j )
          identity &= (fabs( (i == j ? 1.0 : 0.0) - trafo.matrix( i, j ) ) < 1e-10);
      if( !identity )
        DUNE_THROW( DGFException, "SPGrid< ct, " << dim << " > can only handle shifts as periodic face transformations." );

      int numDirs = 0;
      int dir = -1;
      for( int i = 0; i < dim; ++i )
      {
        if( fabs( trafo.shift[ i ] ) < 1e-10 )
          continue;
        dir = i;
        ++numDirs;
      }
      const ct shift = fabs( trafo.shift[ dir ] );
      const ct width = fabs( a[ dir ] - b[ dir ] );
      if( (numDirs != 1) || (fabs( shift - width ) >= 1e-10) )
      {
        std::cerr << "Tranformation '" << trafo
                  << "' does not map boundaries on boundaries." << std::endl;
      }
      else
        periodic |= (1 << dir);
    }

    dgf::SPGridParameterBlock< dim > parameter( input );

    typedef typename Grid::Domain Domain;
    std::vector< typename Domain::Cube > cubes;
    cubes.push_back( typename Domain::Cube( a, b ) );
    Domain domain( cubes, typename Domain::Topology( periodic ) );
    grid_ = new Grid( domain, cells, parameter.overlap(), comm );

    boundaryDomainBlock_ = new BoundaryDomainBlock( input, dimension );
  }



  // DGFGridInfo< SPGrid >
  // ---------------------

  template< class ct, int dim, template< int > class Ref, class Comm >
  struct DGFGridInfo< SPGrid< ct, dim, Ref, Comm > >
  {
    typedef SPGrid< ct, dim, Ref, Comm > Grid;
    typedef typename Grid::RefinementPolicy RefinementPolicy;

    static int refineStepsForHalf ( const RefinementPolicy &policy = RefinementPolicy() )
    {
      const unsigned int weight = policy.weight();
      return (dim + weight - 1) / weight;
    }

    static double refineWeight ( const RefinementPolicy &policy = RefinementPolicy() )
    {
      return 1.0 / double( 1 << policy.weight() );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_DGFPARSER_HH
