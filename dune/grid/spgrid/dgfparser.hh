#ifndef DUNE_SPGRID_DGFPARSER_HH
#define DUNE_SPGRID_DGFPARSER_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/spgrid.hh>

namespace Dune
{

  // DGFGridFactory< SPGrid >
  // ------------------------

  template< class ct, int dim, SPRefinementStrategy strategy >
  struct DGFGridFactory< SPGrid< ct, dim, strategy > >
  {
    typedef SPGrid< ct, dim, strategy > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    static const int dimension = Grid::dimension;
    typedef typename Grid::template Codim< 0 >::Entity Element;
    typedef typename Grid::template Codim< dimension >::Entity Vertex;

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() );

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
      return intersection.boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      DUNE_THROW( InvalidStateException,
                  "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
    }

  private:
    Grid *grid_;
  };


  template< class ct, int dim, SPRefinementStrategy strategy >
  inline DGFGridFactory< SPGrid< ct, dim, strategy > >
    ::DGFGridFactory ( const std::string &filename, MPICommunicatorType )
  {
    typedef SPGrid< ct, dim, strategy > Grid;

    std::ifstream file( filename.c_str() );
    dgf::IntervalBlock intervalBlock( file );

    if( !intervalBlock.isactive() )
      DUNE_THROW( DGFException, "DGF file must contain an interval block to be used with SPGrid< ct, " << dim << " >." );
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
    dgf::PeriodicFaceTransformationBlock trafoBlock( file, dim );

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

    dgf::GridParameterBlock parameter( file );
    const std::string gridName = parameter.name( "SPGrid" );

    grid_ = new Grid( a, b, cells, periodic, gridName );
  }



  // DGFGridInfo< SPGrid >
  // ---------------------

  template< class ct, int dim, SPRefinementStrategy strategy >
  struct DGFGridInfo< SPGrid< ct, dim, strategy > >
  {
    static int refineStepsForHalf ()
    {
      return 1;
    }

    static double refineWeight ()
    {
      return 1.0 / double( 1 << dim );
    }
  };

}

#endif // #ifndef DUNE_SPGRID_DGFPARSER_HH
