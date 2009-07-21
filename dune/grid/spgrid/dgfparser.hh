#ifndef DUNE_SPGRID_DGFPARSER_HH
#define DUNE_SPGRID_DGFPARSER_HH

#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/spgrid.hh>

namespace Dune
{

  // MacroGrid::Impl for SPGrid
  // --------------------------

  template< class ct, int dim >
  struct MacroGrid::Impl< SPGrid< ct, dim > >
  {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    
    static SPGrid< ct, dim > *
    generate ( MacroGrid &macroGrid, const std::string &filename,
               MPICommunicatorType MPICOMM = MPIHelper::getCommunicator() );
  };


  template< class ct, int dim >
  inline SPGrid< ct, dim > *
  MacroGrid::Impl< SPGrid< ct, dim > >
    ::generate ( MacroGrid &macroGrid, const std::string &filename, MPICommunicatorType )
  {
    macroGrid.element = Cube;
    macroGrid.dimgrid = SPGrid< ct, dim >::dimension;
    macroGrid.dimw = SPGrid< ct, dim >::dimensionworld;

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

    bool periodic[ dim ];
    for( int i = 0; i < dim; ++i )
      periodic[ i ] = false;

    const int numTrafos = trafoBlock.numTransformations();
    for( int k = 0; k < numTrafos; ++k )
    {
      const Transformation &trafo = trafoBlock.transformation( k );

      bool identity = true;
      for( int i = 0; i < dim; ++i )
        for( int j = 0; j < dim; ++j )
          identity &= (fabs( (i == j ? 1.0 : 0.0) - trafo.matrix( i, j ) ) < 1e-10);
      if( !identity )
        DUNE_THROW( DGFException, "YaspGrid can only handle shifts as periodic face transformations." );

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
        periodic[ dir ] = true;
    }

    dgf::GridParameterBlock parameter( file );
    std::string gridName = parameter.name( "SPGrid" );

    return new SPGrid< ct, dim >( a, b, cells, gridName );
  }



  // DGFGridInfo for SPGrid
  // ----------------------

  template< class ct, int dim >
  struct DGFGridInfo< SPGrid< ct, dim > >
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
