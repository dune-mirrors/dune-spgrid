#include <config.h>

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/common/mpihelper.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/test/checkgeometryinfather.cc>

static const int dimGrid = DIMGRID;

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  typedef Dune::SPGrid< double, dimGrid > Grid;

  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf file> <max level>" << std::endl;
    return 1;
  }

  std::string dgfFile( argv[ 1 ] );
  const int maxLevel = atoi( argv[ 2 ] );

  Dune::GridPtr< Grid > gridPtr( dgfFile );
  Grid &grid = *gridPtr;

#if 0
  Dune::FieldVector< double, dimGrid > a( 0.0 );
  Dune::FieldVector< double, dimGrid > b( 1.0 );
  int n[ dimGrid ];
  for( int i = 0; i < dimGrid; ++i )
    n[ i ] = 4;
  Grid grid( a, b, n );
#endif

  for( int i = 0; i <= maxLevel; ++i )
  {
    if( i > 0 )
    {
      std::cerr << ">>> Refining grid globally..." << std::endl;
      grid.globalRefine( 1 );
    }
    std::cerr << ">>> Checking grid..." << std::endl;
    gridcheck( grid );
    std::cerr << ">>> Checking intersections..." << std::endl;
    checkIntersectionIterator( grid );
    if( i > 0 )
    {
      std::cerr << ">>> Checking geometry in father..." << std::endl;
      checkGeometryInFather( grid );
    }
  }

  std::cerr << ">>> Writing out grid..." << std::endl;
  grid.writeGrid< Dune::ascii >( "gridcheck.spgrid", 0.0 );

  std::cerr << ">>> Reading back grid..." << std::endl;
  Grid rgrid;
  double time;
  rgrid.readGrid< Dune::ascii >( "gridcheck.spgrid", time );

  std::cerr << ">>> Checking grid..." << std::endl;
  gridcheck( rgrid );
  std::cerr << ">>> Checking intersections..." << std::endl;
  checkIntersectionIterator( rgrid );
  if( rgrid.maxLevel() > 0 )
  {
    std::cerr << ">>> Checking geometry in father..." << std::endl;
    checkGeometryInFather( rgrid );
  }

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
