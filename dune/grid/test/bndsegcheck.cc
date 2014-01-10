#include <config.h>

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

//- C++ includes
#include <cassert>
#include <string>

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

//- dune-spgrid includes
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>

//- local includes
#include "checkbndsegiterator.hh"


static const int dimGrid = DIMGRID;


template< class Grid >
void performCheck ( Grid &grid, const int maxLevel )
{
  std::cerr << ">>> Refining grid globally..." << std::endl;
  for( int level = 0; level < maxLevel; ++level )
    grid.globalRefine( 1 );

  assert( grid.maxLevel() == maxLevel );

  std::cerr << ">>> Checking boundary segment iterator..." << std::endl;
  for( int level = 0; level <= maxLevel; ++level )
    checkBoundarySegmentIterator( grid.levelGridView( level ) );
}


int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf file> <max level>" << std::endl;
    return 1;
  }

  std::string dgfFile( argv[ 1 ] );
  const int maxLevel = atoi( argv[ 2 ] );

  std::cout << "Isotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPIsotropicRefinement > > isoGrid( dgfFile );
  performCheck( *isoGrid, maxLevel );
  std::cout << std::endl;

  std::cout << "Anisotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPAnisotropicRefinement > > anisoGrid( dgfFile );
  performCheck( *anisoGrid, maxLevel );
  std::cout << std::endl;

  std::cout << "Bisection grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPBisectionRefinement > > bisectionGrid( dgfFile );
  performCheck( *bisectionGrid, dimGrid*maxLevel );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
