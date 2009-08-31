#include <config.h>

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/common/mpihelper.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>

#include "checkseiterator.hh"

static const int dimGrid = DIMGRID;

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  typedef Dune::SPGrid< double, dimGrid > Grid;

  if( argc < 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf file>" << std::endl;
    return 1;
  }

  std::string dgfFile( argv[ 1 ] );

  Dune::GridPtr< Grid > gridPtr( dgfFile );
  Grid &grid = *gridPtr;

  Dune::checkSuperEntityIterator( grid.leafView() );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
