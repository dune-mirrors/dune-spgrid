#include <config.h>

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/common/mpihelper.hh>

#include <dune/grid/spgrid.hh>

#include <dune/grid/test/gridcheck.cc>

static const int dimGrid = DIMGRID;

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  typedef Dune::SPGrid< double, dimGrid > Grid;

  Dune::FieldVector< double, dimGrid > a( 0.0 );
  Dune::FieldVector< double, dimGrid > b( 1.0 );
  int n[ dimGrid ];
  for( int i = 0; i < dimGrid; ++i )
    n[ i ] = 4;
  Grid grid( a, b, n );

  const int maxLevel = 2;
  std::cerr << "Checking grid..." << std::endl;
  for( int i = 0; i <= maxLevel; ++i )
  {
    if( i > 0 )
      grid.globalRefine( 1 );
    gridcheck( grid );
  }

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
