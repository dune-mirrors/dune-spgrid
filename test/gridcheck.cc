#include <config.h>

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/grid/spgrid.hh>

#include <dune/grid/test/gridcheck.cc>

static const int dimGrid = DIMGRID;

int main ( int argc, char **argv )
try
{
  typedef Dune::SPGrid< double, dimGrid > Grid;

  Dune::FieldVector< double, dimGrid > a( 0.0 );
  Dune::FieldVector< double, dimGrid > b( 1.0 );
  int n[ dimGrid ];
  for( int i = 0; i < dimGrid; ++i )
    n[ i ] = 4;
  Grid grid( a, b, n );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
