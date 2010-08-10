#include <config.h>

#include <iostream>

#include <dune/grid/spgrid.hh>

typedef Dune::SPGrid< double, 3 > Grid;

typedef Grid::GlobalVector GlobalVector;
typedef Grid::Domain Domain;
typedef Grid::MultiIndex MultiIndex;

int main ( int argc, char **argv )
try
{
  int cells[ 3 ] = { 4, 4, 4 };

  Domain domain( Domain::unitCube() );
  Grid grid( domain, MultiIndex( cells ) );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
