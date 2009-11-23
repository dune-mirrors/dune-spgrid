#include <sstream>

#include <dune/common/iostream.hh>
#include <dune/grid/spgrid/decomposition.hh>

#ifndef DIMGRID
#error "DIMGRID not defined."
#endif

using namespace Dune;

static const int dimGrid = DIMGRID;

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <width> <size> [periodic]" << std::endl;
    return 1;
  }

  SPMultiIndex< dimGrid > width;
  std::istringstream widthStream( argv[ 1 ] );
  widthStream >> width;

  unsigned int size;
  std::istringstream sizeStream( argv[ 2 ] );
  sizeStream >> size;

  unsigned int periodic = 0;
  for( int i = 3; i < argc; ++i )
  {
    unsigned int dir;
    std::istringstream periodicStream;
    periodicStream >> dir;
    periodic |= (1 << dir);
  }

  std::cout << "width =  " << width << ", size = " << size << ", periodic = " << periodic << std::endl;
}
