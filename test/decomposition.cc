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
    std::istringstream periodicStream( argv[ i ] );
    periodicStream >> dir;
    periodic |= (1 << dir);
  }

  std::cout << "width = " << width << ", size = " << size << ", periodic = " << periodic << std::endl;
  SPDecomposition< dimGrid > decomposition( width, size, periodic );

  std::cout << std::endl;
  std::cout << "Interior Partitions:" << std::endl;
  std::cout << "--------------------" << std::endl;
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    SPPartition< dimGrid > partition = decomposition.partition( rank, 0 );
    std::cout << "rank " << rank << ": " << partition << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Overlap Partitions (1 level of overlap):" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    SPPartition< dimGrid > partition = decomposition.partition( rank, 1 );
    std::cout << "rank " << rank << ": " << partition << std::endl;
  }

  return 0;
}
