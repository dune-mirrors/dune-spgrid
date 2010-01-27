#include <limits>
#include <sstream>

#include <dune/common/iostream.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/decomposition.hh>

#ifndef DIMGRID
#error "DIMGRID not defined."
#endif

using namespace Dune;

static const int dimGrid = DIMGRID;

void listPartitions ( const SPDecomposition< dimGrid > &decomposition, const int overlap )
{
  const unsigned int size = decomposition.size();
  SPPartition< dimGrid > globalPartition( decomposition.mesh() );

  int maxload = std::numeric_limits< int >::min();
  int minload = std::numeric_limits< int >::max();
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    SPPartition< dimGrid > localPartition( decomposition.subMesh( rank ) );
    //localPartition = globalPartition.intersect( localPartition.grow( overlap ) );

    const int load = localPartition.volume();
    minload = std::min( minload, load );
    maxload = std::max( maxload, load );
    std::cout << "rank " << rank << ": " << localPartition;
    std::cout << " (load: " << load << ")" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "maximal load: " << maxload << ", minimal load: " << minload
            << ", ratio: " << (double( maxload ) / double( minload )) << std::endl;
}

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    //std::cerr << "Usage: " << argv[ 0 ] << " <width> <size> [periodic]" << std::endl;
    std::cerr << "Usage: " << argv[ 0 ] << " <width> <size>" << std::endl;
    return 1;
  }

  SPMultiIndex< dimGrid > width;
  std::istringstream widthStream( argv[ 1 ] );
  widthStream >> width;

  unsigned int size;
  std::istringstream sizeStream( argv[ 2 ] );
  sizeStream >> size;

#if 0
  unsigned int periodic = 0;
  for( int i = 3; i < argc; ++i )
  {
    unsigned int dir;
    std::istringstream periodicStream( argv[ i ] );
    periodicStream >> dir;
    periodic |= (1 << dir);
  }
#endif

  //std::cout << "width = " << width << ", size = " << size << ", periodic = " << periodic << std::endl;
  std::cout << "width = " << width << ", size = " << size << std::endl;
  SPDecomposition< dimGrid > decomposition( width, size );

  std::cout << std::endl;
  std::cout << "Interior Partitions:" << std::endl;
  std::cout << "--------------------" << std::endl;
  listPartitions( decomposition, 0 );

  std::cout << std::endl;
  std::cout << "Overlap Partitions (1 level of overlap):" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  listPartitions( decomposition, 1 );

  typedef SPGrid< double, dimGrid > Grid;
  FieldVector< double, dimGrid > a( 0.0 ), b( 1.0 );
  SPDomain< double, dimGrid > domain( a, b, width );
  Grid grid( domain );

  std::vector< double > data( grid.size( 0 ) );
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    const SPMesh< dimGrid > &mesh = decomposition.subMesh( rank );
    SPMultiIndex< dimGrid > pos = mesh.begin();
    SPMultiIndex< dimGrid > end = mesh.end();
    for( int d = 0; d < dimGrid; )
    {
      unsigned int index = 0;
      unsigned int factor = 1;
      for( int i = 0; i < dimGrid; ++i )
      {
        index += pos[ i ] * factor;
        factor *= width[ i ];
      }
      data[ index ] = double( rank );

      for( d = 0; d < dimGrid; ++d )
      {
        if( ++pos[ d ] < end[ d ] )
          break;
        pos[ d ] = mesh.begin()[ 0 ];
      }
    }
  }
  VTKWriter< Grid::LeafGridView > vtkWriter( grid.leafView() );
  vtkWriter.addCellData( data, "rank" );
  vtkWriter.write( "decomposition" );

  return 0;
}
