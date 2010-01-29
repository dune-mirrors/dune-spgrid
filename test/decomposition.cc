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

typedef SPMultiIndex< dimGrid > MultiIndex;
typedef SPPartitionPool< dimGrid > PartitionPool;
typedef SPPartitionList< dimGrid > PartitionList;

template< PartitionIteratorType pitype >
void listPartitions ( const SPDecomposition< dimGrid > &decomposition, const MultiIndex &overlap, unsigned int periodic )
{
  const unsigned int size = decomposition.size();

  int maxload = std::numeric_limits< int >::min();
  int minload = std::numeric_limits< int >::max();
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    PartitionPool partitionPool( decomposition.subMesh( rank ), decomposition.mesh(), overlap, periodic );
    const PartitionList &partition = partitionPool.template get< pitype >();

    int load = 0;
    std::cout << "rank " << rank;
    char separator = ':';
    for( typename PartitionList::Iterator it = partition.begin(); it; ++it )
    {
      load += it->volume();
      std::cout << separator << " " << *it;
      separator = ';';
    }
    std::cout << " (load: " << load << ")" << std::endl;

    minload = std::min( minload, load );
    maxload = std::max( maxload, load );
  }
  std::cout << std::endl;
  std::cout << "maximal load: " << maxload << ", minimal load: " << minload
            << ", ratio: " << (double( maxload ) / double( minload )) << std::endl;
}

int main ( int argc, char **argv )
{
  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <width> <size> [overlap] [periodic]" << std::endl;
    return 1;
  }

  MultiIndex width;
  std::istringstream widthStream( argv[ 1 ] );
  widthStream >> width;

  unsigned int size;
  std::istringstream sizeStream( argv[ 2 ] );
  sizeStream >> size;

  MultiIndex overlap = MultiIndex::zero();
  if( argc > 3 )
  {
    std::istringstream overlapStream( argv[ 3 ] );
    overlapStream >> overlap;
  }

  unsigned int periodic = 0;
  for( int i = 4; i < argc; ++i )
  {
    unsigned int dir;
    std::istringstream periodicStream( argv[ i ] );
    periodicStream >> dir;
    periodic |= (1 << dir);
  }

  std::cout << "width = " << width << ", size = " << size << ", periodic = " << periodic << std::endl;
  SPDecomposition< dimGrid > decomposition( width, size );

  std::cout << std::endl;
  std::cout << "Interior Partitions:" << std::endl;
  std::cout << "--------------------" << std::endl;
  listPartitions< Interior_Partition >( decomposition, overlap, periodic );

  std::cout << std::endl;
  std::cout << "InteriorBorder Partitions:" << std::endl;
  std::cout << "--------------------------" << std::endl;
  listPartitions< InteriorBorder_Partition >( decomposition, overlap, periodic );

  std::cout << std::endl;
  std::cout << "Overlap Partitions:" << std::endl;
  std::cout << "-------------------" << std::endl;
  listPartitions< Overlap_Partition >( decomposition, overlap, periodic );

  std::cout << std::endl;
  std::cout << "OverlapFront Partitions:" << std::endl;
  std::cout << "------------------------" << std::endl;
  listPartitions< OverlapFront_Partition >( decomposition, overlap, periodic );

  std::cout << std::endl;
  std::cout << "All Partitions:" << std::endl;
  std::cout << "---------------" << std::endl;
  listPartitions< All_Partition >( decomposition, overlap, periodic );
  
  typedef SPGrid< double, dimGrid > Grid;
  FieldVector< double, dimGrid > a( 0.0 ), b( 1.0 );
  SPDomain< double, dimGrid > domain( a, b, width );
  Grid grid( domain );

  std::vector< double > data( grid.size( 0 ) );
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    const SPMesh< dimGrid > &mesh = decomposition.subMesh( rank );
    MultiIndex pos = mesh.begin();
    MultiIndex end = mesh.end();
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
