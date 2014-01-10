#include <config.h>

#include <limits>
#include <sstream>

#include <dune/common/parallel/mpihelper.hh>
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
typedef SPTopology< dimGrid > Topology;
typedef SPPartitionPool< dimGrid > PartitionPool;
typedef SPPartitionList< dimGrid > PartitionList;


template< PartitionIteratorType pitype >
void listPartitions ( const SPDecomposition< dimGrid > &decomposition, const MultiIndex &overlap, const Topology &topology )
{
  const unsigned int size = decomposition.size();

  int maxload = std::numeric_limits< int >::min();
  int minload = std::numeric_limits< int >::max();
  for( unsigned int rank = 0; rank < size; ++rank )
  {
    PartitionPool partitionPool( decomposition.subMesh( rank ), decomposition.mesh(), overlap, topology );
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

#if 0
    for( typename PartitionList::Iterator it = partition.begin(); it; ++it )
    {
      std::cout << "        - " << it->number();
      separator = ':';
      for( int i = 0; i < 2*dimGrid; ++i )
      {
        std::cout << separator << " " << int( it->neighbor( i ) );
        separator = ',';
      }
      std::cout << std::endl;
    }
#endif

    minload = std::min( minload, load );
    maxload = std::max( maxload, load );
  }
  std::cout << std::endl;
  std::cout << "maximal load: " << maxload << ", minimal load: " << minload
            << ", ratio: " << (double( maxload ) / double( minload )) << std::endl;
}



template< InterfaceType iftype >
void listLinkage ( const SPDecomposition< dimGrid > &decomposition, const MultiIndex &overlap, const Topology &topology )
{
  typedef SPLinkage< dimGrid > Linkage;
  typedef typename Linkage::Interface Interface;
  const int size = decomposition.size();
  for( int rank = 0; rank < size; ++rank )
  {
    PartitionPool localPool( decomposition.subMesh( rank ), decomposition.mesh(), overlap, topology );
    Linkage linkage( rank, localPool, decomposition.subMeshes() );

    std::cout << "rank " << rank << ":" << std::endl;
    const Interface &interface = linkage.interface( iftype );

    const typename Interface::Iterator end = interface.end();
    for( typename Interface::Iterator it = interface.begin(); it != end; ++it )
    {
      const PartitionList &sendList = it->sendList();
      if( !sendList.empty() )
      {
        std::cout << "- snd " << it->rank();
        char separator = ':';
        for( typename PartitionList::Iterator pit = sendList.begin(); pit; ++pit )
        {
          std::cout << separator << " " << *pit;
          separator = ';';
        }
        std::cout << std::endl;
      }

      const PartitionList &receiveList = it->receiveList();
      if( !receiveList.empty() )
      {
        std::cout << "- rcv " << it->rank();
        char separator = ':';
        for( typename PartitionList::Iterator pit = receiveList.begin(); pit; ++pit )
        {
          std::cout << separator << " " << *pit;
          separator = ';';
        }
        std::cout << std::endl;
      }
    }
  }
}


int main ( int argc, char **argv )
{
  Dune::MPIHelper::instance( argc, argv );

  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <width> <size> [overlap] [periodic]" << std::endl;
    return 1;
  }

  MultiIndex width;
  std::string widthArg( argv[ 1 ] );
  std::istringstream widthStream( widthArg );
  widthStream >> width;

  unsigned int size;
  std::string sizeArg( argv[ 2 ] );
  std::istringstream sizeStream( sizeArg );
  sizeStream >> size;

  MultiIndex overlap = MultiIndex::zero();
  if( argc > 3 )
  {
    std::string overlapArg( argv[ 3 ] );
    std::istringstream overlapStream( overlapArg );
    overlapStream >> overlap;
  }

  unsigned int periodic = 0;
  for( int i = 4; i < argc; ++i )
  {
    unsigned int dir;
    std::string periodicArg( argv[ i ] );
    std::istringstream periodicStream( periodicArg );
    periodicStream >> dir;
    periodic |= (1 << dir);
  }
  Topology topology( periodic );

  std::cout << "width = " << width << ", size = " << size << ", periodic = " << periodic << std::endl;
  SPDecomposition< dimGrid > decomposition( width, size );

  std::cout << std::endl;
  std::cout << "Interior Partitions:" << std::endl;
  std::cout << "--------------------" << std::endl;
  listPartitions< Interior_Partition >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "InteriorBorder Partitions:" << std::endl;
  std::cout << "--------------------------" << std::endl;
  listPartitions< InteriorBorder_Partition >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "Overlap Partitions:" << std::endl;
  std::cout << "-------------------" << std::endl;
  listPartitions< Overlap_Partition >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "OverlapFront Partitions:" << std::endl;
  std::cout << "------------------------" << std::endl;
  listPartitions< OverlapFront_Partition >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "All Partitions:" << std::endl;
  std::cout << "---------------" << std::endl;
  listPartitions< All_Partition >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "InteriorBorder InteriorBorder Communication:" << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
  listLinkage< InteriorBorder_InteriorBorder_Interface >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "InteriorBorder All Communication:" << std::endl;
  std::cout << "---------------------------------" << std::endl;
  listLinkage< InteriorBorder_All_Interface >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "Overlap OverlapFront Communication:" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
  listLinkage< Overlap_OverlapFront_Interface >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "Overlap All Communication:" << std::endl;
  std::cout << "--------------------------" << std::endl;
  listLinkage< Overlap_All_Interface >( decomposition, overlap, topology );

  std::cout << std::endl;
  std::cout << "All All Communication:" << std::endl;
  std::cout << "----------------------" << std::endl;
  listLinkage< All_All_Interface >( decomposition, overlap, topology );

  typedef SPGrid< double, dimGrid > Grid;
  FieldVector< double, dimGrid > a( 0.0 ), b( 1.0 );
  SPDomain< double, dimGrid > domain( a, b );
  Grid grid( domain, width );

  std::cout << "grid created." << std::endl;

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
  VTKWriter< Grid::LeafGridView > vtkWriter( grid.leafGridView() );
  vtkWriter.addCellData( data, "rank" );
  vtkWriter.write( "decomposition" );

  return 0;
}
