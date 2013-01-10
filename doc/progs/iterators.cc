#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/spgrid.hh>

typedef Dune::SPGrid< double, 3 > Grid;

typedef Grid::GlobalVector GlobalVector;
typedef Grid::Domain Domain;
typedef Grid::MultiIndex MultiIndex;

typedef Grid::Codim< 0 >::LeafIterator Iterator;
typedef Grid::LeafGridView GridView;
typedef Grid::LeafIndexSet IndexSet;

typedef Dune::VTKWriter< GridView > VTKWriter;

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );

  int cells[ 3 ] = { 4, 4, 4 };

  Domain domain( Domain::unitCube() );
  Grid grid( domain, MultiIndex( cells ) );

  const IndexSet &indexSet = grid.leafIndexSet();
  std::vector< double > values( indexSet.size( 0 ), double( 0 ) );

  for( unsigned int dir = 0; dir < (1u << 3); ++dir )
  {
    const Iterator end = grid.leafend< 0 >( dir );
    unsigned int count = 0;
    for( Iterator it = grid.leafbegin< 0 >( dir ); it != end; ++it, ++count )
    {
      const IndexSet::IndexType index = indexSet.index( *it );
      values[ index ] = double( 1 );

      VTKWriter vtkWriter( grid.leafView(), Dune::VTK::nonconforming );
      vtkWriter.addCellData( values, "Currect Cell" );

      std::ostringstream namestream;
      namestream << "iterator-" << dir << "-" << count;
      std::string name = namestream.str();
      vtkWriter.pwrite( name.c_str(), "vtu", "." );

      values[ index ] = double( 0 );
    }
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
