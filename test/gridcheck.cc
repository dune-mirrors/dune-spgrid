#include <config.h>

#define NEW_SUBENTITY_NUMBERING 1
#define DISABLE_DEPRECATED_METHOD_CHECK 1

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/common/mpihelper.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>

//#include <dune/grid/utility/hierarchicsearch.hh>
#include <dune/grid/utility/sphierarchicsearch.hh>

#include <dune/grid/test/gridcheck.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkiterators.cc>
#include <dune/grid/test/checkpartition.cc>
#include <dune/grid/test/checkcommunicate.cc>

#include "checkidcommunication.hh"

static const int dimGrid = DIMGRID;


template< class GridView >
void checkHierarchicSearch ( const GridView &gridView )
{
  typedef Dune::HierarchicSearch< typename GridView::Grid, typename GridView::IndexSet > HierarchicSearch;

  typedef typename GridView::template Codim< 0 >::EntityPointer EntityPointer;
  typedef typename EntityPointer::Entity::Geometry Geometry;
  typedef typename Geometry::ctype ctype;
  typedef typename Geometry::LocalCoordinate LocalVector;
  typedef typename Geometry::GlobalCoordinate GlobalVector;

  HierarchicSearch hsearch( gridView.grid(), gridView.indexSet() );

  GlobalVector x = gridView.grid().domain().origin();
  x.axpy( 0.49, gridView.grid().domain().width() );

  EntityPointer ep = hsearch.findEntity( x );
  const Geometry &geometry = ep->geometry();
  const LocalVector xl = geometry.local( x );
  if( !Dune::GenericReferenceElements< ctype, dimGrid >::cube().checkInside( xl ) )
  {
    std::cerr << "Cannot find " << x << " in entity returned by hierarchic search." << std::endl;
    std::cerr << "Entity:";
    for( int i = 0; i < geometry.corners(); ++i )
      std::cerr << "  ( " << geometry.corner( i ) << " )";
    std::cerr << std::endl;
  }
}


template< class Grid >
void performCheck ( Grid &grid, const int maxLevel )
{
  for( int i = 0; i <= maxLevel; ++i )
  {
    if( i > 0 )
    {
      std::cerr << ">>> Refining grid globally..." << std::endl;
      grid.globalRefine( 1 );
    }
    std::cerr << ">>> Checking grid..." << std::endl;
    gridcheck( grid );
    checkIterators( grid.leafView() );
    checkPartitionType( grid.leafView() );
    std::cerr << ">>> Checking intersections..." << std::endl;
    checkIntersectionIterator( grid );
    if( i > 0 )
    {
      std::cerr << ">>> Checking geometry in father..." << std::endl;
      checkGeometryInFather( grid );
    }

    std::cerr << ">>> Checking communication..." << std::endl;
    checkIdCommunication( grid.leafView() );
    checkCommunication( grid, -1, std::cout );

    std::cerr << ">>> Checking hierarchic search..." << std::endl;
    checkHierarchicSearch( grid.leafView() );
  }

  std::cerr << ">>> Writing out grid..." << std::endl;
  grid.template writeGrid< Dune::ascii >( "gridcheck.spgrid", 0.0 );

  std::cerr << ">>> Reading back grid..." << std::endl;
  Grid rgrid;
  double time;
  rgrid.template readGrid< Dune::ascii >( "gridcheck.spgrid", time );

  std::cerr << ">>> Checking grid..." << std::endl;
  gridcheck( rgrid );
  std::cerr << ">>> Checking intersections..." << std::endl;
  checkIntersectionIterator( rgrid );
  if( rgrid.maxLevel() > 0 )
  {
    std::cerr << ">>> Checking geometry in father..." << std::endl;
    checkGeometryInFather( rgrid );
  }
}


int main ( int argc, char **argv )
try
{
  Dune::MPIHelper::instance( argc, argv );


  if( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " <dgf file> <max level>" << std::endl;
    return 1;
  }

  std::string dgfFile( argv[ 1 ] );
  const int maxLevel = atoi( argv[ 2 ] );

  std::cout << "Isotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPIsotropicRefinement > > isoGrid( dgfFile );
  performCheck( *isoGrid, maxLevel );

  std::cout << std::endl;
  std::cout << "Anisotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPAnisotropicRefinement > > anisoGrid( dgfFile );
  performCheck( *anisoGrid, maxLevel );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
