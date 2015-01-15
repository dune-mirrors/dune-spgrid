#include <config.h>

#define NEW_SUBENTITY_NUMBERING 1

#ifndef DIMGRID
#error "DIMGRID not defined. Please compile with -DDIMGRID=n"
#endif

#include <dune/common/forloop.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checkintersectionit.hh>
#include <dune/grid/test/checkgeometryinfather.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkpartition.hh>
#include <dune/grid/test/checkcommunicate.hh>

#include "checkidcommunication.hh"
#include "checkseiterator.hh"

static const int dimGrid = DIMGRID;


template< int codim >
struct CheckSubIndex
{
  template< class GridView >
  static void apply ( const GridView &gridView )
  {
    std::cerr << ">>> Checking subIndex() for codim " << codim << "..." << std::endl;
    typedef typename GridView::ctype ctype;

    static const int dimension = GridView::dimension;
    static const int mydimension = dimension - codim;

    const typename GridView::IndexSet &indexSet = gridView.indexSet();

    typedef typename GridView::template Codim< codim >::Iterator Iterator;
    const Iterator end = gridView.template end< codim >();
    for( Iterator it = gridView.template begin< codim >(); it != end; ++it )
    {
      const typename Iterator::Entity &entity = *it;
      const Dune::ReferenceElement< ctype, mydimension > &referenceElement
        = Dune::ReferenceElements< ctype, mydimension >::general( entity.type() );
      for( int i = 0; i < referenceElement.size( mydimension ); ++i )
        indexSet.subIndex( entity, i, dimension );
    }
  }
};

template< class GridView >
void checkSubIndex ( const GridView &gridView )
{
  Dune::ForLoop< CheckSubIndex, 0, GridView::dimension >::apply( gridView );
}


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

  GlobalVector x = gridView.grid().domain().cube().origin();
  x.axpy( 0.49, gridView.grid().domain().cube().width() );

  EntityPointer ep = hsearch.findEntity( x );
  const Geometry &geometry = ep->geometry();
  const LocalVector xl = geometry.local( x );
  if( !Dune::ReferenceElements< ctype, dimGrid >::cube().checkInside( xl ) )
  {
    std::cerr << "Cannot find " << x << " in entity returned by hierarchic search." << std::endl;
    std::cerr << "Entity:";
    for( int i = 0; i < geometry.corners(); ++i )
      std::cerr << "  ( " << geometry.corner( i ) << " )";
    std::cerr << std::endl;
  }
}


template< class Grid >
void performCheck ( Grid &grid, int maxLevel, const typename Grid::RefinementPolicy &policy = typename Grid::RefinementPolicy() )
{
  for( int i = 0; i <= maxLevel; ++i )
  {
    if( i > 0 )
    {
      std::cerr << ">>> Refining grid globally..." << std::endl;
      grid.globalRefine( 1, policy );
    }
    std::cerr << ">>> Checking grid..." << std::endl;
    gridcheck( grid );
    checkIterators( grid.leafGridView() );
    checkPartitionType( grid.leafGridView() );
    std::cerr << ">>> Checking intersections..." << std::endl;
    checkIntersectionIterator( grid );
    if( i > 0 )
    {
      std::cerr << ">>> Checking geometry in father..." << std::endl;
      checkGeometryInFather( grid );
    }

    std::cerr << ">>> Checking communication..." << std::endl;
    checkIdCommunication( grid.leafGridView() );
    checkCommunication( grid, -1, std::cout );

    checkSubIndex( grid.leafGridView() );

    if( grid.comm().size() <= 1 )
    {
      checkSuperEntityIterator( grid.leafGridView() );

      std::cerr << ">>> Checking hierarchic search..." << std::endl;
      checkHierarchicSearch( grid.leafGridView() );
    }
    else
      std::cerr << "WARNING: SuperEntityIterators currently don't work in parallel; test disabled." << std::endl;
  }

  std::ostringstream sFilename;
  sFilename << "gridcheck." << Grid::Refinement::type() << ".spgrid";
  const std::string filename = sFilename.str();

  std::cerr << ">>> Writing out grid..." << std::endl;
  Dune::BackupRestoreFacility< Grid >::backup( grid, filename );

  std::cerr << ">>> Reading back grid..." << std::endl;
  Grid *rgrid = Dune::BackupRestoreFacility< Grid >::restore( filename );

  if( rgrid )
  {
    std::cerr << ">>> Checking grid..." << std::endl;
    gridcheck( *rgrid );
    std::cerr << ">>> Checking intersections..." << std::endl;
    checkIntersectionIterator( *rgrid );
    if( rgrid->maxLevel() > 0 )
    {
      std::cerr << ">>> Checking geometry in father..." << std::endl;
      checkGeometryInFather( *rgrid );
    }

    delete rgrid;
  }
  else
    std::cerr << "Could not read back grid." << std::endl;
}


int main ( int argc, char **argv )
try
{
  const Dune::MPIHelper &mpi = Dune::MPIHelper::instance( argc, argv );

  if( (argc > 1) && (std::string( argv[ 1 ] ) == std::string( "--help" )) )
  {
    if( mpi.rank() == 0 )
      std::cerr << "Usage: " << argv[ 0 ] << " <dgf file> <max level>" << std::endl;
    return 0;
  }

  std::string dgfFile( argc > 1 ? argv[ 1 ] : std::to_string( dimGrid ) + "dcube.dgf" );
  const int maxLevel = (argc > 2 ? atoi( argv[ 2 ] ) : 1);

  std::cout << "Isotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPIsotropicRefinement > > isoGrid( dgfFile );
  performCheck( *isoGrid, maxLevel );

  std::cout << std::endl;
  std::cout << "Anisotropic grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPAnisotropicRefinement > > anisoGrid( dgfFile );
  performCheck( *anisoGrid, maxLevel );

  std::cout << std::endl;
  std::cout << "Bisection grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPBisectionRefinement > > bisectionGrid( dgfFile );
  performCheck( *bisectionGrid, dimGrid*maxLevel );

  std::cout << std::endl;
  std::cout << "Arbitrary grid" << std::endl;
  Dune::GridPtr< Dune::SPGrid< double, dimGrid, Dune::SPArbitraryRefinement > > arbitraryGrid( dgfFile );
  performCheck( *arbitraryGrid, maxLevel, Dune::SPArbitraryRefinementPolicy< dimGrid >( 3 ) );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
