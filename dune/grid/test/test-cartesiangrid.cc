#include <config.h>

#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/mpihelper.hh>

#include <dune/grid/cartesiangrid/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/test/gridcheck.cc>

#include <dune/grid/test/checkgeometryinfather.cc>
#include <dune/grid/test/checkintersectionit.cc>
#include <dune/grid/test/checkcommunicate.cc>

#include <dune/grid/io/visual/grapegriddisplay.hh>

using namespace Dune;

template< class GridType >
void makeNonConfGrid ( GridType &grid, int level, int adapt )
{
  int myrank = grid.comm().rank();
  grid.loadBalance();
  grid.globalRefine(level);
  grid.loadBalance();
  for (int i=0;i<adapt;i++) 
  {
    if (myrank==0) 
    {
      typedef typename GridType :: template Codim<0> :: 
            template Partition<Interior_Partition> :: LeafIterator LeafIterator;
      
      LeafIterator endit = grid.template leafend<0,Interior_Partition>   ();
      int nr = 0;
      int size = grid.size(0);
      for(LeafIterator it    = grid.template leafbegin<0,Interior_Partition> ();
          it != endit ; ++it,nr++ ) 
      {
        grid.mark(1, *it );
        if (nr>size*0.8) break;
      }
    }
    grid.adapt();
    grid.postAdapt();
    grid.loadBalance();
  }
}

template< class GridView >
void writeFile ( const GridView &gridView )
{
  DGFWriter< GridView > writer( gridView );
  writer.write( "dump.dgf" );
}

template< class GridType >
void checkSerial ( GridType &grid, int mxl = 2, const bool display = false )
{
  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

  // be careful, each global refine create 8 x maxlevel elements 
  std::cout << "  CHECKING: Macro" << std::endl;
  gridcheck(grid);
  std::cout << "  CHECKING: Macro-intersections" << std::endl;
  checkIntersectionIterator(grid);

  for(int i=0; i<mxl; i++) 
  {
    grid.globalRefine( 1 );//DGFGridInfo<GridType> :: refineStepsForHalf() );
    std::cout << "  CHECKING: Refined" << std::endl;
    gridcheck(grid);
    std::cout << "  CHECKING: intersections" << std::endl;
    checkIntersectionIterator(grid);
      
    if( display )
    {
      GrapeGridDisplay< GridType > grape( grid );
      grape.display();
    }
  }

  // check also non-conform grids 
  makeNonConfGrid(grid,0,1);

  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

  std::cout << "  CHECKING: non-conform" << std::endl;
  gridcheck(grid);
  std::cout << "  CHECKING: twists " << std::endl;
  // if( checkTwist )
  //  checkTwists( grid.leafView(), NoMapTwist() );
  
  // check the method geometryInFather()
  std::cout << "  CHECKING: geometry in father" << std::endl;
  checkGeometryInFather(grid);
  // check the intersection iterator and the geometries it returns
  std::cout << "  CHECKING: intersections" << std::endl;
  checkIntersectionIterator(grid);

  std::cout << std::endl << std::endl;
}

template< class GridType >
void checkParallel ( GridType &grid, int gref, int mxl = 3, const bool display = false )
{
  if( display )
  {
    GrapeGridDisplay< GridType > grape( grid );
    grape.display();
  }

#if HAVE_MPI
  makeNonConfGrid(grid,gref,mxl);

  // -1 stands for leaf check 
  checkCommunication(grid, -1, std::cout);

  for(int l=0; l<= mxl; ++l)
    checkCommunication(grid, l , Dune::dvverb);
#endif
}


int main ( int argc , char **argv )
{
  
  // this method calls MPI_Init, if MPI is enabled
  MPIHelper & mpihelper = MPIHelper::instance(argc,argv);
  int myrank = mpihelper.rank(); 
  int mysize = mpihelper.size();

  try {
    /* use grid-file appropriate for dimensions */

    if( argc < 2 ) 
    {
      std::cerr << "Usage: " << argv[0] << " <dgf file of hostgrid>" << std::endl;
      return 1;
    }

    const bool display = (argc > 2);

    typedef Dune::GridSelector :: GridType HostGridType;
    typedef CartesianGrid< HostGridType > CartesianGridType;
    std::string filename( argv[1] );
    GridPtr< CartesianGridType > gridPtr( filename );
    CartesianGridType &grid = *gridPtr;

    //grid.loadBalance();

    {
      std::cout << "Check serial grid" << std::endl;
      checkSerial(grid,
                     (mysize == 1) ? 1 : 0, 
                     (mysize == 1) ? display: false);
    }
    
    // perform parallel check only when more then one proc 
    if(mysize > 1)
    {
      if (myrank == 0) std::cout << "Check conform grid" << std::endl;
      checkParallel(grid,1,0, display );
      if (myrank == 0) std::cout << "Check non-conform grid" << std::endl;
      checkParallel(grid,0,2, display );
    }

  } 
  catch (const Dune::Exception &e) 
  {
    std::cerr << e << std::endl;
    return 1;
  } 
  catch (...) 
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
