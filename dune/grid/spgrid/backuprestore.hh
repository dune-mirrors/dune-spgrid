#ifndef DUNE_SPGRID_BACKUPRESTORE_HH
#define DUNE_SPGRID_BACKUPRESTORE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#include <dune/grid/common/backuprestore.hh>

#include <dune/grid/spgrid/communication.hh>
#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/fileio.hh>

namespace Dune
{

  /** \class BackupRestoreFacility
   *  \brief facility for writing and reading a \ref Dune::SPGrid "SPGrid"
   *
   *  The BackupRestoreFacility allows writing hierarchic grids to disk and
   *  reading them back into another program.
   *
   *  It is guaranteed that all index sets and id sets are preserved by the
   *  backup / restore process.
   *  The result of restore is undefined, if the number of processes in a
   *  parallel program differs from the number of processes used on backup.
   *
   *  There are two pairs of backup / restore methods:
   *  - methods writing into one or more dedicated files,
   *  - methods operating on a std::stream.
   *  .
   *  These techniques may not be mixed, i.e., you cannot write the grid into
   *  files and read it back from a stream or vice versa.
   *
   *  \tparam  Grid  type of grid
   */
  template< class ct, int dim, template< int > class Ref, class Comm >
  struct BackupRestoreFacility< SPGrid< ct, dim, Ref, Comm > >
  {
    typedef SPGrid< ct, dim, Ref, Comm > Grid;

    typedef typename Grid::Communication Communication;

    /** \brief write a hierarchic grid to disk
     *
     *  \param[in]  grid        grid to write
     *  \param[in]  path        path to write the file to
     *  \param[in]  fileprefix  prefix of the file name
     */
    static void backup ( const Grid &grid, const std::string &path, const std::string &fileprefix )
    {
      backup( path + "/" + fileprefix + ".spgrid" );
    }

    /** \brief write a hierarchic grid to disk
     *
     *  \param[in]  grid        grid to write
     *  \param[in]  path        filename of the file to create
     */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      SPGridIOData< ct, dim, Ref > ioData;
      backup( grid, ioData );

      int result = 0;
      if( grid.comm().rank() == 0 )
      {
        std::ofstream stream( filename.c_str() );
        if( stream )
          result = int( ioData.write( stream ) );
      }
      grid.comm().broadcast( &result, 1, 0 );
      if( !result )
        DUNE_THROW( IOError, "Unable to write SPGrid to file '" + filename + "'." );
    }

    /** \brief write a hierarchic grid into a stream
     *
     *  \param[in]  grid        grid to write
     *  \param[in]  stream      std::stream to write the grid to
     */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      SPGridIOData< ct, dim, Ref > ioData;
      backup( grid, ioData );

      int result = 0;
      if( grid.comm().rank() == 0 )
        result = int( ioData.write( stream ) );
      grid.comm().broadcast( &result, 1, 0 );
      if( !result )
        DUNE_THROW( IOError, "Unable to write SPGrid to stream." );
    }

    /** \brief read a hierarchic grid from disk
     *
     *  \param[in]  path        path to write the file to
     *  \param[in]  fileprefix  prefix of the file name
     *  \param[in]  comm        collective communication (optional)
     *
     *  \returns a pointer to the grid (allocated by new)
     */
    static Grid *restore ( const std::string &path, const std::string &fileprefix,
                           const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() )
    {
      return restore( path + "/" + fileprefix + ".spgrid" );
    }

    /** \brief read a hierarchic grid from disk
     *
     *  \param[in]  filename    name of the file to read
     *  \param[in]  comm        collective communication (optional)
     *
     *  \returns a pointer to the grid (allocated by new)
     */
    static Grid *restore ( const std::string &filename,
                           const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() )
    {
      std::ifstream stream( filename.c_str() );
      if( !stream )
        DUNE_THROW( IOError, "Unable to open file: " + filename );

      Grid *grid = 0;
      SPGridIOData< ct, dim, Ref > ioData;
      if( ioData.read( stream, filename ) )
        grid = restore( ioData, comm );

      if( !parallelAnd( comm, grid ) )
      {
        delete grid;
        DUNE_THROW( IOError, "Unable to read grid from file '" + filename + "'." );
      }
      return grid;
    }

    /** \brief read a hierarchic grid from a stream
     *
     *  \param[in]  stream      std::stream to read the grid from
     *  \param[in]  comm        collective communication (optional)
     */
    static Grid *restore ( std::istream &stream,
                           const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() )
    {
      Grid *grid = 0;
      SPGridIOData< ct, dim, Ref > ioData;
      if( ioData.read( stream ) )
        grid = restore( ioData, comm );

      if( !parallelAnd( comm, grid ) )
      {
        delete grid;
        DUNE_THROW( IOError, "Unable to read grid from stream." );
      }
      return grid;
    }

  private:
    static void backup ( const Grid &grid, SPGridIOData< ct, dim, Ref > &ioData )
    {
      ioData.time = 0;
      ioData.cubes.push_back( grid.domain().cube() );
      ioData.topology = grid.domain().topology();
      ioData.cells = grid.globalMesh_.width();
      ioData.partitions = grid.comm().size();
      ioData.overlap = grid.overlap_;
      ioData.maxLevel = grid.maxLevel();
      ioData.refinements.resize( ioData.maxLevel );
      for( int level = 0; level < ioData.maxLevel; ++level )
        ioData.refinements[ level ] = grid.gridLevel( level+1 ).refinement().policy();
    }

    static Grid *restore ( const SPGridIOData< ct, dim, Ref > &ioData,
                           const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() )
    {
      if( ioData.partitions != comm.size() )
      {
        std::cerr << "Warning: Reading grid with different number of partitions,"
                  << " index sets will not coincide." << std::endl;
      }

      typename Grid::Domain domain( ioData.cubes, ioData.topology );
      Grid *grid = new Grid( domain, ioData.cells, ioData.overlap, comm );

      for( int level = 0; level < ioData.maxLevel; ++level )
      {
        if( level < int( ioData.refinements.size() ) )
          grid->globalRefine( 1, ioData.refinements[ level ] );
        else
          grid->globalRefine( 1 );
      }
      return grid;
    }

    static bool parallelAnd ( const Communication &comm, bool condition )
    {
      int result( condition );
      result = comm.sum( result );
      return (result == comm.size());
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_BACKUPRESTORE_HH
