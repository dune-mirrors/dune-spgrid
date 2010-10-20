#ifndef DUNE_SPGRID_FILEIO_HH
#define DUNE_SPGRID_FILEIO_HH

#include <fstream>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/topology.hh>
#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>
#include <dune/grid/spgrid/version.hh>

namespace Dune
{

  // SPGridIOData
  // ------------

  template< class ctype, int dim, SPRefinementStrategy strategy >
  struct SPGridIOData
  {
    typedef SPTopology< dim > Topology;
    typedef SPCube< ctype, dim > Cube;
    typedef typename Cube::GlobalVector GlobalVector;
    typedef SPMultiIndex< dim > MultiIndex;
    typedef SPRefinement< dim, strategy > Refinement;
    typedef typename Refinement::Policy RefinementPolicy;

    std::string name;
    ctype time;
    std::vector< Cube > cubes;
    Topology topology;
    MultiIndex cells;
    MultiIndex overlap;
    int partitions;
    int maxLevel;
    std::vector< RefinementPolicy > refinements;

    bool write ( const std::string &filename ) const;
    bool read ( const std::string &filename );

  private:
    static std::string readLine ( std::istream &in, unsigned int *count = 0 );
  };



  // Implementation of SPGridIOData
  // ------------------------------

  template< class ctype, int dim, SPRefinementStrategy strategy >
  inline bool
  SPGridIOData< ctype, dim, strategy >::write ( const std::string &filename ) const
  {
    std::ofstream fileOut( filename.c_str() );
    if( !fileOut )
      return false;

    // write header
    fileOut << "SPGrid";
    fileOut << "  dimension=" << dim;
    const unsigned int versionMajor = SPGridVersion::major;
    const unsigned int versionMinor = SPGridVersion::minor;
    fileOut << "  version=" << versionMajor << "." << versionMinor;
    fileOut << std::endl << std::endl;

    // write generic information
    fileOut << "name " << name << std::endl;
    fileOut << "time " << time << std::endl;
    fileOut << std::endl;

    // write domain information
    fileOut << "domain";
    for( typename std::vector< Cube >::const_iterator it = cubes.begin(); it != cubes.end(); ++it )
      fileOut << " " << *it;
    fileOut << std::endl;

    fileOut << "periodic";
    for( int i = 0; i < dim; ++i )
    {
      if( topology.periodic( i ) )
        fileOut << " " << i;
    }
    fileOut << std::endl << std::endl;

    // write discretization information
    fileOut << "cells " << cells << std::endl;
    fileOut << "partitions " << partitions << std::endl;
    fileOut << "overlap " << overlap << std::endl;
    fileOut << std::endl;

    // write refinement information
    fileOut << "maxLevel " << maxLevel << std::endl;
    fileOut << "refinement " << Refinement::type() << std::endl;
    fileOut << "refinements";
    for( unsigned int i = 0; i < refinements.size(); ++i )
      fileOut << " " << refinements[ i ];
    fileOut << std::endl;

    fileOut.close();
    return true;
  }


  template< class ctype, int dim, SPRefinementStrategy strategy >
  inline bool
  SPGridIOData< ctype, dim, strategy >::read ( const std::string &filename )
  {
    std::ifstream fileIn( filename.c_str() );
    if( !fileIn )
      return false;

    unsigned int lineNr = 0;
    std::string line = readLine( fileIn, &lineNr );
    std::istringstream lineIn( line );
    lineIn >> match( std::string( "SPGrid" ) );
    if( lineIn.fail() )
    {
      std::cerr << filename << "[ " << lineNr << " ]: 'SPGrid' expected." << std::endl;
      return false;
    }

    int fdim = -1;

    while( isGood( lineIn ) )
    {
      std::string tag;
      lineIn >> tag;

      const size_t eq = tag.find( '=' );
      const std::string key = tag.substr( 0, eq );
      const std::string value = (eq+1 < tag.size() ? tag.substr( eq+1 ) : std::string());
      std::istringstream valueIn( value );
      if( key == "version" )
      {
        // ensure that the check passes on read failure
        unsigned int vMajor = SPGridVersion::major, vMinor = SPGridVersion::minor;
        valueIn >> vMajor >> match( '.' ) >> vMinor;
        if( SPGridVersion::later( vMajor, vMinor ) )
        {
          std::cerr << filename << "[ " << lineNr << " ]: File was created by newer version of SPGrid." << std::endl;
          return false;
        }
      }
      else if( key == "dimension" )
        valueIn >> fdim;
      else
      {
        std::cerr << filename << "[ " << lineNr << " ]: Invalid tag: '" << key << "'." << std::endl;
        return false;
      }
      if( !valueIn )
      {
        std::cerr << filename << "[ " << lineNr << " ]: Invalid value for tag '" << key << "'." << std::endl;
        return false;
      }
    }

    if( fdim != dim )
    {
      std::cerr << filename << "[ " << lineNr << " ]: File has wrong grid dimension." << std::endl;
      return false;
    }

    name = "SPGrid";
    partitions = 1;
    overlap = MultiIndex::zero();
    time = ctype( 0 );
    cubes.clear();

    const unsigned int flagDomain = 1;
    const unsigned int flagCells = 2;
    const unsigned int flagMaxLevel = 4;
    const unsigned int flagRefinement = 8;
    const unsigned int flagAll = 15;
    unsigned int flags = 0;

    while( true )
    {
      std::string line = readLine( fileIn, &lineNr );
      if( line.empty() )
        break;
      std::istringstream lineIn( line );

      std::string cmd;
      lineIn >> cmd;

      if( cmd == "name" )
        name = readLine( lineIn );
      else if( cmd == "time" )
        lineIn >> time;
      else if ( cmd == "domain" )
      {
        while( isGood( lineIn ) )
        {
          Cube cube;
          lineIn >> cube;
          if( lineIn )
            cubes.push_back( cube );
        }
        if( lineIn )
          flags |= flagDomain;
      }
      else if( cmd == "periodic" )
      {
        int periodic = 0;
        while( isGood( lineIn ) )
        {
          int axis;
          lineIn >> axis;
          if( (axis < 0) || (axis >= dim) )
          {
            std::cerr << filename << "[ " << lineNr << " ]: Invalid periodic axis: " << axis << "." << std::endl;
            return false;
          }
          periodic |= (1 << axis);
        }
        topology = Topology( periodic );
      }
      else if( cmd == "cells" )
      {
        lineIn >> cells;
        if( lineIn )
          flags |= flagCells;
      }
      else if( cmd == "partitions" )
      {
        lineIn >> partitions;
      }
      else if( cmd == "overlap" )
        lineIn >> overlap;
      else if( cmd == "maxLevel" )
      {
        lineIn >> maxLevel;
        flags |= flagMaxLevel;
      }
      else if( cmd == "refinement" )
      {
        lineIn >> match( Refinement::type() );
        flags |= flagRefinement;
      }
      else if( cmd == "refinements" )
      {
        while( isGood( lineIn ) )
        {
          RefinementPolicy policy;
          lineIn >> policy;
          refinements.push_back( policy );
        }
      }
      else
      {
        std::cerr << filename << "[ " << lineNr << " ]: Invalid statement: '" << cmd << "'." << std::endl;
        return false;
      }
      if( lineIn.fail() )
      {
        std::cerr << filename << "[ " << lineNr << " ]: Invalid arguments for '" << cmd << "'." << std::endl;
        return false;
      }
    }

    if( flags != flagAll )
    {
      std::cerr << filename << ": File misses required field." << std::endl;
      return false;
    }
    return true;
  }


  template< class ctype, int dim, SPRefinementStrategy strategy >
  inline std::string
  SPGridIOData< ctype, dim, strategy >::readLine ( std::istream &in, unsigned int *count )
  {
    std::string line;
    while( line.empty() && !in.eof() )
    {
      std::getline( in, line );
      if( count != 0 )
        ++(*count);

      // remove leading white space
      const size_t first = line.find_first_not_of( " \t" );
      if( first != std::string::npos )
        line = line.substr( first );

      // remove trailing comments
      line = line.substr( 0, line.find_first_of( '#' ) );
    }
    return line;
  }

}

#endif // #ifndef DUNE_SPGRID_FILEIO_HH
