#ifndef DUNE_SPGRID_FILEIO_HH
#define DUNE_SPGRID_FILEIO_HH

#include <fstream>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/cube.hh>
#include <dune/grid/spgrid/topology.hh>
#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPGridIOData
  // ------------

  template< class ctype, int dim, template< int > class Ref >
  struct SPGridIOData
  {
    typedef SPTopology< dim > Topology;
    typedef SPCube< ctype, dim > Cube;
    typedef typename Cube::GlobalVector GlobalVector;
    typedef SPMultiIndex< dim > MultiIndex;
    typedef Ref< dim > Refinement;
    typedef typename Refinement::Policy RefinementPolicy;

    ctype time;
    std::vector< Cube > cubes;
    Topology topology;
    MultiIndex cells;
    MultiIndex overlap;
    int partitions;
    int maxLevel;
    std::vector< RefinementPolicy > refinements;

    bool write ( std::ostream &stream ) const;
    bool write ( const std::string &filename ) const;
    bool read ( std::istream &stream, const std::string &info = "" );
    bool read ( const std::string &filename );

  private:
    static std::pair< unsigned int, unsigned int > version ()
    {
      return std::pair< unsigned int, unsigned int >( DUNE_SPGRID_VERSION_MAJOR, DUNE_SPGRID_VERSION_MINOR );
    }

    static std::string readLine ( std::istream &stream, unsigned int *count = 0 );
  };



  // Implementation of SPGridIOData
  // ------------------------------

  template< class ctype, int dim, template< int > class Ref >
  inline bool SPGridIOData< ctype, dim, Ref >::write ( std::ostream &stream ) const
  {
    // write header
    stream << "SPGrid";
    stream << "  dimension=" << dim;
    stream << "  version=" << version().first << "." << version().second;
    stream << std::endl << std::endl;

    // write generic information
    stream << "time " << time << std::endl;
    stream << std::endl;

    // write domain information
    stream << "domain";
    for( typename std::vector< Cube >::const_iterator it = cubes.begin(); it != cubes.end(); ++it )
      stream << " " << *it;
    stream << std::endl;

    stream << "periodic";
    for( int i = 0; i < dim; ++i )
    {
      if( topology.periodic( i ) )
        stream << " " << i;
    }
    stream << std::endl << std::endl;

    // write discretization information
    stream << "cells " << cells << std::endl;
    stream << "partitions " << partitions << std::endl;
    stream << "overlap " << overlap << std::endl;
    stream << std::endl;

    // write refinement information
    stream << "maxLevel " << maxLevel << std::endl;
    stream << "refinement " << Refinement::type() << std::endl;
    stream << "refinements";
    for( unsigned int i = 0; i < refinements.size(); ++i )
      stream << " " << refinements[ i ];
    stream << std::endl;
    return bool( stream );
  }


  template< class ctype, int dim, template< int > class Ref >
  inline bool SPGridIOData< ctype, dim, Ref >::write ( const std::string &filename ) const
  {
    std::ofstream stream( filename.c_str() );
    return (stream ? write( stream ) : false);
  }


  template< class ctype, int dim, template< int > class Ref >
  inline bool
  SPGridIOData< ctype, dim, Ref >::read ( std::istream &stream, const std::string &info )
  {
    unsigned int lineNr = 0;
    std::string line = readLine( stream, &lineNr );
    std::istringstream lineIn( line );
    lineIn >> match( std::string( "SPGrid" ) );
    if( lineIn.fail() )
    {
      std::cerr << info << "[ " << lineNr << " ]: 'SPGrid' expected." << std::endl;
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
        std::pair< unsigned int, unsigned int > fileVersion;
        valueIn >> fileVersion.first >> match( '.' ) >> fileVersion.second;
        if( fileVersion > version() )
        {
          std::cerr << info << "[ " << lineNr << " ]: File was created by newer version of SPGrid." << std::endl;
          return false;
        }
      }
      else if( key == "dimension" )
        valueIn >> fdim;
      else
      {
        std::cerr << info << "[ " << lineNr << " ]: Invalid tag: '" << key << "'." << std::endl;
        return false;
      }
      if( !valueIn )
      {
        std::cerr << info << "[ " << lineNr << " ]: Invalid value for tag '" << key << "'." << std::endl;
        return false;
      }
    }

    if( fdim != dim )
    {
      std::cerr << info << "[ " << lineNr << " ]: File has wrong grid dimension." << std::endl;
      return false;
    }

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
      std::string line = readLine( stream, &lineNr );
      if( line.empty() )
        break;
      std::istringstream lineIn( line );

      std::string cmd;
      lineIn >> cmd;

      if( cmd == "time" )
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
            std::cerr << info << "[ " << lineNr << " ]: Invalid periodic axis: " << axis << "." << std::endl;
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
        std::cerr << info << "[ " << lineNr << " ]: Invalid statement: '" << cmd << "'." << std::endl;
        return false;
      }
      if( lineIn.fail() )
      {
        std::cerr << info << "[ " << lineNr << " ]: Invalid arguments for '" << cmd << "'." << std::endl;
        return false;
      }
    }

    if( flags != flagAll )
    {
      std::cerr << info << ": File misses required field." << std::endl;
      return false;
    }
    return true;
  }


  template< class ctype, int dim, template< int > class Ref >
  inline bool SPGridIOData< ctype, dim, Ref >::read ( const std::string &filename )
  {
    std::ifstream stream( filename.c_str() );
    return (stream ? read( stream, filename ) : false);
  }


  template< class ctype, int dim, template< int > class Ref >
  inline std::string SPGridIOData< ctype, dim, Ref >::readLine ( std::istream &stream, unsigned int *count )
  {
    std::string line;
    while( line.empty() && !stream.eof() )
    {
      std::getline( stream, line );
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

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_FILEIO_HH
