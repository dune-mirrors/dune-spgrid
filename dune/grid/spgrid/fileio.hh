#ifndef DUNE_SPGRID_FILEIO_HH
#define DUNE_SPGRID_FILEIO_HH

#include <fstream>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPGridIOData
  // ------------

  template< class ctype, int dim, SPRefinementStrategy strategy >
  struct SPGridIOData
  {
    typedef FieldVector< ctype, dim > Vector;
    typedef SPMultiIndex< dim > MultiIndex;
    typedef SPRefinement< dim, strategy > Refinement;

    static const unsigned int versionMajor = DUNE_SPGRID_VERSION_MAJOR;
    static const unsigned int versionMinor = DUNE_SPGRID_VERSION_MINOR;

    std::string name;
    ctype time;
    Vector origin;
    Vector width;
    MultiIndex cells;
    MultiIndex overlap;
    int partitions;
    unsigned int periodic;
    int maxLevel;
    std::vector< Refinement > refinements;

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
    fileOut << "  version=" << versionMajor << "." << versionMinor;
    fileOut << std::endl << std::endl;

    // write generic information
    fileOut << "name " << name << std::endl;
    fileOut << "time " << time << std::endl;
    fileOut << std::endl;

    // write domain information
    fileOut << "origin " << origin << std::endl;
    fileOut << "width " << width << std::endl;

    fileOut << "periodic";
    for( int i = 0; i < dim; ++i )
    {
      if( (periodic & (1 << i)) != 0 )
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

    while( lineIn.good() )
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
        unsigned int vMajor = versionMajor, vMinor = versionMinor;
        valueIn >> vMajor >> match( '.' ) >> vMinor;
        if( (vMajor > versionMajor) || ((vMajor == versionMajor) && (vMinor > versionMinor)) )
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

    const unsigned int flagOrigin = 1;
    const unsigned int flagWidth = 2;
    const unsigned int flagCells = 4;
    const unsigned int flagMaxLevel = 8;
    const unsigned int flagRefinement = 16;
    const unsigned int flagAll = 31;
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
      else if( cmd == "origin" )
      {
        lineIn >> origin;
        if( lineIn )
          flags |= flagOrigin;
      }
      else if( cmd == "width" )
      {
        lineIn >> width;
        if( lineIn )
          flags |= flagWidth;
      }
      else if( cmd == "periodic" )
      {
        while( lineIn.good() )
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
        while( lineIn.good() )
        {
          Refinement refinement;
          lineIn >> refinement;
          refinements.push_back( refinement );
        }
      }
      else
      {
        std::cerr << filename << "[ " << lineNr << " ]: Invalid statement: '" << cmd << "'." << std::endl;
        return false;
      }
      if( !lineIn )
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
