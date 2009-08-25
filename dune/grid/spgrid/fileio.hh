#ifndef DUNE_SPGRID_FILEIO_HH
#define DUNE_SPGRID_FILEIO_HH

#include <fstream>
#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  template< class ctype, int dim >
  struct SPGridIOData
  {
    typedef FieldVector< ctype, dim > Vector;
    typedef SPMultiIndex< dim > MultiIndex;

    std::string name;
    ctype time;
    Vector origin;
    Vector width;
    MultiIndex cells;
    unsigned int periodic;
    int maxLevel;
    std::vector< unsigned int > refDirections;

    void writeAscii ( const std::string &filename ) const;
    void readAscii ( const std::string &filename );

  private:
    static std::string readLine ( std::istream &in, unsigned int *count = 0 );

    template< class T >
    static bool match ( std::istream &in, const T &t );
  };



  template< class ctype, int dim >
  inline void
  SPGridIOData< ctype, dim >::writeAscii ( const std::string &filename ) const
  {
    std::ofstream fileOut( filename.c_str() );

    fileOut << "SPGrid " << dim << std::endl << std::endl;
    fileOut << "name " << name << std::endl;
    fileOut << "time " << time << std::endl;
    fileOut << "origin " << origin << std::endl;
    fileOut << "width " << width << std::endl;

    fileOut << "cells";
    for( int i = 0; i < dim; ++i )
      fileOut << " " << cells[ i ];
    fileOut << std::endl;

    fileOut << "periodic";
    for( int i = 0; i < dim; ++i )
    {
      if( (periodic & (1 << i)) != 0 )
        fileOut << " " << i;
    }
    fileOut << std::endl;

    fileOut << "maxLevel " << maxLevel << std::endl;

    fileOut << "refDirections";
    for( unsigned int i = 0; i < refDirections.size(); ++i )
      fileOut << " " << refDirections[ i ];
    fileOut << std::endl;

    fileOut.close();
  }


  template< class ctype, int dim >
  inline void
  SPGridIOData< ctype, dim >::readAscii ( const std::string &filename )
  {
    std::ifstream fileIn( filename.c_str() );
    if( !fileIn )
      DUNE_THROW( IOError, "Unable to open file: '" << filename << "'." );

    unsigned int lineNr = 0;
    std::string line = readLine( fileIn, &lineNr );
    std::istringstream lineIn( line );
    if( !match( lineIn, std::string( "SPGrid" ) ) || !match( lineIn, dim ) )
      DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: 'SPGrid " << dim << "' expected." );

    name = "SPGrid";
    time = ctype( 0 );

    const unsigned int flagOrigin = 1;
    const unsigned int flagWidth = 2;
    const unsigned int flagCells = 4;
    const unsigned int flagMaxLevel = 8;
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
      else if( cmd == "origin" )
      {
        lineIn >> origin;
        if( !!lineIn )
          flags |= flagOrigin;
      }
      else if( cmd == "width" )
      {
        lineIn >> width;
        if( !!lineIn )
          flags |= flagWidth;
      }
      else if( cmd == "cells" )
      {
        for( int i = 0; i < dim; ++i )
          lineIn >> cells[ i ];
        flags |= flagCells;
      }
      else if( cmd == "periodic" )
      {
        while( !lineIn.eof() )
        {
          int axis;
          lineIn >> axis;
          if( (axis < 0) || (axis >= dim) )
            DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: Invalid periodic axis: " << axis << "." );
          periodic |= (1 << axis);
        }
      }
      else if( cmd == "maxLevel" )
      {
        lineIn >> maxLevel;
        flags |= flagMaxLevel;
      }
      else if( cmd == "refDirections" )
      {
        while( !lineIn.eof() )
        {
          unsigned int dir;
          lineIn >> dir;
          refDirections.push_back( dir );
        }
      }
      else
        DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: Invalid statement: '" << cmd << "'." );
      if( !lineIn )
        DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: Invalid arguments for '" << cmd << "'." );
    }

    if( flags != flagAll )
      DUNE_THROW( IOError, "SPGrid file misses required fields: '" << filename << "'." );

#if 0
    std::cerr << "name=" << name << std::endl;
    std::cerr << "time=" << time << std::endl;
    std::cerr << "origin=" << origin << std::endl;
    std::cerr << "width=" << width << std::endl;
    std::cerr << "cells=" << cells << std::endl;
    std::cerr << "maxLevel=" << maxLevel << std::endl;
    std::cerr << "refDirections=";
    for( size_t i = 0; i < refDirections.size(); ++i )
      std::cerr << (i > 0 ? " ": "") << refDirections[ i ];
    std::cerr << std::endl;
#endif
  }


  template< class ctype, int dim >
  inline std::string
  SPGridIOData< ctype, dim >::readLine ( std::istream &in, unsigned int *count )
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


  template< class ctype, int dim >
  template< class T >
  inline bool
  SPGridIOData< ctype, dim >::match ( std::istream &in, const T &what )
  {
    T t;
    in >> t;
    if( !in )
      return false;
    return (t == what);
  }

}

#endif // #ifndef DUNE_SPGRID_FILEIO_HH
