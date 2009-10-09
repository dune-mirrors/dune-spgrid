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

  template< class ctype, int dim, SPRefinementStrategy strategy >
  struct SPGridIOData
  {
    typedef FieldVector< ctype, dim > Vector;
    typedef SPMultiIndex< dim > MultiIndex;
    typedef SPRefinement< ctype, dim, strategy > Refinement;

    std::string name;
    ctype time;
    Vector origin;
    Vector width;
    MultiIndex cells;
    unsigned int periodic;
    int maxLevel;
    std::vector< Refinement > refinements;

    void writeAscii ( const std::string &filename ) const;
    void readAscii ( const std::string &filename );

  private:
    static std::string readLine ( std::istream &in, unsigned int *count = 0 );

    template< class T >
    static bool match ( std::istream &in, const T &t );
  };



  template< class ctype, int dim, SPRefinementStrategy strategy >
  inline void
  SPGridIOData< ctype, dim, strategy >::writeAscii ( const std::string &filename ) const
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

    fileOut << "refinements";
    for( unsigned int i = 0; i < refinements.size(); ++i )
      fileOut << " " << refinements[ i ];
    fileOut << std::endl;

    fileOut.close();
  }


  template< class ctype, int dim, SPRefinementStrategy strategy >
  inline void
  SPGridIOData< ctype, dim, strategy >::readAscii ( const std::string &filename )
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
      else if( cmd == "refinements" )
      {
        while( !lineIn.eof() )
        {
          Refinement refinement;
          lineIn >> refinement;
          refinements.push_back( refinement );
        }
      }
      else
        DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: Invalid statement: '" << cmd << "'." );
      if( !lineIn )
        DUNE_THROW( IOError, filename << "[ " << lineNr << " ]: Invalid arguments for '" << cmd << "'." );
    }

    if( flags != flagAll )
      DUNE_THROW( IOError, "SPGrid file misses required fields: '" << filename << "'." );
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


  template< class ctype, int dim, SPRefinementStrategy strategy >
  template< class T >
  inline bool
  SPGridIOData< ctype, dim, strategy >::match ( std::istream &in, const T &what )
  {
    T t;
    in >> t;
    if( !in )
      return false;
    return (t == what);
  }

}

#endif // #ifndef DUNE_SPGRID_FILEIO_HH
