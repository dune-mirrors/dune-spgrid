#ifndef DUNE_SPGRID_FILEIO_HH
#define DUNE_SPGRID_FILEIO_HH

#include <fstream>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  template< class ct, int dim >
  struct SPGridIOData
  {
    typedef ct ctype;
    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > Vector;
    typedef SPMultiIndex< dimension > MultiIndex;

    std::string name;
    ctype time;
    Vector origin;
    Vector width;
    MultiIndex cells;
    int maxLevel;
    std::vector< unsigned int > refDirections;

    void writeAscii ( const std::string &filename ) const;
    void readAscii ( const std::string &filename );
  };



  template< class ct, int dim >
  void SPGridIOData< ct, dim >::writeAscii ( const std::string &filename ) const
  {
    std::ofstream fileOut( filename );
    fileOut << "SPGrid " << dimension << std::endl << std::endl;
    fileOut << "name " << name << std::endl;
    fileOut << "time " << time << std::endl;
    fileOut << "origin " << origin << std::endl;
    fileOut << "width " << width << std::endl;
    fileOut << "cells " << cells << std::endl;
    fileOut << "maxLevel " << maxLevel << std::endl;

    fileOut << "refDirections";
    for( unsigned int i = 0; i < refDirections.size(); ++i )
      fileOut << " " << refDirections[ i ];
    fileOut << std::endl;
    fileOut.close();
  }



  template< class ct, int dim >
  void SPGridIOData< ct, dim >::readAscii ( const std::string &filename )
  {
    std::ifstream fileIn( filename );
    if( !fileIn )
      DUNE_THROW( IOError, "Unable to open file: '" << filename << "'." );

    
  }

}

#endif // #ifndef DUNE_SPGRID_FILEIO_HH
