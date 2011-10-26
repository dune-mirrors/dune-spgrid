#ifndef DUNE_CARTESIANGRID_BACKUPRESTORE_HH
#define DUNE_CARTESIANGRID_BACKUPRESTORE_HH

#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/grid/cartesiangrid/capabilities.hh>

namespace Dune
{

  // BackupRestoreFacilities
  // -----------------------

  template< class Grid, bool hasBackupRestoreFacilities = Capabilities::hasBackupRestoreFacilities< Grid > ::v >
  class CartesianGridBackupRestoreFacilities
  {};

  template< class Grid >
  class CartesianGridBackupRestoreFacilities< Grid, true >
  {
    typedef CartesianGridBackupRestoreFacilities< Grid, true > This;

  protected:
    CartesianGridBackupRestoreFacilities ()
    {}

  private:
    CartesianGridBackupRestoreFacilities ( const This & );
    This &operator= ( const This & );

  public:
    template< GrapeIOFileFormatType type >
    bool writeGrid ( const std::string &filename, double time ) const
    {
      return asImp().hostGrid().template writeGrid< type >( filename, time );
    }

    template< GrapeIOFileFormatType type >
    bool readGrid ( const std::string &filename, double &time )
    {
      const bool success
        = asImp().hostGrid().template readGrid< type >( filename, time );
      asImp().update();
      return success;
    }

  protected:
    const Grid &asImp () const
    {
      return static_cast< const Grid & >( *this );
    }

    Grid &asImp ()
    {
      return static_cast< Grid & >( *this );
    }
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_BACKUPRESTORE_HH
