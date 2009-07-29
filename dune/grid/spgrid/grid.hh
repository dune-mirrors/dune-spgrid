#ifndef DUNE_SPGRID_GRID_HH
#define DUNE_SPGRID_GRID_HH

#include <fstream>

#include <dune/common/collectivecommunication.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/genericgeometry/codimtable.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/grid/spgrid/capabilities.hh>
#include <dune/grid/spgrid/gridview.hh>
#include <dune/grid/spgrid/hiterator.hh>
#include <dune/grid/spgrid/idset.hh>
#include <dune/grid/spgrid/indexset.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ct, int dim >
  class SPGrid;



  // SPGridFamily
  // ------------

  template< class ct, int dim >
  struct SPGridFamily
  {
    struct Traits
    {
      typedef SPGrid< ct, dim > Grid;

      typedef SPCube< ct, dim > Cube;
      typedef SPDomain< ct, dim > Domain;

      typedef Dune::CollectiveCommunication< Grid > CollectiveCommunication;

      typedef Dune::Intersection< const Grid, SPIntersection >
        LevelIntersection;
      typedef LevelIntersection LeafIntersection;

      typedef Dune::IntersectionIterator
        < const Grid, SPIntersectionIterator, SPIntersection >
        LevelIntersectionIterator;
      typedef LevelIntersectionIterator LeafIntersectionIterator;

      typedef Dune::HierarchicIterator< const Grid, SPHierarchicIterator >
        HierarchicIterator;

      typedef SPIndexSet< const Grid > LevelIndexSet;
      typedef LevelIndexSet LeafIndexSet;

      typedef SPGlobalIdSet< const Grid > GlobalIdSet;
      typedef SPLocalIdSet< const Grid > LocalIdSet;
      //typedef typename GlobalIdSet::IdType GlobalIdType;
      //typedef typename LocalIdSet::IdType LocalIdType;
      typedef unsigned int GlobalIdType;
      typedef unsigned int LocalIdType;

      template< int codim >
      struct Codim
      {
        typedef SPCube< ct, dim-codim > Cube;

        typedef Dune::Entity< codim, dim, const Grid, SPEntity > Entity;

        typedef SPEntityPointer< codim, const Grid > EntityPointerImpl;
        typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;

        typedef Dune::Geometry< dim - codim, dim, const Grid, SPGeometry > Geometry;
        typedef Dune::Geometry< dim - codim, dim, const Grid, SPLocalGeometry >
          LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::LevelIterator< codim, pitype, const Grid, SPIterator >
            LevelIterator;
          typedef Dune::LeafIterator< codim, pitype, const Grid, SPIterator >
            LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune::GridView< SPLevelGridViewTraits< const Grid, pitype > >
          LevelGridView;
        typedef Dune::GridView< SPLeafGridViewTraits< const Grid, pitype > >
          LeafGridView;
      };
    };
  };



  template< class ct, int dim >
  class SPGrid
  : public GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim > >
  {
    typedef SPGrid< ct, dim > This;
    typedef GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim > > Base;

    friend class SPIntersection< const This >;

  public:
    typedef SPGridFamily< ct, dim > GridFamily;

    typedef typename GridFamily::Traits Traits;

    typedef typename Traits::Cube Cube;
    typedef typename Traits::Domain Domain;

    typedef typename Cube::ctype ctype;

    static const int dimension = Cube::dimension;
    static const int dimensionworld = Cube::dimension;

    typedef typename Cube::GlobalVector GlobalVector;

    typedef typename Traits::GlobalIdSet GlobalIdSet;
    typedef typename Traits::LocalIdSet LocalIdSet;

    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    template< int codim >
    struct Codim
    : Base::template Codim< codim >
    {
      typedef typename Traits::template Codim< codim >::Cube Cube;
    };

    template< PartitionIteratorType pitype >
    struct Partition
    : Base::template Partition< pitype >
    {};
    
    typedef typename Partition< All_Partition >::LevelGridView LevelGridView;
    typedef typename Partition< All_Partition >::LeafGridView LeafGridView;

    typedef typename LevelGridView::IndexSet LevelIndexSet;
    typedef typename LeafGridView::IndexSet LeafIndexSet;

    typedef SPGridLevel< const This > GridLevel;

    static const int numDirections = GridLevel::numDirections;

  private:
    typedef typename LevelGridView::Traits::GridViewImp LevelGridViewImpl;
    typedef typename LeafGridView::Traits::GridViewImp LeafGridViewImpl;

  public:
    SPGrid ()
    : domain_(),
      name_( "SPGrid" ),
      leafLevel_( 0 ),
      leafView_( LeafGridViewImpl() ),
      globalIdSet_(),
      localIdSet_(),
      comm_()
    {
      int cells[ dimension ];
      for( int i = 0; i < dimension; ++i )
        cells[ i ] = 1;
      leafLevel_ = new GridLevel( *this, cells );

      levelViews_.push_back( LevelGridViewImpl( *leafLevel_ ) );
      getRealImplementation( leafView_ ).update( *leafLevel_ );

      createLocalGeometries();
    }

    SPGrid ( const GlobalVector &a, const GlobalVector &b,
             const int (&cells)[ dimension ],
             const std::string &name = "SPGrid" )
    : domain_( a, b ),
      name_( name ),
      leafLevel_( new GridLevel( *this, cells ) ),
      leafView_( LeafGridViewImpl() ),
      globalIdSet_(),
      localIdSet_(),
      comm_()
    {
      levelViews_.push_back( LevelGridViewImpl( *leafLevel_ ) );
      getRealImplementation( leafView_ ).update( *leafLevel_ );

      createLocalGeometries();
    }

    using Base::getRealImplementation;

    const Cube &cube () const
    {
      return cube< 0 >();
    }

    template< int codim >
    const typename Codim< codim >::Cube &cube () const
    {
      Int2Type< codim > codimVariable;
      return cubes_[ codimVariable ];
    }

    const Domain &domain () const
    {
      return domain_;
    }

    const std::string &name () const
    {
      return name_;
    }

    int maxLevel () const
    {
      assert( levelViews_.size() == leafLevel_->level()+1 );
      return leafLevel_->level();
    }

    int size ( const int level, const int codim ) const
    {
      return levelView( level ).size( codim );
    }

    int size ( const int codim ) const
    {
      return leafView().size( codim );
    }

    int size ( const int level, const GeometryType &type ) const
    {
      return levelView( level ).size( type );
    }

    int size ( const GeometryType &type ) const
    {
      return leafView().size( type );
    }

    template< PartitionIteratorType pitype >
    typename Traits::template Partition< pitype >::LevelGridView
    levelView ( const int level ) const
    {
      typedef typename Traits::template Partition< pitype >::LevelGridView GridView;
      typedef typename GridView::Traits::GridViewImp GridViewImpl;
      const LevelGridViewImpl &viewImpl = getRealImplementation( levelViews_[ level ] );
      return GridViewImpl( viewImpl );
    }

    template< PartitionIteratorType pitype >
    typename Traits::template Partition< pitype >::LeafGridView
    leafView () const
    {
      typedef typename Traits::template Partition< pitype >::LeafGridView GridView;
      typedef typename GridView::Traits::GridViewImp GridViewImpl;
      const LeafGridViewImpl &viewImpl = getRealImplementation( leafView_ );
      return GridViewImpl( viewImpl );
    }

    LevelGridView levelView ( const int level ) const
    {
      assert( (level >= 0) && (level <= maxLevel()) );
      return levelViews_[ level ];
    }

    LeafGridView leafView () const
    {
      return leafView_;
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lbegin ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelView( level );
      return getRealImplementation( view ).template begin< codim, pitype >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lend ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelView( level );
      return getRealImplementation( view ).template end< codim, pitype >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lbegin ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelView( level );
      return getRealImplementation( view ).template begin< codim >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lend ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelView( level );
      return getRealImplementation( view ).template end< codim >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafbegin ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafView();
      return getRealImplementation( view ).template begin< codim, pitype >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafend ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafView();
      return getRealImplementation( view ).template end< codim, pitype >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafbegin ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafView();
      return getRealImplementation( view ).template begin< codim >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafend ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafView();
      return getRealImplementation( view ).template end< codim >( sweepDir );
    }

    const GlobalIdSet &globalIdSet () const
    {
      return globalIdSet_;
    }

    const LocalIdSet &localIdSet () const
    {
      return localIdSet_;
    }

    const LevelIndexSet &levelIndexSet ( const int level ) const
    {
      return levelView( level ).indexSet();
    }

    const LeafIndexSet &leafIndexSet () const
    {
      return leafView().indexSet();
    }

    void globalRefine ( const int refCount,
                        const unsigned int refDir = numDirections-1 )
    {
      for( int i = 0; i < refCount; ++i )
      {
        leafLevel_ = new GridLevel( *leafLevel_, refDir );
        levelViews_.push_back( LevelGridViewImpl( *leafLevel_ ) );
      }
      getRealImplementation( leafView_ ).update( *leafLevel_ );
    }

    int overlapSize ( const int level, const int codim ) const
    {
      return levelView( level ).overlapSize( codim );
    }

    int overlapSize ( const int codim ) const
    {
      return leafView().overlapSize( codim );
    }

    int ghostSize ( const int level, const int codim ) const
    {
      return levelView( level ).ghostSize( codim );
    }

    int ghostSize ( const int codim ) const
    {
      return leafView().ghostSize( codim );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType interface, CommunicationDirection dir,
                       int level ) const
    {
      levelView( level ).communicate( data, interface, dir );
    }

    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &data,
                       InterfaceType interface, CommunicationDirection dir ) const
    {
      leafView().communicate( data, interface, dir );
    }

    const CollectiveCommunication &comm () const
    {
      return comm_;
    }

    template< GrapeIOFileFormatType format >
    bool writeGrid ( const std::string &filename, const ctype &time ) const
    {
      if( format == xdr )
      {
        return false;
      }
      else if( format == ascii )
      {
        std::ofstream fileOut( filename );
        fileOut << "SPGrid< " << dimension << " >" << std::endl;
        fileOut << "name: " << name() << std::endl;
        fileOut << "time: " << time << std::endl;
        fileOut << "origin: " << domain().origin() << std::endl;
        fileOut << "width: " << domain().width() << std::endl;
        fileOut << "cells: " << getRealImplemmentation( levelView( 0 ) ).gridLevel().cells() << std::endl;
        fileOut << "maxLevel: " << leafLevel_.level() << std::endl;
      }
      else
        DUNE_THROW( NotImplemented, "SPGrid: Unknwon output format: " << format << "." );
    }

    template< GrapeIOFileFormatType format >
    bool readGrid ( const std::string &filename, ctype &time )
    {
      if( format == xdr )
      {
        return false;
      }
      else if( format == ascii )
      {
        return false;
      }
      else
        DUNE_THROW( NotImplemented, "SPGrid: Unknwon output format: " << format << "." );
    }

  private:
    const typename Codim< 1 >::LocalGeometry &localFaceGeometry ( const int face ) const
    {
      assert( (face >= 0) && (face < Cube::numFaces) );
      return *localFaceGeometry_[ face ];
    }

    void createLocalGeometries ()
    {
      typedef typename Codim< 1 >::LocalGeometry LocalGeo;
      typedef SPLocalGeometry< dimension-1, dimension, const This > LocalGeoImpl;
      const GlobalVector unitH( ctype( 1 ) );
      for( int face = 0; face < Cube::numFaces; ++face )
      {
        const unsigned int direction = ((1 << dimension) - 1) ^ (1 << (face/2));
        GlobalVector origin( ctype( 0 ) );
        origin[ face/2 ] = ctype( face & 1 );
        const SPGeometryCache< ctype, dimension, 1 > cache( unitH, direction );
        localFaceGeometry_[ face ] = new LocalGeo( LocalGeoImpl( cube< 1 >(), cache, origin ) );
      }
    }

    template< int codim >
    struct TheCube
    : public Codim< codim >::Cube
    {};
  
    Domain domain_;
    std::string name_;
    GenericGeometry::CodimTable< TheCube, dimension > cubes_;
    GridLevel *leafLevel_;
    std::vector< LevelGridView > levelViews_;
    LeafGridView leafView_;
    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;
    CollectiveCommunication comm_;
    const typename Codim< 1 >::LocalGeometry *localFaceGeometry_[ Cube::numFaces ];
  };

}

#endif // #ifndef DUNE_SPGRID_GRID_HH
