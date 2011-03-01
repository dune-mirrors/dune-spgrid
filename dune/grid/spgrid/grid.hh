#ifndef DUNE_SPGRID_GRID_HH
#define DUNE_SPGRID_GRID_HH

#include <dune/common/mpicollectivecommunication.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh>
#include <dune/grid/genericgeometry/codimtable.hh>
#include <dune/grid/utility/grapedataioformattypes.hh>

#include <dune/grid/spgrid/capabilities.hh>
#include <dune/grid/spgrid/entityseed.hh>
#include <dune/grid/spgrid/gridview.hh>
#include <dune/grid/spgrid/hierarchiciterator.hh>
#include <dune/grid/spgrid/idset.hh>
#include <dune/grid/spgrid/indexset.hh>
#include <dune/grid/spgrid/hindexset.hh>
#include <dune/grid/spgrid/fileio.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

#if HAVE_MPI
  template< class ct, int dim, SPRefinementStrategy strategy = SPIsotropicRefinement, class Comm = MPI_Comm >
  class SPGrid;
#else
  template< class ct, int dim, SPRefinementStrategy strategy = SPIsotropicRefinement, class Comm = No_Comm >
  class SPGrid;
#endif // #if !HAVE_MPI



  // SPGridFamily
  // ------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  struct SPGridFamily
  {
    struct Traits
    {
      typedef SPGrid< ct, dim, strategy, Comm > Grid;

      typedef SPReferenceCube< ct, dim > ReferenceCube;
      typedef SPDomain< ct, dim > Domain;
      typedef SPMesh< dim > Mesh;
      typedef SPRefinement< dim, strategy > Refinement;
      typedef typename Refinement::Policy RefinementPolicy;

      typedef typename SPCommunicationTraits< Comm >::CollectiveCommunication CollectiveCommunication;

      typedef Dune::Intersection< const Grid, SPIntersection >
        LevelIntersection;
      typedef LevelIntersection LeafIntersection;

      typedef Dune::IntersectionIterator
        < const Grid, SPIntersectionIterator, SPIntersection >
        LevelIntersectionIterator;
      typedef LevelIntersectionIterator LeafIntersectionIterator;

      typedef Dune::EntityIterator< 0, const Grid, SPHierarchicIterator< const Grid > >
        HierarchicIterator;

      typedef SPIndexSet< const Grid > LevelIndexSet;
      typedef LevelIndexSet LeafIndexSet;

      typedef SPGlobalIdSet< const Grid > GlobalIdSet;
      typedef SPLocalIdSet< const Grid > LocalIdSet;
      typedef unsigned long GlobalIdType;
      typedef unsigned long LocalIdType;

      template< int codim >
      struct Codim
      {
        typedef SPReferenceCube< ct, dim-codim > ReferenceCube;

        typedef Dune::Entity< codim, dim, const Grid, SPEntity > Entity;

        typedef SPEntitySeed< codim, const Grid > EntitySeed;

        typedef SPEntityPointer< codim, const Grid > EntityPointerImpl;
        typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;

        typedef Dune::Geometry< dim - codim, dim, const Grid, SPGeometry > Geometry;
        typedef Dune::Geometry< dim - codim, dim, const Grid, SPLocalGeometry >
          LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::EntityIterator< codim, const Grid, SPIterator< codim, pitype, const Grid > >
            LevelIterator;
          typedef Dune::EntityIterator< codim, const Grid, SPIterator< codim, pitype, const Grid > >
            LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune::GridView< SPGridViewTraits< const Grid, pitype > > LevelGridView;
        typedef Dune::GridView< SPGridViewTraits< const Grid, pitype > > LeafGridView;
      };
    };
  };



  /** \class SPGrid
   *  \brief structured, parallel <a href="http://www.dune-project.org">%DUNE</a> grid
   *
   *  \tparam  ct        coordinate type (e.g., double)
   *  \tparam  dim       dimension of the grid
   *  \tparam  strategy  refinement strategy (default is \ref SPRefinementStrategy "SPIsotropicRefinement")
   *  \tparam  Comm      type of communicator (default depends on HAVE_MPI)
   */
  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  class SPGrid
  : public GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim, strategy, Comm > >
  {
    typedef SPGrid< ct, dim, strategy, Comm > This;
    typedef GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim, strategy, Comm > > Base;

    friend class SPIntersection< const This >;
    friend class SPGridLevel< const This >;

  public:
    typedef SPGridFamily< ct, dim, strategy, Comm > GridFamily;

    typedef typename GridFamily::Traits Traits;

    typedef typename Traits::ReferenceCube ReferenceCube;
    typedef typename Traits::Domain Domain;
    typedef typename Traits::Mesh Mesh;
    typedef typename Traits::Refinement Refinement;
    typedef typename Traits::RefinementPolicy RefinementPolicy;

    typedef typename ReferenceCube::ctype ctype;

    static const int dimension = ReferenceCube::dimension;
    static const int dimensionworld = ReferenceCube::dimension;

    typedef typename ReferenceCube::GlobalVector GlobalVector;

    typedef typename Traits::GlobalIdSet GlobalIdSet;
    typedef typename Traits::LocalIdSet LocalIdSet;

    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    template< int codim >
    struct Codim
    : Base::template Codim< codim >
    {
      typedef typename Traits::template Codim< codim >::ReferenceCube ReferenceCube;
    };

    template< PartitionIteratorType pitype >
    struct Partition
    : Base::template Partition< pitype >
    {};
    
    typedef typename Partition< All_Partition >::LevelGridView LevelGridView;
    typedef typename Partition< All_Partition >::LeafGridView LeafGridView;

    typedef typename LevelGridView::IndexSet LevelIndexSet;
    typedef typename LeafGridView::IndexSet LeafIndexSet;

    typedef SPHierarchyIndexSet< const This > HierarchicIndexSet;

    typedef SPGridLevel< const This > GridLevel;

    typedef typename GridLevel::MultiIndex MultiIndex;
    static const int numDirections = GridLevel::numDirections;

  private:
    typedef typename GridLevel::PartitionList PartitionList;

    typedef typename LevelGridView::Traits::GridViewImp LevelGridViewImpl;
    typedef typename LeafGridView::Traits::GridViewImp LeafGridViewImpl;

  public:
    SPGrid ( const CollectiveCommunication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const Domain &domain, const MultiIndex &cells,
             const CollectiveCommunication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const Domain &domain, const MultiIndex &cells, const MultiIndex &overlap,
             const CollectiveCommunication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
             const CollectiveCommunication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
             const MultiIndex &overlap,
             const CollectiveCommunication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    ~SPGrid ()
    {
      clear();
    }

    using Base::getRealImplementation;

    const ReferenceCube &referenceCube () const
    {
      return referenceCube< 0 >();
    }

    SPDecomposition< dimension > decomposition () const 
    {
      return SPDecomposition< dimension > ( globalMesh_, comm().size() );
    }

    template< int codim >
    const typename Codim< codim >::ReferenceCube &referenceCube () const
    {
      integral_constant< int, codim > codimVariable;
      return refCubes_[ codimVariable ];
    }

    const Domain &domain () const
    {
      return domain_;
    }

    const MultiIndex &overlap () const
    {
      return overlap_;
    }

    int maxLevel () const
    {
      return leafLevel().level();
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
    levelView ( const int level ) const;

    template< PartitionIteratorType pitype >
    typename Traits::template Partition< pitype >::LeafGridView
    leafView () const;

    LevelGridView levelView ( const int level ) const;
    LeafGridView leafView () const;

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

    const HierarchicIndexSet &hierarchicIndexSet () const
    {
      return hierarchicIndexSet_;
    }

    bool mark ( const int refCount, const typename Codim< 0 >::Entity &e );
    int getMark ( const typename Codim< 0 >::Entity &e ) const;

    bool preAdapt ();

    bool adapt ();

    template< class DataHandle >
    bool adapt ( AdaptDataHandleInterface< This, DataHandle > &handle );

    void postAdapt ();

    void globalRefine ( const int refCount,
                        const RefinementPolicy &policy = RefinementPolicy() );

    template< class DataHandle >
    void globalRefine ( const int refCount,
                        AdaptDataHandleInterface< This, DataHandle > &handle,
                        const RefinementPolicy &policy = RefinementPolicy() );

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

    const CollectiveCommunication &comm () const;

    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::EntityPointer
    entityPointer ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( gridLevel( seed.level() ), seed.id(), seed.partitionNumber() );
    }

    template< GrapeIOFileFormatType format >
    bool writeGrid ( const std::string &filename, const ctype &time ) const;

    template< GrapeIOFileFormatType format >
    bool readGrid ( const std::string &filename, ctype &time );

    const GridLevel &gridLevel ( const int level ) const;
    const GridLevel &leafLevel () const;

    size_t numBoundarySegments () const;

    typename Codim< 0 >::EntityPointer findEntity ( const GlobalVector &x, const int level ) const;

  private:
    // note: this method ignores the last bit of the macroId
    size_t boundaryIndex ( const MultiIndex &macroId,
                           const unsigned int partitionNumber,
                           const int face ) const;

    const typename Codim< 1 >::LocalGeometry &localFaceGeometry ( const int face ) const;

    void clear ();

    void createLocalGeometries ();
    void setupMacroGrid ();
    void setupBoundaryIndices ();

    static CollectiveCommunication defaultCommunication ();

    template< int codim >
    struct RefCube
    : public Codim< codim >::ReferenceCube
    {};
  
    Domain domain_;
    Mesh globalMesh_;
    MultiIndex overlap_;
    GenericGeometry::CodimTable< RefCube, dimension > refCubes_;
    std::vector< GridLevel * > gridLevels_;
    std::vector< LevelGridView > levelViews_;
    LeafGridView leafView_;
    HierarchicIndexSet hierarchicIndexSet_;
    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;
    CollectiveCommunication comm_;
    size_t boundarySize_;
    std::vector< array< size_t, 2*dimension > > boundaryOffset_;
    const typename Codim< 1 >::LocalGeometry *localFaceGeometry_[ ReferenceCube::numFaces ];
  };



  // Implementation of SPGrid
  // ------------------------

  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline SPGrid< ct, dim, strategy, Comm >::SPGrid ( const CollectiveCommunication &comm )
  : domain_( Domain::unitCube() ),
    globalMesh_( Mesh::unitMesh() ),
    overlap_( MultiIndex::zero() ),
    leafView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline SPGrid< ct, dim, strategy, Comm >
    ::SPGrid ( const Domain &domain, const MultiIndex &cells,
               const CollectiveCommunication &comm )
  : domain_( domain ),
    globalMesh_( cells ),
    overlap_( MultiIndex::zero() ),
    leafView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline SPGrid< ct, dim, strategy, Comm >
    ::SPGrid ( const Domain &domain, const MultiIndex &cells, const MultiIndex &overlap,
               const CollectiveCommunication &comm )
  : domain_( domain ),
    globalMesh_( cells ),
    overlap_( overlap ),
    leafView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline SPGrid< ct, dim, strategy, Comm >
    ::SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const CollectiveCommunication &comm )
  : domain_( a, b ),
    globalMesh_( cells ),
    overlap_( MultiIndex::zero() ),
    leafView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline SPGrid< ct, dim, strategy, Comm >
    ::SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const MultiIndex &overlap, const CollectiveCommunication &comm )
  : domain_( a, b ),
    globalMesh_( cells ),
    overlap_( overlap ),
    leafView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< PartitionIteratorType pitype >
  typename SPGrid< ct, dim, strategy, Comm >::Traits::template Partition< pitype >::LevelGridView
  SPGrid< ct, dim, strategy, Comm >::levelView ( const int level ) const
  {
    typedef typename Traits::template Partition< pitype >::LevelGridView GridView;
    typedef typename GridView::Traits::GridViewImp GridViewImpl;
    assert( (level >= 0) && (level <= maxLevel()) );
    const LevelGridViewImpl &viewImpl = getRealImplementation( levelViews_[ level ] );
    return GridViewImpl( viewImpl );
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< PartitionIteratorType pitype >
  typename SPGrid< ct, dim, strategy, Comm >::Traits::template Partition< pitype >::LeafGridView
  SPGrid< ct, dim, strategy, Comm >::leafView () const
  {
    typedef typename Traits::template Partition< pitype >::LeafGridView GridView;
    typedef typename GridView::Traits::GridViewImp GridViewImpl;
    const LeafGridViewImpl &viewImpl = getRealImplementation( leafView_ );
    return GridViewImpl( viewImpl );
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  typename SPGrid< ct, dim, strategy, Comm >::LevelGridView
  SPGrid< ct, dim, strategy, Comm >::levelView ( const int level ) const
  {
    assert( (level >= 0) && (level <= maxLevel()) );
    return levelViews_[ level ];
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  typename SPGrid< ct, dim, strategy, Comm >::LeafGridView
  SPGrid< ct, dim, strategy, Comm >::leafView () const
  {
    return leafView_;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline bool SPGrid< ct, dim, strategy, Comm >
    ::mark ( const int refCount, const typename Codim< 0 >::Entity &e )
  {
    return false;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline int SPGrid< ct, dim, strategy, Comm >
    ::getMark ( const typename Codim< 0 >::Entity &e ) const
  {
    return 0;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline bool SPGrid< ct, dim, strategy, Comm >::preAdapt ()
  {
    return false;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline bool SPGrid< ct, dim, strategy, Comm >::adapt ()
  {
    return false;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< class DataHandle >
  inline bool SPGrid< ct, dim, strategy, Comm >
    ::adapt ( AdaptDataHandleInterface< This, DataHandle > &handle )
  {
    return false;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >::postAdapt ()
  {}


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >
    ::globalRefine ( const int refCount, const RefinementPolicy &policy )
  {
    for( int i = 0; i < refCount; ++i )
    {
      gridLevels_.push_back( new GridLevel( leafLevel(), policy ) );
      levelViews_.push_back( LevelGridViewImpl( leafLevel() ) );
    }
    getRealImplementation( leafView_ ).update( leafLevel() );
    hierarchicIndexSet_.update();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< class DataHandle >
  inline void SPGrid< ct, dim, strategy, Comm >
    ::globalRefine ( const int refCount,
                     AdaptDataHandleInterface< This, DataHandle > &handle,
                     const RefinementPolicy &policy )
  {
    for( int i = 0; i < refCount; ++i )
    {
      const LevelGridView fatherView = levelView( maxLevel() );

      gridLevels_.push_back( new GridLevel( leafLevel(), policy ) );
      levelViews_.push_back( LevelGridViewImpl( leafLevel() ) );

      hierarchicIndexSet_.update();
      getRealImplementation( leafView_ ).update( leafLevel() );

      handle.preAdapt( leafLevel().size() );
      typedef typename Codim< 0 >::LevelIterator LevelIterator;
      const LevelIterator end = fatherView.template end< 0 >();
      for( LevelIterator it = fatherView.template begin< 0 >(); it != end; ++it )
        handle.postRefinement( *it );
      handle.postAdapt();
    }
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline const typename SPGrid< ct, dim, strategy, Comm >::CollectiveCommunication &
  SPGrid< ct, dim, strategy, Comm >::comm () const
  {
    return comm_;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< GrapeIOFileFormatType format >
  inline bool SPGrid< ct, dim, strategy, Comm >
    ::writeGrid ( const std::string &filename, const ctype &time ) const
  {
    // we ignore the format and always write ascii

    int result = 0;
    if( comm().rank() == 0 )
    {
      SPGridIOData< ctype, dimension, strategy > ioData;

      ioData.time = time;
      ioData.cubes.push_back( domain().cube() );
      ioData.topology = domain().topology();
      ioData.cells = globalMesh_.width();
      ioData.partitions = comm().size();
      ioData.overlap = overlap_;
      ioData.maxLevel = maxLevel();
      ioData.refinements.resize( maxLevel() );
      for( int level = 0; level < maxLevel(); ++level )
        ioData.refinements[ level ] = gridLevel( level+1 ).refinement().policy();

      result = int( ioData.write( filename ) );
    }
    comm().broadcast( &result, 1, 0 );
    return (result != 0);
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  template< GrapeIOFileFormatType format >
  inline bool SPGrid< ct, dim, strategy, Comm >
    ::readGrid ( const std::string &filename, ctype &time )
  {
    // we ignore the format and always read ascii

    int result = 0;
    SPGridIOData< ctype, dimension, strategy > ioData;

    if( ioData.read( filename ) )
    {
      time = ioData.time;

      if( ioData.partitions != comm().size() )
      {
        std::cerr << "Warning: Reading grid with different number of partitions,"
                  << " index sets will not coincide." << std::endl;
      }

      domain_ = Domain( ioData.cubes, ioData.topology );
      globalMesh_ = Mesh( ioData.cells );
      overlap_ = ioData.overlap;
      setupMacroGrid();

      for( int level = 0; level < ioData.maxLevel; ++level )
      {
        if( level < int( ioData.refinements.size() ) )
          globalRefine( 1, ioData.refinements[ level ] );
        else
          globalRefine( 1 );
      }
    }

    result = comm().sum( result );
    return (result == comm().size());
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline const typename SPGrid< ct, dim, strategy, Comm >::GridLevel &
  SPGrid< ct, dim, strategy, Comm >::gridLevel ( const int level ) const
  {
    assert( (level >= 0) && (level < int( gridLevels_.size() )) );
    return *gridLevels_[ level ];
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline const typename SPGrid< ct, dim, strategy, Comm >::GridLevel &
  SPGrid< ct, dim, strategy, Comm >::leafLevel () const
  {
    assert( !gridLevels_.empty() );
    return *gridLevels_.back();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline size_t SPGrid< ct, dim, strategy, Comm >::numBoundarySegments () const
  {
    return boundarySize_;
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline typename SPGrid< ct, dim, strategy, Comm >::template Codim< 0 >::EntityPointer
  SPGrid< ct, dim, strategy, Comm >
    ::findEntity ( const GlobalVector &x, const int level ) const
  {
    typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
    assert( domain().contains( x ) );
    const GridLevel &gLevel = gridLevel( level );
    const PartitionList &partitionList = gLevel.template partition< All_Partition >();

    const GlobalVector y = x - domain().cube().origin();
    GlobalVector z;
    gLevel.template geometryCache< 0 >( (1 << dimension) - 1 ).jacobianInverseTransposed().mv( y, z );

    MultiIndex id;
    for( int i = 0; i < dimension; ++i )
      id[ i ] = 2*int( z[ i ] ) + 1;

    const typename PartitionList::Partition *partition = partitionList.findPartition( id );
    if( partition )
      return EntityPointerImpl( gLevel, id, partition->number() );
    else
      DUNE_THROW( GridError, "Coordinate " << x << " is outside the grid." );
  }


  // note: this method ignores the last bit of the macroId
  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline size_t SPGrid< ct, dim, strategy, Comm >
    ::boundaryIndex ( const MultiIndex &macroId,
                      const unsigned int partitionNumber,
                      const int face ) const
  {
    assert( (face >= 0) && (face < 2*dimension) );

    const LevelGridView &macroView = levelView( 0 );
    const GridLevel &gridLevel = getRealImplementation( macroView ).gridLevel();
    const PartitionList &partitions = gridLevel.template partition< OverlapFront_Partition >();
    const typename PartitionList::Partition &partition = partitions.partition( partitionNumber );

    size_t index = 0;
    size_t factor = 1;
    for( int i = 0; i < dimension; ++i )
    {
      if( i == face/2 )
        continue;
      // note: the OverlapFront_Partition is closed, i.e.,
      //       begin()[ i ] & 1 == end()[ i ] & 1 == 0.
      const int k = (macroId[ i ] - partition.begin()[ i ]) >> 1;
      const int w = (partition.end()[ i ] - partition.begin()[ i ]) >> 1;
      assert( (k >= 0) && (k < w) );
      index += size_t( k ) * factor;
      factor *= size_t( w );
    }
    return index + boundaryOffset_[ partitionNumber - partitions.minNumber() ][ face ];
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline const typename SPGrid< ct, dim, strategy, Comm >::template Codim< 1 >::LocalGeometry &
  SPGrid< ct, dim, strategy, Comm >::localFaceGeometry ( const int face ) const
  {
    assert( (face >= 0) && (face < ReferenceCube::numFaces) );
    return *localFaceGeometry_[ face ];
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >::clear ()
  {
    levelViews_.clear();
    leafView_ = LeafGridView( LeafGridViewImpl() );

    typedef typename std::vector< GridLevel * >::iterator Iterator;
    const Iterator end = gridLevels_.end();
    for( Iterator it = gridLevels_.begin(); it != end; ++it )
      delete *it;
    gridLevels_.clear();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >::createLocalGeometries ()
  {
    typedef typename Codim< 1 >::LocalGeometry LocalGeo;
    typedef SPLocalGeometry< dimension-1, dimension, const This > LocalGeoImpl;
    const GlobalVector unitH( ctype( 1 ) );
    for( int face = 0; face < ReferenceCube::numFaces; ++face )
    {
      const unsigned int direction = ((1 << dimension) - 1) ^ (1 << (face/2));
      GlobalVector origin( ctype( 0 ) );
      origin[ face/2 ] = ctype( face & 1 );
      const SPGeometryCache< ctype, dimension, 1 > cache( unitH, direction );
      localFaceGeometry_[ face ] = new LocalGeo( LocalGeoImpl( referenceCube< 1 >(), cache, origin ) );
    }
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >::setupMacroGrid ()
  {
    clear();

    SPDecomposition< dimension > decomposition( globalMesh_, comm().size() );

    GridLevel *leafLevel = new GridLevel( *this, decomposition );
    gridLevels_.push_back( leafLevel );
    levelViews_.push_back( LevelGridViewImpl( *leafLevel ) );
    getRealImplementation( leafView_ ).update( *leafLevel );
    hierarchicIndexSet_.update();
    setupBoundaryIndices();
  }


  template< class ct, int dim, SPRefinementStrategy strategy, class Comm >
  inline void SPGrid< ct, dim, strategy, Comm >::setupBoundaryIndices ()
  {
    const LevelGridView &macroView = levelView( 0 );
    const GridLevel &gridLevel = getRealImplementation( macroView ).gridLevel();
    const PartitionList &partitions = gridLevel.template partition< OverlapFront_Partition >();

    boundarySize_ = 0;
    boundaryOffset_.resize( partitions.maxNumber() - partitions.minNumber() + 1 );
    for( typename PartitionList::Iterator it = partitions.begin(); it; ++it )
    {
      const int partitionIndex = it->number() - partitions.minNumber();
      for( int i = 0; i < dimension; ++i )
      {
        // note: the OverlapFront_Partition is closed, i.e.,
        //       begin()[ i ] & 1 == end()[ i ] & 1 == 0.
        size_t size = 1;
        for( int j = 0; j < dimension; ++j )
          size *= (i == j ? 1 : size_t( (it->end()[ j ] - it->begin()[ j ]) >> 1 ));

        // offset for lower boundary
        boundaryOffset_[ partitionIndex ][ 2*i ] = boundarySize_;
        boundarySize_ += (it->boundary( 2*i ) ? size : 0);

        // offset for upper boundary
        boundaryOffset_[ partitionIndex ][ 2*i+1 ] = boundarySize_;
        boundarySize_ += (it->boundary( 2*i+1 ) ? size : 0);
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GRID_HH
