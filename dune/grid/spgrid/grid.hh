#ifndef DUNE_SPGRID_GRID_HH
#define DUNE_SPGRID_GRID_HH

#include <cstddef>

#include <array>
#include <memory>
#include <utility>

#include <dune/common/parallel/mpicommunication.hh>

#include <dune/grid/albertagrid/geometryreference.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh>

#include <dune/grid/spgrid/capabilities.hh>
#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/entityseed.hh>
#include <dune/grid/spgrid/gridview.hh>
#include <dune/grid/spgrid/hierarchiciterator.hh>
#include <dune/grid/spgrid/idset.hh>
#include <dune/grid/spgrid/indexset.hh>
#include <dune/grid/spgrid/hindexset.hh>
#include <dune/grid/spgrid/fileio.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  struct BackupRestoreFacility;

  namespace __SPGrid
  {

    template< class, class >
    class TreeIterator;

  } // namespace __SPGrid



  // Internal Forward Declarations
  // -----------------------------

#if HAVE_MPI
  template< class ct, int dim, template< int > class Ref = SPIsotropicRefinement, class Comm = MPI_Comm >
  class SPGrid;
#else
  template< class ct, int dim, template< int > class Ref = SPIsotropicRefinement, class Comm = No_Comm >
  class SPGrid;
#endif // #if !HAVE_MPI



  // SPGridFamily
  // ------------

  template< class ct, int dim, template< int > class Ref, class Comm >
  struct SPGridFamily
  {
    struct Traits
    {
      typedef SPGrid< ct, dim, Ref, Comm > Grid;

      typedef SPReferenceCubeContainer< ct, dim > ReferenceCubeContainer;
      typedef typename ReferenceCubeContainer::ReferenceCube ReferenceCube;
      typedef SPDomain< ct, dim > Domain;
      typedef SPMesh< dim > Mesh;
      typedef Ref< dim > Refinement;
      typedef typename Refinement::Policy RefinementPolicy;

      typedef typename SPCommunicationTraits< Comm >::Communication Communication;
      typedef Communication CollectiveCommunication;

      typedef Dune::Intersection< const Grid, SPIntersection< const Grid > > LevelIntersection;
      typedef LevelIntersection LeafIntersection;

      typedef Dune::IntersectionIterator< const Grid, SPIntersectionIterator< const Grid >, SPIntersection< const Grid > > LevelIntersectionIterator;
      typedef LevelIntersectionIterator LeafIntersectionIterator;

      typedef Dune::EntityIterator< 0, const Grid, SPHierarchicIterator< const Grid, 0 > > HierarchicIterator;

      typedef SPIndexSet< const Grid > LevelIndexSet;
      typedef LevelIndexSet LeafIndexSet;

      typedef SPGlobalIdSet< const Grid > GlobalIdSet;
      typedef SPLocalIdSet< const Grid > LocalIdSet;
      typedef unsigned long GlobalIdType;
      typedef unsigned long LocalIdType;

      template< int codim >
      struct Codim
      {
        typedef typename ReferenceCubeContainer::template Codim< codim >::ReferenceCube ReferenceCube;

        typedef Dune::Entity< codim, dim, const Grid, SPEntity > Entity;

        typedef Dune::EntitySeed< const Grid, SPEntitySeed< codim, const Grid > > EntitySeed;

        typedef SPLocalGeometry< dim - codim, dim, const Grid > LocalGeometryImpl;

        typedef Dune::Geometry< dim - codim, dim, const Grid, SPGeometry > Geometry;
        typedef Dune::Geometry< dim - codim, dim, const Grid, LocalGeometryReference > LocalGeometry;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef Dune::EntityIterator< codim, const Grid, SPPartitionIterator< codim, const Grid > > LevelIterator;
          typedef Dune::EntityIterator< codim, const Grid, SPPartitionIterator< codim, const Grid > > LeafIterator;
        };

        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
      };

      typedef Dune::GridView< SPGridViewTraits< const Grid > > LevelGridView;
      typedef Dune::GridView< SPGridViewTraits< const Grid > > LeafGridView;
    };
  };



  /** \class SPGrid
   *  \brief structured, parallel <a href="http://www.dune-project.org">%DUNE</a> grid
   *
   *  \tparam  ct        coordinate type (e.g., double)
   *  \tparam  dim       dimension of the grid
   *  \tparam  Ref       refinement (default is SPIsotropicRefinement)
   *  \tparam  Comm      type of communicator (default depends on HAVE_MPI)
   */
  template< class ct, int dim, template< int > class Ref, class Comm >
  class SPGrid
    : public GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim, Ref, Comm > >
  {
    typedef SPGrid< ct, dim, Ref, Comm > This;
    typedef GridDefaultImplementation< dim, dim, ct, SPGridFamily< ct, dim, Ref, Comm > > Base;

    friend struct BackupRestoreFacility< This >;
    friend class SPIntersection< const This >;
    friend class SPGridLevel< This >;

    template< class, class > friend class __SPGrid::TreeIterator;

  public:
    typedef SPGridFamily< ct, dim, Ref, Comm > GridFamily;

    typedef typename GridFamily::Traits Traits;

    typedef typename Traits::ReferenceCubeContainer ReferenceCubeContainer;
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

    typedef typename Traits::Communication Communication;
    typedef Communication CollectiveCommunication;

    template< int codim >
    struct Codim
    : Base::template Codim< codim >
    {
      typedef typename Traits::template Codim< codim >::ReferenceCube ReferenceCube;

      typedef typename Traits::template Codim< codim >::LocalGeometryImpl LocalGeometryImpl;
    };

    typedef typename Base::LevelGridView LevelGridView;
    typedef typename Base::LeafGridView LeafGridView;

    typedef typename LevelGridView::IndexSet LevelIndexSet;
    typedef typename LeafGridView::IndexSet LeafIndexSet;

    typedef SPHierarchyIndexSet< const This > HierarchicIndexSet;

    typedef SPGridLevel< This > GridLevel;

    typedef typename GridLevel::MultiIndex MultiIndex;
    static const int numDirections = GridLevel::numDirections;

  private:
    typedef typename GridLevel::PartitionList PartitionList;

    typedef typename LevelGridView::Traits::GridViewImp LevelGridViewImpl;
    typedef typename LeafGridView::Traits::GridViewImp LeafGridViewImpl;

  public:
    SPGrid ( const Domain &domain, const MultiIndex &cells,
             const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const Domain &domain, const MultiIndex &cells, const MultiIndex &overlap,
             const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
             const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
             const MultiIndex &overlap,
             const Communication &comm = SPCommunicationTraits< Comm >::defaultComm() );

    SPGrid ( const This & ) = delete;
    SPGrid ( This &&other );

    const ReferenceCube &referenceCube () const
    {
      return refCubes_.get();
    }

    template< int codim >
    const typename Codim< codim >::ReferenceCube &referenceCube () const
    {
      return refCubes_.template get< codim >();
    }

    const Domain &domain () const { return domain_; }

    const MultiIndex &overlap () const { return overlap_; }

    int maxLevel () const
    {
      return leafLevel().level();
    }

    int size ( const int level, const int codim ) const
    {
      return levelGridView( level ).size( codim );
    }

    int size ( const int codim ) const
    {
      return leafGridView().size( codim );
    }

    int size ( const int level, const GeometryType &type ) const
    {
      return levelGridView( level ).size( type );
    }

    int size ( const GeometryType &type ) const
    {
      return leafGridView().size( type );
    }

    LevelGridView levelGridView ( int level ) const
    {
      assert( (level >= 0) && (level <= maxLevel()) );
      return levelGridViews_[ level ];
    }

    LeafGridView leafGridView () const { return leafGridView_; }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lbegin ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelGridView( level );
      return view.impl().template begin< codim, pitype >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LevelIterator
    lend ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelGridView( level );
      return view.impl().template end< codim, pitype >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lbegin ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelGridView( level );
      return view.impl().template begin< codim >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LevelIterator
    lend ( const int level, const unsigned int sweepDir = 0 ) const
    {
      const LevelGridView &view = levelGridView( level );
      return view.impl().template end< codim >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafbegin ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafGridView();
      return view.impl().template begin< codim, pitype >( sweepDir );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Traits::template Codim< codim >::template Partition< pitype >::LeafIterator
    leafend ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafGridView();
      return view.impl().template end< codim, pitype >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafbegin ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafGridView();
      return view.impl().template begin< codim >( sweepDir );
    }

    template< int codim >
    typename Traits::template Codim< codim >::LeafIterator
    leafend ( const unsigned int sweepDir = 0 ) const
    {
      const LeafGridView &view = leafGridView();
      return view.impl().template end< codim >( sweepDir );
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
      return levelGridView( level ).indexSet();
    }

    const LeafIndexSet &leafIndexSet () const
    {
      return leafGridView().indexSet();
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
      return levelGridView( level ).overlapSize( codim );
    }

    int overlapSize ( const int codim ) const
    {
      return leafGridView().overlapSize( codim );
    }

    int ghostSize ( const int level, const int codim ) const
    {
      return levelGridView( level ).ghostSize( codim );
    }

    int ghostSize ( const int codim ) const
    {
      return leafGridView().ghostSize( codim );
    }

    template< class DataHandle, class Data >
    SPCommunication< This, CommDataHandleIF< DataHandle, Data > >
    communicate ( CommDataHandleIF< DataHandle, Data > &data,
                  InterfaceType interface, CommunicationDirection dir,
                  int level ) const
    {
      LevelGridView view = levelGridView( level );
      return view.impl().communicate( data, interface, dir );
    }

    template< class DataHandle, class Data >
    SPCommunication< This, CommDataHandleIF< DataHandle, Data > >
    communicate ( CommDataHandleIF< DataHandle, Data > &data,
                  InterfaceType interface, CommunicationDirection dir ) const
    {
      LeafGridView view = leafGridView();
      return view.impl().communicate( data, interface, dir );
    }

    const Communication &comm () const;

    template< class Seed >
    typename Traits::template Codim< Seed::codimension >::Entity entity ( const Seed &seed ) const
    {
      typedef typename Traits::template Codim< Seed::codimension >::Entity Entity;
      typedef SPEntity< Seed::codimension, dimension, const This > EntityImpl;
      typename EntityImpl::EntityInfo entityInfo( gridLevel( seed.impl().level() ), seed.impl().id(), seed.impl().partitionNumber() );
      return Entity( EntityImpl( std::move( entityInfo ) ) );
    }

    template< int codim >
    bool hasFather ( const Dune::Entity< codim, dimension, const This, SPEntity > &entity ) const
    {
      return ((entity.level() > 0) && entity.impl().entityInfo().hasFather());
    }

    bool hasFather ( const Dune::Intersection< const This, SPIntersection< const This > > &intersection ) const
    {
      return ((intersection.impl().gridLevel().level() > 0) && intersection.impl().entityInfo().hasFather());
    }

    template< int codim >
    Dune::Entity< codim, dimension, const This, SPEntity >
    father ( const Dune::Entity< codim, dimension, const This, SPEntity > &entity ) const
    {
      assert( hasFather( entity ) );
      Dune::Entity< codim, dimension, const This, SPEntity > father( entity );
      father.impl().entityInfo().up();
      return std::move( father );
    }

    Dune::Intersection< const This, SPIntersection< const This > >
    father ( const Dune::Intersection< const This, SPIntersection< const This > > &intersection ) const
    {
      typedef SPIntersection< const This > IntersectionImpl;
      typedef Dune::Intersection< const This, IntersectionImpl > Intersection;

      assert( hasFather( intersection ) );
      typename IntersectionImpl::EntityInfo fatherInfo( intersection.impl().entityInfo() );
      fatherInfo.up();
      return Intersection( IntersectionImpl( std::move( fatherInfo ), intersection.indexInInside() ) );
    }

    const GridLevel &gridLevel ( const int level ) const;
    const GridLevel &leafLevel () const;

    std::size_t numBoundarySegments () const;

  private:
    // note: this method ignores the last bit of the macroId
    std::size_t boundaryIndex ( const MultiIndex &macroId,
                           const unsigned int partitionNumber,
                           const int face ) const;

    typename Codim< 1 >::LocalGeometry localFaceGeometry ( int face ) const
    {
      assert( (face >= 0) && (face < ReferenceCube::numFaces) );
      return typename Codim< 1 >::LocalGeometry( *localFaceGeometry_[ face ] );
    }

    void createLocalGeometries ();
    void setupMacroGrid ();
    void setupBoundaryIndices ();

    static Communication defaultCommunication ();

    Domain domain_;
    Mesh globalMesh_;
    MultiIndex overlap_;
    ReferenceCubeContainer refCubes_;
    std::vector< std::unique_ptr< GridLevel > > gridLevels_;
    std::vector< LevelGridView > levelGridViews_;
    LeafGridView leafGridView_;
    HierarchicIndexSet hierarchicIndexSet_;
    GlobalIdSet globalIdSet_;
    LocalIdSet localIdSet_;
    Communication comm_;
    std::size_t boundarySize_;
    std::vector< std::array< std::size_t, 2*dimension > > boundaryOffset_;
    std::array< std::unique_ptr< const typename Codim< 1 >::LocalGeometryImpl >, ReferenceCube::numFaces > localFaceGeometry_;
  };



  // Implementation of SPGrid
  // ------------------------

  template< class ct, int dim, template< int > class Ref, class Comm >
  inline SPGrid< ct, dim, Ref, Comm >
    ::SPGrid ( const Domain &domain, const MultiIndex &cells,
               const Communication &comm )
  : domain_( domain ),
    globalMesh_( cells ),
    overlap_( MultiIndex::zero() ),
    leafGridView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline SPGrid< ct, dim, Ref, Comm >
    ::SPGrid ( const Domain &domain, const MultiIndex &cells, const MultiIndex &overlap,
               const Communication &comm )
  : domain_( domain ),
    globalMesh_( cells ),
    overlap_( overlap ),
    leafGridView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline SPGrid< ct, dim, Ref, Comm >
    ::SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const Communication &comm )
  : domain_( a, b ),
    globalMesh_( cells ),
    overlap_( MultiIndex::zero() ),
    leafGridView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline SPGrid< ct, dim, Ref, Comm >
    ::SPGrid ( const GlobalVector &a, const GlobalVector &b, const MultiIndex &cells,
               const MultiIndex &overlap, const Communication &comm )
  : domain_( a, b ),
    globalMesh_( cells ),
    overlap_( overlap ),
    leafGridView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( comm )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline SPGrid< ct, dim, Ref, Comm >::SPGrid ( This &&other )
  : domain_( std::move( other.domain_ ) ),
    globalMesh_( std::move( other.globalMesh_ ) ),
    overlap_( std::move( other.overlap_ ) ),
    leafGridView_( LeafGridViewImpl() ),
    hierarchicIndexSet_( *this ),
    comm_( std::move( other.comm_ ) )
  {
    createLocalGeometries();
    setupMacroGrid();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline bool SPGrid< ct, dim, Ref, Comm >
    ::mark ( const int refCount, const typename Codim< 0 >::Entity &e )
  {
    return false;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline int SPGrid< ct, dim, Ref, Comm >
    ::getMark ( const typename Codim< 0 >::Entity &e ) const
  {
    return 0;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline bool SPGrid< ct, dim, Ref, Comm >::preAdapt ()
  {
    return false;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline bool SPGrid< ct, dim, Ref, Comm >::adapt ()
  {
    return false;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  template< class DataHandle >
  inline bool SPGrid< ct, dim, Ref, Comm >
    ::adapt ( AdaptDataHandleInterface< This, DataHandle > &handle )
  {
    return false;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void SPGrid< ct, dim, Ref, Comm >::postAdapt ()
  {}


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void SPGrid< ct, dim, Ref, Comm >
    ::globalRefine ( const int refCount, const RefinementPolicy &policy )
  {
    for( int i = 0; i < refCount; ++i )
    {
      gridLevels_.emplace_back( new GridLevel( leafLevel(), policy ) );
      levelGridViews_.push_back( LevelGridViewImpl( leafLevel() ) );
    }
    leafGridView_.impl().update( leafLevel() );
    hierarchicIndexSet_.update();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  template< class DataHandle >
  inline void SPGrid< ct, dim, Ref, Comm >
    ::globalRefine ( const int refCount,
                     AdaptDataHandleInterface< This, DataHandle > &handle,
                     const RefinementPolicy &policy )
  {
    for( int i = 0; i < refCount; ++i )
    {
      const LevelGridView fatherView = levelGridView( maxLevel() );

      gridLevels_.emplace_back( new GridLevel( leafLevel(), policy ) );
      levelGridViews_.push_back( LevelGridViewImpl( leafLevel() ) );

      hierarchicIndexSet_.update();
      leafGridView_.impl().update( leafLevel() );

      handle.preAdapt( leafLevel().size() );
      typedef typename Codim< 0 >::LevelIterator LevelIterator;
      const LevelIterator end = fatherView.template end< 0 >();
      for( LevelIterator it = fatherView.template begin< 0 >(); it != end; ++it )
        handle.postRefinement( *it );
      handle.postAdapt();
    }
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline const typename SPGrid< ct, dim, Ref, Comm >::Communication &
  SPGrid< ct, dim, Ref, Comm >::comm () const
  {
    return comm_;
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline const typename SPGrid< ct, dim, Ref, Comm >::GridLevel &
  SPGrid< ct, dim, Ref, Comm >::gridLevel ( const int level ) const
  {
    assert( (level >= 0) && (level < int( gridLevels_.size() )) );
    return *gridLevels_[ level ];
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline const typename SPGrid< ct, dim, Ref, Comm >::GridLevel &
  SPGrid< ct, dim, Ref, Comm >::leafLevel () const
  {
    assert( !gridLevels_.empty() );
    return *gridLevels_.back();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline std::size_t SPGrid< ct, dim, Ref, Comm >::numBoundarySegments () const
  {
    return boundarySize_;
  }


  // note: this method ignores the last bit of the macroId
  template< class ct, int dim, template< int > class Ref, class Comm >
  inline std::size_t SPGrid< ct, dim, Ref, Comm >
    ::boundaryIndex ( const MultiIndex &macroId,
                      const unsigned int partitionNumber,
                      const int face ) const
  {
    assert( (face >= 0) && (face < 2*dimension) );

    const LevelGridView &macroView = levelGridView( 0 );
    const GridLevel &gridLevel = macroView.impl().gridLevel();
    const PartitionList &partitions = gridLevel.template partition< OverlapFront_Partition >();
    const typename PartitionList::Partition &partition = partitions.partition( partitionNumber );

    std::size_t index = 0;
    std::size_t factor = 1;
    for( int i = 0; i < dimension; ++i )
    {
      if( i == face/2 )
        continue;
      // note: the OverlapFront_Partition is closed, i.e.,
      //       begin()[ i ] & 1 == end()[ i ] & 1 == 0.
      const int k = (macroId[ i ] - partition.begin()[ i ]) >> 1;
      const int w = (partition.end()[ i ] - partition.begin()[ i ]) >> 1;
      assert( (k >= 0) && (k < w) );
      index += std::size_t( k ) * factor;
      factor *= std::size_t( w );
    }
    return index + boundaryOffset_[ partitionNumber - partitions.minNumber() ][ face ];
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void SPGrid< ct, dim, Ref, Comm >::createLocalGeometries ()
  {
    typedef typename Codim< 1 >::LocalGeometryImpl LocalGeometryImpl;

    const GlobalVector unitH( ctype( 1 ) );
    for( int face = 0; face < ReferenceCube::numFaces; ++face )
    {
      MultiIndex id;
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 1;
      id += referenceCube().subId( 1, face );
      const SPDirection< dimension > direction( id );
      //const unsigned int direction = ((1 << dimension) - 1) ^ (1 << (face/2));
      GlobalVector origin( ctype( 0 ) );
      origin[ face/2 ] = ctype( face & 1 );
      const SPGeometryCache< ctype, dimension, 1 > cache( unitH, direction );
      localFaceGeometry_[ face ].reset( new LocalGeometryImpl( cache, origin ) );
    }
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void SPGrid< ct, dim, Ref, Comm >::setupMacroGrid ()
  {
    SPDecomposition< dimension > decomposition( globalMesh_, comm().size() );

    GridLevel *leafLevel = new GridLevel( *this, decomposition );
    gridLevels_.emplace_back( leafLevel );
    levelGridViews_.push_back( LevelGridViewImpl( *leafLevel ) );
    leafGridView_.impl().update( *leafLevel );
    hierarchicIndexSet_.update();
    setupBoundaryIndices();
  }


  template< class ct, int dim, template< int > class Ref, class Comm >
  inline void SPGrid< ct, dim, Ref, Comm >::setupBoundaryIndices ()
  {
    const LevelGridView &macroView = levelGridView( 0 );
    const GridLevel &gridLevel = macroView.impl().gridLevel();
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
        std::size_t size = 1;
        for( int j = 0; j < dimension; ++j )
          size *= (i == j ? 1 : std::size_t( (it->end()[ j ] - it->begin()[ j ]) >> 1 ));

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
