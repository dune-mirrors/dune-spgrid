#ifndef DUNE_CARTESIANGRID_GRID_HH
#define DUNE_CARTESIANGRID_GRID_HH

#include <string>

#include <dune/common/version.hh>
#include <dune/common/static_assert.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/alugrid/common/interfaces.hh>

#include <dune/grid/cartesiangrid/capabilities.hh>
#include <dune/grid/cartesiangrid/entity.hh>
#include <dune/grid/cartesiangrid/entitypointer.hh>
#include <dune/grid/cartesiangrid/entityseed.hh>
#include <dune/grid/cartesiangrid/geometry.hh>
#include <dune/grid/cartesiangrid/intersection.hh>
#include <dune/grid/cartesiangrid/intersectioniterator.hh>
#include <dune/grid/cartesiangrid/iterator.hh>
#include <dune/grid/cartesiangrid/idset.hh>
#include <dune/grid/cartesiangrid/indexsets.hh>
#include <dune/grid/cartesiangrid/datahandle.hh>
#include <dune/grid/cartesiangrid/adaptcallback.hh>
#include <dune/grid/cartesiangrid/backuprestore.hh>

#include <dune/grid/geometrygrid/identity.hh>
#include <dune/grid/spgrid/geometricgridlevel.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< class HostGrid >
  class CartesianGrid;



  // CartesianGridExportParams
  // -------------------------

  template< class HG >
  struct CartesianGridExportParams
  {
    typedef HG HostGrid;
  };



  // CartesianGridFamily
  // -------------------
  
  template< class HostGrid >
  struct CartesianGridFamily
  {
    struct Traits
    : public CartesianGridExportParams< HostGrid >
    {
      typedef CartesianGrid< HostGrid > Grid;

      // export host grid type
      typedef HostGrid HostGridType;

      typedef typename HostGrid::ctype ctype;

      // type of data passed to entities, intersections, and iterators 
      typedef const Grid* ExtraData;

      static const int dimension = HostGrid::dimension;
      static const int dimensionworld = HostGrid::dimensionworld;

      typedef SPReferenceCubeContainer< ctype, dimension > ReferenceCubeContainer;
      typedef typename ReferenceCubeContainer::ReferenceCube ReferenceCube;

      typedef Dune::Intersection< const Grid, CartesianGridLeafIntersection > LeafIntersection;
      typedef Dune::Intersection< const Grid, CartesianGridLevelIntersection > LevelIntersection;

      typedef Dune::IntersectionIterator< const Grid, CartesianGridLeafIntersectionIterator, CartesianGridLeafIntersection >
        LeafIntersectionIterator;
      typedef Dune::IntersectionIterator< const Grid, CartesianGridLevelIntersectionIterator, CartesianGridLevelIntersection >
        LevelIntersectionIterator;

      typedef Dune::EntityIterator< 0, const Grid, CartesianGridIterator< CartesianGridHierarchicIteratorTraits< const Grid > > >
        HierarchicIterator;

      template< int codim >
      struct Codim
      {
        typedef typename ReferenceCubeContainer::template Codim< codim >::ReferenceCube ReferenceCube;

        typedef Dune::Geometry< dimension-codim, dimensionworld, const Grid, CartesianGridGeometry > Geometry;
        typedef Dune::Geometry< dimension-codim, dimension, const Grid, SPLocalGeometry > LocalGeometry;

        typedef CartesianGridEntityPointerTraits< codim, const Grid > EntityPointerTraits;
        typedef CartesianGridEntityPointer< EntityPointerTraits > EntityPointerImpl;
        typedef Dune::EntityPointer< const Grid, EntityPointerImpl > EntityPointer;
        typedef typename EntityPointerTraits::Entity Entity;
        typedef CartesianGridEntitySeed< codim, const Grid > EntitySeed;

        template< PartitionIteratorType pitype >
        struct Partition
        {
          typedef CartesianGridLeafIteratorTraits< codim, pitype, const Grid > LeafIteratorTraits;
          typedef Dune::EntityIterator< codim, const Grid, CartesianGridIterator< LeafIteratorTraits > >
            LeafIterator;

          typedef CartesianGridLevelIteratorTraits< codim, pitype, const Grid > LevelIteratorTraits;
          typedef Dune::EntityIterator< codim, const Grid, CartesianGridIterator< LevelIteratorTraits > >
            LevelIterator;
        };

        typedef typename Partition< All_Partition >::LeafIterator LeafIterator;
        typedef typename Partition< All_Partition >::LevelIterator LevelIterator;
      };

      typedef CartesianGridIndexSet< const Grid, typename HostGrid::Traits::LeafIndexSet > LeafIndexSet;
      typedef CartesianGridIndexSet< const Grid, typename HostGrid::Traits::LevelIndexSet > LevelIndexSet;

      typedef CartesianGridIdSet< const Grid, typename HostGrid::Traits::GlobalIdSet > GlobalIdSet;
      typedef CartesianGridIdSet< const Grid, typename HostGrid::Traits::LocalIdSet > LocalIdSet;

      typedef typename HostGrid::Traits::CollectiveCommunication CollectiveCommunication;

      template< PartitionIteratorType pitype >
      struct Partition
      {
        typedef Dune::GridView< DefaultLeafGridViewTraits< const Grid, pitype > >
          LeafGridView;
        typedef Dune::GridView< DefaultLevelGridViewTraits< const Grid, pitype > >
          LevelGridView;
      };
    };
  };

  template <class HostGrid, bool hasObjectStream>
  struct CartesianGridObjectStream 
  {
    // dummy typedef  
    typedef int ObjectStreamType ; 
  };

  template <class HostGrid>
  struct CartesianGridObjectStream<HostGrid, true> : public HasObjectStream 
  {
#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,2,0)
    typedef typename HostGrid :: InStreamType      InStreamType ; 
    typedef typename HostGrid :: OutStreamType     OutStreamType ; 
#else
    // use ObjectStream of HostGrid 
    typedef typename HostGrid :: ObjectStreamType  ObjectStreamType ; 
    typedef ObjectStreamType  InStreamType ;
    typedef ObjectStreamType  OutStreamType ;
#endif
  };



  // CartesianGrid
  // -------------

  /** \class CartesianGrid
   *  \brief identical grid wrapper
   *  \ingroup CartesianGrid
   *
   *  \tparam  HostGrid   DUNE grid to be wrapped (called host grid)
   *
   *  \nosubgrouping
   */
  template< class HostGrid >
  class CartesianGrid
  /** \cond */
  : public GridDefaultImplementation
      < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, CartesianGridFamily< HostGrid > >,
    public CartesianGridExportParams< HostGrid >,
    public CartesianGridBackupRestoreFacilities< CartesianGrid< HostGrid > >,
    public CartesianGridObjectStream< HostGrid, Conversion<HostGrid, HasObjectStream>::exists >
  /** \endcond */
  {
    typedef CartesianGrid< HostGrid > Grid;

    typedef GridDefaultImplementation
      < HostGrid::dimension, HostGrid::dimensionworld, typename HostGrid::ctype, CartesianGridFamily< HostGrid > >
      Base;

    template< int, int, class > friend class CartesianGridEntity;
    template< class > friend class CartesianGridEntityPointer;
    template< class, class > friend class CartesianGridIntersection;
    template< class > friend class CartesianGridIntersectionIterator;
    template< class, class > friend class CartesianGridIdSet;
    template< class, class > friend class CartesianGridIndexSet;
    template< class > friend class HostGridAccess;
    template< class > friend class DGFGridFactory;

    template< class, class > friend class CartesianGridDataHandle;

    template< int, PartitionIteratorType, class > friend class CartesianGridLevelIteratorTraits;
    template< int, PartitionIteratorType, class > friend class CartesianGridLeafIteratorTraits;

    friend class Conversion< Grid, HasObjectStream >;
    friend class Conversion< const Grid, HasObjectStream >;

  public:
    /** \cond */
    typedef CartesianGridFamily< HostGrid > GridFamily;
    /** \endcond */

    /** \name Traits
     *  \{ */

    //! type of the grid traits
    typedef typename GridFamily::Traits Traits;
    
    /** \brief traits structure containing types for a codimension
     *
     *  \tparam codim  codimension
     *
     *  \nosubgrouping
     */
    template< int codim >
    struct Codim;

    /** \} */

    /** \name Iterator Types
     *  \{ */
    
    //! iterator over the grid hierarchy
    typedef typename Traits::HierarchicIterator HierarchicIterator;
    //! iterator over intersections with other entities on the leaf level
    typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
    //! iterator over intersections with other entities on the same level
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */

    /** \name Index and Id Set Types
     *  \{ */

    /** \brief type of leaf index set
     *
     *  The index set assigns consecutive indices to the entities of the
     *  leaf grid. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The leaf index set is a model of Dune::IndexSet.
     */
    typedef typename Traits::LeafIndexSet LeafIndexSet;
    
    /** \brief type of level index set
     *
     *  The index set assigns consecutive indices to the entities of a grid
     *  level. The indices are of integral type and can be used to access
     *  arrays.
     *
     *  The level index set is a model of Dune::IndexSet.
     */
    typedef typename Traits::LevelIndexSet LevelIndexSet;

    /** \brief type of global id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier is unique over all processes sharing this grid.
     *
     *  \note Id's are neither consecutive nor necessarily of an integral
     *        type.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits::GlobalIdSet GlobalIdSet;

    /** \brief type of local id set
     *
     *  The id set assigns a unique identifier to each entity within the
     *  grid. This identifier needs only to be unique over this process.
     *
     *  Though the local id set may be identical to the global id set, it is
     *  often implemented more efficiently.
     *
     *  \note Ids are neither consecutive nor necessarily of an integral
     *        type.
     *  \note Local ids need not be compatible with global ids. Also, no
     *        mapping from local ids to global ones needs to exist.
     *
     *  The global id set is a model of Dune::IdSet.
     */
    typedef typename Traits::LocalIdSet LocalIdSet;

    /** \} */

    /** \name Miscellaneous Types
     * \{ */
    
    //! type of vector coordinates (e.g., double)
    typedef typename Traits::ctype ctype;

    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;

    //! communicator with all other processes having some part of the grid
    typedef typename Traits::CollectiveCommunication CollectiveCommunication;

    /** \} */

  private:
    typedef typename Traits::ReferenceCubeContainer ReferenceCubeContainer;
    typedef typename Traits::ReferenceCube ReferenceCube;
    typedef SPGeometricGridLevel< ctype, dimension > GridLevel;

    enum { MAXL = 32 };

  public:
    /** \name Construction and Destruction
     *  \{ */
    
    /** \brief constructor
     *
     *  The references to host grid and coordinate function are stored in the
     *  grid. Therefore, they must remain valid until the grid is destroyed.
     *
     *  \param[in]  hostGrid       pointer to the grid to wrap
     *
     *  \note The host grid pointer will be deleted upon destruction of the
     *        CartesianGrid.
     */
    explicit CartesianGrid ( HostGrid *hg )
    : hostGrid_( hg ),
      gridLevels_( MAXL ),
      levelIndexSets_( hostGrid().maxLevel()+1, (LevelIndexSet *)0 ),
      leafIndexSet_( 0 ),
      globalIdSet_( 0 ),
      localIdSet_( 0 )
    {
      FieldVector< ctype, dimensionworld > h = gridWidth();

      if( h.infinity_norm() < 1e-10 ) 
      {
        // communicate h in case of empty procs 
        assert( false );
        abort();
      }

      for( int level = 0; level < MAXL; ++level)
      {
        gridLevels_[ level ] = new GridLevel( refCubes_, h );
        h *= 0.5 ;
      }

      buildLocalIntersectionGeometries();
      buildGeometryInFather();
    }
    
    /** \brief destructor
     */
    ~CartesianGrid ()
    {
      // TODO: Delete local geometries

      if( localIdSet_ )
        delete localIdSet_;
      if( globalIdSet_ )
        delete globalIdSet_;

      if( leafIndexSet_ )
        delete leafIndexSet_;
      
      for( unsigned int i = 0; i < levelIndexSets_.size(); ++i )
      {
        if( levelIndexSets_[ i ] )
          delete( levelIndexSets_[ i ] );
      }

      for( unsigned int i = 0; i < gridLevels_.size(); ++i )
      {
        if( gridLevels_[ i ] )
          delete( gridLevels_[ i ] );
      }

      delete hostGrid_;
    }

    /** \} */
    
    /** \name Grid Identification Methods
     *  \{ */
    
    /** \brief obtain a string naming the grid
     *
     *  \returns ''CartesianGrid\< \em host \em grid \em name \>''
     */
    std::string name () const DUNE_DEPRECATED
    {
      return std::string( "CartesianGrid< " ) + hostGrid().name() + std::string( " >" );
    }

    /** \} */
    

    /** \name Size Methods
     *  \{ */

    /** \brief obtain maximal grid level
     *  
     *  Grid levels are numbered 0, ..., L, where L is the value returned by
     *  this method.
     *
     *  \returns maximal grid level
     */
    int maxLevel () const
    {
      return hostGrid().maxLevel();
    }
   
    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of entities of codimension \em codim on grid level
     *           \em level.
     */
    int size ( int level, int codim ) const
    {
      return hostGrid().size( level, codim );        
    }
    
    /** \brief obtain number of leaf entities
     *
     *  \param[in]  codim  codimension to consider
     *
     *  \returns number of leaf entities of codimension \em codim
     */
    int size ( int codim ) const
    {
      return hostGrid().size( codim );
    }

    /** \brief obtain number of entites on a level
     *
     *  \param[in]  level  level to consider
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of entities with a geometry of type \em type on grid
     *           level \em level.
     */
    int size ( int level, GeometryType type ) const
    {
      return hostGrid().size( level, type );
    }
    
    /** \brief returns the number of boundary segments within the macro grid 
     *
     *  \returns number of boundary segments within the macro grid
     */
    int size ( GeometryType type ) const
    {
      return hostGrid().size( type );
    }

    /** \brief obtain number of leaf entities
     *
     *  \param[in]  type   geometry type to consider
     *
     *  \returns number of leaf entities with a geometry of type \em type
     */
    size_t numBoundarySegments () const
    {
      return hostGrid().numBoundarySegments( );
    }
    /** \} */
    
    template< int codim >
    typename Codim< codim >::LevelIterator lbegin ( int level ) const
    {
      typedef CartesianGridLevelIteratorTraits< codim, All_Partition, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(),
                   hostGrid().template lbegin< codim >( level ) );
    }

    template< int codim >
    typename Codim< codim >::LevelIterator lend ( int level ) const
    {
      typedef CartesianGridLevelIteratorTraits< codim, All_Partition, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(), 
                   hostGrid().template lend< codim >( level ) );
    }
    
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LevelIterator
    lbegin ( int level ) const
    {
      typedef CartesianGridLevelIteratorTraits< codim, pitype, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(), 
                   hostGrid().template lbegin< codim, pitype >( level ) );
    }
    
    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LevelIterator
    lend ( int level ) const
    {
      typedef CartesianGridLevelIteratorTraits< codim, pitype, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(),
                   hostGrid().template lend< codim, pitype >( level ) );
    }
    
    template< int codim >
    typename Codim< codim >::LeafIterator leafbegin () const
    {
      typedef CartesianGridLeafIteratorTraits< codim, All_Partition, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(),
                   hostGrid().template leafbegin< codim >() ); 
    }
    
    template< int codim >
    typename Codim< codim >::LeafIterator leafend () const
    {
      typedef CartesianGridLeafIteratorTraits< codim, All_Partition, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(), hostGrid().template leafend< codim >() );
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LeafIterator
    leafbegin () const
    {
      typedef CartesianGridLeafIteratorTraits< codim, pitype, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(), 
                   hostGrid().template leafbegin< codim, pitype >() ); 
    }

    template< int codim, PartitionIteratorType pitype >
    typename Codim< codim >::template Partition< pitype >::LeafIterator
    leafend () const
    {
      typedef CartesianGridLeafIteratorTraits< codim, pitype, const Grid > T;
      typedef CartesianGridIterator< T > Impl;
      return Impl( extraData(), hostGrid().template leafend< codim, pitype >() );
    }
    
    const GlobalIdSet &globalIdSet () const
    {
      if( !globalIdSet_ )
        globalIdSet_ = new GlobalIdSet( hostGrid().globalIdSet() );
      assert( globalIdSet_ );
      return *globalIdSet_;
    }
    
    const LocalIdSet &localIdSet () const
    {
      if( !localIdSet_ )
        localIdSet_ = new LocalIdSet( hostGrid().localIdSet() );
      assert( localIdSet_ );
      return *localIdSet_;
    }
    
    const LevelIndexSet &levelIndexSet ( int level ) const
    {
      assert( levelIndexSets_.size() == (size_t)(maxLevel()+1) );
      if( (level < 0) || (level > maxLevel()) )
      {
        DUNE_THROW( GridError, "LevelIndexSet for nonexisting level " << level
                               << " requested." );
      }

      LevelIndexSet *&levelIndexSet = levelIndexSets_[ level ];
      if( !levelIndexSet )
        levelIndexSet = new LevelIndexSet( hostGrid().levelIndexSet( level ) );
      assert( levelIndexSet );
      return *levelIndexSet;
    }
    
    const LeafIndexSet &leafIndexSet () const
    {
      if( !leafIndexSet_ )
        leafIndexSet_ = new LeafIndexSet( hostGrid().leafIndexSet() );
      assert( leafIndexSet_ );
      return *leafIndexSet_;
    }
    
    void globalRefine ( int refCount )
    {
      hostGrid().globalRefine( refCount );
      // update overall status  
      update();
    }

    bool mark ( int refCount, const typename Codim< 0 >::Entity &entity )
    {
      return hostGrid().mark( refCount, getHostEntity< 0 >( entity ) );
    }

    int getMark ( const typename Codim< 0 >::Entity &entity ) const
    {
      return hostGrid().getMark( getHostEntity< 0 >( entity ) );
    }

    bool preAdapt ()
    {
      return hostGrid().preAdapt();
    }
    
    bool adapt ()
    {
      bool ret = hostGrid().adapt();
      update();
      return ret;
    }

    /** \brief  @copydoc Dune::Grid::adapt() 
        \param handle handler for restriction and prolongation operations 
        which is a Model of the AdaptDataHandleInterface class. 
    */
    template< class GridImp, class DataHandle >
    bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &datahandle )
    { 
      typedef CartesianGridAdaptDataHandle< Grid, 
              AdaptDataHandleInterface< GridImp, DataHandle > > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( extraData(), datahandle );
      const bool ret = hostGrid().adapt( wrappedDataHandle );
      update();
      return ret;
    }

    void postAdapt ()
    {
      hostGrid().postAdapt();
    }

    /** \name Parallel Data Distribution and Communication Methods
     *  \{ */

    /** \brief obtain size of overlap region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int overlapSize ( int codim ) const
    {
      return hostGrid().overlapSize( codim );
    }
    
    /** \brief obtain size of ghost region for the leaf grid
     *
     *  \param[in]  codim  codimension for with the information is desired
     */
    int ghostSize( int codim ) const
    {
      return hostGrid().ghostSize( codim );
    }
    
    /** \brief obtain size of overlap region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int overlapSize ( int level, int codim ) const
    {
      return hostGrid().overlapSize( level, codim );
    }
    
    /** \brief obtain size of ghost region for a grid level
     *
     *  \param[in]  level  grid level (0, ..., maxLevel())
     *  \param[in]  codim  codimension (0, ..., dimension)
     */
    int ghostSize ( int level, int codim ) const
    {
      return hostGrid().ghostSize( level, codim );
    }
        
    /** \brief communicate information on a grid level
     *
     *  \param      datahandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication or BackwardCommunication)
     *  \param[in]  level       grid level to communicate
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &datahandle,
                       InterfaceType interface,
                       CommunicationDirection direction,
                       int level ) const
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef CartesianGridDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( extraData(), datahandle );
      hostGrid().communicate( wrappedDataHandle, interface, direction, level );
    }
    
    /** \brief communicate information on leaf entities
     *
     *  \param      datahandle  communication data handle (user defined)
     *  \param[in]  interface   communication interface (one of
     *                          InteriorBorder_InteriorBorder_Interface,
     *                          InteriorBorder_All_Interface,
     *                          Overlap_OverlapFront_Interface,
     *                          Overlap_All_Interface,
     *                          All_All_Interface)
     *  \param[in]  direction   communication direction (one of
     *                          ForwardCommunication, BackwardCommunication)
     */
    template< class DataHandle, class Data >
    void communicate ( CommDataHandleIF< DataHandle, Data > &datahandle,
                       InterfaceType interface,
                       CommunicationDirection direction ) const
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef CartesianGridDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( extraData(), datahandle );
      hostGrid().communicate( wrappedDataHandle, interface, direction );
    }

    /** \brief obtain CollectiveCommunication object
     *
     *  The CollectiveCommunication object should be used to globally
     *  communicate information between all processes sharing this grid.
     *
     *  \note The CollectiveCommunication object returned is identical to the
     *        one returned by the host grid.
     */
    const CollectiveCommunication &comm () const
    {
      return hostGrid().comm();
    }

    // data handle interface different between geo and interface

    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \returns \b true, if the grid has changed.
     */
    bool loadBalance ()
    {
      const bool gridChanged = hostGrid().loadBalance();
      if( gridChanged ) 
        update();
      return gridChanged;
    }
    
    /** \brief rebalance the load each process has to handle
     *
     *  A parallel grid is redistributed such that each process has about
     *  the same load (e.g., the same number of leaf entites).
     *
     *  The data handle is used to communicate the data associated with
     *  entities that move from one process to another.
     *
     *  \note DUNE does not specify, how the load is measured.
     *
     *  \param  datahandle  communication data handle (user defined)
     *
     *  \returns \b true, if the grid has changed.
     */
    template< class DataHandle, class Data >
    bool loadBalance ( CommDataHandleIF< DataHandle, Data > &datahandle )
    {
      typedef CommDataHandleIF< DataHandle, Data > DataHandleIF;
      typedef CartesianGridDataHandle< Grid, DataHandleIF > WrappedDataHandle;

      WrappedDataHandle wrappedDataHandle( extraData(), datahandle );
      const bool gridChanged = hostGrid().loadBalance( wrappedDataHandle );
      if ( gridChanged ) 
        update();
      return gridChanged;
    }

    template< class DofManager >
    bool loadBalance ( DofManager& dofManager ) 
    {
      typedef CartesianGridWrappedDofManager< Grid, DofManager > WrappedDofManager;

      WrappedDofManager wrappedDofManager( extraData(), dofManager );
      const bool gridChanged = hostGrid().loadBalance( wrappedDofManager );
      if ( gridChanged ) update();
      return gridChanged;
    }

    /** \brief obtain EntityPointer from EntitySeed */
    template< class EntitySeed >
    typename Traits::template Codim< EntitySeed::codimension >::EntityPointer
    entityPointer ( const EntitySeed &seed ) const
    {
      typedef typename Traits::template Codim< EntitySeed::codimension >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( extraData(), hostGrid().entityPointer( seed.hostEntitySeed() ) );
    }

    /** \} */

    /** \name Miscellaneous Methods
     *  \{ */

    const HostGrid &hostGrid () const
    {
      return *hostGrid_;
    }

    HostGrid &hostGrid ()
    {
      return *hostGrid_;
    }

    const ReferenceCube &referenceCube () const
    {
      return refCubes_.get();
    }

    template< int codim >
    const typename Traits::template Codim< codim >::ReferenceCube &
    referenceCube () const
    {
      return refCubes_.template get< codim >();
    }

    //! return geometric information of given level 
    const GridLevel &geometricLevel ( const size_t level ) const 
    { 
      assert( level < gridLevels_.size() );
      assert( gridLevels_[ level ] );
      return *(gridLevels_[ level ]);
    }

    /** \brief update grid caches
     *
     *  This method has to be called whenever the underlying host grid changes.
     *
     *  \note If you adapt the host grid through this geometry grid's
     *        adaptation or load balancing methods, update is automatically
     *        called.
     */
    void update ()
    {
      const int newNumLevels = maxLevel()+1;
      const int oldNumLevels = levelIndexSets_.size();

      for( int i = newNumLevels; i < oldNumLevels; ++i )
      {
        if( levelIndexSets_[ i ] != 0 )
          delete levelIndexSets_[ i ];
      }
      levelIndexSets_.resize( newNumLevels, (LevelIndexSet *)0 );
    }

    /** \} */

  protected:
    FieldVector< ctype, dimensionworld > gridWidth () const 
    {
      typedef typename HostGrid::template Codim< 0 >::Geometry HostGeometry; 
      typedef typename HostGeometry::GlobalCoordinate GlobalCoordinate;

      GlobalCoordinate h( 0 );

      typedef typename HostGrid::template Codim< 0 >::LevelIterator MacroIterator;
      MacroIterator it = hostGrid().template lbegin< 0 >( 0 );
      if( it != hostGrid().template lend< 0 >( 0 ) )
      {
        if( it->type().id() != Capabilities::hasSingleGeometryType< Grid >::topologyId )
          DUNE_THROW( InvalidStateException, "HostGrid has wrong geometry type " << it->type() );

        const HostGeometry &geo = it->geometry();
        const GlobalCoordinate zero = geo.corner( 0 );
        for(int d = 0; d < dimensionworld; ++d ) 
          h[ d ] = std::abs( geo.corner( 1 << d )[ d ] - zero[ d ] );
      }
      return h;
    }

    void buildLocalIntersectionGeometries ()
    {
      typedef typename Codim< 1 >::LocalGeometry LocalGeo;
      typedef typename LocalGeo::GlobalCoordinate Coordinate;
      typedef SPLocalGeometry< dimension-1, dimension, const Grid > LocalGeoImpl;
      for( int face = 0; face < ReferenceCube::numFaces; ++face )
      {
        const unsigned int direction = ((1 << dimension) - 1) ^ (1 << (face/2));
        Coordinate origin( ctype( 0 ) );
        origin[ face/2 ] = ctype( face & 1 );
        const SPGeometryCache< ctype, dimension, 1 > cache( Coordinate( 1.0 ), direction );
        localFaceGeometry_[ face ] = new LocalGeo( LocalGeoImpl( referenceCube< 1 >(), cache, origin ) );

        for( int index = 0; index < (1 << (dimension-1)); ++index )
        {
          // origin[ face/2 ] = ctype( face & 1 );
          for( int i = 0; i < dimension-1; ++i )
          {
            const int j = (i < face/2 ? i : i+1 );
            origin[ j ] = 0.5 * ((index >> i) & 1);
          }
          const SPGeometryCache< ctype, dimension,1 > cache( Coordinate( 0.5 ), direction );
          localRefinedFaceGeometry_[ (1 << (dimension-1))*face + index ]
            = new LocalGeo( LocalGeoImpl( referenceCube< 1 >(), cache, origin ) );
        }
      }
    }

    void buildGeometryInFather ()
    {
      typedef typename Codim< 0 >::LocalGeometry LocalGeo;
      typedef SPLocalGeometry< dimension, dimension, const Grid > LocalGeoImpl;
      const unsigned int numChildren = (1 << dimension);
      const SPGeometryCache< ctype, dimension, 0 > cache( typename LocalGeo::GlobalCoordinate( 0.5 ), (1 << dimension)-1 );
      for( unsigned int index = 0; index < numChildren; ++index )
      {
        typename LocalGeo::GlobalCoordinate origin;
        for( int i = 0; i < dimension; ++i )
          origin[ i ] = 0.5 * ((index >> i) & 1);
        geometryInFather_[ index ] = new LocalGeo( LocalGeoImpl( referenceCube(), cache, origin ) );
      }
    }

    typedef typename Traits :: ExtraData ExtraData;
    //! return extra data passed to entities (e.g. const Grid*)
    ExtraData extraData () const  { return this; }

    using Base::getRealImplementation;

    template< int codim >
    static const typename HostGrid::template Codim< codim >::Entity &
    getHostEntity( const typename Codim< codim >::Entity &entity )
    {
      return getRealImplementation( entity ).hostEntity();
    }

  protected:
    HostGrid *const hostGrid_;
    ReferenceCubeContainer refCubes_;
    std::vector< GridLevel * > gridLevels_;
    mutable std::vector< LevelIndexSet * > levelIndexSets_;
    mutable LeafIndexSet *leafIndexSet_;
    mutable GlobalIdSet *globalIdSet_;
    mutable LocalIdSet *localIdSet_;

    const typename Codim< 0 >::LocalGeometry *geometryInFather_[ 1 << dimension ];
    const typename Codim< 1 >::LocalGeometry *localFaceGeometry_[ ReferenceCube::numFaces ];
    const typename Codim< 1 >::LocalGeometry *localRefinedFaceGeometry_[ (1 << (dimension-1))*ReferenceCube::numFaces ];
  };



  // CartesianGrid::Codim
  // --------------------

  template< class HostGrid >
  template< int codim >
  struct CartesianGrid< HostGrid >::Codim
  : public Base::template Codim< codim >
  {
    /** \name Entity and Entity Pointer Types
     *  \{ */
    
    /** \brief type of entity
     *
     *  The entity is a model of Dune::Entity.
     */
    typedef typename Traits::template Codim< codim >::Entity Entity;
    
    /** \brief type of entity pointer
     *
     *  The entity pointer is a model of Dune::EntityPointer.
     */
    typedef typename Traits::template Codim< codim >::EntityPointer EntityPointer;

    /** \} */

    /** \name Geometry Types
     *  \{ */

    /** \brief type of world geometry
     * 
     *  Models the geomtry mapping of the entity, i.e., the mapping from the
     *  reference element into world coordinates.
     *
     *  The geometry is a model of Dune::Geometry, implemented through the
     *  generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::Geometry Geometry;

    /** \brief type of local geometry
     * 
     *  Models the geomtry mapping into the reference element of dimension
     *  \em dimension.
     *
     *  The local geometry is a model of Dune::Geometry, implemented through
     *  the generic geometries provided by dune-grid.
     */
    typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

    /** \} */

    /** \name Iterator Types
     *  \{ */

    template< PartitionIteratorType pitype >
    struct Partition
    {
      typedef typename Traits::template Codim< codim >
        ::template Partition< pitype >::LeafIterator
        LeafIterator;
      typedef typename Traits::template Codim< codim >
        ::template Partition< pitype >::LevelIterator
        LevelIterator;
    };

    /** \brief type of level iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of a
     *  grid level.
     *
     *  The level iterator is a model of Dune::LevelIterator.
     */
    typedef typename Partition< All_Partition >::LeafIterator LeafIterator;

    /** \brief type of leaf iterator
     *
     *  This iterator enumerates the entites of codimension \em codim of the
     *  leaf grid.
     *
     *  The leaf iterator is a model of Dune::LeafIterator.
     */
    typedef typename Partition< All_Partition >::LevelIterator LevelIterator;

    /** \} */
  };

} // namespace Dune

#include <dune/grid/cartesiangrid/persistentcontainer.hh>
#include <dune/grid/cartesiangrid/twistutility.hh>
#endif // #ifndef DUNE_CARTESIANGRID_GRID_HH
