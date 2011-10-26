#ifndef DUNE_CARTESIANGRID_ENTITY_HH
#define DUNE_CARTESIANGRID_ENTITY_HH

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/cartesiangrid/entityseed.hh>
#include <dune/grid/cartesiangrid/geometry.hh>
#include <dune/grid/cartesiangrid/hostgridinfo.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class > class CartesianGridLeafIntersectionIterator;
  template< class > class CartesianGridLevelIntersectionIterator;

  template< class > class CartesianGridIterator;
  template< class > class CartesianGridHierarchicIteratorTraits;



  // CartesianGridEntity
  // -------------------

  /** \copydoc CartesianGridEntity
   *
   *  \nosubgrouping
   */
  template< int codim, int dim, class Grid >
  class CartesianGridEntity
  {
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    /** \name Attributes
     *  \{ */

    //! codimensioon of the entity
    static const int codimension = codim;
    //! dimension of the grid
    static const int dimension = Traits::dimension;
    //! dimension of the entity
    static const int mydimension = dimension - codimension;
    //! dimension of the world
    static const int dimensionworld = Traits::dimensionworld;

    /** \} */

    /** \name Types Required by DUNE
     *  \{ */

    //! coordinate type of the grid
    typedef typename Traits::ctype ctype;

    //! type of corresponding entity seed 
    typedef typename Grid::template Codim< codimension >::EntitySeed EntitySeed;
    //! type of corresponding geometry
    typedef typename Traits::template Codim< codimension >::Geometry Geometry;

    /** \} */
    
  private:
    typedef typename Traits::HostGrid HostGrid;
    typedef typename Traits::ExtraData ExtraData; 

    // type of host grid info 
    typedef CartesianGridHostGridInfo< HostGrid > HostGridInfo;

  public:
    /** \name Host Types
     *  \{ */

    //! type of corresponding host entity
    typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
    //! type of corresponding host entity pointer
    typedef typename HostGrid::template Codim< codimension >::EntityPointer HostEntityPointer;
    /** \} */

  private:
    typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

    typedef CartesianGridGeometry< mydimension, dimensionworld, Grid > GeometryImpl;

  public:
    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    CartesianGridEntity ( ExtraData data ) 
    : hostEntity_( 0 ),
      data_( data ),
      geo_( GeometryImpl( grid().geometricLevel( 0 ), 
                          HostGridInfo::defaultDirection( mydimension ), 
                          HostGridInfo::defaultOrigin() ) )
    {}

    /** \brief construct an initialized entity
     *
     *  \param[in]  hostEntity  corresponding entity in the host grid
     *
     *  \note The reference to the host entity must remain valid  as long as
     *        this entity is in use.
     */
    explicit CartesianGridEntity ( ExtraData data, 
                                   const HostEntity &hostEntity )
    : hostEntity_( &hostEntity ),
      data_( data ),
      geo_( GeometryImpl( grid().geometricLevel( level() ), 
                          HostGridInfo::direction( hostEntity ), 
                          HostGridInfo::origin( hostEntity ) ) )
    {}

    /** \brief copy constructor */
    CartesianGridEntity ( const CartesianGridEntity &other )
    : hostEntity_( other.hostEntity_ ),
      data_( other.data_ ),
      geo_( Grid::getRealImplementation( other.geo_ ) )
    {}

    const CartesianGridEntity &operator= ( const CartesianGridEntity &other )
    {
      hostEntity_ = other.hostEntity_;
      data_ = other.data_;
      Grid::getRealImplementation( geo_ ) = Grid::getRealImplementation( other.geo_ );
      return *this;
    }

    /** \} */

    operator bool () const { return bool( hostEntity_ ); }

    /** \name Methods Shared by Entities of All Codimensions
     *  \{ */

    /** \brief obtain the name of the corresponding reference element
     *
     *  This type can be used to access the DUNE reference element.
     */
    GeometryType type () const
    {
      return hostEntity().type();
    }

    /** \brief obtain the level of this entity */
    int level () const
    {
      return hostEntity().level();
    }
    
    /** \brief obtain the partition type of this entity */
    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

    /** obtain the geometry of this entity
     *
     *  Each DUNE entity encapsulates a geometry object, representing the map
     *  from the reference element to world coordinates. Wrapping the geometry
     *  is the main objective of the CartesianGrid.
     *
     *  The CartesianGrid provides geometries of order 1, obtained by
     *  interpolation of its corners \f$y_i\f$. There corners are calculated
     *  from the corners \f$x_i\f$ of the host geometry through the
     *  CartesianGrid's coordinate function \f$c\f$, i.e.,
     *  \f$y_i = c( x_i )\f$.
     *
     *  \returns a const reference to the geometry
     */
    const Geometry &geometry () const
    {
      return geo_;
    }

    /** \brief return EntitySeed of host grid entity */
    EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

    /** \} */


    /** \name Methods Supporting the Grid Implementation
     *  \{ */

    const HostEntity &hostEntity () const
    {
      assert( *this );
      return *hostEntity_;
    }

    const Grid& grid() const 
    {
      return *data();
    }

    ExtraData data() const 
    {
      assert( data_ );
      return data_; 
    }

    /** \} */

  private:
    const HostEntity *hostEntity_;
    ExtraData data_;
    mutable Geometry geo_;
  };



  // CartesianGridEntity for codimension 0
  // -------------------------------------

  /** \copydoc CartesianGridEntity
   *
   *  \nosubgrouping
   */
  template< int dim, class Grid >
  class CartesianGridEntity< 0, dim, Grid >
  {
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    /** \name Attributes
     *  \{ */
        
    //! codimensioon of the entity
    static const int codimension = 0;
    //! dimension of the grid
    static const int dimension = Traits::dimension;
    //! dimension of the entity
    static const int mydimension = dimension - codimension;
    //! dimension of the world
    static const int dimensionworld = Traits::dimensionworld;

    /** \} */

    /** \name Types Required by DUNE
     *  \{ */

    //! coordinate type of the grid
    typedef typename Traits::ctype ctype;

    //! type of corresponding entity seed 
    typedef typename Grid::template Codim< codimension >::EntitySeed EntitySeed;
    //! type of corresponding geometry
    typedef typename Traits::template Codim< codimension >::Geometry Geometry;
    //! type of corresponding local geometry
    typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;
    //! type of corresponding entity pointer
    typedef typename Traits::template Codim< codimension >::EntityPointer EntityPointer;

    //! type of hierarchic iterator
    typedef typename Traits::HierarchicIterator HierarchicIterator;
    //! type of leaf intersection iterator
    typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;
    //! type of level intersection iterator
    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;

    /** \} */
    
  private:
    typedef typename Traits::HostGrid   HostGrid;
    typedef typename Traits::ExtraData  ExtraData;

    // type of host grid info 
    typedef CartesianGridHostGridInfo< HostGrid > HostGridInfo;

  public:
    /** \name Host Types
     *  \{ */

    //! type of corresponding host entity
    typedef typename HostGrid::template Codim< codimension >::Entity HostEntity;
    //! type of corresponding host entity pointer
    typedef typename HostGrid::template Codim< codimension >::EntityPointer HostEntityPointer;
    /** \} */

  private:
    typedef typename HostGrid::template Codim< codimension >::Geometry HostGeometry;

    typedef CartesianGridGeometry< mydimension, dimensionworld, Grid > GeometryImpl;

    typedef CartesianGridLeafIntersectionIterator< Grid > LeafIntersectionIteratorImpl;
    typedef CartesianGridLevelIntersectionIterator< Grid > LevelIntersectionIteratorImpl;

  public:
    /** \name Construction, Initialization and Destruction
     *  \{ */

    /** \brief construct a null entity */
    CartesianGridEntity ( ExtraData data )
    : hostEntity_( 0 ),
      data_( data ),
      geo_( GeometryImpl( grid().geometricLevel( 0 ), 
            HostGridInfo::defaultDirection( mydimension ), 
            HostGridInfo::defaultOrigin() ) )
    {}

    /** \brief construct an initialized entity
     *
     *  \param[in]  hostEntity  corresponding entity in the host grid
     *
     *  \note The reference to the host entity must remain valid as long as
     *        this entity is in use.
     */
    explicit CartesianGridEntity ( ExtraData data, const HostEntity &hostEntity )
    : hostEntity_( &hostEntity ),
      data_( data ),
      geo_( GeometryImpl( grid().geometricLevel( level() ), 
                          HostGridInfo::defaultDirection( mydimension ), 
                          HostGridInfo::origin( hostEntity ) ) )
    {}

    /** \brief copy constructor */
    CartesianGridEntity ( const CartesianGridEntity &other )
    : hostEntity_( other.hostEntity_ ),
      data_( other.data_ ),
      geo_( Grid::getRealImplementation( other.geo_ ) )
    {}

    const CartesianGridEntity &operator= ( const CartesianGridEntity &other )
    {
      hostEntity_ = other.hostEntity_;
      Grid::getRealImplementation( geo_ ) = Grid::getRealImplementation( other.geo_ );
      data_ = other.data_ ;
      return *this;
    }

    /** \} */

    operator bool () const { return bool( hostEntity_ ); }

    /** \name Methods Shared by Entities of All Codimensions
     *  \{ */

    /** \brief obtain the name of the corresponding reference element
     *
     *  This type can be used to access the DUNE reference element.
     */
    GeometryType type () const
    {
      assert( hostEntity().type().isCube() );
      return hostEntity().type();
    }

    /** \brief obtain the level of this entity */
    int level () const
    {
      return hostEntity().level();
    }
    
    /** \brief obtain the partition type of this entity */
    PartitionType partitionType () const
    {
      return hostEntity().partitionType();
    }

    /** obtain the geometry of this entity
     *
     *  Each DUNE entity encapsulates a geometry object, representing the map
     *  from the reference element to world coordinates. Wrapping the geometry
     *  is the main objective of the CartesianGrid.
     *
     *  The CartesianGrid provides geometries of order 1, obtained by
     *  interpolation of its corners \f$y_i\f$. There corners are calculated
     *  from the corners \f$x_i\f$ of the host geometry through the
     *  CartesianGrid's coordinate function \f$c\f$, i.e.,
     *  \f$y_i = c( x_i )\f$.
     *
     *  \returns a const reference to the geometry
     */
    const Geometry &geometry () const
    {
      return geo_;
    }

    /** \brief return EntitySeed of host grid entity */
    EntitySeed seed () const { return EntitySeed( hostEntity().seed() ); }

    /** \} */

    template< int codim >
    int count () const
    {
      return hostEntity().template count< codim >();
    }
    
    template< int codim >
    typename Grid::template Codim< codim >::EntityPointer
    subEntity ( int i ) const
    {
      typedef typename Traits::template Codim< codim >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( data(), hostEntity().template subEntity< codim >( i ) );
    }

    LevelIntersectionIterator ilevelbegin () const
    {
      return LevelIntersectionIteratorImpl( data(), hostEntity().ilevelbegin() );
    }
    
    LevelIntersectionIterator ilevelend () const
    {
      return LevelIntersectionIteratorImpl( data(), hostEntity().ilevelend() );
    }
    
    LeafIntersectionIterator ileafbegin () const
    {
      return LeafIntersectionIteratorImpl( data(), hostEntity().ileafbegin() );
    }
    
    LeafIntersectionIterator ileafend () const
    {
      return LeafIntersectionIteratorImpl( data(), hostEntity().ileafend() );
    }

    bool hasBoundaryIntersections () const
    {
      return hostEntity().hasBoundaryIntersections();
    }

    bool isLeaf () const
    {
      return hostEntity().isLeaf();
    }
 
    EntityPointer father () const
    {
      typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;
      return EntityPointerImpl( data(), hostEntity().father() );
    }

    bool hasFather () const
    {
      return hostEntity().hasFather();
    }
      
    const LocalGeometry &geometryInFather () const
    {
      return *(grid().geometryInFather_[ HostGridInfo::childIndex( hostEntity() ) ]);
    }
 
    HierarchicIterator hbegin ( int maxLevel ) const
    {
      typedef CartesianGridHierarchicIteratorTraits< Grid > T;
      return CartesianGridIterator< T >( data(), hostEntity().hbegin( maxLevel ) );
    }
    
    HierarchicIterator hend ( int maxLevel ) const
    {
      typedef CartesianGridHierarchicIteratorTraits< Grid > T;
      return CartesianGridIterator< T >( data(), hostEntity().hend( maxLevel ) );
    }

    bool isRegular () const
    {
      return hostEntity().isRegular();
    }

    bool isNew () const
    {
      return hostEntity().isNew();
    }
          
    bool mightVanish () const
    {
      return hostEntity().mightVanish();
    }

    const HostEntity &hostEntity () const
    {
      assert( *this );
      return *hostEntity_;
    }

    const Grid& grid() const 
    {
      return *data();
    }

    ExtraData data() const 
    {
      assert( data_ );
      return data_ ;
    }

    /** \} */

  private:
    const HostEntity *hostEntity_;
    ExtraData data_;
    mutable Geometry geo_;
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_ENTITY_HH
