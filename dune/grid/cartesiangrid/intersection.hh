#ifndef DUNE_CARTESIANGRID_INTERSECTION_HH
#define DUNE_CARTESIANGRID_INTERSECTION_HH

#include <dune/grid/cartesiangrid/entitypointer.hh>
#include <dune/grid/cartesiangrid/hostgridinfo.hh>

namespace Dune
{

  // External Forward Declataions
  // ----------------------------

  template< class >
  class CartesianGrid;



  // Internal Forward Declarations
  // -----------------------------
  
  template< class Grid >
  class CartesianGridLeafIntersection;
  
  template< class Grid >
  class CartesianGridLevelIntersection;
  


  // CartesianGridIntersection
  // ------------------

  template< class Grid, class HostIntersection >
  class CartesianGridIntersection
  {
    typedef typename HostIntersection::Geometry HostGeometry;

    typedef typename remove_const< Grid >::type::Traits Traits;
  protected:  
    typedef typename Traits :: ExtraData   ExtraData ;
    typedef typename Traits :: HostGrid    HostGrid ;

    typedef CartesianGridHostGridInfo< HostGrid > HostGridInfo;

  public:
    typedef typename Traits::ctype ctype;
    
    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;
    static const int mydimension = dimension-1;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;
    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

    typedef typename Traits::ReferenceCube::NormalVector NormalVector; 

  private:
    typedef typename Traits::template Codim< 0 >::EntityPointerImpl EntityPointerImpl;

    typedef CartesianGridGeometry< mydimension, dimensionworld, Grid > GeometryImpl;
    typedef typename remove_const< Grid >::type::GridLevel GridLevel;

  public:
    explicit CartesianGridIntersection ( ExtraData data )
    : hostIntersection_( 0 ),
      data_( data ),
      insideLevel_( -1 ),
      outsideLevel_( -1 ),
      gridLevel_( &grid().geometricLevel( 0 ) ),
      geo_( GeometryImpl( gridLevel(), 
                          HostGridInfo::direction( 0, dimension ),
                          HostGridInfo::defaultOrigin() ) )
    {}

    explicit CartesianGridIntersection ( ExtraData data, const HostIntersection &hostIntersection )
    : hostIntersection_( &hostIntersection ),
      data_( data ),
      insideLevel_( HostGridInfo::insideLevel( hostIntersection ) ),
      outsideLevel_( HostGridInfo::outsideLevel( hostIntersection, insideLevel_ ) ),
      gridLevel_( &grid().geometricLevel( level() ) ),
      geo_( GeometryImpl( gridLevel(), 
            HostGridInfo::direction( indexInInside(), dimension ), 
            HostGridInfo::originIntersection( hostIntersection ) ) )
    {}

    CartesianGridIntersection ( const CartesianGridIntersection &other )
    : hostIntersection_( other.hostIntersection_ ),
      data_( other.data_ ),
      insideLevel_( other.insideLevel_ ),
      outsideLevel_( other.outsideLevel_ ),
      gridLevel_( other.gridLevel_ ),
      geo_( Grid::getRealImplementation( other.geo_ ) )
    {}

    const CartesianGridIntersection &operator= ( const CartesianGridIntersection &other )
    {
      hostIntersection_ = other.hostIntersection_;
      data_ = other.data_;
      insideLevel_ = other.insideLevel_;
      outsideLevel_ = other.outsideLevel_;
      gridLevel_ = other.gridLevel_;
      Grid::getRealImplementation( geo_ ) = Grid::getRealImplementation( other.geo_ );
      return *this;
    }

    operator bool () const { return bool( hostIntersection_ ); }

    const EntityPointer inside () const
    {
      return EntityPointer( EntityPointerImpl( data(), hostIntersection().inside() ) );
    }
    
    EntityPointer outside () const
    {
      return EntityPointer( EntityPointerImpl( data(), hostIntersection().outside() ) );
    }

    bool boundary () const
    {
      return hostIntersection().boundary();
    }

    bool conforming () const
    {
      assert( !neighbor() || ((insideLevel_ == outsideLevel_) == hostIntersection().conforming()) );
      return hostIntersection().conforming();
    }
        
    bool neighbor () const
    {
      return hostIntersection().neighbor();
    }
        
    int boundaryId () const
    {
      return hostIntersection().boundaryId();
    }
        
    size_t boundarySegmentIndex () const
    {
      return hostIntersection().boundarySegmentIndex();
    }
        
    const LocalGeometry &geometryInInside () const
    {
      return localGeometry( indexInInside(), 
                            HostGridInfo::childIndexInInside( hostIntersection(),
                                                              insideLevel_,
                                                              outsideLevel_ ) );
    }
    
    const LocalGeometry &geometryInOutside () const
    {
      return localGeometry( indexInOutside(), 
                            HostGridInfo::childIndexInOutside( hostIntersection(),
                                                               insideLevel_,
                                                               outsideLevel_ ) );
    }

    const Geometry &geometry () const
    {
      return geo_;
    }

    GeometryType type () const
    {
      return hostIntersection().type();
    }

    int indexInInside () const
    {
      return hostIntersection().indexInInside();
    }
    
    int indexInOutside () const
    {
      assert( hostIntersection().indexInOutside() == (indexInInside()^1) );
      return hostIntersection().indexInOutside();
    }

    NormalVector
    integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return gridLevel().faceVolume( indexInInside() ) * centerUnitOuterNormal();
    }

    NormalVector
    outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return centerUnitOuterNormal();
    }

    NormalVector
    unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
    {
      return centerUnitOuterNormal();
    }

    NormalVector centerUnitOuterNormal () const
    {
      return grid().referenceCube().normal( indexInInside() );
    }

    const HostIntersection &hostIntersection () const
    {
      assert( *this );
      return *hostIntersection_;
    }

    const GridLevel &gridLevel () const
    {
      assert( gridLevel_ );
      return *gridLevel_;
    }

    const Grid &grid () const
    {
      return *data();
    }

    ExtraData data () const 
    { 
      assert( data_ );
      return data_ ;
    }

  private:
    int level () const
    {
      return std::max( insideLevel_, outsideLevel_ );
    }

    const LocalGeometry &localGeometry ( int faceIndex, int childIndex ) const
    {
      if( childIndex < 0 )
        return *(grid().localFaceGeometry_[ faceIndex ]);
      else
      {
        const int index = faceIndex*(1 << mydimension) + childIndex;
        return *(grid().localRefinedFaceGeometry_[ index ]);
      }
    }
   
    const HostIntersection *hostIntersection_;
    ExtraData data_;
    int insideLevel_, outsideLevel_;
    const GridLevel *gridLevel_;
    mutable Geometry geo_;
  };



  // CartesianGridLeafIntersection
  // -----------------------------

  template< class HostGrid >
  class CartesianGridLeafIntersection< const CartesianGrid< HostGrid > >
  : public CartesianGridIntersection< const CartesianGrid< HostGrid >, typename HostGrid::Traits::LeafIntersection >
  {
    typedef CartesianGrid< HostGrid > Grid;
    typedef typename HostGrid::Traits::LeafIntersection HostIntersection;

    typedef CartesianGridIntersection< const Grid, HostIntersection > Base;

    typedef typename Base :: ExtraData  ExtraData;
    typedef typename Base::EntityPointer EntityPointer;

  public:
    explicit CartesianGridLeafIntersection ( ExtraData data )
      : Base( data )
    {}

    CartesianGridLeafIntersection ( ExtraData data , const HostIntersection &hostIntersection )
    : Base( data, hostIntersection )
    {}
  };


  
  // CartesianGridLevelIntersection
  // -----------------------

  template< class HostGrid >
  class CartesianGridLevelIntersection< const CartesianGrid< HostGrid > >
  : public CartesianGridIntersection< const CartesianGrid< HostGrid >, typename HostGrid::Traits::LevelIntersection >
  {
    typedef CartesianGrid< HostGrid > Grid;
    typedef typename HostGrid::Traits::LevelIntersection HostIntersection;

    typedef CartesianGridIntersection< const Grid, HostIntersection > Base;

    typedef typename Base :: ExtraData  ExtraData;
    typedef typename Base::EntityPointer EntityPointer;

  public:
    explicit CartesianGridLevelIntersection ( ExtraData data )
      : Base( data )
    {}

    CartesianGridLevelIntersection ( ExtraData data, const HostIntersection &hostIntersection )
    : Base( data, hostIntersection )
    {}
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_INTERSECTION_HH
