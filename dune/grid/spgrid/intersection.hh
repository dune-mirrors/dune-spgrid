#ifndef DUNE_SPGRID_INTERSECTION_HH
#define DUNE_SPGRID_INTERSECTION_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int, class >
  class SPEntity;

  template< int, class >
  class SPEntityPointer;



  // SPIntersection
  // --------------

  template< class Grid >
  class SPIntersection
  {
    typedef SPIntersection< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::Cube::ctype ctype;

    static const int dimension = Traits::Cube::dimension;
    static const int dimensionworld = Traits::Cube::dimensionworld;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  private:
    typedef SPEntity< 0, dimension, Grid > EntityImpl;
    typedef SPEntityPointer< 0, Grid > EntityPointerImpl;
    typedef SPGeometry< dimension-1, dimension, Grid > GeometryImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    typedef typename EntityInfo::GlobalVector GlobalVector;
    typedef typename GeometryImpl::LocalVector LocalVector;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    SPIntersection ( const EntityImpl &entityImpl, const int face )
    : inside_( &entityImpl ),
      geometry_( GeometryImpl( entityImpl.gridLevel()) )
    {
      setFace( face );
    }

    SPIntersection ( const This &other )
    : inside_( other.inside_ ),
      face_( other.face_ ),
      geometry_( GeometryImpl( Grid::getRealImplementation( other.geometry_ ) ) )
    {}

    This &operator= ( const This &other )
    {
      inside_ = other.inside_;
      face_ = other.face_;
      Grid::getRealImplementation( geometry_ ) = Grid::getRealImplementation( other.geometry_ );
      return *this;
    }

    bool boundary () const
    {
      return !neighbor();
    }

    int boundaryId () const
    {
      return 1;
    }

    bool neighbor () const
    {
      // this should be done much more efficiently
      MultiIndex id = inside_->entityInfo().id();
      id.axpy( 2, gridLevel().cube().subId( 1, face_ ) );
      bool neighbor = true;
      for( int i = 0; i < dimension; ++i )
        neighbor &= ((id[ i ] > 0) && (id[ i ] < 2*gridLevel().cells()[ i ]));
      return neighbor;
    }

    EntityPointer inside () const
    {
      return EntityPointer( *inside_ );
    }

    EntityPointer outside () const
    {
      MultiIndex id = inside_->entityInfo().id();
      id.axpy( 2, gridLevel().cube().subId( 1, face_ ) );
      return EntityPointerImpl( gridLevel(), id );
    }

    bool conforming () const
    {
      return true;
    }

    const LocalGeometry &geometryInInside () const
    {
      return gridLevel().grid().localFaceGeometry ( indexInInside() );
    }

    const LocalGeometry &geometryInOutside () const
    {
      return gridLevel().grid().localFaceGeometry ( indexInOutside() );
    }

    const Geometry &geometry () const
    {
      return geometry_;
    }

    GeometryType type () const
    {
      return GeometryType( GeometryType::cube, dimension-1 );
    }

    int indexInInside () const
    {
      return face_;
    }

    int indexInOutside () const
    {
      return face_ ^ 1;
    }

    GlobalVector outerNormal ( const LocalVector &local ) const
    {
      return integrationOuterNormal( local );
    }

    GlobalVector integrationOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().volumeNormal( face_ );
    }

    GlobalVector unitOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().cube().normal( face_ );
    }

    bool equals ( const This &other ) const
    {
      return (face_ == other.face_) && inside_->equals( *other.inside_ );
    }

    const GridLevel &gridLevel () const
    {
      return inside_->gridLevel();
    }

    void setFace ( const int face )
    {
      assert( face >= 0 );
      face_ = face;
      if( face < GridLevel::Cube::numFaces )
      {
        MultiIndex &id = Grid::getRealImplementation( geometry_ ).entityInfo().id();
        id = inside_->entityInfo().id();
        id += gridLevel().cube().subId( 1, face );
        Grid::getRealImplementation( geometry_ ).entityInfo().update();
      }
    }

  private:
    const EntityImpl *inside_;
    int face_;
    Geometry geometry_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
