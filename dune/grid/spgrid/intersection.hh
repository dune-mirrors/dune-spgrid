#ifndef DUNE_SPGRID_INTERSECTION_HH
#define DUNE_SPGRID_INTERSECTION_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  template< class Grid >
  class SPIntersection
  {
    typedef SPIntersection< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::ctype ctype;

    static const int dimension = Traits::dimension;
    static const int dimensionworld = Traits::dimensionworld;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

    typedef typename Entity::EntityInfo EntityInfo;
    typedef typename Entity::GridLevel GridLevel;

    typedef typename EntityInfo::GlobalVector GlobalVector;
    typedef FieldVector< ctype, dimension-1 > LocalVector;

  private:
    typedef SPGeometry< dimension-1, dimension, Grid > GeometryImpl;

    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    SPIntersection ( const Entity &entity, const unsigned int face )
    : inside_( &entity )
    {
      setFace( face );
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
      MultiIndex &id = inside_->entityInfo().id();
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
      return EntityPointer( gridLevel(), id );
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
      return gridLevel().volumeNormal();
    }

    GlobalVector unitOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().cube().normal( face_ );
    }

    bool equals ( const This &other ) const
    {
      return (*inside_ == *other.inside_) && (face_ == other.face_);
    }

    const GridLevel &gridLevel () const
    {
      return inside_->gridLevel();
    }

    void setFace ( const unsigned int face )
    {
      face_ = face;
      if( face < GridLevel::Cube::numFaces )
      {
        MultiIndex id = inside_->entityInfo().id();
        id += gridLevel().cube().subId( 1, face );
        geometry_ = Geometry( GeometryImpl( gridLevel(), id ) );
      }
    }

  private:
    const Entity *inside_;
    unsigned int face_;
    Geometry geometry_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
