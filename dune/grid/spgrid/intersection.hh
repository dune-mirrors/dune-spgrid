#ifndef DUNE_SPGRID_INTERSECTION_HH
#define DUNE_SPGRID_INTERSECTION_HH

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>
#include <dune/common/typetraits.hh>

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

    typedef SPEntityInfo< ctype, dimension, 0 > EntityInfo;

    typedef typename EntityInfo::GlobalVector GlobalVector;
    typedef FieldVector< ctype, dimension-1 > Localvector;

  private:
    typedef SPGeometry< dimension-1, dimension, Grid > GeometryImpl;

  public:
    bool boundary () const
    {
      // ...
    }

    int boundaryId () const
    {
      // ...
    }

    bool neighbor () const
    {
      // ...
    }

    EntityPointer inside () const
    {
      return EntityPointer( *inside_ );
    }

    EntityPointer outside () const
    {
      MultiIndex id = inside_->entityInfo().id();
      id.axpy( 2, gridLevel().cube().subId( 1, face ) );
      return EntityPointer( gridLevel(), id );
    }

    bool conforming () const
    {
      return true;
    }

    const LocalGeometry &geometryInInside () const
    {
      // ...
    }

    const LocalGeometry &geometryInOutside () const
    {
      // ...
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
      if( face < GridLevel::numFaces )
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
