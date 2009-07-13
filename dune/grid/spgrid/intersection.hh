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

    typedef FieldVector< ctype, dimensionworld > GlobalVector;
    typedef FieldVector< ctype, dimension-1 > Localvector;

    bool boundary () constt
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
      // ...
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
      // ...
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
      return unitOuterNormal( local );
    }

    GlobalVector integrationOuterNormal ( const LocalVector &local ) const
    {
      // ...
    }

    GlobalVector unitOuterNormal ( const LocalVector &local ) const
    {
      // ...
    }

  private:
    const Entity *inside_;
    unsigned int face_;
  };

}

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
