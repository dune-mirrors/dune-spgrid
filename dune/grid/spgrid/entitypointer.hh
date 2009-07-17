#ifndef DUNE_SPGRID_ENTITYPOINTER_HH
#define DUNE_SPGRID_ENTITYPOINTER_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/spgrid/entity.hh>

namespace Dune
{

  // External Forward Declarations

  template< int codim, int dim, class Grid >
  class SPEntity;



  // SPEntityPointer
  // ---------------

  template< int codim, class Grid >
  class SPEntityPointer
  {
    typedef SPEntityPointer< codim, Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;

    typedef typename Entity::EntityInfo EntityInfo;
    typedef typename Entity::GridLevel GridLevel;

  protected:
    typedef SPEntity< codimension, dimension, Grid > EntityImpl;
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    SPEntityPointer ( const EntityInfo &entityInfo )
    : entity_( EntityImpl( entityInfo ) )
    {}

    SPEntityPointer ( const GridLevel &gridLevel, const MultiIndex &id )
    : entity_( EntityImpl( EntityInfo( gridLevel, id ) ) )
    {}

    SPEntityPointer ( const Entity &entity )
    : entity_( entity )
    {}

    SPEntityPointer ( const This &other )
    : entity_( other.dereference() )
    {}

    void compactify ()
    {}

    const Entity &dereference () const
    {
      return entity_;
    }

    bool equals ( const This &other ) const
    {
      return entity_.equals( other.entity_ );
    }

    int level () const
    {
      return entity_.level();
    }

    const GridLevel &gridLevel () const
    {
      return Grid::getRealImplementation( entity_ ).gridLevel();
    }

  protected:
    Entity entity_;
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITYPOINTER_HH
