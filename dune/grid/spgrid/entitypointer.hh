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
    static const int dimension = Traits::Cube::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;

    typedef This EntityPointerImp;

  protected:
    typedef SPEntity< codimension, dimension, Grid > EntityImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

  public:
    SPEntityPointer ( const EntityInfo &entityInfo )
    : entity_( EntityImpl( entityInfo ) )
    {}

    SPEntityPointer ( const GridLevel &gridLevel, const MultiIndex &id )
    : entity_( EntityImpl( EntityInfo( gridLevel, id ) ) )
    {}

    SPEntityPointer ( const EntityImpl &entityImpl )
    : entity_( entityImpl )
    {}

    SPEntityPointer ( const This &other )
    : entity_( EntityImpl( Grid::getRealImplementation( other.dereference() ) ) )
    {}

    This &operator= ( const This &other )
    {
      Grid::getRealImplementation( entity_ )
        = Grid::getRealImplementation( other.dereference() );
      return *this;
    }

    void compactify ()
    {}

    Entity &dereference () const
    {
      return const_cast< Entity & >( entity_ );
    }

    bool equals ( const This &other ) const
    {
      return Grid::getRealImplementation( entity_ ).equals( Grid::getRealImplementation( other.entity_ ) );
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
