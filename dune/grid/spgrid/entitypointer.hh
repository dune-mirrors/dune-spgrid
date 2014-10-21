#ifndef DUNE_SPGRID_ENTITYPOINTER_HH
#define DUNE_SPGRID_ENTITYPOINTER_HH

#include <type_traits>

#include <dune/grid/spgrid/entity.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;



  // SPEntityPointer
  // ---------------

  template< int codim, class Grid >
  class SPEntityPointer
  {
    typedef SPEntityPointer< codim, Grid > This;

  public:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

    static const int dimension = Traits::ReferenceCube::dimension;
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
    SPEntityPointer () : entity_( EntityImpl( EntityInfo() ) ) {}

    SPEntityPointer ( const EntityInfo &entityInfo );
    SPEntityPointer ( const GridLevel &gridLevel, const MultiIndex &id,
                      const unsigned int partitionNumber );
    SPEntityPointer ( const EntityImpl &entityImpl );
    SPEntityPointer ( const This &other );

    const This &operator= ( const This &other );

    const Entity &operator* () const;
    const Entity *operator-> () const;

    bool operator== ( const This &other ) const;
    bool operator!= ( const This &other ) const;

    const Entity &dereference () const { return entity_; }

    bool equals ( const This &other ) const;

    int level () const;
    const GridLevel &gridLevel () const;

    const EntityInfo &entityInfo () const { return Grid::getRealImplementation( entity_ ).entityInfo(); }
    EntityInfo &entityInfo () { return Grid::getRealImplementation( entity_ ).entityInfo(); }

  private:
    Entity entity_;
  };



  // Implementation of SPEntityPointer
  // ---------------------------------

  template< int codim, class Grid >
  inline SPEntityPointer< codim, Grid >
    ::SPEntityPointer ( const EntityInfo &entityInfo )
  : entity_( EntityImpl( entityInfo ) )
  {}


  template< int codim, class Grid >
  inline SPEntityPointer< codim, Grid >
    ::SPEntityPointer ( const GridLevel &gridLevel, const MultiIndex &id,
                        const unsigned int partitionNumber )
  : entity_( EntityImpl( EntityInfo( gridLevel, id, partitionNumber ) ) )
  {}


  template< int codim, class Grid >
  inline SPEntityPointer< codim, Grid >
    ::SPEntityPointer ( const EntityImpl &entityImpl )
  : entity_( entityImpl )
  {}


  template< int codim, class Grid >
  inline SPEntityPointer< codim, Grid >::SPEntityPointer ( const This &other )
  : entity_( EntityImpl( Grid::getRealImplementation( other.dereference() ) ) )
  {}


  template< int codim, class Grid >
  inline const typename SPEntityPointer< codim, Grid >::This &
  SPEntityPointer< codim, Grid >::operator= ( const This &other )
  {
    Grid::getRealImplementation( entity_ )
      = Grid::getRealImplementation( other.dereference() );
    return *this;
  }


  template< int codim, class Grid >
  inline const typename SPEntityPointer< codim, Grid >::Entity &
  SPEntityPointer< codim, Grid >::operator* () const
  {
    return entity_;
  }


  template< int codim, class Grid >
  inline const typename SPEntityPointer< codim, Grid >::Entity *
  SPEntityPointer< codim, Grid >::operator-> () const
  {
    return &entity_;
  }


  template< int codim, class Grid >
  inline bool
  SPEntityPointer< codim, Grid >::operator== ( const This &other ) const
  {
    return equals( other );
  }

  template< int codim, class Grid >
  inline bool
  SPEntityPointer< codim, Grid >::operator!= ( const This &other ) const
  {
    return !equals( other );
  }


  template< int codim, class Grid >
  inline bool
  SPEntityPointer< codim, Grid >::equals ( const This &other ) const
  {
    return Grid::getRealImplementation( entity_ ).equals( Grid::getRealImplementation( other.entity_ ) );
  }


  template< int codim, class Grid >
  inline int SPEntityPointer< codim, Grid >::level () const
  {
    return entity_.level();
  }


  template< int codim, class Grid >
  inline const typename SPEntityPointer< codim, Grid >::GridLevel &
  SPEntityPointer< codim, Grid >::gridLevel () const
  {
    return Grid::getRealImplementation( entity_ ).gridLevel();
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYPOINTER_HH
