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
    SPEntityPointer () = default;
    SPEntityPointer ( const EntityInfo &entityInfo ) : entityInfo_( entityInfo ) {}

    SPEntityPointer ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
      : entityInfo_( gridLevel, id, partitionNumber )
    {}

    SPEntityPointer ( const EntityImpl &entityImpl ) : entityInfo_( entityImpl.entityInfo() ) {}

    Entity operator* () const { return dereference(); }

    int level () const { return gridLevel().level(); }

    Entity dereference () const { return EntityImpl( entityInfo() ); }

    bool equals ( const This &other ) const { return entityInfo().equals( other.entityInfo() ); }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

  protected:
    EntityInfo entityInfo_;
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITYPOINTER_HH
