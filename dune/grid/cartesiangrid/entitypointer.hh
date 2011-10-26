#ifndef DUNE_CARTESIANGRID_ENTITYPOINTER_HH
#define DUNE_CARTESIANGRID_ENTITYPOINTER_HH

#include <dune/grid/common/grid.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< class >
  class CartesianGrid;

  template< class >
  struct CartesianGridExportParams;

  template< int, int, class >
  class CartesianGridEntity;



  // EntityPointerTraits
  // -------------------
  
  template< int codim, class Grid >
  struct CartesianGridEntityPointerTraits;

  /** \cond */
  template< int codim, class Grid >
  struct CartesianGridEntityPointerTraits< codim, const Grid >
  : public CartesianGridEntityPointerTraits< codim, Grid >
  {};
  /** \endcond */

  template< int codim, class HostGrid >
  struct CartesianGridEntityPointerTraits< codim, CartesianGrid< HostGrid > >
  : public CartesianGridExportParams< HostGrid >
  {
    typedef Dune::CartesianGrid< HostGrid > Grid;

    typedef typename HostGrid::ctype ctype;

    static const int dimension = HostGrid::dimension;
    static const int codimension = codim;

    typedef Dune::Entity< codimension, dimension, const Grid, CartesianGridEntity > Entity;

    typedef typename HostGrid::template Codim< codim >::Entity HostEntity;
    typedef typename HostGrid::template Codim< codim >::EntityPointer HostEntityPointer;
    typedef HostEntityPointer HostIterator;
  };



  // CartesianGridEntityPointer
  // -------------------

  template< class Traits >
  class CartesianGridEntityPointer
  {
    typedef CartesianGridEntityPointer< Traits > This;

    typedef typename Traits::Grid Grid;
  
    typedef CartesianGridEntityPointerTraits< Traits::codimension, const Grid > BaseTraits;
    friend class CartesianGridEntityPointer< BaseTraits >;

  public:
    static const int dimension = Traits::dimension;
    static const int codimension = Traits::codimension;
 
    typedef typename Traits::Entity Entity;

    typedef CartesianGridEntityPointer< BaseTraits > EntityPointerImp;

  protected:
    typedef typename Traits::HostEntityPointer HostEntityPointer;
    typedef typename Traits::HostIterator HostIterator;
    typedef typename Grid::Traits::ExtraData  ExtraData; 

  private:
    typedef CartesianGridEntity< codimension, dimension, const Grid > EntityImpl;

  public:
    CartesianGridEntityPointer ( ExtraData data, 
                                const HostIterator &hostIterator )
    : entity_( EntityImpl( data ) ),
      hostIterator_( hostIterator )
    {}

    CartesianGridEntityPointer ( const EntityImpl &entity )
    : entity_( EntityImpl( entity.data() ) ),
      hostIterator_( entity.hostEntity() )
    {}

    CartesianGridEntityPointer ( const This &other )
    : entity_( EntityImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}

    template< class T >
    explicit CartesianGridEntityPointer ( const CartesianGridEntityPointer< T > &other )
    : entity_( EntityImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}
    
    This &operator= ( const This &other )
    {
      Grid::getRealImplementation( entity_ ) = EntityImpl( other.data() );
      hostIterator_ = other.hostIterator_;
      return *this;
    }

    operator const EntityPointerImp & () const
    {
      return reinterpret_cast< const EntityPointerImp & >( *this );
    }

    template< class T >
    bool equals ( const CartesianGridEntityPointer< T > &other ) const
    {
      return (hostIterator() == other.hostIterator());
    }
    
    Entity &dereference () const
    {
      if( ! entityImpl() )
        entityImpl() = EntityImpl( data(), *hostIterator() );
      return entity_;
    }
    
    int level () const
    {
      return hostIterator().level();
    }

    void compactify ()
    {
      releaseEntity();
      hostIterator_.compactify();
    }

    const HostIterator &hostIterator() const
    {
      return hostIterator_;
    }

  protected:
    void releaseEntity ()
    {
      entityImpl() = EntityImpl( data() );
    }

    ExtraData data() const { return entityImpl().data(); }

  private:
    EntityImpl &entityImpl () const
    {
      return Grid::getRealImplementation( entity_ );
    }

    mutable Entity entity_;

  protected:
    HostIterator hostIterator_;
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_ENTITYPOINTER_HH
