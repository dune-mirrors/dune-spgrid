#ifndef DUNE_SPGRID_ENTITYPOINTER_HH
#define DUNE_SPGRID_ENTITYPOINTER_HH

namespace Dune
{

  template< int codim, class Grid >
  class SPEntityPointer
  {
    typedef SPEntityPointer< codim, Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::Codim< codimension >::Entity Entity;

    SPEntityPointer ( const This &other )
    : entity_( other.dereference() )
    {}

    SPEntityPointer ( const Entity &entity )
    : entity_( entity )
    {}

    void compactify ()
    {}

    const Entity &dereference () const
    {
      assert( !!entity_ );
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

  private:
    Entity entity_;
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITYPOINTER_HH
