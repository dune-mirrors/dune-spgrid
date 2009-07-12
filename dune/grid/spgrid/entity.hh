#ifndef DUNE_SPGRID_ENTITY_HH
#define DUNE_SPGRID_ENTITY_HH

namespace Dune
{

  template< int codim, int dim, class Grid >
  class SPEntity;


  template< int codim, class Grid >
  class SPBasicEntity
  {
    typedef SPBasicEntity< codim, Grid > This;

  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;

    typedef typename Traits::Codim< 0 >::Geometry Geometry;

    typedef SPEntityInfo< ctype, dimension, codimension > EntityInfo;

    int level () const
    {
      return entityInfo().gridLevel().level();
    }

    bool isLeaf () const
    {
      return entityInfo().gridLevel().isLeaf();
    }

    PartitionType partitionType () const
    {
      return entityInfo().partitionType();
    }

    GeometryType type () const
    {
      return geometry().type();
    }

    const Geometry &geometry () const
    {
      return geometry_;
    }
  
    const EntityInfo &entityInfo () const
    {
      return Grid::getRealImplementation( geometry_ ).entityInfo_;
    }

  private:
    Geometry geometry_;
  };



  template< int codim, int dim, class Grid >
  class SPEntity
  : public SPBasicEntity< codim, Grid >
  {
    typedef SPEntity< codim, dim, Grid > This;
    typedef SPBasicEntity< codim, Grid > Base;
  };



  template< int dim, class Grid >
  class SPEntity
  : public SPBasicEntity< 0, Grid >
  {
    typedef SPEntity< 0, dim, Grid > This;
    typedef SPBasicEntity< 0, Grid > Base;

  protected:
    using Base::Traits;
    using Base::entityInfo;

  public:
    using Base::isLeaf;

  public:
    template< int codim >
    int count () const
    {
      // ...
    }

    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( const int i ) const
    {
      // ...
    }

    LeafIntersectionIterator ileafbegin () const
    {
      return (isLeaf() ? ilevelbegin() : ilevelend());
    }

    LeafIntersectionIterator ileafend () const
    {
      return ilevelend();
    }

    LevelIntersectionIterator ilevelbegin () const
    {
      // :..
    }

    LevelIntersectionIterator ilevelend () const
    {
      // ...
    }

    EntityPointer father () const
    {
      // ...
    }

    const LocalGeometry &geometryInFather () const
    {
      // ...
    }

    HierarchicIterator hbegin ( int maxlevel ) const
    {
      // ...
    }

    HierarchicIterator hend ( int maxlevel ) const
    {
      // ...
    }

    bool isRegular () const
    {
      return true;
    }

    bool isNew () const
    {
      return false;
    }

    bool mightVanish () const
    {
      return false;
    }

    bool hasBoundaryIntersections () const
    {
      // ...
    }
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITY_HH
