#ifndef DUNE_SPGRID_ENTITY_HH
#define DUNE_SPGRID_ENTITY_HH

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;



  // External Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid >
  class SPGeometry;



  // SPBasicEntity
  // -------------

  template< int codim, class Grid >
  class SPBasicEntity
  {
    typedef SPBasicEntity< codim, Grid > This;

  protected:
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    static const int dimension = Traits::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::Codim< 0 >::Geometry Geometry;

    typedef SPEntityInfo< ctype, dimension, codimension > EntityInfo;

  private:
    typedef SPGeometry< mydimension, dimension, Grid > GeometryImpl;

  public:
    explicit SPBasicEntity ( const EntityInfo &entityInfo )
    : geometry_( Geometry( GeometryImpl( entityInfo ) ) )
    {}

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

    bool equals ( const This &other ) const
    {
      return entityInfo().equals( other.entityInfo() );
    }
  
    const EntityInfo &entityInfo () const
    {
      return Grid::getRealImplementation( geometry_ ).entityInfo_;
    }

  private:
    Geometry geometry_;
  };



  // SPEntity
  // --------

  template< int codim, int dim, class Grid >
  class SPEntity
  : public SPBasicEntity< codim, Grid >
  {
    typedef SPEntity< codim, dim, Grid > This;
    typedef SPBasicEntity< codim, Grid > Base;

  public:
    typedef typename Base::EntityInfo EntityInfo;

    explicit Entity ( const EntityInfo &entityInfo )
    : Base( entityInfo )
    {}
  };



  // SPEntity (for codimension 0)
  // ----------------------------

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
    typedef typename Base::EntityInfo EntityInfo;

    explicit Entity ( const EntityInfo &entityInfo )
    : Base( entityInfo )
    {}

    using Base::isLeaf;

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

    bool hasBoundaryIntersections () const
    {
      const MultiIndex &midx = entityInfo().multiIndex();
      const MultiIndex &n = entityInfo().gridLevel().n();

      bool hasBoundaryIntersections = false;
      for( int i = 0; i < dimension; ++i )
        hasBoundaryIntersections |= ((midx[ i ] == 0) || (midx[ i ] == n[ i ]-1));
      return hasBoundaryIntersections;
    }

    EntityPointer father () const
    {
      return EntityPointer( entityInfo.father() );
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
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITY_HH
