#ifndef DUNE_SPGRID_ENTITY_HH
#define DUNE_SPGRID_ENTITY_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/entityseed.hh>
#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;



  // External Forward Declarations
  // -----------------------------

  template< int, class >
  class SPEntityPointer;

  template< class >
  class SPHierarchicIterator;

  template< class >
  class SPIntersectionIterator;



  // SPBasicEntity
  // -------------

  template< int codim, class Grid >
  class SPBasicEntity
  {
    typedef SPBasicEntity< codim, Grid > This;

  public:
    typedef SPEntityInfo< Grid, codim > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;
    typedef typename EntityInfo::Traits Traits;

    static const int dimension = EntityInfo::dimension;
    static const int codimension = EntityInfo::codimension;
    static const int mydimension = EntityInfo::mydimension;

    typedef typename Traits::template Codim< codim >::EntitySeed EntitySeed;
    typedef typename Traits::template Codim< codim >::Geometry Geometry;
    typedef typename Traits::template Codim< codim >::LocalGeometry LocalGeometry;

  protected:
    typedef SPGeometry< mydimension, dimension, Grid > GeometryImpl;
    typedef SPEntitySeed< codimension, Grid > EntitySeedImpl;

  public:
    explicit SPBasicEntity ( const EntityInfo &entityInfo )
    : entityInfo_( entityInfo )
    {}

    int level () const
    {
      return entityInfo().gridLevel().level();
    }

    bool isLeaf () const
    {
      return (level() == entityInfo().gridLevel().grid().maxLevel());
    }

    PartitionType partitionType () const
    {
      return entityInfo().partitionType();
    }

    GeometryType type () const
    {
      typedef typename GenericGeometry::CubeTopology< mydimension >::type Topology;
      return GeometryType( Topology() );
    }

    Geometry geometry () const
    {
      return Geometry( GeometryImpl( entityInfo_ ) );
    }

    bool equals ( const This &other ) const
    {
      return entityInfo().equals( other.entityInfo() );
    }

    EntitySeed seed () const
    {
      return EntitySeed( EntitySeedImpl( gridLevel().level(), entityInfo().id(), entityInfo().partitionNumber() ) );
    }
  
    const EntityInfo &entityInfo () const
    {
      return entityInfo_;
    }

    EntityInfo &entityInfo ()
    {
      return entityInfo_;
    }

    const GridLevel &gridLevel () const
    {
      return entityInfo().gridLevel();
    }

    const Grid &grid () const
    {
      return gridLevel().grid();
    }

  private:
    EntityInfo entityInfo_;
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

    explicit SPEntity ( const EntityInfo &entityInfo )
    : Base( entityInfo )
    {}
  };



  // SPEntity (for codimension 0)
  // ----------------------------

  template< int dim, class Grid >
  class SPEntity< 0, dim, Grid >
  : public SPBasicEntity< 0, Grid >
  {
    typedef SPEntity< 0, dim, Grid > This;
    typedef SPBasicEntity< 0, Grid > Base;

  public:
    typedef typename Base::EntityInfo EntityInfo;
    typedef typename Base::GridLevel GridLevel;
    typedef typename Base::Traits Traits;

    static const int dimension = Base::dimension;

    typedef typename Base::Geometry Geometry;
    typedef typename Base::LocalGeometry LocalGeometry;

    typedef typename Traits::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename Traits::LeafIntersectionIterator LeafIntersectionIterator;

    typedef typename Traits::HierarchicIterator HierarchicIterator;

    template< int codim >
    struct Codim
    {
      typedef typename Traits::template Codim< codim >::EntityPointer EntityPointer;
    };

    typedef typename Codim< 0 >::EntityPointer EntityPointer;

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    typedef typename GridLevel::Mesh Mesh;

    static const int numFaces = GridLevel::ReferenceCube::numFaces;

    typedef SPHierarchicIterator< Grid > HierarchicIteratorImpl;
    typedef SPIntersectionIterator< Grid > IntersectionIteratorImpl;

  public:
    explicit SPEntity ( const EntityInfo &entityInfo )
    : Base( entityInfo )
    {}

    using Base::entityInfo;
    using Base::gridLevel;
    using Base::isLeaf;
    using Base::level;

    template< int codim >
    int count () const
    {
      return gridLevel().referenceCube().count( codim );
    }

    unsigned int subEntities ( unsigned int codim ) const
    {
      return gridLevel().referenceCube().count( codim );
    }

    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( const int i ) const;

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
      return IntersectionIteratorImpl( entityInfo(), 0 );
    }

    LevelIntersectionIterator ilevelend () const
    {
      return IntersectionIteratorImpl( entityInfo(), numFaces );
    }

    bool hasBoundaryIntersections () const;

    bool hasFather () const
    {
      return (level() > 0);
    }

    EntityPointer father () const;

    LocalGeometry geometryInFather () const
    {
      return gridLevel().geometryInFather( entityInfo().id() );
    }

    HierarchicIterator hbegin ( int maxlevel ) const
    {
      assert( maxlevel >= level() );
      return HierarchicIteratorImpl( *this, maxlevel );
    }

    HierarchicIterator hend ( int maxlevel ) const
    {
      assert( maxlevel >= level() );
      return HierarchicIteratorImpl( *this, level() );
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



  // Implementation of SPEntity (for codimension 0)
  // ----------------------------------------------


  template< int dim, class Grid >
  template< int codim >
  inline typename SPEntity< 0, dim, Grid >::template Codim< codim >::EntityPointer
  SPEntity< 0, dim, Grid >::subEntity ( const int i ) const
  {
    // warning: this is only true for closed partitions
    const unsigned int partitionNumber = entityInfo().partitionNumber();
    MultiIndex id = entityInfo().id();
    id += gridLevel().referenceCube().subId( codim, i );
    return SPEntityPointer< codim, Grid >( gridLevel(), id, partitionNumber );
  }


  template< int dim, class Grid >
  inline bool SPEntity< 0, dim, Grid >::hasBoundaryIntersections () const
  {
    const Mesh &globalMesh = gridLevel().globalMesh();
    const MultiIndex &id = entityInfo().id();

    bool hasBoundaryIntersections = false;
    for( int i = 0; i < dimension; ++i )
    {
      hasBoundaryIntersections |= (id[ i ] == 2*globalMesh.begin()[ i ] + 1);
      hasBoundaryIntersections |= (id[ i ] == 2*globalMesh.end()[ i ] - 1);
    }
    return hasBoundaryIntersections;
  }


  template< int dim, class Grid >
  inline typename SPEntity< 0, dim, Grid >::EntityPointer
  SPEntity< 0, dim, Grid >::father () const
  {
    EntityInfo fatherInfo( entityInfo() );
    fatherInfo.up();
    fatherInfo.update();
    return EntityPointer( fatherInfo );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITY_HH
