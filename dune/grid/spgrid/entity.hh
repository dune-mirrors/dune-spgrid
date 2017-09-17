#ifndef DUNE_SPGRID_ENTITY_HH
#define DUNE_SPGRID_ENTITY_HH

#include <utility>

#include <dune/geometry/type.hh>

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

  template< class, int >
  class SPHierarchicIterator;



  // SPBasicEntity
  // -------------

  template< int codim, class Grid >
  class SPBasicEntity
  {
    typedef SPBasicEntity< codim, Grid > This;

  public:
    typedef __SPGrid::EntityInfo< Grid, codim > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;
    typedef typename EntityInfo::Traits Traits;

    static const int dimension = EntityInfo::dimension;
    static const int codimension = EntityInfo::codimension;
    static const int mydimension = EntityInfo::mydimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;
    typedef typename Traits::template Codim< codimension >::EntitySeed EntitySeed;
    typedef typename Traits::template Codim< codimension >::Geometry Geometry;
    typedef typename Traits::template Codim< codimension >::LocalGeometry LocalGeometry;

    typedef typename Traits::HierarchicIterator HierarchicIterator;

  protected:
    typedef SPGeometry< mydimension, dimension, Grid > GeometryImpl;
    typedef SPEntitySeed< codimension, Grid > EntitySeedImpl;

    typedef SPHierarchicIterator< Grid, codimension > HierarchicIteratorImpl;

  public:
    SPBasicEntity () {}

    explicit SPBasicEntity ( const EntityInfo &entityInfo )
    : entityInfo_( entityInfo )
    {}

    int level () const
    {
      return gridLevel().level();
    }

    bool isLeaf () const
    {
      return (level() == entityInfo().gridLevel().grid().maxLevel());
    }

    PartitionType partitionType () const
    {
      return entityInfo().partitionType();
    }

    GeometryType type () const { return GeometryTypes::cube( mydimension ); }

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
      return EntitySeed( EntitySeedImpl( level(), entityInfo().id(), entityInfo().partitionNumber() ) );
    }

    unsigned int subEntities ( unsigned int cd ) const
    {
      return gridLevel().template referenceCube< codim >().count( cd-codim );
    }

    bool hasFather () const { return ((level() > 0) && entityInfo().hasFather()); }

    Entity father () const
    {
      SPEntity< codimension, dimension, Grid > father( entityInfo() );
      father.entityInfo().up();
      return Entity( std::move( father ) );
    }

    HierarchicIterator hbegin ( int maxlevel ) const
    {
      assert( maxlevel >= level() );
      return HierarchicIteratorImpl( entityInfo(), maxlevel );
    }

    HierarchicIterator hend ( int maxlevel ) const
    {
      assert( maxlevel >= level() );
      return HierarchicIteratorImpl( entityInfo(), level() );
    }

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

    const Grid &grid () const { return gridLevel().grid(); }

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
    typedef typename Base::GridLevel GridLevel;

    typedef typename GridLevel::MultiIndex MultiIndex;

    SPEntity () {}

    explicit SPEntity ( const EntityInfo &entityInfo )
     : Base( entityInfo )
    {}

    SPEntity ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
      : Base( EntityInfo( gridLevel, id, partitionNumber ) )
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

    template< int codim >
    struct Codim
    {
      typedef typename Traits::template Codim< codim >::Entity Entity;
    };

    typedef typename GridLevel::MultiIndex MultiIndex;

  protected:
    typedef typename GridLevel::Mesh Mesh;

    static const int numFaces = GridLevel::ReferenceCube::numFaces;

  public:
    SPEntity () {}

    explicit SPEntity ( const EntityInfo &entityInfo )
     : Base( entityInfo )
    {}

    SPEntity ( const GridLevel &gridLevel, const MultiIndex &id, unsigned int partitionNumber )
      : Base( EntityInfo( gridLevel, id, partitionNumber ) )
    {}

    using Base::entityInfo;
    using Base::gridLevel;

    template< int codim >
    typename Codim< codim >::Entity subEntity ( int i ) const;

    bool hasBoundaryIntersections () const;

    LocalGeometry geometryInFather () const
    {
      return gridLevel().geometryInFather( entityInfo().id() );
    }

    bool isRegular () const { return true; }
    bool isNew () const { return false; }
    bool mightVanish () const { return false; }
  };



  // Implementation of SPEntity (for codimension 0)
  // ----------------------------------------------


  template< int dim, class Grid >
  template< int codim >
  inline typename SPEntity< 0, dim, Grid >::template Codim< codim >::Entity
  SPEntity< 0, dim, Grid >::subEntity ( int i ) const
  {
    typedef typename SPEntity< 0, dim, Grid >::template Codim< codim >::Entity SubEntity;
    typedef SPEntity< codim, dimension, Grid > SubEntityImpl;

    // warning: this is only true for closed partitions
    const unsigned int partitionNumber = entityInfo().partitionNumber();
    MultiIndex id = entityInfo().id();
    id += gridLevel().referenceCube().subId( codim, i );
    __SPGrid::EntityInfo< Grid, codim > subInfo( gridLevel(), id, partitionNumber );
    return SubEntity( SubEntityImpl( std::move( subInfo ) ) );
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

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_ENTITY_HH
