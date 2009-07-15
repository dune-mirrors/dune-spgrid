#ifndef DUNE_SPGRID_ENTITY_HH
#define DUNE_SPGRID_ENTITY_HH

#include <dune/grid/common/gridenums.hh>

#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int codim, int dim, class Grid >
  class SPEntity;



  // External Forward Declarations
  // -----------------------------

  template< class Grid >
  class SPIntersectionIterator;



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

    typedef typename Traits::template Codim< 0 >::Geometry Geometry;

    typedef SPEntityInfo< typename Traits::ctype, dimension, codimension > EntityInfo;
    typedef typename EntityInfo::GridLevel GridLevel;

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

    const GridLevel &gridLevel () const
    {
      return entityInfo().gridLevel();
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

  protected:
    typedef typename Base::Traits Traits;

  public:
    typedef typename Base::EntityInfo EntityInfo;
    typedef typename Base::GridLevel GridLevel;

    static const int dimension = Base::dimension;

    typedef typename Traits::template Codim< 0 >::Geometry Geometry;
    typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

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

    static const int numFaces = GridLevel::Cube::numFaces;

    typedef SPIntersectionIterator< Grid > IntersectionIteratorImpl;

  public:
    explicit SPEntity ( const EntityInfo &entityInfo )
    : Base( entityInfo )
    {}

  protected:
    using Base::entityInfo;
    using Base::gridLevel;

  public:
    using Base::isLeaf;

    template< int codim >
    int count () const
    {
      return gridLevel().cube().count( codim );
    }

    template< int codim >
    typename Codim< codim >::EntityPointer subEntity ( const int i ) const
    {
      MultiIndex id = entityInfo().id();
      id += gridLevel().cube().subId( codim, i );
      return Codim< codim >::EntityPointer( gridLevel(), id );
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
      return IntersectionIteratorImpl( *this, 0 );
    }

    LevelIntersectionIterator ilevelend () const
    {
      return IntersectionIteratorImpl( *this, numFaces );
    }

    bool hasBoundaryIntersections () const
    {
      const MultiIndex &id = entityInfo().id();
      const MultiIndex &cells = gridLevel().cells();

      bool hasBoundaryIntersections = false;
      for( int i = 0; i < dimension; ++i )
        hasBoundaryIntersections |= ((id[ i ] == 1) || (id[ i ] == 2*cells[ i ]-1));
      return hasBoundaryIntersections;
    }

    EntityPointer father () const
    {
      const MultiIndex &id = entityInfo().id();
      const unsigned int refDir = gridLevel().refinementDirection();
      MultiIndex fatherId;
      for( int i = 0; i < dimension; ++i )
        fatherId[ i ] = ((refDir >> i) & 1 ? (id[ i ] >> 1) | 1 : id[ i ]);
      return EntityPointer( EntityInfo( gridLevel().fatherLevel(), fatherId ) );
    }

    const LocalGeometry &geometryInFather () const
    {
      // ...
    }

    HierarchicIterator hbegin ( int maxlevel ) const
    {
      assert( maxlevel > level() );
      return HierarchicIteratorImpl( *this, maxlevel );
    }

    HierarchicIterator hend ( int maxlevel ) const
    {
      assert( maxlevel > level() );
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

}

#endif // #ifndef DUNE_SPGRID_ENTITY_HH
