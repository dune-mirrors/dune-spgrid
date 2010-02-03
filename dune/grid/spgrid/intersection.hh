#ifndef DUNE_SPGRID_INTERSECTION_HH
#define DUNE_SPGRID_INTERSECTION_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/grid/spgrid/geometry.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int, class >
  class SPEntity;

  template< int, class >
  class SPEntityPointer;



  // SPIntersection
  // --------------

  template< class Grid >
  class SPIntersection
  {
    typedef SPIntersection< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    typedef typename Traits::Cube::ctype ctype;

    static const int dimension = Traits::Cube::dimension;
    static const int dimensionworld = Traits::Cube::dimensionworld;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  private:
    typedef SPEntity< 0, dimension, Grid > EntityImpl;
    typedef SPEntityPointer< 0, Grid > EntityPointerImpl;
    typedef SPGeometry< dimension-1, dimension, Grid > GeometryImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    typedef typename EntityInfo::GlobalVector GlobalVector;
    typedef typename GeometryImpl::LocalVector LocalVector;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;
    typedef typename GridLevel::Domain Domain;

    typedef typename GridLevel::PartitionList PartitionList;
    typedef typename PartitionList::Partition Partition;

  public:
    SPIntersection ( const EntityImpl &entityImpl, const int face )
    : inside_( &entityImpl ),
      geometry_( GeometryImpl( entityImpl.gridLevel()) )
    {
      setFace( face );
    }

    SPIntersection ( const This &other )
    : inside_( other.inside_ ),
      face_( other.face_ ),
      geometry_( GeometryImpl( Grid::getRealImplementation( other.geometry_ ) ) )
    {}

    This &operator= ( const This &other )
    {
      inside_ = other.inside_;
      face_ = other.face_;
      Grid::getRealImplementation( geometry_ ) = Grid::getRealImplementation( other.geometry_ );
      return *this;
    }

    bool boundary () const;

    int boundaryId () const
    {
      return 1;
    }

    size_t boundarySegmentIndex () const
    {
      assert( boundary() );
      return gridLevel().boundaryIndex( inside_->entityInfo().id(), face_ );
    }

    bool neighbor () const;

    EntityPointer inside () const
    {
      return EntityPointer( *inside_ );
    }

    EntityPointer outside () const;

    bool conforming () const
    {
      return true;
    }

    const LocalGeometry &geometryInInside () const
    {
      return gridLevel().grid().localFaceGeometry( indexInInside() );
    }

    const LocalGeometry &geometryInOutside () const
    {
      return gridLevel().grid().localFaceGeometry( indexInOutside() );
    }

    const Geometry &geometry () const
    {
      return geometry_;
    }

    GeometryType type () const
    {
      return GeometryType( GeometryType::cube, dimension-1 );
    }

    int indexInInside () const
    {
      return face_;
    }

    int indexInOutside () const
    {
      return face_ ^ 1;
    }

    GlobalVector outerNormal ( const LocalVector &local ) const
    {
      return integrationOuterNormal( local );
    }

    GlobalVector integrationOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().volumeNormal( face_ );
    }

    GlobalVector centerUnitOuterNormal () const
    {
      return gridLevel().cube().normal( face_ );
    }

    GlobalVector unitOuterNormal ( const LocalVector &local ) const
    {
      return centerUnitOuterNormal();
    }

    bool equals ( const This &other ) const
    {
      return (face_ == other.face_) && inside_->equals( *other.inside_ );
    }

    const GridLevel &gridLevel () const
    {
      return inside_->gridLevel();
    }

    void setFace ( const int face )
    {
      assert( face >= 0 );
      face_ = face;
      if( face < GridLevel::Cube::numFaces )
      {
        MultiIndex &id = Grid::getRealImplementation( geometry_ ).entityInfo().id();
        id = inside_->entityInfo().id();
        id += gridLevel().cube().subId( 1, face );
        Grid::getRealImplementation( geometry_ ).entityInfo().update();
      }
    }

  private:
    const EntityImpl *inside_;
    int face_;
    Geometry geometry_;
  };



  // Implementation of SPIntersection
  // --------------------------------

  template< class Grid >
  inline bool SPIntersection< Grid >::boundary () const
  {
    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( inside_->entityInfo().partitionNumber() );

    const MultiIndex &id = inside_->entityInfo().id();
    const int i = face_ >> 1;
    const int j = 2*(face_ & 1) - 1;
    return partition.boundary( face_ ) && (id[ i ] + j == partition.bound( face_ & 1 )[ i ]);
  }


  template< class Grid >
  inline bool SPIntersection< Grid >::neighbor () const
  {
    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( inside_->entityInfo().partitionNumber() );

    const MultiIndex &id = inside_->entityInfo().id();
    const int i = face_ >> 1;
    const int j = 2*(face_ & 1) - 1;
    return (partition.hasNeighbor( face_ ) || (j*id[ i ] + 2 <= j*partition.bound( face_ & 1 )[ i ]));
  }


  template< class Grid >
  inline typename SPIntersection< Grid >::EntityPointer
  SPIntersection< Grid >::outside () const
  {
    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( inside_->entityInfo().partitionNumber() );

    MultiIndex id = inside_->entityInfo().id();
    const int i = face_ >> 1;
    const int j = 2*(face_ & 1) - 1;
    id[ i ] += 2*j;
    if( (j*id[ i ] <= j*partition.bound( face_ & 1 )[ i ]) )
      return EntityPointerImpl( gridLevel(), id, partition.number() );

    assert( partition.hasNeighbor( face_ ) );
    const int bound = allPartition.partition( partition.neighbor( face_ ) ).bound( 1 - (face_ & 1) )[ i ];
    id[ i ] = bound + j*(1-(bound & 1));
    return EntityPointerImpl( gridLevel(), id, partition.neighbor( face_ ) );
  }

}

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
