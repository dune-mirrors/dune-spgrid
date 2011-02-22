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

  /** \class SPIntersection
   *
   *  \tparam  Grid  type of the \ref Dune::SPGrid "SPGrid" this intersection
   *                 belongs to
   *
   *  \note The SPIntersection stores a pointer to the inside entity. It is
   *        assumed that this pointer remains valid during the entire life
   *        time of the intersection.
   */
  template< class Grid >
  class SPIntersection
  {
    typedef SPIntersection< Grid > This;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ReferenceCube ReferenceCube;

  public:
    typedef typename ReferenceCube::ctype ctype;

    static const int dimension = ReferenceCube::dimension;
    static const int mydimension = dimension-1;
    static const int dimensionworld = dimension;

    typedef typename Traits::template Codim< 0 >::Entity Entity;
    typedef typename Traits::template Codim< 0 >::EntityPointer EntityPointer;

    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  private:
    typedef SPEntity< 0, dimension, Grid > EntityImpl;
    typedef SPEntityPointer< 0, Grid > EntityPointerImpl;
    typedef SPGeometry< mydimension, dimension, Grid > GeometryImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    typedef typename GeometryImpl::LocalVector LocalVector;
    typedef typename ReferenceCube::NormalVector NormalVector;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    typedef typename GridLevel::Mesh Mesh;
    typedef typename GridLevel::PartitionList PartitionList;
    typedef typename PartitionList::Partition Partition;

  public:
    SPIntersection ( const EntityImpl &entityImpl, const int face )
    : inside_( &entityImpl ),
      geometry_( GeometryImpl( entityImpl.gridLevel() ) )
    {
      setFace( face );
    }

    SPIntersection ( const This &other )
    : inside_( other.inside_ ),
      face_( other.face_ ),
      geometry_( GeometryImpl( Grid::getRealImplementation( other.geometry_ ) ) )
    {}

    const This &operator= ( const This &other )
    {
      inside_ = other.inside_;
      face_ = other.face_;
      Grid::getRealImplementation( geometry_ ) = Grid::getRealImplementation( other.geometry_ );
      return *this;
    }

    bool boundary () const;

    int boundaryId () const
    {
      return (boundary() ? (indexInInside()+1) : 0);
    }

    size_t boundarySegmentIndex () const;

    bool neighbor () const;

    EntityPointer inside () const
    {
      //return EntityPointer( *inside_ );
      return EntityPointer( entityInfo() );
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
      typedef typename GenericGeometry::CubeTopology< mydimension >::type Topology;
      return GeometryType( Topology() );
    }

    int indexInInside () const
    {
      return face_;
    }

    int indexInOutside () const
    {
      return face_ ^ 1;
    }

    NormalVector outerNormal ( const LocalVector &local ) const
    {
      return integrationOuterNormal( local );
    }

    NormalVector integrationOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().faceVolume( face_ ) * centerUnitOuterNormal();
    }

    NormalVector centerUnitOuterNormal () const
    {
      return gridLevel().referenceCube().normal( face_ );
    }

    NormalVector unitOuterNormal ( const LocalVector &local ) const
    {
      return centerUnitOuterNormal();
    }

    bool equals ( const This &other ) const
    {
      return (face_ == other.face_) && entityInfo().equals( other.entityInfo() );
    }

    const GridLevel &gridLevel () const
    {
      return entityInfo().gridLevel();
    }

    void setFace ( const int face )
    {
      assert( face >= 0 );
      face_ = face;
      if( face < ReferenceCube::numFaces )
      {
        const unsigned int partitionNumber = entityInfo().partitionNumber();
        MultiIndex &id = Grid::getRealImplementation( geometry_ ).entityInfo().id();
        id = entityInfo().id();
        id += gridLevel().referenceCube().subId( 1, face );
        Grid::getRealImplementation( geometry_ ).entityInfo().update( partitionNumber );
      }
    }

  private:
    const EntityInfo &entityInfo () const
    {
      return inside_->entityInfo();
    }

    const EntityImpl *inside_;
    int face_;
    Geometry geometry_;
  };



  // Implementation of SPIntersection
  // --------------------------------

  template< class Grid >
  inline bool SPIntersection< Grid >::boundary () const
  {
    const Mesh &globalMesh = gridLevel().globalMesh();
    const MultiIndex &id = entityInfo().id();
    const int i = face_ >> 1;
    const int j = 2*(face_ & 1) - 1;
    return (id[ i ] + j == 2*globalMesh.bound( face_ & 1 )[ i ]);
  }


  template< class Grid >
  inline size_t SPIntersection< Grid >::boundarySegmentIndex () const
  {
    assert( boundary() );
    const MultiIndex &id = entityInfo().id();
    const unsigned int partitionNumber = entityInfo().partitionNumber();
    return gridLevel().boundaryIndex( id, partitionNumber, face_ );
  }


  template< class Grid >
  inline bool SPIntersection< Grid >::neighbor () const
  {
    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( entityInfo().partitionNumber() );

    const MultiIndex &id = entityInfo().id();
    const int i = face_ >> 1;
    const int j = 2*(face_ & 1) - 1;
    return (partition.hasNeighbor( face_ ) || (j*id[ i ] + 2 <= j*partition.bound( face_ & 1 )[ i ]));
  }


  template< class Grid >
  inline typename SPIntersection< Grid >::EntityPointer
  SPIntersection< Grid >::outside () const
  {
    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( entityInfo().partitionNumber() );

    MultiIndex id = entityInfo().id();
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

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
