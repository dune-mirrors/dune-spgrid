#ifndef DUNE_SPGRID_INTERSECTION_HH
#define DUNE_SPGRID_INTERSECTION_HH

#include <type_traits>

#include <dune/common/typetraits.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/intersection.hh>

#include <dune/grid/spgrid/geometry.hh>
#include <dune/grid/spgrid/normal.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int, int, class >
  class SPEntity;



  // SPIntersection
  // --------------

  /**
   * \class SPIntersection
   *
   * \tparam  Grid  type of the \ref Dune::SPGrid "SPGrid" this intersection
   *                belongs to
   */
  template< class Grid >
  class SPIntersection
  {
    typedef SPIntersection< Grid > This;

    typedef typename std::remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::ReferenceCube ReferenceCube;

  public:
    typedef typename ReferenceCube::ctype ctype;

    static const int dimension = ReferenceCube::dimension;
    static const int mydimension = dimension-1;
    static const int dimensionworld = dimension;

    typedef typename Traits::template Codim< 0 >::Entity Entity;

    typedef typename Traits::template Codim< 1 >::Geometry Geometry;
    typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

  private:
    typedef SPEntity< 0, dimension, Grid > EntityImpl;
    typedef SPGeometry< mydimension, dimension, Grid > GeometryImpl;

  public:
    typedef typename EntityImpl::EntityInfo ElementInfo;
    typedef __SPGrid::EntityInfo< Grid, 1 > EntityInfo;

    typedef typename EntityImpl::GridLevel GridLevel;

    typedef typename GeometryImpl::LocalVector LocalVector;
    typedef SPNormalVector< ctype, dimensionworld > NormalVector;

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    typedef typename GridLevel::Mesh Mesh;
    typedef typename GridLevel::PartitionList PartitionList;
    typedef typename PartitionList::Partition Partition;

  public:
    SPIntersection () = default;

    SPIntersection ( const ElementInfo &insideInfo, int face )
      : normalId_( face ),
        insideInfo_( insideInfo )
    {}

    SPIntersection ( const EntityInfo &entityInfo, int face )
      : normalId_( face ),
        insideInfo_( entityInfo.gridLevel(), entityInfo.id() - normalId_, entityInfo.partitionNumber() )
    {}

    bool boundary () const
    {
      return (insideInfo_.id()[ normalId_.axis() ] + normalId_.sign() == 2*gridLevel().globalMesh().bound( normalId_ ));
    }

    int boundaryId () const
    {
      return (boundary() ? (indexInInside()+1) : 0);
    }

    std::size_t boundarySegmentIndex () const
    {
      assert( boundary() );
      return gridLevel().boundaryIndex( insideInfo_.id(), insideInfo_.partitionNumber(), normalId_.face() );
    }

    bool neighbor () const
    {
      const Partition &partition = gridLevel().template partition< All_Partition >().partition( insideInfo_.partitionNumber() );
      return ((insideInfo_.id()[ normalId_.axis() ] + normalId_.sign() != partition.bound( normalId_ )) || partition.hasNeighbor( indexInInside() ));
    }

    Entity inside () const { return Entity( EntityImpl( insideInfo_ ) ); }

    Entity outside () const;

    bool conforming () const
    {
      return true;
    }

    LocalGeometry geometryInInside () const
    {
      return gridLevel().grid().localFaceGeometry( indexInInside() );
    }

    LocalGeometry geometryInOutside () const
    {
      return gridLevel().grid().localFaceGeometry( indexInOutside() );
    }

    Geometry geometry () const
    {
      return Geometry( GeometryImpl( entityInfo() ) );
    }

    GeometryType type () const { return GeometryTypes::cube( mydimension ); }

    int indexInInside () const { return normalId_.face(); }
    int indexInOutside () const { return (-normalId_).face(); }

    NormalVector outerNormal ( const LocalVector &local ) const
    {
      return unitOuterNormal( local );
    }

    NormalVector integrationOuterNormal ( const LocalVector &local ) const
    {
      return gridLevel().faceVolume( indexInInside() ) * centerUnitOuterNormal();
    }

    NormalVector centerUnitOuterNormal () const { return normalId_; }

    NormalVector unitOuterNormal ( const LocalVector &local ) const { return normalId_; }

    bool equals ( const This &other ) const
    {
      return (indexInInside() == other.indexInInside()) && insideInfo_.equals( other.insideInfo_ );
    }

    const GridLevel &gridLevel () const { return insideInfo_.gridLevel(); }

    void setInside ( const ElementInfo &insideInfo ) { insideInfo_ = insideInfo; }

    void setEntityInfo ( const EntityInfo &entityInfo )
    {
      insideInfo_ = ElementInfo( entityInfo.gridLevel(), entityInfo.id() - normalId_, entityInfo.partitionNumber() );
    }

    EntityInfo entityInfo () const
    {
      return EntityInfo( gridLevel(), insideInfo_.id() + normalId_, insideInfo_.partitionNumber() );
    }

  private:
    SPNormalId< dimension > normalId_;
    ElementInfo insideInfo_;
  };



  // Implementation of SPIntersection
  // --------------------------------

  template< class Grid >
  inline typename SPIntersection< Grid >::Entity
  SPIntersection< Grid >::outside () const
  {
    MultiIndex id = insideInfo_.id() + normalId_ + normalId_;
    if( !boundary() )
      return Entity( EntityImpl( ElementInfo( gridLevel(), id, insideInfo_.partitionNumber() ) ) );

    const PartitionList &allPartition = gridLevel().template partition< All_Partition >();
    const Partition &partition = allPartition.partition( insideInfo_.partitionNumber() );

    assert( partition.hasNeighbor( indexInInside() ) );
    const Partition &nbPartition = allPartition.partition( partition.neighbor( indexInInside() ) );

    // manipulate id in case of periodic (i.e., transformed) neighbors
    const int bound = nbPartition.bound( -normalId_ );
    id[ normalId_.axis() ] = bound + normalId_.sign()*(1 - (bound & 1));
    return Entity( EntityImpl( ElementInfo( gridLevel(), id, nbPartition.number() ) ) );
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_INTERSECTION_HH
