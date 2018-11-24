#ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
#define DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH

//- C++ includes
#include <cassert>
#include <type_traits>

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>

//- local includes
#include <dune/grid/spgrid/intersection.hh>
#include <dune/grid/spgrid/iterator.hh>


namespace Dune
{

  // SPBoundarySegmentIterator
  // -------------------------

  template< class Grid >
  class SPBoundarySegmentIterator
  {
    // this type
    typedef SPBoundarySegmentIterator< Grid > This;
    // grid traits
    typedef typename std::remove_const< Grid >::type::Traits Traits;

    // intersection implementation
    typedef SPIntersection< Grid > IntersectionImpl;
    // partition iterator for faces
    typedef SPPartitionIterator< 1, Grid > PartitionIterator;

  public:
    //! \brief intersection type
    typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

    //! \brief entity type
    typedef typename Intersection::Entity Entity;

    //! \brief grid dimension
    static const int dimension = IntersectionImpl::dimension;
    //! \brief single coordinate type
    typedef typename IntersectionImpl::ctype ctype;
    //! \brief entity info type
    typedef typename IntersectionImpl::EntityInfo EntityInfo;
    //! \brief grid level
    typedef typename IntersectionImpl::GridLevel GridLevel;
    //! \brief partition list
    typedef typename GridLevel::PartitionList PartitionList;

    typedef typename PartitionIterator::Begin Begin;
    typedef typename PartitionIterator::End End;

  private:
    // number of faces
    static const int numFaces = GridLevel::ReferenceCube::numFaces;

    // partition
    typedef typename PartitionList::Partition Partition;
    // multi index
    typedef typename PartitionList::MultiIndex MultiIndex;

  public:
    //! \brief constructor
    SPBoundarySegmentIterator ( const GridLevel &gridLevel, int face, Begin );
    SPBoundarySegmentIterator ( const GridLevel &gridLevel, int face, End );

    //! \brief copy constructor
    SPBoundarySegmentIterator ( const This &other );
    //! \brief assignment operator
    const This &operator= ( const This &other );

    //! \brief check for equality
    bool equals ( const This &other ) const;

    //! \brief iterator increment
    void increment ();

    //! \brief dereference intersection
    const Intersection &dereference () const;

  protected:
    int face () const;
    const GridLevel &gridLevel () const;

    void init ( int face );

  private:
    // return intersection implementation
    IntersectionImpl &intersectionImpl () { return intersection_.impl(); }
    // return intersection implementation
    const IntersectionImpl &intersectionImpl () const { return intersection_.impl(); }

    Intersection intersection_;
    PartitionIterator pit_;
  };



  // Implementation of SPBoundarySegmentIterator
  // -------------------------------------------

  template< class Grid >
  inline SPBoundarySegmentIterator< Grid >
    ::SPBoundarySegmentIterator ( const GridLevel &gridLevel, int face, Begin )
  : intersection_( IntersectionImpl( EntityInfo( gridLevel ), face ) ),
    pit_( gridLevel, gridLevel.boundaryPartition( face ), Begin() )
  {
    init( face );
  }


  template< class Grid >
  inline SPBoundarySegmentIterator< Grid >
    ::SPBoundarySegmentIterator ( const GridLevel &gridLevel, int face, End )
  : intersection_( IntersectionImpl( EntityInfo( gridLevel ), face+1 ) ),
    pit_( gridLevel, gridLevel.boundaryPartition( face+1 ), Begin() )
  {}


  template< class Grid >
  inline SPBoundarySegmentIterator< Grid >::SPBoundarySegmentIterator ( const This &other )
  : intersection_( other.intersectionImpl() ),
    pit_( other.pit_ )
  {}


  template< class Grid >
  inline const typename SPBoundarySegmentIterator< Grid >::This &
  SPBoundarySegmentIterator< Grid >::operator= ( const This &other )
  {
    intersectionImpl() = other.intersectionImpl();
    pit_ = other.pit_;
    return *this;
  }


  template< class Grid >
  inline const typename SPBoundarySegmentIterator< Grid >::Intersection &
  SPBoundarySegmentIterator< Grid >::dereference () const
  {
    return intersection_;
  }


  template< class Grid >
  inline bool SPBoundarySegmentIterator< Grid >::equals ( const This &other ) const
  {
    return (pit_ == other.pit_);
  }


  template< class Grid >
  inline void SPBoundarySegmentIterator< Grid >::increment ()
  {
    // get current face
    int face = This::face();

    // try to increment internal iterator
    ++pit_;
    while( !pit_ && (face < numFaces) )
      pit_ = PartitionIterator( gridLevel(), gridLevel().boundaryPartition( ++face ), Begin() );

    init( face );
  }


  template< class Grid >
  inline int SPBoundarySegmentIterator< Grid >::face () const
  {
    return intersection_.indexInInside();
  }


  template< class Grid >
  inline const typename SPBoundarySegmentIterator< Grid >::GridLevel &
  SPBoundarySegmentIterator< Grid >::gridLevel () const
  {
    return intersectionImpl().gridLevel();
  }


  template< class Grid >
  inline void SPBoundarySegmentIterator< Grid >::init ( int face )
  {
    if( !pit_ )
      return;

    // compute id
    const SPEntity< 1, dimension, Grid > &entityImpl = ( *pit_ ).impl();
    MultiIndex id = entityImpl.entityInfo().id();
    const int i = face >> 1;
    const int j = 2*(face & 1) - 1;
    id[ i ] -= j;

    // update intersection
    const unsigned int number = entityImpl.entityInfo().partitionNumber();
    EntityInfo entityInfo = EntityInfo( gridLevel(), id, number );
    intersectionImpl() = IntersectionImpl( entityInfo, face );
  }

} // end namespace Dune

#endif // #ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
