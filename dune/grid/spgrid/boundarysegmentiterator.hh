#ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
#define DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH

//- C++ includes
#include <cassert>

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
    typedef typename remove_const< Grid >::type::Traits Traits;

  public:
    //! \brief intersection type
    typedef Dune::Intersection< Grid, SPIntersection > Intersection;
    //! \brief entity type
    typedef typename Intersection::Entity Entity;
    
  private:
    // intersection implementation
    typedef SPIntersection< Grid > IntersectionImpl;
    // partition iterator for faces
    typedef SPPartitionIterator< 1, Grid > PartitionIterator;

  public:
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
    template< class BeginEnd >
    SPBoundarySegmentIterator ( const GridLevel &gridLevel, const BeginEnd &be, 
                                int firstFace = 0, int lastFace = numFaces );

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
    //! \brief return face
    int face () const;
    //! \brief return grid level
    const GridLevel &gridLevel () const;

  private:
    // initialize internal iterators
    template< class BeginEnd >
    void initialize ( const BeginEnd &be, const int face );

    // get partition list
    PartitionList partitionList ( const int face );

    // return intersection implementation
    IntersectionImpl &intersectionImpl () { return Grid::getRealImplementation( intersection_ ); }
    // return intersection implementation
    const IntersectionImpl &intersectionImpl () const { return Grid::getRealImplementation( intersection_ ); }

    Intersection intersection_;
    int lastFace_;
    PartitionIterator pit_, pend_;
  };



  // Implementation of SPBoundarySegmentIterator
  // -------------------------------------------

  template< class Grid >
  template< class BeginEnd >
  inline SPBoundarySegmentIterator< Grid >
    ::SPBoundarySegmentIterator ( const GridLevel &gridLevel, const BeginEnd &be, 
                                  int firstFace, int lastFace )
  : intersection_( IntersectionImpl( EntityInfo( gridLevel ), face ) ),
    lastFace_( lastFace )
  {  
    initialize( be, face );
  }


  template< class Grid >
  inline SPBoundarySegmentIterator< Grid >::SPBoundarySegmentIterator ( const This &other )
  : intersection_( other.intersectionImpl() ),
    lastFace_( other.lastFace_ ),
    pit_( other.pit_ ),
    pend_( other.pend_ )
  { }


  template< class Grid >
  inline const typename SPBoundarySegmentIterator< Grid >::This &
  SPBoundarySegmentIterator< Grid >::operator= ( const This &other ) 
  {
    intersectionImpl() = other.intersectionImpl();
    lastFace_ = other.lastFace_;
    pit_ = other.pit_;
    pend_ = other.pend_;
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
    assert( lastFace_ == other.lastFace_ );
    return intersectionImpl().equals( other.intersectionImpl() );
  }


  template< class Grid >
  inline void SPBoundarySegmentIterator< Grid >::increment ()
  {
    // get current face
    int face = This::face(); 

    // try to increment internal iterator
    ++pit_;
    if( pit_ == pend_ )
    {
      ++face;
      if( face > lastFace_ )
        return;
      initialize( face );
    }

    // compute id
    MultiIndex id = pit_->entityInfo().id();
    const int i = face >> 1;
    const int j = 2*(face & 1) - 1;
    id[ i ] -= j;

    // update intersection
    EntityInfo entityInfo = EntityInfo( gridLevel(), id );
    intersectionImpl() = IntersectionImp( entityInfo, face );
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
  template< class BeginEnd >
  inline void SPBoundarySegmentIterator< Grid >::initialize ( const BeginEnd &be, const int face )
  {
    PartitionList plist = partitionList( face );
    pit_ = PartitionIterator( gridLevel(), partitionList, be );
    pend_ = PartitionIterator( gridLevel(), partitionList, End() );
  }


  template< class Grid >
  inline typename SPBoundarySegmentIterator< Grid >::PartitionList
  SPBoundarySegmentIterator< Grid >::partitionList ( const int face )
  {
    // return value
    PartitionList boundaryPartitions;

    const MacroCube &globalCube = gridLevel().globalCube();
    const int i = face >> 1;

    // iterate over all partitions in grid level
    typedef typename PartitionList::Iterator Iterator;
    const PartitionList &plist = gridLevel().template partition< All_Partition >();
    const Iterator end = plist.end();
    for( Iterator it = plist.begin(); it != end; ++it )
    {
      // get partition
      const Partition partition = *it;

      // get partition bounds
      MultiIndex bound[ 2 ];
      bound[ 0 ] = partition.begin();
      bound[ 1 ] = partition.end();

      // shrink partition bounds to face bounds
      int bnd = (face & 1)*bound[ 0 ][ i ] + (1 - (face & 1))*bound[ 1 ][ i ];
      bound[ 0 ][ i ] = bnd;
      bound[ 1 ][ i ] = bnd;

      // insert partition iff it is part of the global boundary (see also intersection.hh)
      if( bnd == 2*globalCube.bound( face & 1 )[ i ] )
        boundaryPartitions += Partition( bound[ 0 ], bound[ 1 ], partition.number() );
    }
    return boundaryPartitions;
  }

} // end namespace Dune

#endif // #ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
