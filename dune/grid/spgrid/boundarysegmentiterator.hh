#ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
#define DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH

//- C++ includes
#include <cassert>

//- dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>

//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>

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

    struct Begin
    : public PartitionIterator::Begin
    { };

    struct End
    : public PartitionIterator::End
    { };

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
    IntersectionImpl &intersectionImp () { return Grid::getRealImplementation( intersection_ ); }
    // return intersection implementation
    const IntersectionImpl &intersectionImp () const { return Grid::getRealImplementation( intersection_ ); }

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
  : intersection_( Grid::getRealImplementation( other.intersection_ ) ),
    lastFace_( other.lastFace_ ),
    pit_( other.pit_ ),
    pend_( other.pend_ )
  { }


  template< class Grid >
  inline const typename SPBoundarySegmentIterator< Grid >::This &
  SPBoundarySegmentIterator< Grid >::operator= ( const This &other ) 
  {
    Grid::getRealImplementation( intersection_ ) = Grid::getRealImplementation( other.intersection_ );
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
    return Grid::getRealImplementation( intersection_ ).equals( Grid::getRealImplementation( other.intersection_ ) );
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
    intersectionImp() = IntersectionImp( entityInfo, face );
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
    return intersectionImp().gridLevel();
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

    // iterate over all partitions in grid level
    typedef typename PartitionList::Iterator Iterator;
    const PartitionList &plist = gridLevel().template partition< All_Partition >();
    const Iterator end = plist.end();
    for( Iterator it = plist.begin(); it != end; ++it )
    {
      // get partition
      const Partition partition = *it;

      // get vertex indices of left and right corner of face
      int subEntity[ 2 ];
      const GenericReferenceElement< ctype, dimension > &referenceElement 
        = GenericReferenceElements< ctype, dimension >::cube();
      for( int i = 0; i < 2; ++i )
        subEntity[ i ] = referenceElement.subEntity( face, 1, i*( (1<<(dimension-1))-1 ), dimension );

      // compute boundary partition
      MultiIndex bound[ 2 ];
      const MultiIndex direction = partition.end() - partition.begin();
      for( int i = 0; i < 2; ++i )
      {
        bound[ i ] = partition.begin();
        MultiIndex subId = gridLevel().referenceCube().subId( dimension, subEntity[ i ] );
        for( int j = 0; j < dimension; ++j )
          bound[ i ][ j ] += subId[ j ]*direction[ i ][ j ];
      }
      boundaryPartitions += Partition( bound[ 0 ], bound[ 1 ], partition.number() );
    }
    return boundaryPartitions;
  }

} // end namespace Dune

#endif // #ifndef DUNE_SPGRID_BOUNDARYSEGMENTITERATOR_HH
