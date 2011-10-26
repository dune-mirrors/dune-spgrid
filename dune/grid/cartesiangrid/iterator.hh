#ifndef DUNE_CARTESIANGRID_ITERATOR_HH
#define DUNE_CARTESIANGRID_ITERATOR_HH

#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/grid/cartesiangrid/entitypointer.hh>

namespace Dune
{

  // Iterator
  // --------
  
  template< class Traits >
  class CartesianGridIterator
  : public CartesianGridEntityPointer< Traits >
  {
    typedef CartesianGridEntityPointer< Traits > Base;
    typedef typename Traits::Grid Grid;

  protected:
    using Base::hostIterator_;
    using Base::releaseEntity;

    typedef typename Base::ExtraData ExtraData;
  public:
    typedef typename Base::HostIterator HostIterator;

    CartesianGridIterator( ExtraData data, 
                           const HostIterator &hostIterator )
    : Base( data, hostIterator )
    {}
    
    void increment ()
    {
      ++hostIterator_;
      releaseEntity();
    }
  };


  // CartesianGridLeafIteratorTraits
  // ------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct CartesianGridLeafIteratorTraits
  : public CartesianGridEntityPointerTraits< codim, Grid >
  {
    typedef typename remove_const< Grid >::type::Traits::HostGrid HostGrid;

    typedef typename HostGrid::template Codim< codim >
      ::template Partition< pitype >::LeafIterator
      HostIterator;

    static const PartitionIteratorType partitionType = pitype ;
  };


  // CartesianGridLevelIteratorTraits
  // -------------------------

  template< int codim, PartitionIteratorType pitype, class Grid >
  struct CartesianGridLevelIteratorTraits
  : public CartesianGridEntityPointerTraits< codim, Grid >
  {
    typedef typename remove_const< Grid >::type::Traits::HostGrid HostGrid;

    typedef typename HostGrid::template Codim< codim >
      ::template Partition< pitype >::LevelIterator
      HostIterator;

    static const PartitionIteratorType partitionType = pitype ;
  };
 

  // CartesianGridHierarchicIteratorTraits
  // ------------------------------

  template< class Grid >
  struct CartesianGridHierarchicIteratorTraits
  : public CartesianGridEntityPointerTraits< 0, Grid >
  {
    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::HostGrid::Traits::HierarchicIterator HostIterator;
  };

} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_ITERATOR_HH
