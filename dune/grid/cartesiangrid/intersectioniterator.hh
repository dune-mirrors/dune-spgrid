#ifndef DUNE_CARTESIANGRID_INTERSECTIONITERATOR_HH
#define DUNE_CARTESIANGRID_INTERSECTIONITERATOR_HH

#include <dune/grid/cartesiangrid/entitypointer.hh>
#include <dune/grid/cartesiangrid/intersection.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------
  
  template< class Grid >
  class CartesianGridLeafIntersectionIterator;

  template< class Grid >
  class CartesianGridLevelIntersectionIterator;



  // CartesianGridIntersectionIterator
  // --------------------------

  template< class Traits >
  class CartesianGridIntersectionIterator
  {
    typedef typename Traits::HostIntersectionIterator HostIntersectionIterator;

  public:
    typedef typename Traits::Intersection Intersection;
    typedef typename Traits::GridTraits::Grid Grid;

    typedef typename Grid::template Codim< 0 >::EntityPointer EntityPointer;
   
  private:
    typedef typename Traits::GridTraits::ExtraData  ExtraDataType ;
    typedef typename Traits::IntersectionImpl IntersectionImpl;

  public:
    CartesianGridIntersectionIterator ( ExtraDataType data, 
                                       const HostIntersectionIterator &hostIterator )
    : intersection_( IntersectionImpl( data ) ),
      hostIterator_( hostIterator )
    {}

    CartesianGridIntersectionIterator ( const CartesianGridIntersectionIterator &other )
    : intersection_( IntersectionImpl( other.data() ) ),
      hostIterator_( other.hostIterator_ )
    {}

    const CartesianGridIntersectionIterator &
    operator= ( const CartesianGridIntersectionIterator &other )
    {
      intersectionImpl() = IntersectionImpl( data() );
      hostIterator_ = other.hostIterator_;
      return *this;
    }

    bool equals ( const CartesianGridIntersectionIterator &other ) const
    {
      return (hostIterator_ == other.hostIterator_);
    }
    
    void increment ()
    {
      ++hostIterator_;
      intersectionImpl() = IntersectionImpl( data() );
    }

    const Intersection &dereference () const
    {
      if( !intersectionImpl() )
        intersectionImpl() = IntersectionImpl( data(), *hostIterator_ );
      return intersection_;
    }

  private:
    ExtraDataType data () const { return intersectionImpl().data(); }

    IntersectionImpl &intersectionImpl () const
    {
      return Grid::getRealImplementation( intersection_ );
    }

    mutable Intersection intersection_;
    HostIntersectionIterator hostIterator_;
  };

  

  // CartesianGridLeafIntersectionIteratorTraits
  // ------------------------------------

  template< class Grid >
  struct CartesianGridLeafIntersectionIteratorTraits
  {
    typedef typename remove_const< Grid >::type::Traits GridTraits;

    typedef typename GridTraits::LeafIntersection Intersection;
    typedef CartesianGridLeafIntersection< const Grid > IntersectionImpl;

    typedef typename GridTraits::HostGrid::Traits::LeafIntersectionIterator
      HostIntersectionIterator;
  };
  

  
  // LeafIntersectionIterator
  // ------------------------

  template< class Grid >
  class CartesianGridLeafIntersectionIterator
  : public CartesianGridIntersectionIterator< CartesianGridLeafIntersectionIteratorTraits< Grid > >
  {
    typedef CartesianGridLeafIntersectionIteratorTraits< Grid > Traits;
    typedef CartesianGridIntersectionIterator< Traits > Base;

    typedef typename Traits::HostIntersectionIterator HostIntersectionIterator;
    typedef typename Traits::GridTraits::ExtraData  ExtraDataType ;
    
  public:
    typedef typename Traits::Intersection Intersection;

  public:
    CartesianGridLeafIntersectionIterator ( ExtraDataType data,
                                           const HostIntersectionIterator &hostIterator )
    : Base( data, hostIterator )
    {}
  };



  // CartesianGridLevelIntersectionIteratorTraits
  // -------------------------------------

  template< class Grid >
  struct CartesianGridLevelIntersectionIteratorTraits
  {
    typedef typename remove_const< Grid >::type::Traits GridTraits;

    typedef typename GridTraits::LevelIntersection Intersection;
    typedef CartesianGridLevelIntersection< const Grid > IntersectionImpl;

    typedef typename GridTraits::HostGrid::Traits::LevelIntersectionIterator
      HostIntersectionIterator;
  };
  

  
  // CartesianGridLevelIntersectionIterator
  // -------------------------------

  template< class Grid >
  class CartesianGridLevelIntersectionIterator
  : public CartesianGridIntersectionIterator< CartesianGridLevelIntersectionIteratorTraits< Grid > >
  {
    typedef CartesianGridLevelIntersectionIteratorTraits< Grid > Traits;
    typedef CartesianGridIntersectionIterator< Traits > Base;

    typedef typename Traits::HostIntersectionIterator HostIntersectionIterator;
    typedef typename Traits::GridTraits::ExtraData  ExtraDataType ;
    
  public:
    typedef typename Traits::Intersection Intersection;

  public:
    CartesianGridLevelIntersectionIterator ( ExtraDataType data,
                                            const HostIntersectionIterator &hostIterator )
    : Base( data, hostIterator )
    {}
  };

}

#endif // #ifndef DUNE_CARTESIANGRID_INTERSECTIONITERATOR_HH
