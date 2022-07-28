#ifndef DUNE_GRID_SPGRID_TREE_HH
#define DUNE_GRID_SPGRID_TREE_HH

#include <type_traits>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>

#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/entityinfo.hh>
#include <dune/grid/spgrid/entity.hh>
#include <dune/grid/spgrid/intersection.hh>

namespace Dune
{

  // Forward Declarations
  // --------------------

  template< int codim, class Grid, class IsLeaf >
  class EntityTree;

  template< class Grid, class IsLeaf >
  class IntersectionTree;



  namespace __SPGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class, class >
    class TreeIterator;



    // TreeIterator for Entity
    // -----------------------

    template< int codim, int dim, class Grid, template< int, int, class > class EImpl, class IsLeaf >
    class TreeIterator< Dune::Entity< codim, dim, Grid, EImpl >, IsLeaf >
    {
      typedef TreeIterator< Dune::Entity< codim, dim, Grid, EImpl >, IsLeaf > This;

    public:
      typedef typename std::remove_const< Grid >::type::Traits Traits;

      static const int dimension = Traits::ReferenceCube::dimension;
      static const int codimension = codim;
      static const int mydimension = dimension - codimension;

      typedef typename Traits::template Codim< codimension >::Entity Entity;

    private:
      typedef EImpl< codimension, dimension, Grid > EntityImpl;

    public:
      typedef typename EntityImpl::EntityInfo EntityInfo;
      typedef typename EntityImpl::GridLevel GridLevel;

      TreeIterator () = default;

      explicit TreeIterator ( const IsLeaf &isLeaf ) : isLeaf_( isLeaf ) {}

      explicit TreeIterator ( const Entity &entity, const IsLeaf &isLeaf )
        : entityInfo_( entity.impl().entityInfo() ),
          rootLevel_( &gridLevel() ),
          isLeaf_( isLeaf )
      {}

      Entity dereference () const { return EntityImpl( entityInfo() ); }

      bool equals ( const This &other ) const { return entityInfo().equals( other.entityInfo() ); }

      void increment ()
      {
        if( isLeaf_( dereference() ) || (&gridLevel() == &leafLevel()) )
        {
          while( !isDone() )
          {
            if( entityInfo().nextChild() )
              return;
            entityInfo().up();
          }
          entityInfo() = EntityInfo();
        }
        else
          entityInfo().down();
      }

      const EntityInfo &entityInfo () const { return entityInfo_; }
      EntityInfo &entityInfo () { return entityInfo_; }

      const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

    private:
      const GridLevel &leafLevel () const { return gridLevel().grid().leafLevel(); }

      bool isDone () const { return (&gridLevel() == rootLevel_); }

      EntityInfo entityInfo_;
      const GridLevel *rootLevel_ = nullptr;
      IsLeaf isLeaf_;
    };



    // TreeIterator for Intersection
    // -----------------------------

    template< class Grid, class IntersectionImpl, class IsLeaf >
    class TreeIterator< Dune::Intersection< Grid, IntersectionImpl >, IsLeaf >
    {
      typedef TreeIterator< Dune::Intersection< Grid, IntersectionImpl >, IsLeaf > This;

    public:
      typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;

      static const int dimension = Intersection::dimension;
      static const int codimension = Intersection::codimension;
      static const int mydimension = Intersection::mydimension;

      typedef typename IntersectionImpl::EntityInfo EntityInfo;
      typedef typename IntersectionImpl::ElementInfo ElementInfo;
      typedef typename IntersectionImpl::GridLevel GridLevel;

      TreeIterator () = default;

      explicit TreeIterator ( int face, const IsLeaf &isLeaf )
        : intersection_( IntersectionImpl( ElementInfo(), face ) )
      {}

      explicit TreeIterator ( const Intersection &intersection, const IsLeaf &isLeaf )
        : intersection_( intersection ),
          rootLevel_( &intersection.impl().gridLevel() ),
          isLeaf_( isLeaf )
      {}

      TreeIterator ( const This & ) = default;
      TreeIterator ( This && ) = default;

      This &operator= ( const This & ) = default;
      This &operator= ( This && ) = default;

      const Intersection &dereference () const { return intersection_; }

      bool equals ( const This &other ) const
      {
        return intersection_.impl().equals( other.intersection_.impl() );
      }

      void increment ()
      {
        if( isLeaf_( intersection_ ) || (&gridLevel() == &leafLevel()) )
        {
          EntityInfo info = entityInfo();
          while( !isDone( info ) )
          {
            if( info.nextChild() )
            {
              intersection_.impl().setEntityInfo( info );
              return;
            }
            info.up();
          }
          intersection_.impl().setInside( ElementInfo() );
        }
        else
        {
          EntityInfo info = entityInfo();
          info.down();
          intersection_.impl().setEntityInfo( info );
        }
      }

    private:
      EntityInfo entityInfo () const { return intersection_.impl().entityInfo(); }

      const GridLevel &gridLevel () const  { return entityInfo().gridLevel(); }
      const GridLevel &leafLevel () const { return gridLevel().grid().leafLevel(); }

      bool isDone ( const EntityInfo &entityInfo ) const { return (&entityInfo.gridLevel() == rootLevel_); }

      Intersection intersection_;
      const GridLevel *rootLevel_;
      IsLeaf isLeaf_;
    };

  } // namespace __SPGrid



  // EntityTree for SPGrid
  // ---------------------

  template< int codim, class ct, int dim, template< int > class Ref, class Comm, class IsLeaf >
  class EntityTree< codim, SPGrid< ct, dim, Ref, Comm >, IsLeaf >
  {
  public:
    typedef SPGrid< ct, dim, Ref, Comm > Grid;

    typedef Dune::Entity< codim, dim, const Grid, SPEntity > Entity;
    typedef Dune::EntityIterator< codim, Grid, __SPGrid::TreeIterator< Entity, IsLeaf > > Iterator;

    EntityTree ( const Grid &grid, const Entity &entity, const IsLeaf &isLeaf )
      : entity_( entity ), isLeaf_( isLeaf )
    {}

    Iterator begin () const { return __SPGrid::TreeIterator< Entity, IsLeaf >( entity_, isLeaf_ ); }
    Iterator end () const { return __SPGrid::TreeIterator< Entity, IsLeaf >( isLeaf_ ); }

    bool empty () const { return false; }

  private:
    Entity entity_;
    IsLeaf isLeaf_;
  };



  // IntersectionTree for SPGrid
  // ---------------------------

  template< class ct, int dim, template< int > class Ref, class Comm, class IsLeaf >
  class IntersectionTree< SPGrid< ct, dim, Ref, Comm >, IsLeaf >
  {
  public:
    typedef SPGrid< ct, dim, Ref, Comm > Grid;

    typedef Dune::Intersection< const Grid, SPIntersection< const Grid > > Intersection;
    typedef Dune::IntersectionIterator< const Grid, __SPGrid::TreeIterator< Intersection, IsLeaf >, SPIntersection< const Grid > > Iterator;

    IntersectionTree ( const Grid &grid, const Intersection &intersection, const IsLeaf &isLeaf )
      : intersection_( intersection ), isLeaf_( isLeaf )
    {}

    Iterator begin () const { return __SPGrid::TreeIterator< Intersection, IsLeaf >( intersection_, isLeaf_ ); }
    Iterator end () const { return __SPGrid::TreeIterator< Intersection, IsLeaf >( intersection_.indexInInside(), isLeaf_ ); }

    bool empty () const { return false; }

  private:
    Intersection intersection_;
    IsLeaf isLeaf_;
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_TREE_HH
