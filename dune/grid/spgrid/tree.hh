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
#include <dune/grid/spgrid/entitypointer.hh>
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

    template< int codim, int dim, class Grid, template< int, int, class > class EntityImpl, class IsLeaf >
    class TreeIterator< Dune::Entity< codim, dim, Grid, EntityImpl >, IsLeaf >
      : public SPEntityPointer< codim, Grid >
    {
      typedef TreeIterator< Dune::Entity< codim, dim, Grid, EntityImpl >, IsLeaf > This;
      typedef SPEntityPointer< codim, Grid > Base;

    public:
      typedef typename Base::Entity Entity;
      typedef typename Base::EntityInfo EntityInfo;
      typedef typename Base::GridLevel GridLevel;

    private:
      struct IsDone
      {
        IsDone () : rootLevel_( nullptr ) {}
        IsDone ( const GridLevel &rootLevel ) : rootLevel_( &rootLevel ) {}

        bool operator () ( const Entity &entity ) const
        {
          return (&Grid::getRealImplementation( entity ).gridLevel() == rootLevel_);
        }

      private:
        const GridLevel *rootLevel_;
      };

    public:
      using Base::dereference;
      using Base::entityInfo;
      using Base::gridLevel;

      TreeIterator () = default;

      explicit TreeIterator ( const IsLeaf &isLeaf ) : isLeaf_( isLeaf ) {}

      explicit TreeIterator ( const Entity &entity, const IsLeaf &isLeaf )
        : Base( Grid::getRealImplementation( entity ) ),
          isDone_( gridLevel() ),
          isLeaf_( isLeaf )
      {}

      TreeIterator ( const This & ) = default;
      TreeIterator ( This && ) = default;

      This &operator= ( const This & ) = default;
      This &operator= ( This && ) = default;

      void increment ()
      {
        if( isLeaf_( dereference() ) || (&gridLevel() == &leafLevel()) )
        {
          while( !isDone_( dereference() ) )
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

    private:
      const GridLevel &leafLevel () const { return gridLevel().grid().leafLevel(); }

      IsDone isDone_;
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

    private:
      struct IsDone
      {
        IsDone () : rootLevel_( nullptr ) {}
        IsDone ( const GridLevel &rootLevel ) : rootLevel_( &rootLevel ) {}

        bool operator () ( const EntityInfo &entityInfo ) const
        {
          return (&entityInfo.gridLevel() == rootLevel_);
        }

      private:
        const GridLevel *rootLevel_;
      };

    public:
      TreeIterator () = default;

      explicit TreeIterator ( int face, const IsLeaf &isLeaf )
        : intersection_( IntersectionImpl( ElementInfo(), face ) )
      {}

      explicit TreeIterator ( const Intersection &intersection, const IsLeaf &isLeaf )
        : intersection_( intersection ),
          isDone_( Grid::getRealImplementation( intersection ).gridLevel() ),
          isLeaf_( isLeaf )
      {}

      TreeIterator ( const This & ) = default;
      TreeIterator ( This && ) = default;

      This &operator= ( const This & ) = default;
      This &operator= ( This && ) = default;

      const Intersection &dereference () const { return intersection_; }

      bool equals ( const This &other ) const
      {
        return Grid::getRealImplementation( intersection_ ).equals( Grid::getRealImplementation( other.intersection_ ) );
      }

      void increment ()
      {
        EntityInfo entityInfo = Grid::getRealImplementation( intersection_ ).entityInfo();
        if( isLeaf_( entityInfo ) || (&entityInfo.gridLevel() == &entityInfo.gridLevel().grid().leafLevel()) )
        {
          while( !isDone_( entityInfo ) )
          {
            if( entityInfo.nextChild() )
            {
              Grid::getRealImplementation( intersection_ ).setEntityInfo( entityInfo );
              return;
            }
            entityInfo.up();
          }
          Grid::getRealImplementation( intersection_ ).setInside( ElementInfo() );
        }
        else
        {
          entityInfo.down();
          Grid::getRealImplementation( intersection_ ).setEntityInfo( entityInfo );
        }
      }

    private:
      Intersection intersection_;
      IsDone isDone_;
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
