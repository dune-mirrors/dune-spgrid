#ifndef DUNE_GRID_SPGRID_TREE_HH
#define DUNE_GRID_SPGRID_TREE_HH

#include <type_traits>

#include <dune/grid/common/entity.hh>
#include <dune/grid/common/intersection.hh>

#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/entityinfo.hh>
#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  namespace __SPGrid
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class, class >
    class TreeIterator;

    template< class, class >
    struct Tree;



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

        bool operator () ( const EntityInfo &entityInfo ) const
        {
          return (&entityInfo.gridLevel() == rootLevel_);
        }

      private:
        const GridLevel *rootLevel_;
      };

    public:
      using Base::entityInfo;
      using Base::gridLevel;

      TreeIterator () = default;

      explicit TreeIterator ( const IsLeaf &isLeaf ) : isLeaf_( isLeaf ) {}

      explicit TreeIterator ( const EntityInfo &entityInfo, const IsLeaf &isLeaf )
        : Base( entityInfo ),
          isDone_( entityInfo.gridLevel() ),
          isLeaf_( isLeaf )
      {}

      TreeIterator ( const This & ) = default;
      TreeIterator ( This && ) = default;

      This &operator= ( const This & ) = default;
      This &operator= ( This && ) = default;

      void increment ()
      {
        if( isLeaf( entityInfo() ) || (&gridLevel() == &gridLevel().grid().leafLevel()) )
        {
          while( !isDone( entityInfo() ) )
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
        if( isLeaf( entityInfo ) || (&entityInfo.gridLevel() == &entityInfo.gridLevel().grid().leafLevel()) )
        {
          while( !isDone( entityInfo ) )
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



    // Tree for Entity
    // ---------------

    template< int codim, int dim, class Grid, template< int, int, class > class EntityImpl, class IsLeaf >
    struct Tree< Dune::Entity< codim, dim, Grid, EntityImpl >, IsLeaf >
    {
      typedef Dune::Entity< codim, dim, Grid, EntityImpl > Entity;
      typedef Dune::EntityIterator< codim, Grid, TreeIterator< Entity, IsLeaf > > Iterator;

      Tree ( const Entity &entity, const IsLeaf &isLeaf )
        : entity_( entity ), isLeaf_( isLeaf )
      {}

      Iterator begin () const { return TreeIterator< Entity, IsLeaf >( entity_, isLeaf_ ); }
      Iterator end () const { return TreeIterator< Entity, IsLeaf >( isLeaf_ ); }

      bool empty () const { return false; }

    private:
      Entity entity_;
      IsLeaf isLeaf_;
    };



    // Tree for Intersection
    // ---------------------

    template< class Grid, class IntersectionImpl, class IsLeaf >
    struct Tree< Dune::Intersection< Grid, IntersectionImpl >, IsLeaf >
    {
      typedef Dune::Intersection< Grid, IntersectionImpl > Intersection;
      typedef Dune::IntersectionIterator< Grid, TreeIterator< Intersection, IsLeaf >, IntersectionImpl > Iterator;

      Tree ( const Intersection &intersection, const IsLeaf &isLeaf )
        : intersection_( intersection ), isLeaf_( isLeaf )
      {}

      Iterator begin () const { return TreeIterator< Intersection, IsLeaf >( intersection_, isLeaf_ ); }
      Iterator end () const { return TreeIterator< Intersection, IsLeaf >( intersection_.indexInInside(), isLeaf_ ); }

      bool empty () const { return false; }

    private:
      Intersection intersection_;
      IsLeaf isLeaf_;
    };

  } // namespace __SPGrid



  // tree
  // ----

  template< class ct, int dim, template< int > class Ref, class Comm, class Entity, class IsLeaf >
  inline __SPGrid::Tree< Entity, IsLeaf >
  tree ( const SPGrid< ct, dim, Ref, Comm > &grid, const Entity &entity, const IsLeaf &isLeaf )
  {
    __SPGrid::Tree< Entity, IsLeaf >( entity, isLeaf );
  }

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_TREE_HH
