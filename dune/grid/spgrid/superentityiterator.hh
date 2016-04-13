#ifndef DUNE_SPGRID_SUPERENTITYITERATOR_HH
#define DUNE_SPGRID_SUPERENTITYITERATOR_HH

#include <type_traits>

#include <dune/grid/extensions/superentityiterator.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/entity.hh>
#include <dune/grid/spgrid/referencecube.hh>

namespace Dune
{

  // SPSuperEntityIterator
  // ---------------------

  template< class Grid >
  class SPSuperEntityIterator
  {
    typedef SPSuperEntityIterator< Grid > This;

  public:
    typedef typename std::remove_const< Grid >::type::Traits Traits;

    static const int dimension = Traits::ReferenceCube::dimension;
    static const int codimension = 0;
    static const int mydimension = dimension - codimension;

    typedef typename Traits::template Codim< codimension >::Entity Entity;

  private:
    typedef SPEntity< codimension, dimension, Grid > EntityImpl;

  public:
    typedef typename EntityImpl::EntityInfo EntityInfo;
    typedef typename EntityImpl::GridLevel GridLevel;

    struct Begin {};
    struct End {};

  private:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    struct Sequence
    {
      const Sequence *next;
      MultiIndex idAdd;
      int index;
      unsigned int fBoundary;
    };

    struct SequenceProvider;

  public:
    SPSuperEntityIterator () = default;

    template< class SubInfo, class BeginEnd >
    SPSuperEntityIterator ( const SubInfo &subInfo, const BeginEnd &be )
      : entityInfo_( subInfo.gridLevel() ),
        fBoundary_( 0 )
    {
      sequence_ = SequenceProvider::sequence( subInfo.direction().bits(), be );

      MultiIndex &id = entityInfo().id();

      id = subInfo.id();
      const typename GridLevel::Mesh &globalMesh = entityInfo().gridLevel().globalMesh();
      for( int i = 0; i < dimension; ++i )
      {
        const bool bndLow = (id[ i ] == 2*globalMesh.begin()[ i ]);
        const bool bndHigh = (id[ i ] == 2*globalMesh.end()[ i ]);
        fBoundary_ |= (int( bndLow ) | 2*int( bndHigh )) << (2*i);
      }

      if( next( id ) )
        entityInfo().update( subInfo.partitionNumber() );
    }

    SPSuperEntityIterator ( const This & ) = default;
    SPSuperEntityIterator ( This && ) = default;

    This &operator= ( const This & ) = default;
    This &operator= ( This && ) = default;

    Entity operator* () const { return dereference(); }

    bool operator== ( const This &other ) const { return equals( other ); }
    bool operator!= ( const This &other ) const { return !equals( other ); }

    Entity dereference () const { return EntityImpl( entityInfo() ); }

    bool equals ( const This &other ) const { return entityInfo().equals( other.entityInfo() ); }

    void increment ()
    {
      if( next( entityInfo().id() ) )
        entityInfo().update();
    }

    int index () const { return index_; }

    const EntityInfo &entityInfo () const { return entityInfo_; }
    EntityInfo &entityInfo () { return entityInfo_; }

    const GridLevel &gridLevel () const { return entityInfo().gridLevel(); }

  private:
    bool next ( MultiIndex &id )
    {
      bool skip;
      do {
        assert( sequence_ != 0 );
        id += sequence_->idAdd;
        index_ = sequence_->index;
        skip = ((fBoundary_ & sequence_->fBoundary) != 0);
        sequence_ = sequence_->next;
      } while( skip );
      return (sequence_ != 0);
    }

    EntityInfo entityInfo_;
    const Sequence *sequence_;
    int index_;
    unsigned int fBoundary_;
  };



  // SPSuperEntityIterator::SequenceProvider
  // ---------------------------------------

  template< class Grid >
  struct SPSuperEntityIterator< Grid >::SequenceProvider
  {
    static const unsigned int numDirections = 1 << dimension;

    static const Sequence *
    sequence ( const unsigned int direction, const Begin &b )
    {
      assert( direction < numDirections );
      return instance().begin_[ direction ];
    }

    static const Sequence *
    sequence ( const unsigned int direction, const End &e )
    {
      assert( direction < numDirections );
      return instance().end_[ direction ];
    }

  private:
    SequenceProvider ();
    ~SequenceProvider ();

    static const SequenceProvider &instance ()
    {
      static SequenceProvider instance;
      return instance;
    }

    const Sequence *begin_[ numDirections ];
    const Sequence *end_[ numDirections ];
  };



  template< class Grid >
  SPSuperEntityIterator< Grid >::SequenceProvider::SequenceProvider ()
  {
    SPReferenceCube< typename Grid::ctype, dimension > refCube;

    for( unsigned int dir = 0; dir < numDirections; ++dir )
    {
      const int codim = SPDirection< dimension >( dir ).codimension();

      Sequence *head = new Sequence;

      Sequence *last = head;
      head->idAdd.clear();
      for( unsigned int d = 0; d < numDirections; ++d )
      {
        if( (d & dir) != 0 )
          continue;

        Sequence *next = new Sequence;
        next->fBoundary = 0;
        for( int i = 0; i < dimension; ++i )
        {
          const int dirbit = int( (dir >> i) & 1);
          const int dbit = int( (d >> i) & 1 );
          next->idAdd[ i ] = (1 - dirbit) * (2*dbit - 1);
          next->fBoundary |= (1 - dirbit) * (1 << (2*i + dbit));
        }
        next->index = -1;
        for( int k = 0; k < refCube.count( codim ); ++k )
        {
          const MultiIndex &subId = refCube.subId( codim, k );
          bool found = true;
          for( int i = 0; i < dimension; ++i )
            found &= (next->idAdd[ i ] == -subId[ i ]);
          if( found )
            next->index = k;
        }
        assert( next->index != -1 );
        next->idAdd -= head->idAdd;
        head->idAdd += next->idAdd;
        
        last->next = next;
        last = next;
      }
      begin_[ dir ] = head->next;

      Sequence *end = new Sequence;
      end->next = 0;
      end->index = -1;
      end->fBoundary = 0;
      for( int i = 0; i < dimension; ++i )
        end->idAdd[ i ] = 3 * (1 - int( (dir >> i) & 1 ));
      end_[ dir ] = end;

      head->next = 0;
      head->index = end->index;
      head->fBoundary = end->fBoundary;
      head->idAdd -= end->idAdd;
      head->idAdd *= -1;
      last->next = head;
    }
  }


  template< class Grid >
  SPSuperEntityIterator< Grid >::SequenceProvider::~SequenceProvider ()
  {
    for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
    {
      delete end_[ dir ];
      while( begin_[ dir ] != 0 )
      {
        const Sequence *tmp = begin_[ dir ];
        begin_[ dir ] = begin_[ dir ]->next;
        delete tmp;
      }
    }
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_SUPERENTITYITERATOR_HH
