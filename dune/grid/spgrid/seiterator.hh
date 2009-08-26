#ifndef DUNE_SPGRID_SEITERATOR_HH
#define DUNE_SPGRID_SEITERATOR_HH

#include <dune/grid/common/superelementiterator.hh>

#include <dune/grid/spgrid/entitypointer.hh>

namespace Dune
{

  // SPSuperElementIterator
  // ----------------------

  template< class Grid >
  class SPSuperElementIterator
  : public SPEntityPointer< 0, Grid >
  {
    typedef SPSuperElementIterator< Grid > This;
    typedef SPEntityPointer< 0, Grid > Base;

  public:
    typedef typename Base::Traits Traits;

    typedef typename Base::Entity Entity;
    typedef typename Base::EntityInfo EntityInfo;

    static const int dimension = Base::dimension;

  protected:
    typedef typename EntityInfo::MultiIndex MultiIndex;

    struct Begin {};
    struct End {};

    template< int codim >
    struct Codim
    {
      typedef SPEntity< codim, dimension, Grid > EntityImpl;
    };

  private:
    struct Sequence
    {
      const Sequence *next;
      MultiIndex idAdd;
    };

    struct SequenceProvider;

  protected:
    template< int codim, class BeginEnd >
    SPSuperElementIterator ( const typename Codim< codim >::EntityImpl &entityImpl,
                             const BeginEnd &be )
    : Base( entityImpl.gridLevel() )
    {
      const unsigned int direction = entityImpl.entityInfo().direction();
      sequence_ = SequenceProvider::sequence( direction, be );

      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      entityInfo.id() = entityImpl().entityInfo.id();
      entityInfo.id() += sequence_->idAdd;
      entityInfo.update();
    }

  public:
    void increment ()
    {
      sequence_ = sequence_->next;
      assert( sequence_ != 0 );

      EntityInfo &entityInfo = Grid::getRealImplementation( entity_ ).entityInfo();
      entityInfo.id() += sequence_->idAdd;
      entityInfo.update();
    }

  protected:
    using Base::entity_;

  private:
    const Sequence *sequence_;
  };



  // SPSuperElementIterator::SequenceProvider
  // ----------------------------------------

  template< class Grid >
  struct SPSuperElementIterator< Grid >::SequenceProvider
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
  SPSuperElementIterator< Grid >::SequenceProvider::SequenceProvider ()
  {
    for( unsigned int dir = 0; dir < numDirections; ++dir )
    {
      Sequence *head = new Sequence;

      Sequence *last = head;
      head->idAdd.clear();
      for( unsigned int d = 0; d < numDirections; ++d )
      {
        if( (d & dir) != 0 )
          continue;

        Sequence *next = new Sequence;
        for( int i = 0; i < dimension; ++i )
          next->idAdd[ i ] = (1 - int( (dir >> i) & 1 )) * (2 * int( (d >> i) & 1 ) - 1);
        next->idAdd -= head->idAdd;
        head->idAdd += next->idAdd;
        
        last->next = next;
        last = next;
      }
      begin_[ dir ] = head->next;

      Sequence *end = new Sequence;
      end->next = 0;
      for( int i = 0; i < dimension; ++i )
        end->idAdd[ i ] = 3 * (1 - int( (dir >> i) & 1 ));
      end_[ dir ] = end;

      head->next = 0;
      head->idAdd -= end->idAdd;
      head->idAdd *= -1;
      last->next = head;
    }
  }


  template< class Grid >
  SPSuperElementIterator< Grid >::SequenceProvider::~SequenceProvider ()
  {
    for( unsigned int dir = 0; dir < (1 << dimension); ++dir )
    {
      delete end_[ dir ];
      while( begin_[ dir ] != 0 )
      {
        Sequence *tmp = begin_[ dir ];
        begin_[ dir ] = begin_[ dir ]->next;
        delete tmp;
      }
    }
  }

}

#endif // #ifndef DUNE_SPGRID_HITERATOR_HH
