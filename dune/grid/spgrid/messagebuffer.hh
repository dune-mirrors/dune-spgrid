#ifndef DUNE_GRID_SPGRID_MESSAGEBUFFER_HH
#define DUNE_GRID_SPGRID_MESSAGEBUFFER_HH

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <utility>

#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpitraits.hh>

namespace Dune
{

  // SPSimpleSerialize
  // -----------------

  struct SPSimpleSerialize
  {
    SPSimpleSerialize () {}

    template< class T >
    std::size_t maxSize () const
    {
      return sizeof( T );
    }

    template< class T >
    typename std::enable_if< std::is_trivial< T >::value, std::size_t >::type
    pack ( const T &value, void *buffer, std::size_t position, std::size_t capacity ) const
    {
      assert( position + sizeof( T ) <= capacity );
      std::memcpy( static_cast< char * >( buffer ) + position, &value, sizeof( T ) );
      return position + sizeof( T );
    }

    template< class T >
    typename std::enable_if< std::is_trivial< T >::value, std::size_t >::type
    unpack ( void *buffer, std::size_t position, std::size_t size, T &value ) const
    {
      assert( position + sizeof( T ) <= size );
      std::memcpy( &value, static_cast< char * >( buffer ) + position, sizeof( T ) );
      return position + sizeof( T );
    }
  };



  // SPMPISerialize
  // --------------

#if HAVE_MPI
  struct SPMPISerialize
  {
    SPMPISerialize ( MPI_Comm comm ) : comm_( comm ) {}

    template< class T >
    std::size_t maxSize () const
    {
      MPI_Datatype mpitype = MPITraits< T >::getType();
      int size;
      MPI_Pack_size( 1, mpitype, comm(), &size );
      return size;
    }

    template< class T >
    std::size_t pack ( const T &value, void *buffer, std::size_t position, std::size_t capacity ) const
    {
      assert( position + maxSize< T >() <= capacity );
      MPI_Datatype mpitype = MPITraits< T >::getType();
      int pos = position;
      MPI_Pack( const_cast< T * >( &value ), 1, mpitype, buffer, capacity, &pos, comm() );
      return pos;
    }

    template< class T >
    std::size_t unpack ( void *buffer, std::size_t position, std::size_t size, T &value ) const
    {
      MPI_Datatype mpitype = MPITraits< T >::getType();
      int pos = position;
      MPI_Unpack( buffer, size, &pos, &value, 1, mpitype, comm() );
      return pos;
    }

    MPI_Comm comm () const { return comm_; }

  protected:
    MPI_Comm comm_;
  };
#endif // #if HAVE_MPI



  // SPBasicPackedMessageWriteBuffer
  // -------------------------------

  template< class Serialize >
  class SPBasicPackedMessageWriteBuffer
  {
    typedef SPBasicPackedMessageWriteBuffer< Serialize > This;

  public:
    template< class... Args >
    explicit SPBasicPackedMessageWriteBuffer ( Args &&... args )
      : serialize_( std::forward< Args >( args )... )
    {
      initialize();
    }

    SPBasicPackedMessageWriteBuffer ( const This & ) = delete;

    SPBasicPackedMessageWriteBuffer ( This &&other )
      : serialize_( std::move( other.serialize_ ) ),
        buffer_( other.buffer_ ),
        position_( other.position_ ), capacity_( other.capacity_ )
    {
      other.initialize();
    }

    ~SPBasicPackedMessageWriteBuffer () { std::free( buffer_ ); }

    This &operator= ( const This & ) = delete;

    This &operator= ( This &&other )
    {
      serialize_ = std::move( other.serialize_ );
      buffer_ = other.buffer_;
      position_ = other.position_;
      capacity_ = other.capacity_;
      other.initialize();
    }

    template< class T >
    void write ( const T &value )
    {
      const std::size_t size = serialize().template maxSize< T >();
      reserve( position_ + size );
      position_ = serialize().pack( value, buffer_, position_, capacity_ );
    }

    const Serialize &serialize () const { return serialize_; }

  protected:
    void initialize () { buffer_ = nullptr; position_ = 0; capacity_ = 0; }

    void reserve ( std::size_t size )
    {
      if( size <= capacity_ )
        return;

      std::size_t capacity = std::max( size, 2*capacity_ );
      void *buffer = std::realloc( buffer_, capacity );
      if( !buffer )
      {
        capacity = capacity_ + size;
        buffer = std::realloc( buffer_, capacity );
        if( !buffer )
          DUNE_THROW( OutOfMemoryError, "Cannot allocate sufficiently large buffer." );
      }
      buffer_ = buffer;
      capacity_ = capacity;
    }

    Serialize serialize_;
    void *buffer_;
    std::size_t position_, capacity_;
  };



  // SPPackedMessageWriteBuffer
  // --------------------------

  template< class CollectiveCommunication >
  class SPPackedMessageWriteBuffer;

  template< class C >
  class SPPackedMessageWriteBuffer< CollectiveCommunication< C > >
    : public SPBasicPackedMessageWriteBuffer< SPSimpleSerialize >
  {
    typedef SPPackedMessageWriteBuffer< CollectiveCommunication< C > > This;
    typedef SPBasicPackedMessageWriteBuffer< SPSimpleSerialize > Base;

  public:
    explicit SPPackedMessageWriteBuffer ( const CollectiveCommunication< C > &comm ) {}

    void send ( int rank, int tag ) {}
    void wait () {}
  };

#if HAVE_MPI
  template<>
  class SPPackedMessageWriteBuffer< CollectiveCommunication< MPI_Comm > >
    : public SPBasicPackedMessageWriteBuffer< SPMPISerialize >
  {
    typedef SPPackedMessageWriteBuffer< CollectiveCommunication< MPI_Comm > > This;
    typedef SPBasicPackedMessageWriteBuffer< SPMPISerialize > Base;

  public:
    explicit SPPackedMessageWriteBuffer ( const CollectiveCommunication< MPI_Comm > &comm ) : Base( comm ) {}

    void send ( int rank, int tag )
    {
      MPI_Isend( buffer_, position_, MPI_PACKED, rank, tag, serialize().comm(), &request_ );
    }

    void wait () { MPI_Wait( &request_, MPI_STATUS_IGNORE ); }

  protected:
    MPI_Request request_;
  };
#endif // #if HAVE_MPI



  // SPBasicPackedMessageReadBuffer
  // ------------------------------

  template< class Serialize >
  class SPBasicPackedMessageReadBuffer
  {
    typedef SPBasicPackedMessageReadBuffer< Serialize > This;

  public:
    template< class... Args >
    SPBasicPackedMessageReadBuffer ( Args &&... args )
      : serialize_( std::forward< Args >( args )... )
    {
      initialize();
    }

    SPBasicPackedMessageReadBuffer ( const This & ) = delete;

    SPBasicPackedMessageReadBuffer ( This &&other )
      : serialize_( std::move( other.serialize_ ) ),
        buffer_( other.buffer_ ),
        position_( other.position_ ), size_( other.size_ )
    {
      other.initialize();
    }

    ~SPBasicPackedMessageReadBuffer () { std::free( buffer_ ); }

    This &operator= ( const This & ) = delete;

    This &operator= ( This &&other )
    {
      serialize_ = std::move( other.serialize_ );
      buffer_ = other.buffer_;
      position_ = other.position_;
      size_ = other.size_;
      other.initialize();
    }

    template< class T >
    void read ( T &value )
    {
      if( position_ < size_ )
        position_ = serialize().unpack( buffer_, position_, size_, value );
      else
        DUNE_THROW( IOError, "Cannot read beyond the buffer's end." );
    }

    const Serialize &serialize () const { return serialize_; }

  protected:
    void initialize () { buffer_ = nullptr; position_ = 0; size_ = 0; }

    void reset ( std::size_t size )
    {
      std::free( buffer_ );
      initialize();
      if( size == 0 )
        return;
      buffer_ = std::malloc( size );
      if( !buffer_ )
        DUNE_THROW( OutOfMemoryError, "Cannot allocate sufficiently large buffer." );
      size_ = size;
    }

    Serialize serialize_;
    void *buffer_;
    std::size_t position_, size_;
  };



  // SPPackedMessageReadBuffer
  // -------------------------

  template< class CollectiveCommunication >
  class SPPackedMessageReadBuffer;

  template< class C >
  class SPPackedMessageReadBuffer< CollectiveCommunication< C > >
    : SPBasicPackedMessageReadBuffer< SPSimpleSerialize >
  {
    typedef SPPackedMessageReadBuffer< CollectiveCommunication< C > > This;
    typedef SPBasicPackedMessageReadBuffer< SPSimpleSerialize > Base;

  public:
    explicit SPPackedMessageReadBuffer ( const CollectiveCommunication< C > &comm ) {}

    void receive ( int rank, int rag, std::size_t size )
    {
      DUNE_THROW( IOError, "Nothing to receive in a serial communication." );
    }

    void receive ( int rank, int tag ) { receive( rank, tag, 0 ); }
    void receive ( int tag ) { receive( 0, tag, 0 ); }

    int rank () const { return 0 ; }

    void wait () {}
  };

#if HAVE_MPI
  template<>
  class SPPackedMessageReadBuffer< CollectiveCommunication< MPI_Comm > >
    : public SPBasicPackedMessageReadBuffer< SPMPISerialize >
  {
    typedef SPPackedMessageReadBuffer< CollectiveCommunication< MPI_Comm > > This;
    typedef SPBasicPackedMessageReadBuffer< SPMPISerialize > Base;

  public:
    SPPackedMessageReadBuffer ( const CollectiveCommunication< MPI_Comm > &comm ) : Base( comm ) {}

    void receive ( int rank, int tag, std::size_t size )
    {
      rank_ = rank;
      reset( size );
      MPI_Irecv( buffer_, size_, mpitype(), rank, tag, comm(), &request_ );
    }

    void receive ( int rank, int tag )
    {
      MPI_Status status;
      MPI_Probe( rank, tag, comm(), &status );
      int count;
      MPI_Get_count( &status, mpitype(), &count );
      receive( status.MPI_SOURCE, tag, count );
    }

    void receive ( int tag ) { receive( MPI_ANY_SOURCE, tag ); }

    int rank () const { return rank_; }

    void wait () { MPI_Wait( &request_, MPI_STATUS_IGNORE ); }

  protected:
    static constexpr MPI_Datatype mpitype () { return MPI_PACKED; }

    MPI_Comm comm () const { return serialize().comm(); }

    int rank_;
    MPI_Request request_;
  };
#endif // #if HAVE_MPI

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_MESSAGEBUFFER_HH
