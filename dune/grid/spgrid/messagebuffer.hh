#ifndef DUNE_GRID_SPGRID_MESSAGEBUFFER_HH
#define DUNE_GRID_SPGRID_MESSAGEBUFFER_HH

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>


#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpicommunication.hh>

namespace Dune
{

  // SPBasicPackedMessageWriteBuffer
  // -------------------------------

  class SPBasicPackedMessageWriteBuffer
  {
    typedef SPBasicPackedMessageWriteBuffer This;

  public:
    SPBasicPackedMessageWriteBuffer () { initialize(); }

    SPBasicPackedMessageWriteBuffer ( const This & ) = delete;

    SPBasicPackedMessageWriteBuffer ( This &&other )
      : buffer_( other.buffer_ ),
        position_( other.position_ ), capacity_( other.capacity_ )
    {
      other.initialize();
    }

    ~SPBasicPackedMessageWriteBuffer () { std::free( buffer_ ); }

    This &operator= ( const This & ) = delete;

    This &operator= ( This &&other )
    {
      buffer_ = other.buffer_;
      position_ = other.position_;
      capacity_ = other.capacity_;
      other.initialize();
      return *this;
    }

    template< class T >
    void write ( const T &value )
    {
      reserve( position_ + sizeof( T ) );
      std::memcpy( static_cast< char * >( buffer_ ) + position_, &value, sizeof( T ) );
      position_ += sizeof( T );
    }

    std::size_t position () const { return position_; }

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

    void *buffer_;
    std::size_t position_, capacity_;
  };



  // SPPackedMessageWriteBuffer
  // --------------------------

  template< class Communication >
  class SPPackedMessageWriteBuffer;

  template< class C >
  class SPPackedMessageWriteBuffer< Communication< C > >
    : public SPBasicPackedMessageWriteBuffer
  {
    typedef SPPackedMessageWriteBuffer< Communication< C > > This;
    typedef SPBasicPackedMessageWriteBuffer Base;

  public:
    explicit SPPackedMessageWriteBuffer ( const Communication< C > &comm ) {}

    void send ( int rank, int tag ) {}
    void wait () {}
  };

#if HAVE_MPI
  template<>
  class SPPackedMessageWriteBuffer< Communication< MPI_Comm > >
    : public SPBasicPackedMessageWriteBuffer
  {
    typedef SPPackedMessageWriteBuffer< Communication< MPI_Comm > > This;
    typedef SPBasicPackedMessageWriteBuffer Base;

  public:
    explicit SPPackedMessageWriteBuffer ( const Communication< MPI_Comm > &comm ) : comm_( comm ) {}

    void send ( int rank, int tag )
    {
      MPI_Isend( buffer_, position_, MPI_PACKED, rank, tag, comm_, &request_ );
    }

    void wait () { MPI_Wait( &request_, MPI_STATUS_IGNORE ); }

  protected:
    MPI_Comm comm_;
    MPI_Request request_;
  };
#endif // #if HAVE_MPI



  // SPBasicPackedMessageReadBuffer
  // ------------------------------

  class SPBasicPackedMessageReadBuffer
  {
    typedef SPBasicPackedMessageReadBuffer This;

  public:
    SPBasicPackedMessageReadBuffer () { initialize(); }

    SPBasicPackedMessageReadBuffer ( const This & ) = delete;

    SPBasicPackedMessageReadBuffer ( This &&other )
      : buffer_( other.buffer_ ),
        position_( other.position_ ), size_( other.size_ )
    {
      other.initialize();
    }

    ~SPBasicPackedMessageReadBuffer () { std::free( buffer_ ); }

    This &operator= ( const This & ) = delete;

    This &operator= ( This &&other )
    {
      buffer_ = other.buffer_;
      position_ = other.position_;
      size_ = other.size_;
      other.initialize();
      return *this;
    }

    template< class T >
    void read ( T &value )
    {
      if( position_ + sizeof( T ) <= size_ )
      {
        std::memcpy( static_cast< void * >( &value ), static_cast< char * >( buffer_ ) + position_, sizeof( T ) );
        position_ += sizeof( T );
      }
      else
        DUNE_THROW( IOError, "Cannot read beyond the buffer's end." );
    }

    std::size_t position () const { return position_; }

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

    void *buffer_;
    std::size_t position_, size_;
  };



  // SPPackedMessageReadBuffer
  // -------------------------

  template< class Communication >
  class SPPackedMessageReadBuffer;

  template< class C >
  class SPPackedMessageReadBuffer< Communication< C > >
    : public SPBasicPackedMessageReadBuffer
  {
    typedef SPPackedMessageReadBuffer< Communication< C > > This;
    typedef SPBasicPackedMessageReadBuffer Base;

  public:
    explicit SPPackedMessageReadBuffer ( const Communication< C > &comm ) {}

    void receive ( int rank, int rag, std::size_t size )
    {
      DUNE_THROW( IOError, "Nothing to receive in a serial communication." );
    }

    void receive ( int rank, int tag ) { receive( rank, tag, 0 ); }
    void receive ( int tag ) { receive( 0, tag, 0 ); }

    int rank () const { return 0 ; }

    void wait () {}

    friend inline typename std::vector< This >::iterator waitAny ( std::vector< This > &readBuffers )
    {
      return readBuffers.end();
    }
  };

#if HAVE_MPI
  template<>
  class SPPackedMessageReadBuffer< Communication< MPI_Comm > >
    : public SPBasicPackedMessageReadBuffer
  {
    typedef SPPackedMessageReadBuffer< Communication< MPI_Comm > > This;
    typedef SPBasicPackedMessageReadBuffer Base;

  public:
    SPPackedMessageReadBuffer ( const Communication< MPI_Comm > &comm ) : comm_( comm ) {}

    void receive ( int rank, int tag, std::size_t size )
    {
      rank_ = rank;
      reset( size );
      MPI_Irecv( buffer_, size_, MPI_BYTE, rank, tag, comm_, &request_ );
    }

    void receive ( int rank, int tag )
    {
      MPI_Status status;
      MPI_Probe( rank, tag, comm_, &status );
      int count;
      MPI_Get_count( &status, MPI_BYTE, &count );
      receive( status.MPI_SOURCE, tag, count );
    }

    void receive ( int tag ) { receive( MPI_ANY_SOURCE, tag ); }

    int rank () const { return rank_; }

    void wait () { MPI_Wait( &request_, MPI_STATUS_IGNORE ); }

    friend inline typename std::vector< This >::iterator waitAny ( std::vector< This > &readBuffers )
    {
      const std::size_t numBuffers = readBuffers.size();
      std::vector< MPI_Request > requests( numBuffers );
      for( std::size_t i = 0; i < numBuffers; ++i )
        requests[ i ] = readBuffers[ i ].request_;

      int index = MPI_UNDEFINED;
      MPI_Waitany( numBuffers, requests.data(), &index, MPI_STATUS_IGNORE );
      if( index == MPI_UNDEFINED )
        return readBuffers.end();

      readBuffers[ index ].request_ = requests[ index ];
      return readBuffers.begin() + index;
    }

  protected:
    int rank_;
    MPI_Comm comm_;
    MPI_Request request_;
  };
#endif // #if HAVE_MPI

} // namespace Dune

#endif // #ifndef DUNE_GRID_SPGRID_MESSAGEBUFFER_HH
