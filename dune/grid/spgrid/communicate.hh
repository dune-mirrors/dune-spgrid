#ifndef DUNE_SPGRID_COMMUNICATE_HH
#define DUNE_SPGRID_COMMUNICATE_HH

#include <dune/common/forloop.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpicollectivecommunication.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/grid/spgrid/iterator.hh>

namespace Dune
{

  // SPCommunication
  // ---------------

  template< class Grid, class DataHandle >
  class SPCommunication
  {
    template< class T >
    class WriteBuffer;

    template< class T >
    class ReadBuffer;

    template< int codim >
    struct Codim;

  public:
    typedef SPGridLevel< Grid > GridLevel;
    typedef typename GridLevel::PartitionList PartitionList;

    static const unsigned int dimension = GridLevel::dimension;

    typedef typename DataHandle::DataType DataType;

  private:
    typedef std::pair< WriteBuffer< int >, WriteBuffer< DataType > > Packet;

  public:
    SPCommunication ( const GridLevel &gridLevel, DataHandle &dataHandle )
    : gridLevel_( gridLevel ),
      dataHandle_( dataHandle )
    {}

    ~SPCommunication ()
    {
      typedef typename std::vector< Packet * >::iterator Iterator;
      const Iterator end = packets_.end();
      for( Iterator it = packets_.begin(); it != end; ++it )
      {
        it->first.wait( gridLevel_.grid().comm() );
        it->second.wait( gridLevel_.grid().comm() );
        delete *it;
      }
    }

    void gather ( const unsigned int rank, const PartitionList &partitionList )
    {
      Packet *packet = new Packet;
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, partitionList, packet->first, packet->second );
      packet->first.send( rank, 1, gridLevel_.grid().comm() );
      packet->second.send( rank, 2, gridLevel_.grid().comm() );
      packets_.push_back( packet );
    }

    void scatter ( const unsigned int rank, const PartitionList &partitionList )
    {
      ReadBuffer sizes( rank, 1, gridLevel_.grid().comm() );
      ReadBuffer buffer( rank, 2, gridLevel_.grid().comm() );
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, partitionList, sizes, buffer );
    }

  private:
    SPCommunication ( const SPCommunication &other );

    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
    std::vector< Packet * > packets_;
  };


  template< class Grid, class DataHandle >
  template< class T >
  struct SPCommunication< Grid, DataHandle >::WriteBuffer
  {
    void write ( const T &value )
    {
      buffer_.push_back( value );
    }

    template< class C >
    int send ( int rank, int tag, const CollectiveCommunication< C > &comm )
    {}

#if HAVE_MPI
    int send ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm )
    {
      MPI_Datatype mpiDataType = Generic_MPI_Datatype< T >::get();
      return MPI_Isend( &(buffer_[ 0 ]), buffer_.size(), mpiDataType, rank, tag, comm, &request_ );
    }
#endif // #if HAVE_MPI

    template< class C >
    int wait ( const CollectiveCommunication &comm )
    {}

#if HAVE_MPI
    int wait ( const CollectiveCommunication< MPI_Comm > &comm )
    {
      MPI_Status status;
      return MPI_Wait( &request_, &status );
    }
#endif // #if HAVE_MPI

  private:
    std::vector< T > buffer_;
#if HAVE_MPI
    MPI_Request request_;
#endif // #if HAVE_MPI
  };


  template< class Grid, class DataHandle >
  template< class T >
  struct SPCommunication< Grid, DataHandle >::ReadBuffer
  {
    template< class C >
    ReadBuffer ( int rank, int tag, const CollectiveCommunication< C > &comm )
    : read_( buffer.begin() )
    {}

#if HAVE_MPI
    ReadBuffer ( int rank, int tag, const CollectiveCommunication< MPI_Comm > &comm )
    : read_( buffer_.begin() )
    {
      MPI_Status status;
      MPI_Probe( rank, tag, comm, &status );

      MPI_Datatype mpitype = Generic_MPI_Datatype< T >::get();

      int count;
      MPI_Get_count( &status, mpitype, &count );
      buffer_.resize( count );

      MPI_Recv( &(buffer_[ 0 ]), buffer_.size(), mpitype, rank, tag, comm, &status );

      read_ = buffer_.begin();
    }
#endif // #if HAVE_MPI

    void read ( T &value )
    {
      if( read_ != buffer_.end() )
      {
        value = *read_;
        ++read_;
      }
      else
        DUNE_THROW( IOError, "Cannot read beyond the buffer's end." );
    }

  private:
    std::vector< T > buffer_;
    typename std::vector< T >::const_iterator read_;
  };


  template< class Grid, class DataHandle >
  template< int codim >
  struct SPCommunication< Grid, DataHandle >::Codim
  {
    typedef SPPartitionIterator< codim, Grid > Iterator;

    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const PartitionList &partitionList,
            WriteBuffer< int > &sizes, WriteBuffer< DataType > &buffer )
    {
      if( dataHandle.contains( dimension, codim ) )
      {
        const Iterator end( gridLevel, partitionList, typename Iterator::End() );
        for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
        {
          sizes.write( dataHandle.size( *it ) );
          dataHandle.gather( buffer, *it ); 
        }
      }
    }

    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const PartitionList &partitionList,
            ReadBuffer< int > &sizes, ReadBuffer< DataType > &data )
    {
      if( dataHandle.contains( dimension, codim ) )
      {
        const Iterator end( gridLevel, partitionList, typename Iterator::End() );
        for( Iterator it( gridLevel, partitionList, typename Iterator::Begin() ); it != end; ++it )
        {
          int n;
          sizes.read( n );
          dataHandle.scatter( data, *it, n ); 
        }
      }
    }
  };

}

#endif // #ifndef DUNE_SPGRID_COMMUNICATE_HH
