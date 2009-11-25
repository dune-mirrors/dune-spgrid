#ifndef DUNE_SPGRID_COMMUNICATE_HH
#define DUNE_SPGRID_COMMUNICATE_HH

#if HAVE_MPI
#include <mpi.h>
#endif // #if HAVE_MPI

#include <dune/common/forloop.hh>

#include <dune/grid/spgrid/iterator.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class Grid, class DataHandle >
  struct SPPartitionSend;

  template< class Grid, class DataHandle >
  struct SPPartitionReceive;



  // SPPartitionSend
  // ---------------

#if HAVE_MPI
  template< class Grid, class DataHandle >
  struct SPPartitionSend
  {
    typedef SPGridLevel< Grid > GridLevel;

    static const unsigned int dimension = GridLevel::dimension;

    void operator() ( const unsigned int rank, const Partition *partition )
    {
      std::vector< DataHandle::DataType > buffer;
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, partition, buffer );
    }

  private:
    template< int codim >
    struct Codim;

    const GridLevel &gridLevel_;
    const DataHandle &dataHandle_;
  };


  template< class Grid, class DataHandle >
  template< int codim >
  struct SPPartitionSend< Grid, DataHandle >::Codim
  {
    template< class MessageBuffer >
    static void
    apply ( const GridLevel &gridLevel, const DataHandle &dataHandle,
            const Partition *&partition, MessageBuffer &buffer )
    {
      typedef SPPartitionIterator< codim, Grid > Iterator;

      if( dataHandle.contains( dimension, codim ) )
      {
        const Iterator end( gridLevel, 0 );
        for( Iterator it( gridLevel, partition ); it != end; ++it )
          dataHandle.gather( buffer, *it ); 
      }
    };
  };
#endif // #if HAVE_MPI



  // SPPartitionReceive
  // ------------------

#if HAVE_MPI
  template< class Grid, class DataHandle >
  struct SPPartitionReceive
  {
    typedef SPGridLevel< Grid > GridLevel;

    static const unsigned int dimension = GridLevel::dimension;

    void operator() ( const unsigned int rank, const Partition *partition )
    {
      std::vector< DataHandle::DataType > buffer;
      ForLoop< Codim, 0, dimension >::apply( gridLevel_, dataHandle_, partition, buffer );
    }

  private:
    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
  };


  template< class Grid, class DataHandle >
  template< int codim >
  struct SPPartitionReceive< Grid, DataHandle >::Codim
  {
    template< class MessageBuffer >
    static void
    apply ( const GridLevel &gridLevel, DataHandle &dataHandle,
            const Partition *&partition, MessageBuffer &buffer )
    {
      typedef SPPartitionIterator< codim, Grid > Iterator;

      if( dataHandle.contains( dimension, codim ) )
      {
        const Iterator end( gridLevel, 0 );
        for( Iterator it( gridLevel, partition ); it != end; ++it )
          dataHandle.scatter( buffer, *it, dataHandle.size( *it ) ); 
      }
    };
  };
#endif // #if HAVE_MPI

}

#endif // #ifndef DUNE_SPGRID_COMMUNICATE_HH
