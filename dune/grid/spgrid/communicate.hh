#ifndef DUNE_SPGRID_COMMUNICATE_HH
#define DUNE_SPGRID_COMMUNICATE_HH

#if HAVE_MPI
#include <mpi.h>
#endif // #if HAVE_MPI

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

    void operator() ( const unsigned int rank, const Partition &partition )
    {
      // Here, we need a partition iterator
    }

  private:
    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
  };
#endif // #if HAVE_MPI



  // SPPartitionReceive
  // ------------------

#if HAVE_MPI
  template< class Grid, class DataHandle >
  struct SPPartitionReceive
  {
    typedef SPGridLevel< Grid > GridLevel;

    void operator() ( const unsigned int rank, const Partition &partition )
    {
      // Here, we need a partition iterator
    }

  private:
    const GridLevel &gridLevel_;
    DataHandle &dataHandle_;
  };
#endif // #if HAVE_MPI

}

#endif // #ifndef DUNE_SPGRID_COMMUNICATE_HH
