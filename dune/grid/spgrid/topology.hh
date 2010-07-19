#ifndef DUNE_SPGRID_TOPOLOGY_HH
#define DUNE_SPGRID_TOPOLOGY_HH

#include <limits>

namespace Dune
{

  // SPTopology
  // --------

  /** \class SPTopology
   *  \brief description of the grid's topology
   *
   *  \tparam  dim  dimension of the grid
   */
  template< int dim >
  class SPTopology
  {
    typedef SPTopology< dim > This;

  public:
    /** \brief dimension of the domain */
    static const int dimension = dim;

    static const int numFaces = 2*dimension;

    /** \brief constructor
     *
     *  \param[in]  periodic  bit field specifying which directions should be
     *                        periodic (defaults to 0)
     */
    SPTopology ( const unsigned int periodic = 0 );

    SPTopology ( const This &other );

    ~SPTopology ();

    const This &operator= ( const This &other );

    unsigned int numNodes () const { return data_[ 0 ]; }

    unsigned int neighbor ( const unsigned int node, const int face ) const;

    /** \brief determine whether a direction is periodic
     *
     *  \param[in]  i  direction (0 <= i < dimension)
     *
     *  \returns true, if direction i is periodic
     */
    bool periodic ( const int i ) const;

    /** \brief obtain the periodicity bit field
     *
     *  \returns the bitfield specifying which directions are periodic
     */
    unsigned int periodic () const;

  private:
    unsigned int &refCount () const { return data_[ 1 ]; }
    unsigned int &nb ( const unsigned int node, const int face );

    unsigned int *data_;
  };



  // Implementation of SPTopology
  // ----------------------------

  template< int dim >
  inline SPTopology< dim >::SPTopology ( const unsigned int periodic )
  {
    const unsigned int numNodes = 1;
    data_ = new unsigned int[ 2 + numFaces*numNodes ];
    data_[ 0 ] = numNodes;
    refCount() = 1;

    const unsigned int max = std::numeric_limits< unsigned int >::max();
    for( int i = 0; i < dimension; ++i )
    {
      const unsigned int b = (periodic >> i) & 1;
      nb( 0, 2*i ) = nb( 0, 2*i+1 ) = (1-b)*max;
    }
  }


  template< int dim >
  inline SPTopology< dim >::SPTopology ( const This &other )
  : data_( other.data_ )
  {
    ++refCount();
  }


  template< int dim >
  inline SPTopology< dim >::~SPTopology ()
  {
    if( --refCount() == 0 )
      delete[] data_;
  }


  template< int dim >
  inline const typename SPTopology< dim >::This &
  SPTopology< dim >::operator= ( const This &other )
  {
    ++other.refCount();
    if( --refCount() == 0 )
      delete[] data_;
    data_ = other.data_;
    return *this;
  }


  template< int dim >
  inline unsigned int SPTopology< dim >::neighbor ( const unsigned int node, const int face ) const
  {
    assert( node < numNodes() );
    assert( (face >= 0) && (face < numFaces) );
    return data_[ node * numFaces + face + 2 ];
  }


  template< int dim >
  inline bool SPTopology< dim >::periodic ( const int i ) const
  {
    assert( numNodes() == 1 );
    assert( (i >= 0) && (i < dimension) );
    return (neighbor( 0, 2*i ) == 0);
  }


  template< int dim >
  inline unsigned int SPTopology< dim >::periodic () const
  {
    unsigned int p = 0;
    for( int i = 0; i < dimension; ++i )
      p |= (1u << i) * (unsigned int)periodic( i );
    return p;
  }


  template< int dim >
  inline unsigned int &
  SPTopology< dim >::nb ( const unsigned int node, const int face )
  {
    assert( node < numNodes() );
    assert( (face >= 0) && (face < numFaces) );
    return data_[ node * numFaces + face + 2 ];
  }

}

#endif // #ifndef DUNE_SPGRID_TOPOLOGY_HH
