#ifndef DUNE_SPGRID_TOPOLOGY_HH
#define DUNE_SPGRID_TOPOLOGY_HH

#include <limits>

namespace Dune
{

  // SPTopology
  // ----------

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
    explicit SPTopology ( const unsigned int periodic = 0 );

    SPTopology ( const This &other );

    ~SPTopology ();

    const This &operator= ( const This &other );

    unsigned int numNodes () const { return data_[ 0 ]; }

    bool hasNeighbor ( const unsigned int node, const int face ) const;
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
  inline bool SPTopology< dim >::hasNeighbor ( const unsigned int node, const int face ) const
  {
    return (neighbor( node, face ) < numNodes());
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



  // Auxilliary Functions for SPTopology
  // -----------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out,
               const SPTopology< dim > &topology )
  {
    typedef SPTopology< dim > Topology;

    const unsigned int numNodes = topology().numNodes();
    out << numNodes;
    for( unsigned int node = 0; node < numNodes; ++node )
    {
      out << " [";
      for( int face = 0; face < Topology::numFaces; ++face )
      {
        if( topology.hasNeighbor( node, face ) )
          out << " " << topology.neighbor( node, face );
        else
          out << " *";
      }
      out << " ]";
    }
    return out;
  }


  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator>> ( std::basic_istream< char_type, traits > &in,
               SPTopology< dim > &topology )
  {
    //unsigned int numNodes = 0;
    //in >> numNodes;
    return in;
  }

}

#endif // #ifndef DUNE_SPGRID_TOPOLOGY_HH
