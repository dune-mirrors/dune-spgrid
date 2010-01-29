#ifndef DUNE_SPGRID_PARTITION_HH
#define DUNE_SPGRID_PARTITION_HH

#include <dune/common/iostream.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/mesh.hh>

namespace Dune
{

  // SPPartition
  // -----------

  template< int dim >
  class SPPartition
  {
    typedef SPPartition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;
    typedef SPMesh< dimension > Mesh;

    explicit SPPartition ( const Mesh &mesh, const unsigned int number )
    : begin_( 2*mesh.begin() ),
      end_( 2*mesh.end() ),
      number_( number )
    {}

    SPPartition ( const MultiIndex &begin, const MultiIndex &end, const unsigned int number )
    : begin_( begin ),
      end_( end ),
      number_( number )
    {}

    const MultiIndex &begin () const;
    const MultiIndex &end () const;

    unsigned int number () const;

    bool contains ( const MultiIndex &id ) const;

    int volume () const;
    MultiIndex width () const;

    template< class char_type, class traits >
    void print ( std::basic_ostream< char_type, traits > &out ) const;

    template< class char_type, class traits >
    void print ( std::basic_ostream< char_type, traits > &out, const int i ) const;

  private:
    MultiIndex begin_, end_;
    unsigned int number_;
  };



  // Implementation of SPPartition
  // -----------------------------

  template< int dim >
  inline const typename SPPartition< dim >::MultiIndex &
  SPPartition< dim >::begin () const
  {
    return begin_;
  }


  template< int dim >
  inline const typename SPPartition< dim >::MultiIndex &
  SPPartition< dim >::end () const
  {
    return end_;
  }


  template< int dim >
  inline unsigned int SPPartition< dim >::number () const
  {
    return number_;
  }


  template< int dim >
  inline bool SPPartition< dim >::contains ( const MultiIndex &id ) const
  {
    bool contains = true;
    for( int i = 0; i < dimension; ++i )
      contains &= (id[ i ] >= begin()[ i ]) && (id[ i ] <= end()[ i ]);
    return contains;
  }


  template< int dim >
  inline int SPPartition< dim >::volume () const
  {
    const MultiIndex &w = width();
    int volume = 1;
    for( int i = 0; i < dimension; ++i )
      volume *= w[ i ];
    return volume;
  }


  template< int dim >
  inline typename SPPartition< dim >::MultiIndex SPPartition< dim >::width () const
  {
    MultiIndex width;
    for( int i = 0; i < dimension; ++i )
      width[ i ] = std::max( (end()[ i ]+1)/2 - begin()[ i ]/2, 0 );
    return width;
  }


  template< int dim >
  template< class char_type, class traits >
  inline void SPPartition< dim >
    ::print ( std::basic_ostream< char_type, traits > &out ) const
  {
    print( out, 0 );
    for( int i = 1; i < dimension; ++i )
    {
      out << " x ";
      print( out, i );
    }
  }


  template< int dim >
  template< class char_type, class traits >
  inline void SPPartition< dim >
    ::print ( std::basic_ostream< char_type, traits > &out, const int i ) const
  {
    const int b = begin()[ i ];
    const int e= end()[ i ];

    const char left = '[' + (b & 1)*(']'-'[');
    const char right = ']' + (e & 1)*('['-']');

    out << left << ' ' << (b/2) << ", " << ((e+1)/2) << ' ' << right;
  }



  // Auxilliary functions for SPPartition
  // ------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPPartition< dim > &partition )
  {
    partition.print( out );
    return out;
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITION_HH
