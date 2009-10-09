#ifndef DUNE_SPGRID_REFINEMENT_HH
#define DUNE_SPGRID_REFINEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< class ct, int dim >
  struct SPRefinement;

  template< class ct, int dim >
  std::ostream &operator<< ( std::ostream &, const SPRefinement< ct, dim > & );
  template< class ct, int dim >
  std::istream &operator>> ( std::istream &, SPRefinement< ct, dim > & );



  // SPRefinement
  // ------------

  template< class ct, int dim >
  class SPRefinement
  {
    friend std::ostream &operator<<<> ( std::ostream &, const SPRefinement< ct, dim > & );
    friend std::istream &operator>><> ( std::istream &, SPRefinement< ct, dim > & );

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    SPRefinement ()
    : refDir_( 0 )
    {}

    explicit SPRefinement ( const unsigned int refDir )
    : refDir_( refDir )
    {}

    unsigned int factor ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return ((refDir_ >> i) & 1)+1;
    }

    unsigned int numChildren () const
    {
      return (1 << bitCount( refDir_ ));
    }

    void child ( MultiIndex &id, unsigned int index ) const
    {
      assert( index < numChildren() );
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int b = (refDir_ >> i) & 1;
        id[ i ] = ((id[ i ] + (index & b)) << b) - b;
        index = index >> b;
      }
    }

    unsigned int childIndex ( const MultiIndex &id ) const
    {
      unsigned int index = 0;
      for( int i = dimension-1; i >= 0; --i )
      {
        const unsigned int b = (refDir_ >> i) & 1;
        index = (index << b) | ((id[ i ] >> 1) & b);
      }
      assert( index < numChildren() );
      return index;
    }

    void firstChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = (((refDir_ >> i) & 1) != 0 ? 2*id[ i ] - 1 : id[ i ]);
    }

    bool nextChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        if( ((refDir_ >> i) & 1) == 0 )
          continue;
        id[ i ] ^= 2;
        if( (id[ i ] & 2) != 0 )
          return true;
      }
      return false;
    }

    void father ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = (((refDir_ >> i) & 1) != 0 ? (id[ i ] >> 1) | 1 : id[ i ]);
    }

    GlobalVector hInFather () const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = ctype( 1 ) / ctype( factor( i ) );
      return h;
    }

    GlobalVector originInFather ( unsigned int index ) const
    {
      GlobalVector origin;
      for( int i = 0; i < dimension; ++i )
      {
        unsigned int b = (refDir_ >> i) & 1;
        origin[ i ] = ctype( index & b ) / ctype( b+1 );
        index = index >> b;
      }
      return origin;
    }

  private:
    unsigned int refDir_;
  };



  template< class ct, int dim >
  inline std::ostream &
  operator<< ( std::ostream &out, const SPRefinement< ct, dim > &refinement )
  {
    return out << refinement.refDir_;
  }

  template< class ct, int dim >
  inline std::istream &
  operator>> ( std::istream &in, SPRefinement< ct, dim > &refinement )
  {
    return in >> refinement.refDir_;
  }

}

#endif // #ifndef DUNE_SPGRID_REFINEMENT_HH
