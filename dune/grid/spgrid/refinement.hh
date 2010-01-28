#ifndef DUNE_SPGRID_REFINEMENT_HH
#define DUNE_SPGRID_REFINEMENT_HH

#include <dune/common/fvector.hh>
#include <dune/common/iostream.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPRefinementStrategy
  // --------------------

  enum SPRefinementStrategy
  {
    SPIsotropicRefinement,
    SPAnisotropicRefinement
  };



  // Internal Forward Declarations
  // -----------------------------

  template< int dim, SPRefinementStrategy strategy >
  struct SPRefinement;



  // SPRefinement
  // ------------

  template< int dim >
  class SPRefinement< dim, SPIsotropicRefinement >
  {
    typedef SPRefinement< dim, SPIsotropicRefinement > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    unsigned int factor ( const int i ) const
    {
      return 2;
    }

    unsigned int numChildren () const
    {
      return (1 << dimension);
    }

    void child ( MultiIndex &id, unsigned int index ) const
    {
      assert( index < numChildren() );
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*(id[ i ] + ((index >> i) & 1)) - 1;
    }

    unsigned int childIndex ( const MultiIndex &id ) const
    {
      unsigned int index = 0;
      for( int i = 0; i < dimension; ++i )
        index |= (((id[ i ] >> 1) & 1) << i);
      assert( index < numChildren() );
      return index;
    }

    void firstChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = 2*id[ i ] - 1;
    }

    bool nextChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        id[ i ] ^= 2;
        if( (id[ i ] & 2) != 0 )
          return true;
      }
      return false;
    }

    void father ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = (id[ i ] >> 1) | 1;
    }

    template< class ctype >
    FieldVector< ctype, dimension > hInFather () const
    {
      FieldVector< ctype, dimension > h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = ctype( 1 ) / ctype( factor( i ) );
      return h;
    }

    template< class ctype >
    FieldVector< ctype, dimension > originInFather ( const unsigned int index ) const
    {
      FieldVector< ctype, dimension > origin;
      for( int i = 0; i < dimension; ++i )
        origin[ i ] = ctype( (index >> i) & 1 ) / ctype( 2 );
      return origin;
    }

    static std::string type ()
    {
      return "isotropic";
    }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &refinement )
    {
      return out << ((1 << dimension)-1);
    }
    
    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &refinement )
    {
      in >> match( (1 << dimension)-1 );
      return in;
    }
  };



  template< int dim >
  class SPRefinement< dim, SPAnisotropicRefinement >
  {
    typedef SPRefinement< dim, SPAnisotropicRefinement > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    SPRefinement ()
    : refDir_( (1 << dimension)-1 )
    {}

    explicit SPRefinement ( const unsigned int refDir )
    : refDir_( refDir )
    {
      if( refDir >= (1 << dimension) )
        DUNE_THROW( GridError, "Trying to create anisotropic refinement from invalid value " << refDir << "." );
    }

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

    template< class ctype >
    FieldVector< ctype, dimension > hInFather () const
    {
      FieldVector< ctype, dimension > h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = ctype( 1 ) / ctype( factor( i ) );
      return h;
    }

    template< class ctype >
    FieldVector< ctype, dimension > originInFather ( unsigned int index ) const
    {
      FieldVector< ctype, dimension > origin;
      for( int i = 0; i < dimension; ++i )
      {
        unsigned int b = (refDir_ >> i) & 1;
        origin[ i ] = ctype( index & b ) / ctype( b+1 );
        index = index >> b;
      }
      return origin;
    }

    static std::string type ()
    {
      return "anisotropic";
    }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &refinement )
    {
      return out << refinement.refDir_;
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &refinement )
    {
      unsigned int refDir;
      in >> refDir;
      if( !in.fail() )
        refinement = This( refDir );
      return in;
    }

  private:
    unsigned int refDir_;
  };



  // Auxilliary Functions for SPRefinement
  // -------------------------------------

  template< int dim, SPRefinementStrategy strategy >
  inline SPMultiIndex< dim >
  operator* ( const SPMultiIndex< dim > &width,
              const SPRefinement< dim, strategy > &refinement )
  {
    SPMultiIndex< dim > result;
    for( int i = 0; i < dim; ++i )
      result[ i ] = width[ i ] * refinement.factor( i );
    return result;
  }

}

#endif // #ifndef DUNE_SPGRID_REFINEMENT_HH
