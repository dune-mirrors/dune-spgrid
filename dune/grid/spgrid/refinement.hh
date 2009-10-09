#ifndef DUNE_SPGRID_REFINEMENT_HH
#define DUNE_SPGRID_REFINEMENT_HH

#include <dune/common/fvector.hh>

#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPRefinementStragety

  enum SPRefinementStrategy
  {
    SPIsotropicRefinement,
    SPAnisotropicRefinement
  };



  // Internal Forward Declarations
  // -----------------------------

  template< class ct, int dim, SPRefinementStrategy strategy >
  struct SPRefinement;

  template< class ct, int dim >
  std::ostream &operator<< ( std::ostream &, const SPRefinement< ct, dim, SPAnisotropicRefinement > & );



  // SPRefinement
  // ------------

  template< class ct, int dim >
  class SPRefinement< ct, dim, SPIsotropicRefinement >
  {
    typedef SPRefinement< ct, dim, SPIsotropicRefinement > This;

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
    typedef SPMultiIndex< dimension > MultiIndex;

    SPRefinement ()
    {}

    explicit SPRefinement ( const unsigned int refDir )
    {
      if( refDir >= (1 << dimension) )
        DUNE_THROW( GridError, "Trying to create isotropic refinement from invalid value " << refDir << "." );
      if( refDir < (1 << dimension)-1 )
        DUNE_THROW( GridError, "Trying to create isotropic refinement from anisotropic value " << refDir << "." );
    }

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

    GlobalVector hInFather () const
    {
      GlobalVector h;
      for( int i = 0; i < dimension; ++i )
        h[ i ] = ctype( 1 ) / ctype( factor( i ) );
      return h;
    }

    GlobalVector originInFather ( const unsigned int index ) const
    {
      GlobalVector origin;
      for( int i = 0; i < dimension; ++i )
        origin[ i ] = ctype( (index >> i) & 1 ) / ctype( 2 );
      return origin;
    }

  private:
    unsigned int refDir_;
  };



  template< class ct, int dim >
  class SPRefinement< ct, dim, SPAnisotropicRefinement >
  {
    typedef SPRefinement< ct, dim, SPAnisotropicRefinement > This;

    friend std::ostream &operator<<<> ( std::ostream &, const This & );

  public:
    typedef ct ctype;

    static const int dimension = dim;

    typedef FieldVector< ctype, dimension > GlobalVector;
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
  operator<< ( std::ostream &out, const SPRefinement< ct, dim, SPIsotropicRefinement > &refinement )
  {
    const unsigned int refDir = (1 << dim)-1;
    return out << refDir;
  }

  template< class ct, int dim >
  inline std::ostream &
  operator<< ( std::ostream &out, const SPRefinement< ct, dim, SPAnisotropicRefinement > &refinement )
  {
    return out << refinement.refDir_;
  }

  template< class ct, int dim, SPRefinementStrategy strategy >
  inline std::istream &
  operator>> ( std::istream &in, SPRefinement< ct, dim, strategy > &refinement )
  {
    unsigned int refDir;
    in >> refDir;
    if( in )
      refinement = SPRefinement< ct, dim, strategy >( refDir );
    return in;
  }

}

#endif // #ifndef DUNE_SPGRID_REFINEMENT_HH
