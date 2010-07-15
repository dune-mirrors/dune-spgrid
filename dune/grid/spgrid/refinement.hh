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
  struct SPBasicRefinement;



  // SPBasicRefinement (for SPIsotropicRefinement)
  // ---------------------------------------------

  template< int dim >
  class SPBasicRefinement< dim, SPIsotropicRefinement >
  {
    typedef SPBasicRefinement< dim, SPIsotropicRefinement > This;

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

    template< class ctype >
    FieldVector< ctype, dimension > originInFather ( const unsigned int index ) const
    {
      FieldVector< ctype, dimension > origin;
      for( int i = 0; i < dimension; ++i )
        origin[ i ] = ctype( (index >> i) & 1 ) / ctype( 2 );
      return origin;
    }

    static std::string type ();

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



  // SPBasicRefinement (for SPAnisotropicRefinement)
  // -----------------------------------------------

  template< int dim >
  class SPBasicRefinement< dim, SPAnisotropicRefinement >
  {
    typedef SPBasicRefinement< dim, SPAnisotropicRefinement > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    SPBasicRefinement ()
    : refDir_( (1 << dimension)-1 )
    {}

    explicit SPBasicRefinement ( const unsigned int refDir )
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

    static std::string type ();

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



  // SPRefinement
  // ------------

  template< int dim, SPRefinementStrategy strategy >
  class SPRefinement
  : public SPBasicRefinement< dim, strategy >
  {
    typedef SPRefinement< dim, strategy > This;
    typedef SPBasicRefinement< dim, strategy > Base;

  public:
    static const int dimension = Base::dimension;

    typedef typename Base::MultiIndex MultiIndex;
    //typedef typename Base::Policy Policy;

    using Base::factor;

    template< class Mesh >
    Mesh operator() ( const Mesh &mesh ) const;

    template< class Mesh >
    std::vector< Mesh > operator() ( const std::vector< Mesh > &coarseMeshes ) const;

    void father ( MultiIndex &id ) const;

    template< class ctype >
    FieldVector< ctype, dim > hInFather () const;
  };


  // Implementation of SPRefinement (for SPIsotropicRefinement)
  // ----------------------------------------------------------
  
  template< int dim >
  inline std::string SPBasicRefinement< dim, SPIsotropicRefinement >::type ()
  {
    return "isotropic";
  }



  // Implementation of SPRefinement (for SPAnisotropicRefinement)
  // ------------------------------------------------------------
  
  template< int dim >
  inline std::string SPBasicRefinement< dim, SPAnisotropicRefinement >::type ()
  {
    return "anisotropic";
  }



  // Implementation of SPRefinement
  // ------------------------------

  template< int dim, SPRefinementStrategy strategy >
  template< class Mesh >
  inline Mesh SPRefinement< dim, strategy >::operator() ( const Mesh &mesh ) const
  {
    return mesh.refine( *this );
  }


  template< int dim, SPRefinementStrategy strategy >
  template< class Mesh >
  inline std::vector< Mesh > SPRefinement< dim, strategy >
    ::operator() ( const std::vector< Mesh > &coarseMeshes ) const
  {
    std::vector< Mesh > fineMeshes;
    const size_t size = coarseMeshes.size();
    fineMeshes.reserve( size );
    for( size_t i = 0; i < size; ++i )
      fineMeshes.push_back( coarseMeshes[ i ].refine( *this ) );
    return fineMeshes;
  }


  template< int dim, SPRefinementStrategy strategy >
  inline void SPRefinement< dim, strategy >::father ( MultiIndex &id ) const
  {
    for( int i = 0; i < dimension; ++i )
    {
      assert( (id[ i ] & 1) != 0 );
      id[ i ] = (id[ i ] / factor( i )) | 1;
    }
  }


  template< int dim, SPRefinementStrategy strategy >
  template< class ctype >
  inline FieldVector< ctype, dim >
  SPRefinement< dim, strategy >::hInFather () const
  {
    FieldVector< ctype, dimension > h;
    for( int i = 0; i < dimension; ++i )
      h[ i ] = ctype( 1 ) / ctype( factor( i ) );
    return h;
  }



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
