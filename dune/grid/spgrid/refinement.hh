#ifndef DUNE_SPGRID_REFINEMENT_HH
#define DUNE_SPGRID_REFINEMENT_HH

#include <dune/common/fvector.hh>
#include <dune/common/iostream.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< int dim, SPRefinementStrategy strategy >
  class SPRefinementPolicy;

  template< int dim, SPRefinementStrategy strategy >
  class SPBasicRefinement;



  // SPRefinementPolicy (for SPIsotropicRefinement)
  // ----------------------------------------------

  template< int dim >
  class SPRefinementPolicy< dim, SPIsotropicRefinement >
  {
    typedef SPRefinementPolicy< dim, SPIsotropicRefinement > This;

  public:
    static const int dimension = dim;

    unsigned int weight () const
    {
      return dimension;
    }

    unsigned int factor ( const int i ) const
    {
      return 2;
    }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      return out << ((1 << dimension)-1);
    }
    
    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      in >> match( (1 << dimension)-1 );
      return in;
    }
  };



  // SPRefinementPolicy (for SPAnisotropicRefinement)
  // ------------------------------------------------

  template< int dim >
  class SPRefinementPolicy< dim, SPAnisotropicRefinement >
  {
    typedef SPRefinementPolicy< dim, SPAnisotropicRefinement > This;

    friend class SPBasicRefinement< dim, SPAnisotropicRefinement >;

  public:
    static const int dimension = dim;

    SPRefinementPolicy ()
    : refDir_( (1 << dimension)-1 )
    {}

    explicit SPRefinementPolicy ( const unsigned int refDir )
    : refDir_( refDir )
    {
      if( refDir >= (1 << dimension) )
        DUNE_THROW( GridError, "Trying to create anisotropic refinement policy from invalid value " << refDir << "." );
    }

    unsigned int weight () const
    {
      return bitCount( refDir_ );
    }

    unsigned int factor ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return ((refDir_ >> i) & 1)+1;
    }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      return out << policy.refDir_;
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      unsigned int refDir;
      in >> refDir;
      if( !in.fail() )
        policy = This( refDir );
      return in;
    }

  private:
    unsigned int refDir_;
  };



  // SPRefinementPolicy (for SPBisectionRefinement)
  // ----------------------------------------------

  template< int dim >
  class SPRefinementPolicy< dim, SPBisectionRefinement >
  {
    typedef SPRefinementPolicy< dim, SPBisectionRefinement > This;

    friend class SPBasicRefinement< dim, SPBisectionRefinement >;

  public:
    static const int dimension = dim;

    SPRefinementPolicy ()
    : dir_( -1 )
    {}

    explicit SPRefinementPolicy ( const int dir )
    : dir_( dir )
    {
      if( (dir < 0) || (dir >= dimension) )
        DUNE_THROW( GridError, "Trying to create bisection refinement policy for invalid direction " << dir << "." );
    }

  private:
    SPRefinementPolicy ( const This &father, const This &policy )
    : dir_( policy.dir_ < 0 ? (father.dir_ + 1) % dimension : policy.dir_ )
    {}

  public:
    unsigned int weight () const
    {
      return 1;
    }

    unsigned int factor ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return (i == dir_ ? 2 : 1);
    }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      assert( (policy.dir_ >= 0) && (policy.dir_ < dimension) );
      return out << (1 << policy.dir_);
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      unsigned int refDir;
      in >> refDir;
      if( !in.fail() )
      {
        if( bitCount( refDir ) != 1 )
          DUNE_THROW( GridError, "Trying to create bisection refinement with multiple refined directions." );
        int dir = -1;
        for( ; refDir != 0; ++dir )
          refDir = refDir >> 1;
        policy = This( dir );
      }
      return in;
    }

  private:
    int dir_;
  };



  // SPDefaultRefinement
  // -------------------

  template< int dim, SPRefinementStrategy strategy >
  struct SPDefaultRefinement
  {
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    typedef SPRefinementPolicy< dim, strategy > Policy;

  protected:
    explicit SPDefaultRefinement ( const Policy &policy )
    : policy_( policy )
    {}

  public:
    unsigned int factor ( const int i ) const
    {
      return policy().factor( i );
    }

    unsigned int numChildren () const
    {
      unsigned int numChildren = 1;
      for( int i = 0; i < dimension; ++i )
        numChildren *= factor( i );
      return numChildren;
    }

    void father ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
        id[ i ] = ((id[ i ] / factor( i )) & ~1) | (id[ i ] & 1);
    }

    void child ( MultiIndex &id, unsigned int index ) const
    {
      assert( index < numChildren() );
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int alpha = factor( i );
        id[ i ] = ((id[ i ] + (index % alpha)) * alpha) - (alpha - 1);
        index /= alpha;
      }
    }

    unsigned int childIndex ( const MultiIndex &id ) const
    {
      unsigned int index = 0;
      for( int i = dimension-1; i >= 0; --i )
      {
        const unsigned int alpha = factor( i );
        index = (index * alpha) + ((id[ i ] >> 1) % alpha);
      }
      assert( index < numChildren() );
      return index;
    }

    void firstChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int alpha = factor( i );
        id[ i ] = alpha*id[ i ] - (alpha / 2);
      }
    }

    bool nextChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        if( factor( i ) < 2 )
          continue;
        id[ i ] ^= 2;
        if( (id[ i ] & 2) != 0 )
          return true;
      }
      return false;
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
        const unsigned int alpha = factor( i );
        origin[ i ] = ctype( index % alpha ) / ctype( alpha );
        index /= alpha;
      }
      return origin;
    }

    bool isCopy ( const MultiIndex id ) const
    {
      bool copy = true;
      for( int i = 0; i < dimension; ++i )
        copy &= (factor( i ) == 1) | ((id[ i ] & 3) == 0);
      return copy;
    }

    const Policy &policy () const
    {
      return policy_;
    }

  private:
    Policy policy_;
  };



  // SPBasicRefinement (for SPIsotropicRefinement)
  // ---------------------------------------------

  template< int dim >
  class SPBasicRefinement< dim, SPIsotropicRefinement >
  : public SPDefaultRefinement< dim, SPIsotropicRefinement >
  {
    typedef SPBasicRefinement< dim, SPIsotropicRefinement > This;
    typedef SPDefaultRefinement< dim, SPIsotropicRefinement > Base;

  public:
    using Base::dimension;
    
    typedef typename Base::MultiIndex MultiIndex;
    typedef typename Base::Policy Policy;

  protected:
    SPBasicRefinement ()
    : Base( Policy() )
    {}

    explicit SPBasicRefinement ( const This &father, const Policy &policy )
    : Base( policy )
    {}

  public:
    unsigned int numChildren () const
    {
      return (1 << dimension);
    }

    static std::string type () { return "isotropic"; }
  };



  // SPBasicRefinement (for SPAnisotropicRefinement)
  // -----------------------------------------------

  template< int dim >
  class SPBasicRefinement< dim, SPAnisotropicRefinement >
  : public SPDefaultRefinement< dim, SPAnisotropicRefinement >
  {
    typedef SPBasicRefinement< dim, SPAnisotropicRefinement > This;
    typedef SPDefaultRefinement< dim, SPAnisotropicRefinement > Base;

  public:
    using Base::dimension;

    typedef typename Base::MultiIndex MultiIndex;
    typedef typename Base::Policy Policy;

  protected:
    SPBasicRefinement ()
    : Base( Policy() )
    {}

    explicit SPBasicRefinement ( const This &father, const Policy &policy )
    : Base( policy )
    {}

  public:
    using Base::policy;

    unsigned int numChildren () const
    {
      return (1 << bitCount( policy().refDir_ ));
    }

    static std::string type () { return "anisotropic"; }
  };



  // SPBasicRefinement (for SPBisectionRefinement)
  // ---------------------------------------------

  template< int dim >
  class SPBasicRefinement< dim, SPBisectionRefinement >
  : public SPDefaultRefinement< dim, SPBisectionRefinement >
  {
    typedef SPBasicRefinement< dim, SPBisectionRefinement > This;
    typedef SPDefaultRefinement< dim, SPBisectionRefinement > Base;

  public:
    using Base::dimension;

    typedef typename Base::MultiIndex MultiIndex;
    typedef typename Base::Policy Policy;

  protected:
    SPBasicRefinement ()
    : Base( Policy() )
    {}

    explicit SPBasicRefinement ( const This &father, const Policy &policy )
    : Base( Policy( father.policy(), policy ) )
    {}

  public:
    using Base::policy;

    unsigned int numChildren () const
    {
      return 2;
    }

    void father ( MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] = ((id[ dir ] / 2) & ~1) | (id[ dir ] & 1);
    }

    void child ( MultiIndex &id, unsigned int index ) const
    {
      assert( index < numChildren() );
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] = 2*(id[ dir ] + index) - 1;
    }

    unsigned int childIndex ( const MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      const unsigned int index = (id[ dir ] >> 1) % 2;
      assert( index < numChildren() );
      return index;
    }

    void firstChild ( MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] = 2*id[ dir ] - 1;
    }

    bool nextChild ( MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] ^= 2;
      return ((id[ dir ] & 2) != 0);
    }

    bool isCopy ( const MultiIndex id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      return ((id[ dir ] & 3) == 0);
    }

    static std::string type () { return "bisection"; }
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
    typedef typename Base::Policy Policy;

    SPRefinement ()
    {}

    explicit SPRefinement ( const This &father, const Policy &policy )
    : Base( father, policy )
    {}

    template< class Mesh >
    Mesh operator() ( const Mesh &mesh ) const;

    template< class Mesh >
    std::vector< Mesh > operator() ( const std::vector< Mesh > &coarseMeshes ) const;
  };



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
