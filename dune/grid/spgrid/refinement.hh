#ifndef DUNE_SPGRID_REFINEMENT_HH
#define DUNE_SPGRID_REFINEMENT_HH

#include <algorithm>
#include <bitset>

#include <dune/common/fvector.hh>
#include <dune/common/iostream.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid/spgrid/declaration.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/multiindex.hh>

namespace Dune
{

  // SPIsotropicRefinementPolicy
  // ---------------------------

  template< int dim >
  class SPIsotropicRefinementPolicy
  {
    typedef SPIsotropicRefinementPolicy< dim > This;

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

    static std::string type () { return "isotropic"; }

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



  // SPAnisotropicRefinementPolicy
  // -----------------------------

  template< int dim >
  class SPAnisotropicRefinementPolicy
  {
    typedef SPAnisotropicRefinementPolicy< dim > This;

    friend class SPAnisotropicRefinement< dim >;

  public:
    static const int dimension = dim;

    SPAnisotropicRefinementPolicy ()
      : refDir_( (1u << dimension) - 1u )
    {}

    SPAnisotropicRefinementPolicy ( std::bitset< dimension > refDir )
      : refDir_( refDir )
    {}

    unsigned int weight () const
    {
      return refDir_.count();
    }

    unsigned int factor ( const int i ) const
    {
      assert( (i >= 0) && (i < dimension) );
      return refDir_[ i ]+1;
    }

    static std::string type () { return "anisotropic"; }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      return out << policy.refDir_.to_ulong();
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      unsigned long refDir;
      in >> refDir;
      if( !in.fail() )
        policy.refDir_ = refDir;
      return in;
    }

  private:
    std::bitset< dimension > refDir_;
  };



  // SPBisectionRefinementPolicy
  // ---------------------------

  template< int dim >
  class SPBisectionRefinementPolicy
  {
    typedef SPBisectionRefinementPolicy< dim > This;

    friend class SPBisectionRefinement< dim >;

  public:
    static const int dimension = dim;

    SPBisectionRefinementPolicy ()
      : dir_( -1 )
    {}

    explicit SPBisectionRefinementPolicy ( int dir )
      : dir_( dir )
    {
      if( (dir < 0) || (dir >= dimension) )
        DUNE_THROW( GridError, "Trying to create bisection refinement policy for invalid direction " << dir << "." );
    }

  private:
    SPBisectionRefinementPolicy ( const This &father, const This &policy )
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

    static std::string type () { return "bisection"; }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      assert( (policy.dir_ >= 0) && (policy.dir_ < dimension) );
      return out << (1u << policy.dir_);
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      unsigned long rd( 0 );
      in >> rd;
      std::bitset< dimension > refDir( rd );
      if( !in.fail() )
      {
        if( refDir.count() != 1 )
          DUNE_THROW( GridError, "Trying to create bisection refinement with multiple refined directions." );
        for( policy.dir_ = 0; !refDir[ policy.dir_ ]; ++policy.dir_ )
          continue;
      }
      return in;
    }

  private:
    int dir_;
  };



  // SPArbitratyRefinementPolicy
  // ---------------------------

  template< int dim >
  class SPArbitraryRefinementPolicy
  {
    typedef SPArbitraryRefinementPolicy< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    explicit SPArbitraryRefinementPolicy ( int factor = 2 )
    {
      if( factor <= 0 )
        DUNE_THROW( GridError, "Trying to create arbitrary refinement policy with non-positive factor." );
      for( int i = 0; i < dimension; ++i )
        factor_[ i ] = factor;
    }

    explicit SPArbitraryRefinementPolicy ( const MultiIndex &factor )
      : factor_( factor )
    {
      for( int i = 0; i < dimension; ++i )
      {
        if( factor_[ i ] <= 0 )
          DUNE_THROW( GridError, "Trying to create arbitrary refinement policy with non-positive factor." );
      }
    }

    constexpr unsigned int weight () const { return dimension; }

    unsigned int factor ( int i ) const { return factor_[ i ]; }

    static std::string type () { return "arbitrary"; }

    template< class char_type, class traits >
    friend std::basic_ostream< char_type, traits > &
    operator<< ( std::basic_ostream< char_type, traits > &out, const This &policy )
    {
      return out << policy.factor_;
    }

    template< class char_type, class traits >
    friend std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, This &policy )
    {
      return in >> policy.factor_;
    }

  private:
    MultiIndex factor_;
  };



  // SPDefaultRefinement
  // -------------------

  template< class P >
  class SPDefaultRefinement
  {
    typedef SPDefaultRefinement< P > This;

  public:
    typedef P Policy;

    static const int dimension = Policy::dimension;

    typedef SPMultiIndex< dimension > MultiIndex;

    explicit SPDefaultRefinement ( const Policy &policy = Policy() ) : policy_( policy ) {}

    SPDefaultRefinement ( const This &father, const Policy &policy ) : policy_( policy ) {}

    static std::string type () { return Policy::type(); }

    unsigned int factor ( int i ) const { return policy().factor( i ); }

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
        id[ i ] = (id[ i ] / factor( i )) | (id[ i ] & 1);
    }

    bool hasFather ( const MultiIndex &id ) const
    {
      bool hasFather = true;
      for( int i = 0; i < dimension; ++i )
        hasFather &= (id[ i ] & 1) | ((id[ i ] % (2*factor( i ))) == 0);
      return hasFather;
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
        id[ i ] = (factor( i ) * (id[ i ] & ~1)) | (id[ i ] & 1);
    }

    bool nextChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int alpha = factor( i );
        const int step = 2*(id[ i ] & 1);
        id[ i ] += step;
        if( ((id[ i ] % (2*alpha)) & ~1) != 0 )
          return true;
        id[ i ] -= alpha*step;
      }
      return false;
    }

    bool isCopy ( const MultiIndex id ) const
    {
      bool copy = true;
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int alpha = factor( i );
        copy &= (alpha == 1) | ((id[ i ] % (2*alpha)) == 0);
      }
      return copy;
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

    const Policy &policy () const { return policy_; }

  private:
    Policy policy_;
  };



  // SPBinaryRefinement
  // ------------------

  template< class P >
  class SPBinaryRefinement
    : public SPDefaultRefinement< P >
  {
    typedef SPBinaryRefinement< P > This;
    typedef SPDefaultRefinement< P > Base;

  public:
    using Base::dimension;

    typedef typename Base::Policy Policy;
    typedef typename Base::MultiIndex MultiIndex;

    explicit SPBinaryRefinement ( const Policy &policy = Policy() ) : Base ( policy ) {}

    SPBinaryRefinement ( const This &father, const Policy &policy ) : Base( father, policy ) {}

    using Base::factor;

    bool hasFather ( const MultiIndex &id ) const
    {
      bool hasFather = true;
      for( int i = 0; i < dimension; ++i )
        hasFather &= ((id[ i ] & 3) != 2);
      return hasFather;
    }

    void firstChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        const unsigned int alpha = factor( i );
        id[ i ] = alpha*id[ i ] - (alpha / 2)*(id[ i ] & 1);
      }
    }

    bool nextChild ( MultiIndex &id ) const
    {
      for( int i = 0; i < dimension; ++i )
      {
        if( factor( i ) < 2 )
          continue;
        id[ i ] ^= 2*(id[ i ] & 1);
        if( (id[ i ] & 2) != 0 )
          return true;
      }
      return false;
    }

    bool isCopy ( const MultiIndex id ) const
    {
      bool copy = true;
      for( int i = 0; i < dimension; ++i )
        copy &= (factor( i ) == 1) | ((id[ i ] & 3) == 0);
      return copy;
    }
  };



  // SPIsotropicRefinement
  // ---------------------

  /**
   * \class SPIsotropicRefinement
   * \brief each element is split into 2<sup>dim</sup> children.
   *
   * \note This is the default refinement technique.
   */
  template< int dim >
  class SPIsotropicRefinement
    : public SPBinaryRefinement< SPIsotropicRefinementPolicy< dim > >
  {
    typedef SPIsotropicRefinement< dim > This;
    typedef SPBinaryRefinement< SPIsotropicRefinementPolicy< dim > > Base;

  public:
    using Base::dimension;

    typedef typename Base::Policy Policy;

    SPIsotropicRefinement () : Base( Policy() ) {}

    explicit SPIsotropicRefinement ( const This &father, const Policy &policy ) : Base( policy ) {}

    constexpr unsigned int numChildren () const { return (1 << dimension); }
  };



  // SPAnisotropicRefinement
  // -----------------------

  /**
   * \class SPAnisotropicRefinement
   * \brief the user may choose freely along which axes the grid should be refined
   *
   * \note By default, this coincides with SPIsotropicRefinement.
   */
  template< int dim >
  class SPAnisotropicRefinement
    : public SPBinaryRefinement< SPAnisotropicRefinementPolicy< dim > >
  {
    typedef SPAnisotropicRefinement< dim > This;
    typedef SPBinaryRefinement< SPAnisotropicRefinementPolicy< dim > > Base;

  public:
    using Base::dimension;

    typedef typename Base::Policy Policy;

    SPAnisotropicRefinement () : Base( Policy() ) {}

    explicit SPAnisotropicRefinement ( const This &father, const Policy &policy ) : Base( policy ) {}

    using Base::policy;

    unsigned int numChildren () const { return (1u << policy().refDir_.count()); }
  };



  // SPBisectionRefinement
  // ---------------------

  /**
   * \class SPBisectionRefinement
   * \brief  each element is split into 2 children.
   *
   * \note The axis along which the elements are split may be chosen by the
   *       user.
   *       By default, the axes are cycled periodically.
   */
  template< int dim >
  class SPBisectionRefinement
    : public SPDefaultRefinement< SPBisectionRefinementPolicy< dim > >
  {
    typedef SPBisectionRefinement< dim > This;
    typedef SPDefaultRefinement< SPBisectionRefinementPolicy< dim > > Base;

  public:
    using Base::dimension;

    typedef typename Base::MultiIndex MultiIndex;
    typedef typename Base::Policy Policy;

    SPBisectionRefinement () : Base( Policy() ) {}

    SPBisectionRefinement ( const This &father, const Policy &policy )
      : Base( Policy( father.policy(), policy ) )
    {}

    using Base::policy;

    unsigned int numChildren () const
    {
      return 2;
    }

    void father ( MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] = (id[ dir ] / 2) | (id[ dir ] & 1);
    }

    bool hasFather ( const MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      return ((id[ dir ] & 3) != 2);
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
      id[ dir ] = 2*id[ dir ] - (id[ dir ] & 1);
    }

    bool nextChild ( MultiIndex &id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      id[ dir ] ^= 2*(id[ dir ] & 1);
      return ((id[ dir ] & 2) != 0);
    }

    bool isCopy ( const MultiIndex id ) const
    {
      const int dir = policy().dir_;
      assert( dir >= 0 );
      return ((id[ dir ] & 3) == 0);
    }
  };



  // SPArbitratyRefinement
  // ---------------------

  template< int dim >
  class SPArbitraryRefinement
    : public SPDefaultRefinement< SPArbitraryRefinementPolicy< dim > >
  {
    typedef SPArbitraryRefinement< dim > This;
    typedef SPDefaultRefinement< SPArbitraryRefinementPolicy< dim > > Base;

  public:
    typedef typename Base::Policy Policy;

    SPArbitraryRefinement () : Base( Policy() ) {}

    explicit SPArbitraryRefinement ( const This &father, const Policy &policy ) : Base( policy ) {}
  };

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_REFINEMENT_HH
