#ifndef DUNE_COMMON_IOSTREAM_HH
#define DUNE_COMMON_IOSTREAM_HH

#include <iostream>

namespace Dune
{

  namespace iostream
  {

    template< class T >
    struct MatchTraits
    {
      typedef T Type;
    };

    template< class T >
    struct MatchTraits< const T >
    {
      typedef const typename MatchTraits< T >::Type Type;
    };

    template< int n >
    struct MatchTraits< char[ n ] >
    {
      typedef std::string Type;
    };


    template< class T >
    struct Match
    {
      explicit Match ( const T &value )
      : value_( value )
      {}

      template< class U >
      Match ( const Match< U > &other )
      : value_( other.value_ )
      {}

      bool operator() ( const T &value ) const
      {
        return (value_ == value);
      }

    private:
      T value_;
    };



    template< class char_type, class traits, class T >
    inline std::basic_istream< char_type, traits > &
    operator>> ( std::basic_istream< char_type, traits > &in, const Match< T > &match )
    {
      T value;
      in >> value;
      if( !match( value ) )
        in.clear( std::ios_base::failbit );
      return in;
    }

  } // namespace iostream


  template< class char_type, class traits >
  inline bool isGood ( std::basic_istream< char_type, traits > &in )
  {
    bool good = in.good();
    if( good )
    {
      char_type c;
      in >> c;
      good = !in.fail();
      if( good )
        in.unget();
      in.clear();
    }
    return good;
  }


  template< class T >
  inline iostream::Match< typename iostream::MatchTraits< T >::Type >
  match ( const T &value )
  {
    return iostream::Match< typename iostream::MatchTraits< T >::Type >( value );
  }

} // namespace Dune

#endif // #ifndef DUNE_COMMON_IOSTREAM_HH
