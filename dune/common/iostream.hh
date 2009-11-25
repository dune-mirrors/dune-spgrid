#ifndef DUNE_COMMON_IOSTREAM_HH
#define DUNE_COMMON_IOSTREAM_HH

#include <iostream>

namespace Dune
{

  namespace iostream
  {

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

  }


  template< class T >
  iostream::Match< T > match ( const T &value )
  {
    return iostream::Match< T >( value );
  }


  template< class char_type, class traits, class T >
  inline std::basic_istream< char_type, traits > &
  operator>> ( std::basic_istream< char_type, traits > &in, const iostream::Match< T > &match )
  {
    T value;
    in >> value;
    if( !match( value ) )
      in.clear( std::ios_base::failbit );
    return in;
  }

}

#endif // #ifndef DUNE_COMMON_IOSTREAM_HH
