#ifndef DUNE_COMMON_IOSTREAM_HH
#define DUNE_COMMON_IOSTREAM_HH

#include <iostream>

namespace Dune
{

  namespace iostream
  {

    template< class char_type >
    struct Match
    {
      explicit Match ( const char_type &c )
      : c_( c )
      {}

      template< class T >
      bool operator() ( const T &c ) const
      {
        return (c_ == c);
      }

    private:
      char_type c_;
    };

  }


  template< class char_type >
  iostream::Match< char_type > match ( const char_type &c )
  {
    return iostream::Match< char_type >( c );
  }


  template< class char_type, class Traits >
  inline std::basic_istream< char_type, Traits > &
  operator>> ( std::basic_istream< char_type, Traits > &in, const iostream::Match &match )
  {
    char c;
    in >> c;
    if( !match( c ) )
      in.clear( std::ios_base::badbit );
  }

}

#endif // #ifndef DUNE_COMMON_IOSTREAM_HH
