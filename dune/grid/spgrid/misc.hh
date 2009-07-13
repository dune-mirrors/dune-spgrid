#ifndef DUNE_SPGRID_MISC_HH
#define DUNE_SPGRID_MISC_HH

namespace Dune
{

  inline unsigned int bitcount ( const unsigned int i )
  {
    unsigned int bitcount = 0;
    for( ; i != 0; i /= 2 )
      bitcount += (i & 1);
    return bitcount;
  }

}

#endif // #ifndef DUNE_SPGRID_MISC_HH
