#ifndef DUNE_SPGRID_MISC_HH
#define DUNE_SPGRID_MISC_HH

namespace Dune
{

  inline unsigned int bitCount ( unsigned int i )
  {
    unsigned int bitCount = 0;
    for( ; i != 0; i /= 2 )
      bitCount += (i & 1);
    return bitCount;
  }

}

#endif // #ifndef DUNE_SPGRID_MISC_HH
