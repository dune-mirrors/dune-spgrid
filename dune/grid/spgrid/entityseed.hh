#ifndef DUNE_SPGRID_ENTITYSEED_HH
#define DUNE_SPGRID_ENTITYSEED_HH

namespace Dune
{

  template< int codim, class G >
  struct SPEntitySeed
  {
    typedef typename remove_const< G >::type Grid;

    typedef typename Grid::Traits::ReferenceCube ReferenceCube;

    static const int dimension = ReferenceCube::dimension;
    static const int codimension = codim;
    static const int mydimension = dimension - codimension;
    static const int dimensionworld = dimension;

    typedef typename Grid::Traits::template Codim< codimension >::Entity Entity;

    typedef typename ReferenceCube::MultiIndex MultiIndex;

    SPEntitySeed ( const int level, const MultiIndex &id, const unsigned int partitionNumber )
    : level_( level ), id_( id ), partitionNumber_( partitionNumber )
    {}

    int level () const { return level_; }
    MultiIndex id () const { return id_; }
    unsigned int partitionNumber () const { return partitionNumber_; }

  private:
    int level_;
    MultiIndex id_;
    unsigned int partitionNumber_;
  };

}

#endif // #ifndef DUNE_SPGRID_ENTITYSEED_HH
