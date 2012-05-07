#ifndef DUNE_SPGRID_GRIDLEVEL_HH
#define DUNE_SPGRID_GRIDLEVEL_HH

#include <cassert>
#include <vector>

#include <dune/grid/spgrid/geometricgridlevel.hh>
#include <dune/grid/spgrid/misc.hh>
#include <dune/grid/spgrid/refinement.hh>
#include <dune/grid/spgrid/domain.hh>
#include <dune/grid/spgrid/mesh.hh>
#include <dune/grid/spgrid/partitionpool.hh>
#include <dune/grid/spgrid/linkage.hh>
#include <dune/grid/spgrid/decomposition.hh>

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template< int mydim, int cdim, class Grid >
  class SPLocalGeometry;



  // SPGridLevel
  // -----------

  template< class Grid >
  class SPGridLevel
  : public SPGeometricGridLevel< typename remove_const< Grid >::type::ctype, remove_const< Grid >::type::dimension >
  {
    typedef SPGridLevel< Grid > This;
    typedef SPGeometricGridLevel< typename remove_const< Grid >::type::ctype, remove_const< Grid >::type::dimension > Base;

  public:
    typedef SPGridLevel< Grid > GridLevel;

    typedef typename Base::ReferenceCube ReferenceCube;
    typedef typename Base::ctype ctype;

    static const int dimension = Base::dimension;
    static const unsigned int numDirections = Base::numDirections;

    typedef typename remove_const< Grid >::type::Traits Traits;

    typedef typename Traits::Domain Domain;
    typedef typename Traits::Refinement Refinement;
    typedef typename Traits::RefinementPolicy RefinementPolicy;

    typedef typename ReferenceCube::GlobalVector GlobalVector;
    typedef typename ReferenceCube::MultiIndex MultiIndex;

    typedef SPDecomposition< dimension > Decomposition;
    typedef SPPartitionPool< dimension > PartitionPool;
    typedef SPLinkage< dimension > Linkage;

    typedef typename Decomposition::Mesh Mesh;

    typedef typename PartitionPool::PartitionList PartitionList;

    typedef typename Linkage::Interface CommInterface;

    typedef typename Traits::template Codim< 0 >::LocalGeometry LocalGeometry;

  private:
    typedef SPLocalGeometry< dimension, dimension, Grid > LocalGeometryImpl;

  public:
    SPGridLevel ( const Grid &grid, const Decomposition &decomposition );
    SPGridLevel ( const GridLevel &father, const RefinementPolicy &policy );
    SPGridLevel ( const This &other );

    ~SPGridLevel ();

    using Base::referenceCube;

    const Grid &grid () const { return *grid_; }
    int level () const { return level_; }
    const Domain &domain () const { return domain_; }
    const Refinement &refinement () const { return refinement_; }

    const Mesh &globalMesh () const;
    const Mesh &localMesh () const;

    template< PartitionIteratorType pitype >
    const PartitionList &partition () const;

    template< int codim >
    PartitionType
    partitionType ( const MultiIndex &id, const unsigned int partitionNumber ) const;

    const CommInterface &commInterface ( const InterfaceType iftype ) const;

    MultiIndex macroId ( const MultiIndex &id ) const;

    size_t boundaryIndex ( const MultiIndex &id,
                           const unsigned int partitionNumber,
                           const int face ) const;

    LocalGeometry geometryInFather ( const MultiIndex &id ) const;

    int size () const;

  private:
    void buildLocalGeometry ();

    static MultiIndex coarseMacroFactor ();
    static GlobalVector meshWidth ( const Domain &domain, const Mesh &mesh );

    MultiIndex overlap () const;

    const Grid *grid_;
    int level_;

    const Refinement refinement_;
    MultiIndex macroFactor_;
    Domain domain_;
    std::vector< Mesh > decomposition_;
    Mesh localMesh_;
    PartitionPool partitionPool_;
    Linkage linkage_;

    LocalGeometryImpl **geometryInFather_;
  };



  // Implementation of SPGridLevel
  // -----------------------------

  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( const Grid &grid, const Decomposition &decomposition )
  : Base( grid.refCubes_, meshWidth( grid.domain(), decomposition.mesh() ) ),
    grid_( &grid ),
    level_( 0 ),
    refinement_(),
    macroFactor_( coarseMacroFactor() ),
    domain_( grid.domain() ),
    decomposition_( decomposition.subMeshes() ),
    localMesh_( decomposition_[ grid.comm().rank() ] ),
    partitionPool_( localMesh_, decomposition.mesh(), overlap(), domain_.topology() ),
    linkage_( grid.comm().rank(), partitionPool_, decomposition_ )
  {
    buildLocalGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >
    ::SPGridLevel ( const GridLevel &father, const RefinementPolicy &policy )
  : Base( father.grid().refCubes_, meshWidth( father.domain(), Refinement( father.refinement(), policy )( father.globalMesh() ) ) ),
    grid_( father.grid_ ),
    level_( father.level() + 1 ),
    refinement_( father.refinement(), policy ),
    macroFactor_( father.macroFactor_ * refinement_ ),
    domain_( father.domain() ),
    decomposition_( refinement_( father.decomposition_ ) ),
    localMesh_( refinement_( father.localMesh() ) ),
    partitionPool_( localMesh_, refinement_( father.globalMesh() ), overlap(), domain_.topology() ),
    linkage_( father.grid().comm().rank(), partitionPool_, decomposition_ )
  {
    buildLocalGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >::SPGridLevel ( const This &other )
  : Base( other ),
    grid_( other.grid_ ),
    refinement_( other.refinement_ ),
    macroFactor_( other.macroFactor_ ),
    domain_( other.domain_ ),
    decomposition_( other.decomposition_ ),
    localMesh_( other.localMesh_ ),
    partitionPool_( other.partitionPool_ ),
    linkage_( other.linkage_ )
  {
    buildLocalGeometry();
  }


  template< class Grid >
  inline SPGridLevel< Grid >::~SPGridLevel ()
  {
    if( geometryInFather_ )
    {
      unsigned int numChildren = refinement().numChildren();
      for( unsigned int index = 0; index < numChildren; ++index )
        delete geometryInFather_[ index ];
      delete geometryInFather_;
    }
  }

  template< class Grid >
  inline const typename SPGridLevel< Grid >::Mesh &
  SPGridLevel< Grid >::globalMesh () const
  {
    return partitionPool_.globalMesh();
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::Mesh &
  SPGridLevel< Grid >::localMesh () const
  {
    return localMesh_;
  }


  template< class Grid >
  template< PartitionIteratorType pitype >
  inline const typename SPGridLevel< Grid >::PartitionList &
  SPGridLevel< Grid >::partition () const
  {
    return partitionPool_.template get< pitype >();
  }


  template< class Grid >
  template< int codim >
  inline PartitionType SPGridLevel< Grid >
    ::partitionType ( const MultiIndex &id, const unsigned int partitionNumber ) const
  {
    return partitionPool_.template partitionType< codim >( id, partitionNumber );
  }


  template< class Grid >
  inline const typename SPGridLevel< Grid >::CommInterface &
  SPGridLevel< Grid >::commInterface ( const InterfaceType iftype ) const
  {
    return linkage_.interface( iftype );
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::macroId ( const MultiIndex &id ) const
  {
    MultiIndex macroId;
    for( int i = 0; i < dimension; ++i )
      macroId[ i ] = (((id[ i ] >> 1) / macroFactor_[ i ]) << 1) | (id[ i ] & 1);
    return macroId;
  }


  template< class Grid >
  inline size_t SPGridLevel< Grid >
    ::boundaryIndex ( const MultiIndex &id,
                      const unsigned int partitionNumber,
                      const int face ) const
  {
    // note: boundaryIndex ignores the last bit of macroId,
    //       hence we can use this fast computation
    MultiIndex macroId;
    for( int i = 0; i < dimension; ++i )
      macroId[ i ] = id[ i ] / macroFactor_[ i ];
    return grid().boundaryIndex( macroId, partitionNumber, face );
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::LocalGeometry
  SPGridLevel< Grid >::geometryInFather ( const MultiIndex &id ) const
  {
    assert( (level() > 0) && (geometryInFather_ != 0) );
    return LocalGeometry( *(geometryInFather_[ refinement().childIndex( id ) ]) );
  }


  template< class Grid >
  inline int SPGridLevel< Grid >::size () const
  {
    return globalMesh().volume();
  }


  template< class Grid >
  inline void SPGridLevel< Grid >::buildLocalGeometry ()
  {
    geometryInFather_ = 0;
    if( level() > 0 )
    {
      const unsigned int numChildren = refinement().numChildren();
      geometryInFather_ = new LocalGeometryImpl *[ numChildren ];
      const GlobalVector hInFather = refinement().template hInFather< ctype >();
      const typename Base::template Codim< 0 >::GeometryCache cacheInFather( hInFather, numDirections-1 );
      for( unsigned int index = 0; index < numChildren; ++index )
      {
        const GlobalVector origin = refinement().template originInFather< ctype >( index );
        geometryInFather_[ index ] = new LocalGeometryImpl( cacheInFather, origin );
      }
    }
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::coarseMacroFactor ()
  {
    MultiIndex macroFactor;
    for( int i = 0; i < dimension; ++i )
      macroFactor[ i ] = 1;
    return macroFactor;
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::GlobalVector
  SPGridLevel< Grid >::meshWidth ( const Domain &domain, const Mesh &mesh )
  {
    GlobalVector h = domain.cube().width();
    const MultiIndex meshWidth = mesh.width();
    for( int i = 0; i < dimension; ++i )
      h[ i ] /= ctype( meshWidth[ i ] );
    return h;
  }


  template< class Grid >
  inline typename SPGridLevel< Grid >::MultiIndex
  SPGridLevel< Grid >::overlap () const
  {
    MultiIndex overlap;
    for( int i = 0; i < dimension; ++i )
      overlap[ i ] = macroFactor_[ i ] * grid().overlap()[ i ];
    return overlap;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_GRIDLEVEL_HH
