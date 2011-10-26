#ifndef DUNE_CARTESIANGRID_DGFPARSER_HH
#define DUNE_CARTESIANGRID_DGFPARSER_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridfactory.hh>

#if HAVE_ALUGRID
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif // #if HAVE_ALUGRID

#include <dune/grid/spgrid/dgfparser.hh>

#include <dune/grid/cartesiangrid/grid.hh>
#include <dune/grid/cartesiangrid/hostgridaccess.hh>

namespace Dune
{

  // DGFGridFactory for CartesianGrid
  // --------------------------------

  template< class HostGrid >
  struct DGFGridFactory< CartesianGrid< HostGrid > >
  {
    typedef CartesianGrid< HostGrid > Grid;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;

    static const int dimension = Grid::dimension;
    typedef typename Grid::template Codim< 0 >::Entity Element;
    typedef typename Grid::template Codim< dimension >::Entity Vertex;

    typedef FieldVector< double, dimension > Point;
    typedef dgf::BoundaryDomBlock BoundaryDomainBlock;

    explicit DGFGridFactory ( std::istream &input,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      generate( input, comm );
    }

    explicit DGFGridFactory ( const std::string &filename,
                              MPICommunicatorType comm = MPIHelper::getCommunicator() )
    {
      std::ifstream input( filename.c_str() );
      if( !input )
        DUNE_THROW( DGFException, "Unable to open file: " << filename << "." );
      generate( input, comm, filename );
      input.close();
    }

    ~DGFGridFactory ()
    {
      delete boundaryDomainBlock_;
    }

    Grid *grid () const
    {
      return grid_;
    }

    template< class Intersection >
    bool wasInserted ( const Intersection &intersection ) const
    {
      return false;
    }

    template< class Intersection >
    int boundaryId ( const Intersection &intersection ) const
    {
      return Grid::getRealImplementation( intersection ).boundaryId();
    }

    template< int codim >
    int numParameters () const
    {
      return 0;
    }

    template< class Entity >
    std::vector< double > &parameter ( const Entity &entity )
    {
      DUNE_THROW( InvalidStateException,
                  "Calling DGFGridFactory::parameter is only allowed if there are parameters." );
    }

    bool haveBoundaryParameters () const
    {
      return boundaryDomainBlock_->hasParameter();
    }

    template< class Intersection >
    const typename DGFBoundaryParameter::type &
    boundaryParameter ( const Intersection &intersection ) const
    {
      std::vector< Point > corners;
      getCorners( intersection.geometry(), corners );
      const dgf::DomainData *data = boundaryDomainBlock_->contains( corners );
      if( data )
        return data->parameter();
      else
        return DGFBoundaryParameter::defaultValue();
    }

  private:
    void generate ( std::istream &input, MPICommunicatorType comm, const std::string &name = "" );

    template< class Geometry >
    static void getCorners ( const Geometry &geometry, std::vector< Point > &corners )
    {
      corners.resize( geometry.corners() );
      for( int i = 0; i < geometry.corners(); ++i )
      {
        const typename Geometry::GlobalCoordinate corner = geometry.corner( i );
        for( int j = 0; j < dimension; ++j )
          corners[ i ][ j ] = corner[ j ];
      }
    }

    Grid *grid_;
    BoundaryDomainBlock *boundaryDomainBlock_;
  };



  // FactoryWrapper
  // --------------

  template <class Grid>
  struct FactoryWrapper
  {
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef GridFactory< Grid > Factory ;
    typedef unsigned int VertexId ;

    Factory factory_ ;
    VertexId counter_; 

    FactoryWrapper(MPICommunicatorType comm) 
      : factory_( ),
        counter_( 0 )
    {} 

    template <class VertexType>
    VertexId insertVertex( const VertexType& vx, const size_t globalId ) 
    {
      factory_.insertVertex( vx );
      return counter_++;
    }

    void insertElement( const GeometryType type, const std::vector< VertexId >& indices) 
    {
      factory_.insertElement( type, indices );
    }

    void insertBoundary( const int element, const int face, const int id ) 
    {
      factory_.insertBoundary( element, face, id );
    }

    void insertProcessBorder( const int element, const int face ) 
    {
      assert( false );
      abort();
    }

    Grid* createGrid(const std::string& name) { return factory_.createGrid(); }
  };



  // FactoryWrapper for ALUCubeGrid
  // ------------------------------

#if ENABLE_ALUGRID
  template<>
  struct FactoryWrapper< ALUCubeGrid< 3, 3 > >
  {
    typedef ALUCubeGrid< 3, 3 > Grid;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef GridFactory< Grid > Factory;
    typedef Factory::VertexId VertexId;

    Factory factory_ ;

    FactoryWrapper(MPICommunicatorType comm) : factory_( comm ) 
    {} 

    template <class VertexType>
    VertexId insertVertex( const VertexType& vx, const size_t globalId ) 
    {
      return factory_.insertVertex( vx, globalId );
    }

    void insertElement( const GeometryType type, const std::vector< VertexId >& indices) 
    {
      factory_.insertElement( type, indices );
    }

    void insertBoundary( const int element, const int face, const int id ) 
    {
      factory_.insertBoundary( element, face, id );
    }

    void insertProcessBorder( const int element, const int face ) 
    {
      factory_.insertProcessBorder( element, face );
    }

    Grid* createGrid(const std::string& name) { return factory_.createGrid( true, true, name ); }
  };
#endif // #if ENABLE_ALUGRID



  // FactoryWrapper for UGGrid
  // -------------------------

  template< int dim >
  class UGGrid;

  template< int dim >
  struct FactoryWrapper< UGGrid< dim > >
  {
    typedef UGGrid< dim > Grid;
    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef GridFactory< Grid > Factory;
    typedef unsigned int VertexId;

    FactoryWrapper ( MPICommunicatorType comm )
    : factory_(),
      counter_( 0 )
    {}

    template< class VertexType >
    VertexId insertVertex( const VertexType &vx, const size_t globalId )
    {
      factory_.insertVertex( vx );
      return counter_++;
    }

    void insertElement( const GeometryType type, const std::vector< VertexId > &indices )
    {
      factory_.insertElement( type, indices );
    }

    void insertBoundary( const int element, const int face, const int id ) 
    {}

    void insertProcessBorder( const int element, const int face ) 
    {
      assert( false );
      abort();
    }

    Grid *createGrid ( const std::string &name )
    {
      Grid *grid = factory_.createGrid();
      grid->setClosureType( Grid::NONE );
      std::cerr << "created UGGrid from '" << name << "'." << std::endl;
      return grid;
    }

  private:
    Factory factory_;
    VertexId counter_; 
  };



  // Implementation of DGFGridFactory for CartesianGrid
  // --------------------------------------------------

  template< class HostGrid >
  inline void DGFGridFactory< CartesianGrid< HostGrid > >
    ::generate ( std::istream &input, MPICommunicatorType comm, const std::string &name )
  {
    static const int dim = Grid::dimension;
    typedef Dune::SPGrid< typename Grid::ctype, dim > SPGrid;
    typedef typename SPGrid::LevelGridView GridView;
    typedef typename SPGrid::GlobalIdSet IdSet;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim< dim >::template Partition< InteriorBorder_Partition >::Iterator VertexIterator;
    typedef typename GridView::template Codim< 0 >::template Partition< InteriorBorder_Partition >::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef FactoryWrapper< HostGrid > Factory;
    typedef typename Factory::VertexId VertexId;

    GridPtr< SPGrid > spGrid( input, comm );
    Factory factory( comm );

    const GridView spView = spGrid->levelView( 0 );
    const IndexSet &indexSet = spView.indexSet();
    const IdSet &idSet = spGrid->globalIdSet();

    std::vector< VertexId > vertexId( indexSet.size( dim ) );

    const VertexIterator vend = spView.template end< dim, InteriorBorder_Partition >();
    for( VertexIterator vit = spView.template begin< dim, InteriorBorder_Partition >(); vit != vend; ++vit )
    {
      vertexId[ indexSet.index( *vit ) ]
        = factory.insertVertex( (*vit).geometry().center(), idSet.id( *vit ) );
    }

    int elIndex = 0;
    const int numVertices = (1 << dim);
    std::vector< VertexId > vertices( numVertices );
    const ElementIterator end = spView.template end< 0, InteriorBorder_Partition >();
    for( ElementIterator it = spView.template begin< 0, InteriorBorder_Partition >(); it != end; ++it )
    {
      const typename ElementIterator::Entity &entity = *it;
      assert( numVertices == entity.template count< dim >() );
      for( int i = 0; i < numVertices; ++i )
        vertices[ i ] = vertexId[ indexSet.subIndex( entity, i, dim ) ];
      factory.insertElement( entity.type(), vertices );

      const IntersectionIterator iend = spView.iend( entity );
      for( IntersectionIterator iit = spView.ibegin( entity ); iit != iend; ++iit )
      {
        const int face = iit->indexInInside();
        if( iit->boundary() )
          factory.insertBoundary( elIndex, face, iit->boundaryId() );
        if( (*(entity.template subEntity< 1 >( face ))).partitionType() == BorderEntity )
          factory.insertProcessBorder( elIndex, face );
      }

      ++elIndex;
    }

    grid_ = new Grid( factory.createGrid( name ) );
    boundaryDomainBlock_ = new BoundaryDomainBlock( input, dimension );
  }



  // DGFGridInfo for CartesianGrid
  // -----------------------------

  template< class HostGrid >
  struct DGFGridInfo< CartesianGrid< HostGrid > >
  {
    static int refineStepsForHalf ()
    {
      return DGFGridInfo< HostGrid >::refineStepsForHalf();
    }

    static double refineWeight ()
    {
      return DGFGridInfo< HostGrid >::refineWeight();
    }
  };
  
} // namespace Dune

#endif // #ifndef DUNE_CARTESIANGRID_DGFPARSER_HH
