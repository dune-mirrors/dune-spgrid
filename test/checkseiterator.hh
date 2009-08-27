#ifndef DUNE_SPGRID_CHECKSEITERATOR_HH
#define DUNE_SPGRID_CHECKSEITERATOR_HH

#include <dune/grid/common/gridview.hh>

#include <dune/grid/extensions/superentityiterator.hh>

namespace Dune
{

  template< int codim, class VT >
  void checkSuperEntityIterator ( const Dune::GridView< VT > &gridView )
  {
    typedef Dune::GridView< VT > GridView;
    typedef typename GridView::IndexSet IndexSet;

    const IndexSet &indexSet = gridView.indexSet();
    std::vector< int > count( indexSet.size( codim ), 0 );

    typedef typename GridView::template Codim< 0 >::Iterator ElementIterator;

    const ElementIterator elEnd = gridView.template end< 0 >();
    for( ElementIterator elIt = gridView.template begin< 0 >(); elIt != elEnd; ++elIt )
    {
      const typename ElementIterator::Entity &entity =*elIt;
      for( int i = 0; i < entity.template count< codim >(); ++i )
        ++count[ indexSet.subIndex( entity, i, codim ) ];
    }

    typedef SuperEntityIteratorExtension< GridView > ExtGridView;
    ExtGridView extGridView( gridView );

    typedef typename GridView::template Codim< codim >::Iterator CodimIterator;
    typedef typename ExtGridView::template Codim< codim >::SuperEntityIterator SEIterator;

    const CodimIterator codimEnd = gridView.template end< codim >();
    for( CodimIterator codimIt = gridView.template begin< codim >(); codimIt != codimEnd; ++codimIt )
    {
      const typename CodimIterator::Entity &entity = *codimIt;

      int cnt = 0;
      const SEIterator seEnd = extGridView.superEntityEnd( entity );
      for( SEIterator seIt = extGridView.superEntityBegin( entity ); seIt != seEnd; ++seIt )
      {
        const typename SEIterator::Entity &element = *seIt;

        ++cnt;

        bool found = false;
        for( int i = 0; i < element.template count< codim >(); ++i )
          found |= (element.template subEntity< codim >( i ) == codimIt);
        if( !found )
        {
          std::cout << "Entity " << indexSet.index( entity )
                    << "of codimension " << codim
                    << "is not a subentity of its super entity iterator."
                    << std::endl;
        }
      }

      if( cnt != count[ indexSet.index( entity ) ] )
      {
        std::cout << "Number of super entities differs for entity "
                  << indexSet.index( entity )
                  << " of codimension " << codim
                  << " (" << cnt << " != " << count[ indexSet.index( entity ) ] << ")."
                  << std::endl;
      }
    }
  }

}

#endif // #ifndef DUNE_SPGRID_CHECKSEITERATOR_HH
