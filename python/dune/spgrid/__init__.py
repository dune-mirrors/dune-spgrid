from __future__ import absolute_import, division, print_function, unicode_literals

from dune.common.checkconfiguration import assertHave, ConfigurationError

def spBisectionGrid(domain, dimgrid=None, ctype="double"):
    from ..grid.grid_generator import module, getDimgrid

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPBisectionRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(domain))


def spIsotropicGrid(domain, dimgrid=None, ctype="double"):
    from ..grid.grid_generator import module, getDimgrid

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPIsotropicRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    gridModule = module(includes, typeName)

    return gridModule.LeafGrid(gridModule.reader(domain))

def spAnisotropicGrid(domain, dimgrid=None, ctype="double"):
    from ..grid.grid_generator import module, getDimgrid
    from dune.generator import Method

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPAnisotropicRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    globalRefineMethod = Method('globalRefine', '''[]( DuneType &self, int level, std::vector<int> refDir )
                            {
                              static const int dim = ''' +str(dimgrid)+''';
                              std::bitset< dim > s(0);
                              const int size = (int(refDir.size()) >= dim) ? dim : refDir.size();
                              for( int i=0; i<size; ++i )
                              {
                                if( refDir[i] ) s[ i ] = 1;
                              }
                              self.globalRefine(level, Dune::SPAnisotropicRefinementPolicy< dim >( s ) );
                              return;
                            }''' )
    gridModule = module(includes, typeName, globalRefineMethod)

    return gridModule.LeafGrid(gridModule.reader(domain))

#############################################################
##  Registry
#############################################################

registry = {}
registry["grid"] = {
        "SP" : spIsotropicGrid,
        "SPBisection": spBisectionGrid,
        "SPIsotropic": spIsotropicGrid,
        "SPAnisotropic": spAnisotropicGrid
    }


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
