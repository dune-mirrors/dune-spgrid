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

#-------------------------------------------------------------------
# grid module loading
def _checkModule(includes, typeName, typeTag):
    from importlib import import_module

    # check if pre-compiled module exists and if so load it
    try:
        gridModule = import_module("dune.spgrid._spgrid._spgrid_" + typeTag)
        return gridModule
    except ImportError:
        from dune.grid.grid_generator import module
        # otherwise proceed with generate, compile, load
        gridModule = module(includes, typeName)
        return gridModule



def spIsotropicGrid(domain, dimgrid=None, ctype="double", **kwargs):
    from ..grid.grid_generator import module, getDimgrid

    if dimgrid is None:
        dimgrid = getDimgrid(domain)

    typeName = "Dune::SPGrid< " + ctype + ", " + str(dimgrid) + ", Dune::SPIsotropicRefinement >"
    includes = ["dune/grid/spgrid.hh", "dune/grid/spgrid/dgfparser.hh"]
    typeTag = str(dimgrid) + "_" + ctype
    gridModule = _checkModule(includes, typeName, typeTag)

    return gridModule.LeafGrid(gridModule.reader(domain))

def spAnisotropicGrid(domain, dimgrid=None, ctype="double", **kwargs):
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

def spGrid(*args, **kwargs):
    return spIsotropicGrid(*args, **kwargs)

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
