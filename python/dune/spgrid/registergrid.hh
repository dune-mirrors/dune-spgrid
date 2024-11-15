#pragma once

#include <string>
#include <sstream>
#include <dune/common/tupleutility.hh>

#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>

template <int dim>
void registerSPGrid(pybind11::module module)
#ifdef INCLUDE_REGALUGRID_INLINE
{
  // add commonly used P4estGrid variants
  using pybind11::operator""_a;
  std::string modname = std::string("_spgrid_" + std::to_string(dim) + "_double");
  std::string descr("Precompiled ");
  descr += modname;
  pybind11::module cls0 = module.def_submodule( modname.c_str(), descr.c_str());
  {
    using DuneType = Dune::SPGrid< double, dim, Dune::SPIsotropicRefinement >;
    std::string gridTypeName;
    {
      std::stringstream gridStr;
      gridStr << "Dune::SPGrid< double, " << dim << ", Dune::SPIsotropicRefinement >";
      gridTypeName = gridStr.str();
    }

    auto cls = Dune::Python::insertClass< DuneType, std::shared_ptr<DuneType> >( cls0, "HierarchicalGrid",pybind11::dynamic_attr(),
        Dune::Python::GenerateTypeName(gridTypeName),
        Dune::Python::IncludeFiles{"dune/grid/spgrid.hh","dune/grid/spgrid/dgfparser.hh","dune/python/grid/hierarchical.hh"}).first;
    Dune::Python::registerHierarchicalGrid( cls0, cls );
  }
}
#else
;
#endif
