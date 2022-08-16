# ContGridMod

*ContGridMod* is a package proving tools to construct and simulate continuous models of large power systems.

**Contributors**: Laurent Pagnier and Julian Fritzsch

Please cite our work as 
L. Pagnier, J. Fritzsch, P. Jacquod and M. Chertkov, "Toward Model Reduction for Power System Transients With Physics-Informed PDE,"
in IEEE Access, vol. 10, pp. 65118-65125, 2022, doi: 10.1109/ACCESS.2022.3183336.

# How to install

Being in the ContGridMod folder:
```julia
using Pkg
Pkg.add(path=".")
```

# Adding a dependency

Being in the ContGridMod folder:
```julia
using Pkg
Pkg.activate(".")
Pkg.add(<some_package>)
```

# Updating the package

If we modify the package, it seems that we have to push (or at least commit) first


