# NeutrinoScalarNSI.jl

Scalar non-standard neutrino interactions for [GLoBES.jl](https://github.com/JuliaNeutrino/GLoBES.jl). A light scalar coupled to neutrinos and matter corrects the neutrino mass matrix, M → M + δM, rather than adding a potential like vector NSI. The effect therefore depends on the absolute neutrino masses, and it scales as 1/E like the vacuum term, whereas a matter potential matters more as the energy grows.

Physically δM is proportional to the density of the matter the neutrinos cross. The package takes η as a fixed parameter, the value for the matter along the baseline, so give each experiment the η of its density. With one fixed η, scalar NSI are exactly standard oscillations with shifted mixing angles, phase and mass splittings, at every energy and density: they show up as a disagreement between experiments in different environments, not within one experiment.

```julia
using Pkg
Pkg.develop(path = "/path/to/NeutrinoScalarNSI.jl")   # not registered yet; Julia 1.11 or later
using GLoBES, NeutrinoScalarNSI
```

`ScalarNSI` is a GLoBES.jl propagation block, with parameters `η_αβ`, the phases `δ_αβ` and the lightest mass `m_lightest`. δM = √|Δm²31| η by default, the convention of Ge and Parke (arXiv:1812.08376).

```julia
model = NewPhysics(ScalarNSI(:eτ))       # or ScalarNSI(:ee, :eμ; scale = 0.05) for δM in units of 0.05 eV
use_model!(model)                        # the probability engine of all experiments
p = OscParams(model, nufit_params("6.1"); η_eτ = 0.1, m_lightest = 0.01)
probability(p, :μ, :e; E = 2.5, L = 1300, ρ = 2.848)
```

It combines with the other blocks, e.g. `NewPhysics(ScalarNSI(:eτ), Decoherence())` with [NeutrinoDecoherence.jl](https://github.com/JuliaNeutrino/NeutrinoDecoherence.jl). Rates, χ² and fits in GLoBES.jl then use the new parameters like the standard ones.

The probabilities agree with the C engine of the SNSI-decoherence DUNE study to 10⁻¹⁰ (`test/snsi_c_reference.txt`). The notebook [`notebooks/01_scalar_NSI.ipynb`](notebooks/01_scalar_NSI.ipynb) goes through the package: the dependence on the absolute masses, the equivalence to shifted standard parameters, the C-engine cross-check and the DUNE sensitivity to η_ee ([open in Colab](https://colab.research.google.com/github/JuliaNeutrino/NeutrinoScalarNSI.jl/blob/main/notebooks/01_scalar_NSI.ipynb); choose **Runtime → Change runtime type → Julia** first).

Tests: `julia --project=. test/runtests.jl`.
