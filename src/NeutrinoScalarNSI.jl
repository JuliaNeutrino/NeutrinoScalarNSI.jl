"""
Scalar non-standard neutrino interactions for GLoBES.jl: a correction to
the neutrino mass matrix from a light scalar, as a propagation block that
combines with the other blocks in `NewPhysics`.
"""
module NeutrinoScalarNSI

using GLoBES
using GLoBES:
    PropagationBlock,
    flavour_pairs,
    pair_label,
    coupling_parameters,
    coupling_aliases,
    add_couplings!
import GLoBES: block_parameters, parameter_aliases, uses_absolute_masses, mass_correction!

export ScalarNSI

"""
    ScalarNSI(couplings...; scale=:Δm31)

Scalar non-standard interactions: a correction δM to the neutrino mass
matrix, M → M + δM with δM_αβ = s |η_αβ| e^{iδ_αβ} (Hermitian), so that the
vacuum Hamiltonian is (M + δM)(M + δM)†/2E. The scale s is √|Δm²31| (the
convention of Ge and Parke, arXiv:1812.08376, and Medhi et al.) or a fixed
scale in eV.

Couplings are flavour pairs (`:ee`, `:eμ`, `:eτ`, `:μμ`, `:μτ`, `:ττ`, or
ASCII `:emu`, ...; all six by default). Parameters: `η_αβ` and, off the
diagonal, the phase `δ_αβ` (ASCII `eta_etau`, `delta_etau`, ...).

Scalar NSI depend on the absolute neutrino masses: the parameter
`m_lightest` (eV, default 0) is the lightest mass.

```julia
ScalarNSI(:eτ)   # one coupling: parameters η_eτ and δ_eτ
ScalarNSI(:ee, :eμ; scale = 0.05)   # couplings; δM in units of 0.05 eV instead of √|Δm²31|
ScalarNSI()   # all six couplings
```
"""
struct ScalarNSI <: PropagationBlock
    pairs::Vector{Tuple{Int, Int}}
    scale::Union{Symbol, Float64}
    function ScalarNSI(couplings...; scale = :Δm31)
        scale === :Δm31 ||
            scale isa Real ||
            throw(ArgumentError("scale must be :Δm31 or a mass in eV, not $scale"))
        return new(flavour_pairs(couplings), scale isa Symbol ? scale : Float64(scale))
    end
end

block_parameters(b::ScalarNSI) = coupling_parameters(b.pairs, :η, :δ)
parameter_aliases(b::ScalarNSI) = coupling_aliases(b.pairs, :η, :δ, :eta, :delta)
uses_absolute_masses(::ScalarNSI) = true
mass_correction!(δM, b::ScalarNSI, p, x) =
    add_couplings!(δM, b.pairs, x, :η, :δ, b.scale === :Δm31 ? sqrt(abs(p.Δm31)) : b.scale)

Base.show(io::IO, b::ScalarNSI) = print(
    io,
    "ScalarNSI(",
    join((pair_label(p...) for p in b.pairs), ", "),
    b.scale === :Δm31 ? "" : "; scale = $(b.scale)",
    ")",
)

end
