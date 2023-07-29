import DFTK
using AtomsBase

abstract type AbstractCalculator end
abstract type AbstractState end
abstract type AbstractParameters end

struct DftkCalculator <: AbstractCalculator end

Base.@kwdef struct DftkParameters <: AbstractParameters
    system::AbstractSystem
    functionals::Vector{Symbol}    = [:gga_x_pbe, :gga_c_pbe]  # default to model_PBE
    model_kwargs = (; )
    basis_kwargs = (; )
    scf_kwargs   = (; )
end

struct DftkState <: AbstractState
    params::DftkParameters
    model::DFTK.Model
    basis::DFTK.PlaneWaveBasis
    scfres
end
function DftkState(params::DftkParameters)
    model = DFTK.model_DFT(params.system, params.functionals; params.model_kwargs...)
    basis = DFTK.PlaneWaveBasis(model; params.basis_kwargs...)
    ψ = nothing
    ρ = DFTK.guess_density(basis, params.system)
    DftkState(params, model, basis, (; basis, ψ, ρ))
end

function calculate(calc::DftkCalculator, params::DftkParameters)
    calculate(calc, DftkState(params))
end

function calculate(calc::DftkCalculator, state::DftkState)
    scfres = state.scfres
    scfres = DFTK.self_consistent_field(state.basis;
                                        scfres.ψ, scfres.ρ, state.params.scf_kwargs...)
    DftkState(state.params, state.model, state.basis, scfres)
end

energy(state::DftkState) = state.scfres.energies.total
