module ElectronicStructure

const PROJECT_ROOT = pkgdir(ElectronicStructure)

export AbstractCalculator, AbstractState, AbstractParameters

using AtomsBase
using LinearAlgebra

include("interface.jl")
export DftkCalculator, DftkParameters, DftkState, calculate, energy
include("dftk.jl")

export AbinitParameters, AbinitState
include("abinit.jl")
export QeParameters, QeState
include("quantum_espresso.jl")

end
