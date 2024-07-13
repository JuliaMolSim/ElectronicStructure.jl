using ASEconvert
using AtomsCalculators
using PythonCall
using QuantumEspresso_jll
using MPI
using Unitful
using UnitfulAtomic


### !!!
### TODO: As of now, this relies on ASE, but we want to have this to only depend on QuantumEspresso_jll
### !!!


Base.@kwdef struct QeParameters <: AbstractParameters
    # TODO Keywords currently based on ASE
    system::AbstractSystem
    ecutwfc     = 40
    conv_thr    = 1e-11
    tstress     = true
    tprnfor     = true
    smearing    = "gaussian"
    mixing_mode = "plain"
    mixing_beta = 0.7
    mixing_ndim = 10
    kpts        = (1, 1, 1)
    occupations = "smearing"
    degauss     = 0.01
    input_dft   = "pbe"
    electron_maxstep  = 100
    pseudopotentials  = Dict{String,String}()
    extra_parameter   = Dict{Symbol,Any}()
    working_directory = mktempdir(pwd())
    n_mpi_procs       = MPI.Comm_size(MPI.COMM_WORLD)
    n_threads         = BLAS.get_num_threads()
end



function convert(::Type{QeParameters}, params::DftkParameters)
    # Convert DFTK parameters to QE parameters
    #
    # keep in mind the unit conversion
    error("TODO")
end


struct QeState <: AbstractState
    params::QeParameters
    ase_atoms::Py
end

function QeState(params::QeParameters)
    ase_atoms = convert_ase(params.system)
    ase_atoms.calc = pyimport("ase.calculators.espresso").Espresso(;
        label="espresso",
        params.input_dft,
        params.pseudopotentials,
        params.kpts,
        params.ecutwfc,
        params.tstress,
        params.tprnfor,
        params.mixing_mode,
        params.mixing_beta,
        params.conv_thr,
        params.occupations,
        params.smearing,
        params.degauss,
        params.electron_maxstep,
        params.mixing_ndim,
        params.extra_parameter...
    )
    QeState(params, ase_atoms)
end


struct QeCalculator <: AbstractCalculator
    state::QeState

    QeCalculator(params::QeParameters) = new(QeState(params))
end


function calculate(calc::QeCalculator, params::QeParameters)
    calculate(calc, QeState(params))
end

function calculate(::QeCalculator, state::QeState)
    n_mpi_procs = state.params.n_mpi_procs
    MPI.mpiexec() do mpirun
        QuantumEspresso_jll.pwscf() do pwscf
            qe_command = "$mpirun -np $n_mpi_procs $pwscf -in PREFIX.pwi > PREFIX.pwo"
            state.ase_atoms.calc.command = qe_command
            withenv("OMP_NUM_THREADS" => state.params.n_threads) do
                state.ase_atoms.calc.get_potential_energy()
            end
        end
    end
end

function energy(state::QeState)
    austrip(state.ase_calculator.get_potential_energy() * u"eV")
end

# # AtomsCalculators interface

# This might need some work reorganizing how the system and the state are wrapped
# Ignoring `sys` for now`for all methods as we wrap it

AtomsCalculators.energy_unit(calc::QeCalculator) = u"eV"
AtomsCalculators.length_unit(calc::QeCalculator) = u"Å"

AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(sys, calc::QeCalculator; kwargs...)
    unit = energy_unit(calc)
    # This might need some work reorganizing how the system and the state are wrapped
    # Ignoring `sys` for now`
    energy = calc.state.ase_calculator.get_potential_energy()
    return pyconvert(Float64, energy) * unit
end

AtomsCalculators.@generate_interface function AtomsCalculators.forces(sys, calc::QeCalculator; kwargs...)
    unit = energy_unit(calc) / energy_unit(calc)
    forces = calc.state.ase_calculator.get_forces()
    return pyconvert(Array, forces) * unit
end

AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, calc::QeCalculator; kwargs...)
    state = calc.state
    ase_system = state.params.system

    unit = energy_unit(calc)
    stress = state.ase_calculator.get_stress(system)
    cons = ase.constraints
    virial = cons.voigt_6_to_full_3x3_stress(stress) * ( -ase_system.get_volume() )
    return pyconvert(Array, stress) * unit
end