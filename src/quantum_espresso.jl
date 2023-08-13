using ASEconvert
using PythonCall
using QuantumEspresso_jll
using MPI
using Unitful
using UnitfulAtomic


Base.@kwdef struct QeCalculator <: AbstractCalculator
    # TODO QE setup parameters
end

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
