using ASEconvert
using PythonCall
using ABINIT_jll
using Unitful
using UnitfulAtomic

Base.@kwdef struct AbinitCalculator <: AbstractCalculator
    # TODO
end

Base.@kwdef struct AbinitParameters <: AbstractParameters
    system::AbstractSystem
    ecut        = ustrip(u"eV", 40u"Eh_au")
    kpts        = [1, 1, 1]
    tolwfr      = 1e-22
    xc          = "LDA"
    n_threads   = BLAS.get_num_threads()
end



function convert(::Type{AbinitParameters}, params::DftkParameters)
    error("TODO")
end


struct AbinitState <: AbstractState
    params::AbinitParameters
    ase_atoms::Py
end

function AbinitState(params::AbinitParameters)
    ase_atoms = convert_ase(params.system)
    ase_atoms.calc = pyimport("ase.calculators.abinit").Abinit(;
        label="abinit",
        params.ecut,
        params.kpts,
        params.tolwfr,
        params.xc,
        nsym=1,
        ixc="-001007",
        v8_legacy_format=false,
    )
    AbinitState(params, ase_atoms)
end

function calculate(calc::AbinitCalculator, params::AbinitParameters)
    calculate(calc, AbinitState(params))
end

function calculate(::AbinitCalculator, state::AbinitState)
    ABINIT_jll.abinit() do abinit
        abinit_command = "$abinit $(state.ase_atoms.calc.label).in"
        state.ase_atoms.calc.command = abinit_command
        withenv("OMP_NUM_THREADS" => state.params.n_threads) do
            state.ase_atoms.get_potential_energy()
        end
    end
end

function energy(state::AbinitState)
    energy = mktempdir() do tmpdir
        cdir = pwd()
        try
            cd(tmpdir)
            energy = ABINIT_jll.abinit() do abinit
                # `withenv` does not work…
                state.ase_atoms.calc.command = "$abinit abinit.in > abinit.out"
                @show state.ase_atoms.get_forces()
                state.ase_atoms.get_potential_energy()
            end
            @info open("$tmpdir/abinit.in") do f
                while !eof(f)
                    @info readline(f)
                end
            end
            energy
        finally
            cd(cdir)
            energy
        end
    end
    austrip(pyconvert(Float64, energy) * u"eV")
end
