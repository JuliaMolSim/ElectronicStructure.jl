using ASEconvert
using PythonCall
using ABINIT_jll
using Unitful
using UnitfulAtomic
using LazyArtifacts

Base.@kwdef struct AbinitCalculator <: AbstractCalculator
    # TODO
end

Base.@kwdef struct AbinitParameters <: AbstractParameters
    system::AbstractSystem
    ecut        = ustrip(u"eV", 10u"Eh_au")
    kpts        = [1, 1, 1]
    tolwfr      = 1e-12
    xc          = "PBE"
    pps         = "hgh"
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
        #ixc = "-001012", # explicit [:lda_x, :lda_c_pw]
        # Let's not use ABINIT_PP_PATH.
        pp_paths=["$PROJECT_ROOT/data/psp"],
        # Cannot work because of duplicate…
        # pseudos="tt",
        params.pps,
        nsym=1,
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
    pd_pbe_family = artifact"pd_nc_sr_pbe_standard_0.4.1_upf"
    energy = mktempdir() do tmpdir
        cdir = pwd()
        try
            cd(tmpdir)
            energy = ABINIT_jll.abinit() do abinit
                # `withenv` does not work…
                state.ase_atoms.calc.command = "$abinit abinit.in > abinit.out"
                #= Hack if other pseudo…
                delete_first_pseudo_line = `sed -i '0,/^pseudos/{//d}' abinit.in`
                add_pseudo_line = `sed -i "1 i\\pseudos \"$(pd_pbe_family)/Si.upf\"" abinit.in`
                state.ase_atoms.calc.command = "$(delete_first_pseudo_line); $(add_pseudo_line); " *
                                               pyconvert(String, state.ase_atoms.calc.command)
                =#
                state.ase_atoms.get_potential_energy()
            end
        finally
            @info open("$tmpdir/abinit.in") do f
                while !eof(f)
                    @info readline(f)
                end
            end
            cd(cdir)
            energy
        end
    end
    austrip(pyconvert(Float64, energy) * u"eV")
end
