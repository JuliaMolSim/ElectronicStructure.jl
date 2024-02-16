using ASEconvert
using PythonCall
using ABINIT_jll
using Unitful
using UnitfulAtomic
using LazyArtifacts
using NCDatasets

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

function AbinitState(params::AbinitParameters; extra_parameters...)
    ase_atoms = convert_ase(params.system)
    ase_atoms.calc = pyimport("ase.calculators.abinit").Abinit(;
        label="abinit",
        params.ecut,
        #params.kpts,
        params.tolwfr,
        params.xc,
        #ixc = "-001012", # explicit [:lda_x, :lda_c_pw]
        # Let's not use ABINIT_PP_PATH.
        pp_paths=["$PROJECT_ROOT/data/psp"],
        # Cannot work because of duplicate…
        # pseudos="tt",
        params.pps,
        nsym=1,
        # Otherwise, cannot have non-symmetric kgrid because of ASE that uses kptopt=1…
        kptopt=1,
        nshiftk=1,
        ngkpt="1 2 1",
        shiftk="0.0 0.0 0.0",
        v8_legacy_format=false,
        extra_parameters...,
    )
    AbinitState(params, ase_atoms)
end

function calculate(calc::AbinitCalculator, params::AbinitParameters)
    calculate(calc, AbinitState(params))
end

function calculate(::AbinitCalculator, state::AbinitState)
    mktempdir() do tmpdir
        cdir = pwd()
        data = (; )
        try
            cd(tmpdir)
            potential_energy = ABINIT_jll.abinit() do abinit
                ase_workaround = ";ln -sf abinito_DS1_EIG abinito_EIG"
                state.ase_atoms.calc.command = "$abinit $(state.ase_atoms.calc.label).in >/dev/null" *
                                               ase_workaround
                austrip(pyconvert(Float64, state.ase_atoms.get_potential_energy()) * u"eV")
            end
            dynmats = NCDataset("$(state.ase_atoms.calc.label)o_DS2_DDB.nc", "r") do ds
                copy(ds["second_derivative_of_energy"])
            end
            data = merge(data, (; potential_energy, dynmats))
        finally
            cd(cdir)
            data
        end
    end
end

function energy(state::AbinitState)
    # pd_pbe_family = artifact"pd_nc_sr_pbe_standard_0.4.1_upf"
    potential_energy = mktempdir() do tmpdir
        cdir = pwd()
        potential_energy = nothing
        try
            cd(tmpdir)
            potential_energy = ABINIT_jll.abinit() do abinit
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
            @debug open("$tmpdir/abinit.in") do f
                while !eof(f)
                    @info readline(f)
                end
            end
            cd(cdir)
            potential_energy
        end
    end
    austrip(pyconvert(Float64, potential_energy) * u"eV")
end
