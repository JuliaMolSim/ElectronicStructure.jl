using ElectronicStructure
using Test
# Not sure what happens…
ENV["HDF5_DISABLE_VERSION_CHECK"] = "2"

@testset "ElectronicStructure.jl" begin
    # Write your tests here.
    @testset "Quantum Espresso" begin
        using ASEconvert
        using DFTK

        system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
        state = QeState(QeParameters(; system))
        @test state isa AbstractState
    end
end

@testset "Abinit: Phonon" begin
    using ASEconvert
    using DFTK
    using LinearAlgebra

    system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
    system = attach_psp(system; Si="hgh/lda/si-q4")

    energy_abinit = let
        phonon_parameters = (;
                             ndtset=2,
                             tolwfr1=1e-22,
                             tolwfr2=1e-6,
                             rfphon2=1,
                             rfatpol2="1 2",
                             rfdir2="1 1 1",
                             nqpt2="1",
                             qpt2="0.0 0.0 0.0",
                             getwfk2=-1,
                            )
        state_abinit = AbinitState(AbinitParameters(; system); phonon_parameters...)
        @show result = calculate(AbinitCalculator(), state_abinit)
    end
    #energy_dftk = let
    #    state_dftk = DftkState(DftkParameters(; system,
    #                                          # functionals=[:lda_x, :lda_c_pw],
    #                                          basis_kwargs = (; Ecut=10, kgrid=[1, 1, 1]),
    #                                          scf_kwargs = (; tol=1e-6),
    #                                         ))
    #    state_dftk = calculate(DftkCalculator(), state_dftk)
    #    energy(state_dftk)
    #end
    # Not equivalent for the moment.
end

@testset "Abinit" begin
    using ASEconvert
    using DFTK
    using LinearAlgebra
    using LazyArtifacts

    system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
    system = attach_psp(system; Si="hgh/lda/si-q4")
    #= Systematic error…
    pd_pbe_family = artifact"pd_nc_sr_pbe_standard_0.4.1_upf"
    system = attach_psp(system; Si=joinpath(pd_pbe_family, "Si.upf"))
    =#

    energy_abinit = let
        state_abinit = AbinitState(AbinitParameters(; system))
        energy(state_abinit)
    end
    energy_dftk = let
        state_dftk = DftkState(DftkParameters(; system,
                                              # functionals=[:lda_x, :lda_c_pw],
                                              basis_kwargs = (; Ecut=10, kgrid=[1, 1, 1]),
                                              scf_kwargs = (; tol=1e-6),
                                             ))
        state_dftk = calculate(DftkCalculator(), state_dftk)
        energy(state_dftk)
    end
    @test norm(energy_dftk - energy_abinit) < 1e-6
end
