using ElectronicStructure
using Test

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

@testset "Abinit" begin
    using ASEconvert
    using DFTK
    using LinearAlgebra

    system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
    system = attach_psp(system; Si="hgh/lda/si-q4")

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
    # Not equivalent for the moment.
    @test norm(energy_dftk - energy_abinit) < 1e-6
end

