using ASEconvert
using AtomsBase
using AtomsCalculators
using DFTK
using ElectronicStructure


system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
system = attach_psp(system; Si="hgh/lda/Si-q4")

params = DftkParameters(; system, functionals=[:lda_x, :lda_c_pw],
                        basis_kwargs=(; Ecut=15, kgrid=(4, 4, 4)),
                        scf_kwargs=(; tol=1e-8))
calc = DftkCalculator(params)

state = calculate(calc, calc.state.params)
@show energy(state)

@show potential_energy(system, calc)