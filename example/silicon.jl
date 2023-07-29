using ASEconvert
using AtomsBase
using DFTK
using ElectronicStructure


system = pyconvert(AbstractSystem, ase.build.bulk("Si"))
system = attach_psp(system; Si="hgh/lda/Si-q4")

params = DftkParameters(; system, functionals=[:lda_x, :lda_c_pw],
                        basis_kwargs=(; Ecut=15, kgrid=(4, 4, 4)),
                        scf_kwargs=(; tol=1e-8))

state = calculate(DftkCalculator(), params)
@show energy(state)
