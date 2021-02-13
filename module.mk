$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
# module_units_cpp-h := 
module_units_f := binomial_coeff clebsch_gordan I_S_module orthonormalization_matrix outer_multiplicity \
    racah su2racah transformation_coeff wigner_canonical_extremal wigner_canonical wigner_physical wigner_physical_wrap Z_coeff
# module_programs_cpp :=
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
