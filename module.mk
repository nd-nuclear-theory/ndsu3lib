$(eval $(begin-module))

################################################################
# unit definitions
################################################################

# module_units_h := 
module_units_cpp-h := c++wrappers
module_units_f := ndsu3lib_recoupling ndsu3lib_wigner_canonical ndsu3lib_wigner_su3so3
module_programs_cpp := ndsu3lib_example_cpp
module_programs_f := ndsu3lib_example 
# module_programs_f_test := # DOES NOT EXIST YET
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
