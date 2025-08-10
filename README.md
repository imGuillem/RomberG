Prerequisite: Have performed all the calculations in Gaussian. Secondly, have all the fchk files into a folder of the molecule name - all data extraction will be performed from there. Downloading the "f90getopt.F90" code is mandatory as it is a flag parser that is used in the code.

Compilation: gfortran f90getopt.F90 -ffree-line-length-none -o ROMBERG.exe ROMBERG.f95

Execution: ./ROMBERG.exe -i {input_NLOP} -o {output_NLOP} -F {number_of_total_fields_scoped} {-P} (name_of_the_molecule)
