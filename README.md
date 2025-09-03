+  Several prerequisites to execute RomberG: 
    -  Have performed all the calculations in Gaussian, the NLOP are read from their respective .fchk file.
    -  Have all the fchk files into a folder under the molecule's name - all data extraction will be performed from there.
    -  IMPORTANT! Downloading the "f90getopt.F90" code is mandatory as it is a flag parser that is used in the code and a key module for the script.
    -  IMPORTANT! All fchk files must have the following name structure: sp(NAME)(x/X/y/Y/z/Z)(FieldStrength·10^-4).com, except for the null-field one - sp(NAME)0.com. In this repository it is provided a simple Python script that generates the files, download the 'Documentation.pdf' for full details.

+ Execution: ./ROMBERG.exe -i (input_NLOP) -o (output_NLOP) -F (number_of_total_fields_scoped) {-P} {--isotropic} { (name_of_the_molecule)
