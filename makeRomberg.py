    # Execution of makeromberg.py
    # python makeromberg.py (name_of_the_molecule)
    # Works under the assumption that a Field=X+000 calculation has been done under with the name sp(name_of_the_molecule)0.com
import sys
import os
molname = sys.argv[1]
coordinates=("X","Y","Z")
for j in coordinates:
    for i in range(1,7):
        ffield = 2**i
        field = round(float(ffield)/(10**4),4)
        val="%i" % ffield
        val_str="%.4f" % field
        os.system('cp sp%s0.com sp%s%s+%s.com' % (molname,molname,j,val_str))
        os.system('sed -i "s/Field=X+000/Field=%s+%s/g" sp%s%s+%s.com' % (j,val,molname,j,val_str))
        os.system('sed -i "s/chk=sp%s0/chk=sp%s%s+%s/g" sp%s%s+%s.com' % (molname,molname,j,val_str,molname,j,val_str))
        os.system('cp sp%s0.com sp%s%s-%s.com' % (molname,molname,j,val_str))
        os.system('sed -i "s/Field=X+000/Field=%s-%s/g" sp%s%s-%s.com' % (j,val,molname,j,val_str))
        os.system('sed -i "s/chk=sp%s0/chk=sp%s%s-%s/g" sp%s%s-%s.com' % (molname,molname,j,val_str,molname,j,val_str))
