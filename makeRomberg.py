import sys
import os
molname = sys.argv[1]
for i in range(1,7):
    ffield = 2**i
    field = round(float(ffield)/(10**4),4)
    val="%i" % ffield
    val_str="%.4f" % field
    os.system('cp sp%s0.com sp%sX+%s.com' % (molname,molname,val_str))
    os.system('sed -i "s/Field=X+000/Field=X+%s/g" sp%sX+%s.com' % (val,molname,val_str))
    os.system('sed -i "s/chk=sp%s0/chk=sp%sX+%s/g" sp%sX+%s.com' % (molname,molname,val_str,molname,val_str))
    os.system('cp sp%s0.com sp%sX-%s.com' % (molname,molname,val_str))
    os.system('sed -i "s/Field=X+000/Field=X-%s/g" sp%sX-%s.com' % (val,molname,val_str))
    os.system('sed -i "s/chk=sp%s0/chk=sp%sX-%s/g" sp%sX-%s.com' % (molname,molname,val_str,molname,val_str))
