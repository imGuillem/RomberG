#!/bin/bash

molname=$1
if [[ -f "fullDipole.out" ]]; then
    rm fullDipole.out
fi

for i in A B G; do
	if [ "$i" = "A" ]; then
		Prop="Alpha"
	fi	
	if [ "$i" = "B" ]; then
		Prop="Beta"
	fi	
	if [ "$i" = "G" ]; then
		Prop="Gamma"
	fi	
	./ROMBERG.exe -i M -o $i -F 13 -T $molname > ${molname}_Dipole${Prop}.out
    cat ${molname}_Dipole${Prop}.out >> fullDipole.out
    echo >> fullDipole.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullDipole.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullDipole.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullDipole.out
    echo >> fullDipole.out
done

