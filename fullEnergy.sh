#!/bin/bash

molname=$1
if [[ -e "fullEnergy.out" ]]; then
    rm fullEnergy.out
fi

for i in M A B G; do
	if [ "$i" = "M" ]; then
		Prop="Dipole"
	fi	
	if [ "$i" = "A" ]; then
		Prop="Alpha"
	fi	
	if [ "$i" = "B" ]; then
		Prop="Beta"
	fi	
	if [ "$i" = "G" ]; then
		Prop="Gamma"
	fi	
	./ROMBERG.exe -i E -o $i -F 13 -T $molname > ${molname}_Energy${Prop}.out
    cat ${molname}_Energy${Prop}.out >> fullEnergy.out
    echo >> fullEnergy.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullEnergy.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullEnergy.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullEnergy.out
    echo >> fullEnergy.out
done
