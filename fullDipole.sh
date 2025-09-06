#!/bin/bash

molname=$1

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
	/home/gpey/ROMBERG/ROMBERG.exe -i M -o $i -F 13 -T $molname > ${molname}_Dipole${Prop}.out
done
