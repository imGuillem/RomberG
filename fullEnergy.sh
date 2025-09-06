#!/bin/bash

molname=$1

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
	/home/gpey/ROMBERG/ROMBERG.exe -i E -o $i -F 13 -T $molname > ${molname}_Energy${Prop}.out
done
