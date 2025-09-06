#!/bin/bash

molname=$1

for i in B G; do
	if [ "$i" = "B" ]; then
		Prop="Beta"
	fi	
	if [ "$i" = "G" ]; then
		Prop="Gamma"
	fi	
	/home/gpey/ROMBERG/ROMBERG.exe -i A -o $i -F 13 -T $molname > ${molname}_Alpha${Prop}.out
done
