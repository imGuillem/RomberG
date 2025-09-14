#!/bin/bash

molname=$1
if [[ -e "fullAlpha.out" ]]; then
    rm fullAlpha.out
fi

for i in B G; do
	if [ "$i" = "B" ]; then
		Prop="Beta"
	fi	
	if [ "$i" = "G" ]; then
		Prop="Gamma"
	fi	
	./ROMBERG.exe -i A -o $i -F 13 -T $molname > ${molname}_Alpha${Prop}.out
    cat ${molname}_Alpha${Prop}.out >> fullAlpha.out
    echo >> fullAlpha.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullAlpha.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullAlpha.out
    echo "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" >> fullAlpha.out
    echo >> fullAlpha.out
done
