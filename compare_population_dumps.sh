#!/bin/bash

echo "Insert the number of population dumps to compare"
read number_population_dumps
dump_dir=./

echo "Begin comparing files"
for (( i = 0; i <= $number_population_dumps; i++ )) 
do
    printf -v serial_dump_name "${dump_dir}ser_Population%04d.dmp" "$i"
    printf -v parallel_dump_name "${dump_dir}par_Population%04d.dmp" "$i"
    diff "$serial_dump_name" "$parallel_dump_name" > /dev/null 2>&1
    result=$?
    if [ $result -eq 1 ]
    then
        echo "First difference found at step ${i} in files ${serial_dump_name} and ${parallel_dump_name}"
        exit 1
    fi
done
echo "Population dumps are equal at every step."
    

