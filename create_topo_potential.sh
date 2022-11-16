#!/bin/bash

if [ "$#" -ne 3 ]; then
	echo "ERROR! Correct usage: $0 topo_potential_file_name grid_step grid_max"
	exit 1
fi

topo_filename=$1
grid_step=$2
grid_max=$3

rm -f ${topo_filename}

# create topo potential file
for x in $(seq -${grid_max} ${grid_step} ${grid_max}); do
	awk -v x=${x} 'BEGIN{ printf "%.5lf %.18lf\n", x, -2.0*i  }' >> ${topo_filename}
done
