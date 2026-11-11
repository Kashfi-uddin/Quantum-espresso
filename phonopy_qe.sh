#create header.in
phonopy-init --qe -d --dim="2 1 1" -c pwscf.in 
chmod +x run.sh
./run.sh
phonopy-init --qe -f pwscf-{001..005}.out
phonopy --config band.conf -p
phonopy-bandplot --gnuplot band.yaml > phonon_dispersion.dat

