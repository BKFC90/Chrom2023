log.lammps	: in.lmp poly.init interactions
	export OMP_NUM_THREADS=;
	mpiexec --oversubscribe -np 16 lmp -v Px 4 -v Py 2 -v Pz 2 -in in.lmp;
poly.init	: ConfigPoly.dat
	PolyInit
clean	:
	rm -f atom.dat *.mp4 *.jpg *~ log*
realclean		: clean
	rm poly.init
