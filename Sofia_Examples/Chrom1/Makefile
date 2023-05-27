log.lammps	: in.lmp Body Header interactions softinteractions variables
	lmprun in.lmp
Body Header	: ConfigChrom.dat
	DemoInit_2.0; cat Header Body > chrom.init
realclean	:	clean
	rm -f Bconfig.dat Body Header Atom*
clean		:
	rm -f atom*.dat *.jpg *.mp4 DropletSize* */DropletSize* restart_state.* log.lammps *~ slurm* beadtrack* *.pdf
