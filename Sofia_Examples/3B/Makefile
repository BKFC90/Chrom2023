.SILENT		: distclean
.SILENT		: clean
log.lammps	: in.lmp
	lmprun in.lmp
clean		:
	rm -f atom*.dat *.jpg *.mp4 restart_state.* log.lammps *~
distclean	:
	for dir in 2B 3B LJ1 LJ2 LJ3 Poly1 Poly2 Poly3;	\
	do						\
	cd $$dir;					\
	make clean;					\
	cd ..;						\
	done;
