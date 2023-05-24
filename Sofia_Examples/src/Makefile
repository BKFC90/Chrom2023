#
# Makefile for C++ library files
#
CPPTOP=		$(HOME)/g++
CXX=		g++ -std=c++11
INCLUDEPATH=	$(CPPTOP)/include
LIBPATH=	$(CPPTOP)/lib
LIBFLAGS=	 -g -fopenmp -L$(LIBPATH) -lconfig++
COMPILEFLAGS=	-g -c -fopenmp -I$(INCLUDEPATH)
MEMFLAGS=

clean:	
	rm -f *.o *.man *~ *#  *0

PolyInit%:	PolyInit%.o
	$(CXX) -O3 -o $@ $(MEMFLAGS) $< $(EXTRA_LIB_FLAGS) $(LIBFLAGS) -lm

.cc.o:
	$(CXX) $(MEMFLAGS) $(COMPILEFLAGS) $*.cc
#
#       Reevaluate dependencies
#
depend:
	@ echo making dependencies...
	for i in *.cc; do \
	$(CXX) -MM $(COMPILEFLAGS) $$i >> makedep; done
	@ echo '/^# DO NOT DELETE THIS LINE/+2,$$d' >eddep
	@ echo '$$r makedep' >>eddep
	@ echo 'w' >>eddep
	@ ed - Makefile < eddep
	@ rm eddep makedep 
	@ echo '# DEPENDENCIES MUST END AT END OF FILE' >> Makefile
	@ echo '# IF YOU PUT STUFF HERE IT WILL GO AWAY' >> Makefile
	@ echo '# see make depend above' >> Makefile

#-----------------------------------------------------------------
# DO NOT DELETE THIS LINE -- make depend uses it

PolyInit_1.0.o: PolyInit_1.0.cc
# DEPENDENCIES MUST END AT END OF FILE
# IF YOU PUT STUFF HERE IT WILL GO AWAY
# see make depend above
