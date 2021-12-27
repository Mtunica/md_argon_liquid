

SRCDIR = src
BUILDDIR = build


fc = gfortran
OPT= -O3

program.exe: main0.o main1.o main2.o function_module.o statistic_module.o
	$(fc) -o program.exe $(OPT) $(BUILDDIR)/function_module.o $(BUILDDIR)/statistic_module.o $(BUILDDIR)/main0.o -J$(BUILDDIR)
	$(fc) -o program1.exe $(OPT) $(BUILDDIR)/function_module.o $(BUILDDIR)/statistic_module.o $(BUILDDIR)/main1.o -J$(BUILDDIR)
	$(fc) -o program2.exe $(OPT) $(BUILDDIR)/function_module.o $(BUILDDIR)/statistic_module.o $(BUILDDIR)/main2.o -J$(BUILDDIR)
	
function_module.mod: function_module.o $(SRCDIR)/function_module.f90
	$(fc) -o $(BUILDDIR)/function_module.o -c $(OPT) $(SRCDIR)/function_module.f90 -J$(BUILDDIR)

function_module.o: $(SRCDIR)/function_module.f90
	$(fc) -o $(BUILDDIR)/function_module.o -c $(OPT) $(SRCDIR)/function_module.f90 -J$(BUILDDIR)

statistic_module.mod: statistic_module.o $(SRCDIR)/statistic_module.f90
	$(fc) -o $(BUILDDIR)/statistic_module.o -c $(OPT) $(SRCDIR)/statistic_module.f90 -J$(BUILDDIR)

statistic_module.o: $(SRCDIR)/statistic_module.f90
	$(fc) -o $(BUILDDIR)/statistic_module.o -c $(OPT) $(SRCDIR)/statistic_module.f90 -J$(BUILDDIR)

main0.o: function_module.mod statistic_module.mod $(SRCDIR)/main0.f90
	gfortran -o $(BUILDDIR)/main0.o -c $(OPT) $(SRCDIR)/main0.f90 -J$(BUILDDIR)

main1.o: function_module.mod statistic_module.mod $(SRCDIR)/main1.f90
	gfortran -o $(BUILDDIR)/main1.o -c $(OPT) $(SRCDIR)/main1.f90 -J$(BUILDDIR)

main2.o: function_module.mod statistic_module.mod $(SRCDIR)/main2.f90
	gfortran -o $(BUILDDIR)/main2.o -c $(OPT) $(SRCDIR)/main2.f90 -J$(BUILDDIR)
	
clean:
	rm -f $(BUILDDIR)/main0.o $(BUILDDIR)/main1.o $(BUILDDIR)/main2.o $(BUILDDIR)/function_module.o $(BUILDDIR)/function_module.mod $(BUILDDIR)/statistic_module.o $(BUILDDIR)/statistic_module.mod program.exe program1.exe program2.exe

build:
	mkdir $(BUILDDIR)

run:
	./program.exe
	./program1.exe
	./program2.exe
	
plot:

	gnuplot gnuplot/plot.gp
	gnuplot gnuplot/plot2.gp
