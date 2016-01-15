LIBS= -L/usr/X11R6/lib -lX11 -lm
CPGPLOT_LIB= -lcpgplot -lpgplot -lpng -lz


%: %.c
	gcc -c $<
	gcc -W -Wall -o $* $*.o $(LIBS) $(CPGPLOT_LIB)
