
FLAG = -I. -L./omd -DFDTD_CONSOLE
INSTALL_DIR=$(HOME)
MAKEHOST =

all:
	make -C drift -f Makefile$(MAKEHOST)
	make -C drift -f Makefile$(MAKEHOST) shlib
	#make -C frontend -f Makefile$(MAKEHOST)

document:
	doxygen doc/dox.cfg

install: all
	sh install$(MAKEHOST)

clean:
#	make -C omd clean
	make -C drift clean
#	make -C frontend clean
