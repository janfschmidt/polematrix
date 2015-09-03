CC=g++
CCFLAGS=-std=c++11 -O3 -Wall -Wl,--no-as-needed #-g -O0 #last 2 for valgrind
INSTALL_PATH=/usr/local/
PROG_NAME=polematrix
ALL_O=main.o TrackingTask.o Tracking.o Configuration.o


all: $(PROG_NAME)
.PHONY: all

$(PROG_NAME): $(ALL_O)
	$(CC) $(CCFLAGS) -o $@ $(ALL_O) -lopenblas -llapack -larmadillo -pthread -lboost_program_options -lpalattice -lgsl -lgslcblas -lm

main.o: main.cpp TrackingTask.hpp Configuration.hpp
	$(CC) $(CCFLAGS) -c $<
TrackingTask.o: TrackingTask.cpp TrackingTask.hpp Configuration.hpp
	$(CC) $(CCFLAGS) -c $<
Tracking.o: Tracking.cpp Tracking.hpp TrackingTask.hpp Configuration.hpp
	$(CC) $(CCFLAGS) -c $<
Configuration.o: Configuration.cpp Configuration.hpp
	$(CC) $(CCFLAGS) -c $<


clean: 
	rm $(PROG_NAME) $(ALL_O)

install:
	install -m 755 -p -D -v $(PROG_NAME) $(INSTALL_PATH)/bin/$(PROG_NAME)
uninstall:
	rm -f $(INSTALL_PATH)/bin/$(PROG_NAME)

