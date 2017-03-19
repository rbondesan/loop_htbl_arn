all: rmxeq  mkxeq
.PHONY: all

# main programs and modules to be compiled 

MAIN = loop_htbl_arn

MOD = comp_matrix hamiltonian loop_htbl_states parameters TLopen TLperiodic alterHS diag_exact  alterCS 

MODULES = $(MOD)

MDIR = modules/

INCPATH = -I include -I arp++_include/

VPATH = .:$(MDIR)

CC = g++ 

#TODO: are they all needed??
LIBS = /usr/local/lib/libsuperlu.a  -larpack -llapack -lblas -lgsl -lgslcblas

OFLAGS = -O2 -Wall #-DDEBUG #-g -DSHOW_EIGENVECTORS #-g -DDEBUG #-pg 

PGMS = $(MAIN) $(MODULES) 

$(addsuffix .o,$(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(OFLAGS) $(INCPATH) 

$(MAIN): %: %.o $(addsuffix .o,$(MODULES)) Makefile
	$(CC) $< $(addsuffix .o,$(MODULES)) $(OFLAGS) -o $(MAIN) $(LIBS) 

# produce executables

mkxeq: $(MAIN)

# remove old executables

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		

# clean directory 

clean:
	@ -rm -rf *.o *~  $(MAIN)

# Once the following line is specified, `make clean' will run the commands 
#regardless of whether there is a file named `clean'.
.PHONY: clean

# create a compressed package and move it in the upper directory
tar:
	@ -tar -cvzf tball.tgz *
	@ -mv tball.tgz ../.		
.PHONY: tar
