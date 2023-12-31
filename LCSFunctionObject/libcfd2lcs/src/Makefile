#
#Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
#libcfd2lcs is distributed is under the terms of the GNU General Public License
#
#Include any platform specific settings
include ./../Makefile.in

#Define the headers for single & double prec.:
FHEADER_DP = cfd2lcs_inc_double.f90
CHEADER_DP = cfd2lcs_inc_double.h
FHEADER_SP = cfd2lcs_inc_single.f90
CHEADER_SP = cfd2lcs_inc_single.h

#Set the correct io_m.f90 source file
ifeq ($(HDF5_SUPPORT),"TRUE")
IO = io_hdf5_m.f90
else
IO = io_nohdf5_m.f90
endif

OBJS = \
	data_m.o \
	geometry_m.o \
	structured_m.o \
	unstructured_m.o \
	lp_m.o \
	comms_m.o \
	gradient_m.o \
	io_m.o \
	sgrid_m.o \
	lp_tracking_m.o \
	lp_motion_m.o \
	flowmap_m.o \
      invariants_m.o \
	lcs_m.o \
	cfd2lcs_c.o \
	cfd2lcs.o

.SUFFIXES:
.SUFFIXES: .f90 .o
.SUFFIXES: .c .o

default:
	$(info -----------------------------------------------------------)
	$(info ERROR: Compile with one of the following:				  )
	$(info 1.	make sp:  single precision          				  )
	$(info 2.	make dp:  double precision          				  )
	$(info 2.	make sp_shared:  shared lib single precision          )
	$(info 2.	make dp_shared:  shared lib double precision          )
	$(info -----------------------------------------------------------)

sp: headers_sp libbspline_sp.a libcfd2lcs_sp.a
	mkdir -p ./../lib
	cp libcfd2lcs_sp.a ./../lib
	$(info *** LIBCFD2LCS SINGLE PRECISION COMPILATION SUCCESSFUL *****)

sp_shared: headers_sp libbspline_sp.so libcfd2lcs_sp.so
	mkdir -p ./../lib
	cp libcfd2lcs_sp.so ./../lib
	cp ./../lib/libcfd2lcs_sp.so $(FOAM_USER_LIBBIN)/
	$(info *** DYNAMIC LIBCFD2LCS SINGLE PRECISION COMPILATION SUCCESSFUL *****)

dp: headers_dp libbspline_dp.a libcfd2lcs_dp.a
	mkdir -p ./../lib
	cp libcfd2lcs_dp.a ./../lib
	$(info *** LIBCFD2LCS DOUBLE PRECISION COMPILATION SUCCESSFUL *****)

dp_shared: headers_dp libbspline_dp.so libcfd2lcs_dp.so
	mkdir -p ./../lib
	cp libcfd2lcs_dp.so ./../lib
	cp ./../lib/libcfd2lcs_dp.so $(FOAM_USER_LIBBIN)/
	$(info *** DYNAMIC LIBCFD2LCS DOUBLE PRECISION COMPILATION SUCCESSFUL *****)


libcfd2lcs_sp.so libcfd2lcs_dp.so:	$(OBJS)
	$(CC) -fPIC -shared -o $@  $(OBJS)

sp_dep: sp libcfd2lcs_sp_dep.a

headers_sp:
	mkdir -p ./../include
	$(info *******Compiling libcfd2lcs for SINGLE precision************)
	cp $(FHEADER_SP) ./../include/cfd2lcs_inc_sp.f90
	cp $(CHEADER_SP) ./../include/cfd2lcs_inc_sp.h
	ln -fs $(FHEADER_SP) cfd2lcs_inc.f90
	ln -fs $(CHEADER_SP) cfd2lcs_inc.h
	ln -fs $(IO) io_m.f90

headers_dp:
	mkdir -p ./../include
	$(info *******Compiling libcfd2lcs for DOUBLE precision************)
	cp $(FHEADER_DP) ./../include/cfd2lcs_inc_dp.f90
	cp $(CHEADER_DP) ./../include/cfd2lcs_inc_dp.h
	ln -fs $(FHEADER_DP) cfd2lcs_inc.f90
	ln -fs $(CHEADER_DP) cfd2lcs_inc.h
	ln -fs $(IO) io_m.f90

libbspline_sp.a:
	cd bspline; make sp

libbspline_sp.so:
	cd bspline; make sp_shared

libbspline_dp.a:
	cd bspline; make dp

libbspline_dp.so:
	cd bspline; make dp_shared

libcfd2lcs_sp.a libcfd2lcs_dp.a:	$(OBJS)
	$(AR) $@ $(OBJS)

.c.o:
	$(CC) $(CFLAGS) -$(INC). -$(INC)./bspline/ $(MPI_INC) $(HDF5_INC) $*.c

.f90.o:
	$(F90) $(FFLAGS) -$(INC). -$(INC)./bspline/ $(MPI_INC) $(HDF5_INC) $*.f90

clean:
	cd bspline; make clean
	rm -f *.o *.mod *.a
	rm -f cfd2lcs_inc.f90
	rm -f cfd2lcs_inc.h

libclean:
	cd bspline; make clean
	rm -f *.o *.mod *.a
	cd ../include/; rm -f *.f90 *.h
	cd ../lib/; rm -f *.a
	rm -f cfd2lcs_inc.f90
	rm -f cfd2lcs_inc.h
