#
#

DEBUG    := TRUE

SST_OBJS := test_sst.o sst_data.o precision.o pmgrid.o rgrid.o \
            commap.o mpishorthand.o bnddyi.o calendr.o wrap_nf.o ioFileMod.o \
            endrun.o shr_sys_mod.o comsrf.o ppgrid.o physconst.o phys_grid.o \
            error_messages.o dyn_grid.o

test_sst: $(SST_OBJS)
	$(FC) -o $@ $(LDFLAGS) $(SST_OBJS) $(LDFLAGS)

-include ../../../../bld/Makefile
