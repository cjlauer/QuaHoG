include ../Makefile.in

ifndef V
       V = 0
endif

ifeq ($(V), 1) 
	E = @echo > /dev/null
	C = 
else
	E = @echo
	C = @
endif

all: libqhg.a

CFLAGS += -I./
SOURCES=\
	qhg_lattice			\
	qhg_comms			\
	qhg_alloc			\
	qhg_stop_watch			\
	qhg_gauge_field			\
	qhg_spinor_field		\
	qhg_correlator			\
	qhg_xchange_gauge		\
	qhg_xchange_spinor		\
	qhg_import_gauge_field		\
	qhg_gauge_field_io		\
	qhg_point_spinor_field		\
	qhg_plaquette			\
	qhg_mom_list			\
	qhg_mesons			\
	qhg_write_mesons		\
	qhg_mesons_udsc			\
	qhg_write_mesons_udsc		\
	qhg_nucleons			\
	qhg_write_nucleons		\
	qhg_ape_smear			\
	qhg_su3_project			\
	qhg_gauss_smear			\
	qhg_io_utils			\
	qhg_spinor_field_io		\
	qhg_spinor_field_bc		\
	qhg_spinor_field_linalg		\
	qhg_spinor_field_rms		\
	qhg_nn_thrp			\
	qhg_write_nn_thrp		\
	qhg_write_mesons_thrp		\
	qhg_prop_field_linalg		\
	qhg_ft				\
	qhg_nn_sequential_sink		\
	qhg_conjugate_thrp 		\
	qhg_correlator_shift

include $(wildcard *.d)

libqhg.a: ${addsuffix .o, $(SOURCES)}
	$(E) LD $@
	$(C)$(AR) cr $@  ${addsuffix .o, $(SOURCES)}
	$(C)$(RANLIB) $@

%.o: %.c
	$(E) CC $<
	$(C)$(CC) $(CFLAGS) -c $< -o $@
	$(C)$(CC) -M $(CFLAGS) -c $< > $*.d

clean:
	$(E) CLEAN in $(PWD)
	$(C)$(RM) *.o *.d

cleanall: clean
	$(E) CLEANALL in $(PWD)
	$(C)$(RM) libqhg.a
