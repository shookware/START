#
# This Makefile was automatically generated by Code::Blocks IDE.
#

EXE = START
FC = gfortran #ifort
IDIR =

#CFLAGS = -O2 -g -module $(OBJS_DIR) $(IDIR) -CB -ftrapuv -init=snan -traceback  -cpp
#CFLAGS = -O0 -g -module $(OBJS_DIR) $(IDIR) -CB -traceback -cpp -D DEBUG
#CFLAGS = -O2 -g  $(IDIR) -J$(OBJS_DIR) -cpp -I/opt/intel/compilers_and_libraries_2017/mac/mkl/include/
CFLAGS = -O0 -g -J$(OBJS_DIR) $(IDIR) -cpp   -g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Wunused-parameter -I/opt/intel/compilers_and_libraries_2017/mac/mkl/include/

#LFLAGS = -mkl
LFLAGS = -L/opt/intel/compilers_and_libraries_2017/mac/mkl/lib/ -lmkl_rt
LIBS =

OBJS_DIR = obj/
EXE_DIR = bin/

SRC_DIR_f90_MAIN = src/
SRC_DIR_f90_BeFor64 = src/BeFoR64/
SRC_DIR_f90_cfgio = src/cfgio/
SRC_DIR_F90_PENF = src/PENF/
SRC_DIR_F90_StringiFor = src/StringiFor/

VPATH = $(SRC_DIR_f90_MAIN):$(OBJS_DIR):$(SRC_DIR_f90_BeFor64):$(OBJS_DIR):$(SRC_DIR_f90_cfgio):$(OBJS_DIR):$(SRC_DIR_F90_PENF):$(OBJS_DIR):$(SRC_DIR_F90_StringiFor):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90_MAIN) $(OBJS_f90_BeFor64) $(OBJS_f90_cfgio) $(OBJS_F90_PENF) $(OBJS_F90_StringiFor))

SRCS_f90_MAIN = \
module_debug.f90 \
module_spectral_method.f90 \
module_nonlinear_term.f90 \
module_npse.f90 \
module_npse_eqn.f90 \
module_lpse_apse.f90 \
module_apse_lpse.f90 \
module_lst.f90 \
module_apse.f90 \
module_lst_eqn_IR.f90 \
module_lpse.f90 \
module_lpse_eqn.f90 \
module_lame.f90 \
module_baseflow_normal.f90 \
module_coordinate.f90 \
module_lns_op_normal.f90 \
module_baseflow_org.f90 \
module_parameter.f90 \
module_sparsematrix.f90 \
module_lst_dis_op_point.f90 \
module_baseflow.f90 \
module_lpse_bc.f90 \
module_context.f90 \
module_vector.f90 \
module_lst_dis_op_normal.f90 \
module_gas.f90 \
module_dis_flux.f90 \
module_op_mat_cmplex.f90 \
main.f90 \
module_op_mat.f90 \
module_apse_eqn.f90 \
module_dis_normal.f90 \
module_dis_shape.f90 \
module_tensor.f90 \
module_vector_cmplx.f90 \
module_lns_op_point.f90 \
module_pardiso_adapter.f90 \
module_sparsekit.f90 \
module_pardiso.f90 \
module_grid.f90 \
module_print.f90 \
module_basis.f90 \
module_cfgio_adapter.f90 \
module_generic_bc.f90 \
module_lpse_dis_normal.f90 \
module_lpse_dis_op_point.f90 \
module_dis_wavenum.f90 \
module_difference.f90 \
module_dis.f90 \
module_point.f90 \
module_inverse_rayleigh.f90 \
module_lpse_dis_op_normal.f90 \
module_solver.f90 \
module_gmres_adapter.f90 \
module_mkl_gmres.f90 \
module_fft.f90

SRCS_f90_BeFor64 = \
befor64_pack_data_m.F90 \
befor64.F90

SRCS_f90_cfgio = \
string_conv_mod.f90 \
cfgio_mod.f90

SRCS_F90_PENF = \
penf_global_parameters_variables.F90 \
penf.F90 \
penf_b_size.F90 \
penf_stringify.F90

SRCS_F90_StringiFor = \
stringifor_string_t.F90 \
stringifor.F90

OBJS_f90_MAIN = \
module_debug.o \
module_spectral_method.o \
module_nonlinear_term.o \
module_npse.o \
module_npse_eqn.o \
module_lpse_apse.o \
module_apse_lpse.o \
module_lst.o \
module_apse.o \
module_lst_eqn_IR.o \
module_lpse_eqn.o \
module_lame.o \
module_baseflow_normal.o \
module_coordinate.o \
module_lns_op_normal.o \
module_baseflow_org.o \
module_parameter.o \
module_sparsematrix.o \
module_lst_dis_op_point.o \
module_baseflow.o \
module_lpse_bc.o \
module_context.o \
module_vector.o \
module_lst_dis_op_normal.o \
module_gas.o \
module_dis_flux.o \
module_op_mat_cmplex.o \
main.o \
module_op_mat.o \
module_apse_eqn.o \
module_dis_normal.o \
module_dis_shape.o \
module_tensor.o \
module_vector_cmplx.o \
module_lns_op_point.o \
module_pardiso_adapter.o \
module_sparsekit.o \
module_pardiso.o \
module_grid.o \
module_print.o \
module_basis.o \
module_cfgio_adapter.o \
module_generic_bc.o \
module_lpse_dis_normal.o \
module_lpse_dis_op_point.o \
module_dis_wavenum.o \
module_difference.o \
module_dis.o \
module_point.o \
module_inverse_rayleigh.o \
module_lpse_dis_op_normal.o \
module_solver.o \
module_gmres_adapter.o \
module_mkl_gmres.o \
module_fft.o \
module_lpse.o

OBJS_f90_BeFor64 = \
befor64_pack_data_m.o \
befor64.o

OBJS_f90_cfgio = \
string_conv_mod.o \
cfgio_mod.o

OBJS_F90_PENF = \
penf_global_parameters_variables.o \
penf.o \
penf_b_size.o \
penf_stringify.o

OBJS_F90_StringiFor = \
stringifor_string_t.o \
stringifor.o


all : $(EXE)

$(EXE) : $(OBJS_f90_MAIN) $(OBJS_f90_BeFor64) $(OBJS_f90_cfgio) $(OBJS_F90_PENF) $(OBJS_F90_StringiFor)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90_MAIN):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90_MAIN)$(@:.o=.f90)  -o $(OBJS_DIR)$@

$(OBJS_f90_BeFor64):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90_BeFor64)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_f90_cfgio):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90_cfgio)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_F90_PENF):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90_PENF)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_F90_StringiFor):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_F90_StringiFor)$(@:.o=.F90) -o $(OBJS_DIR)$@

clean :
	rm -fr $(OBJS_DIR)
	rm  *.mod
	rm -f obj/*.o

# Dependencies of files
module_lpse.o:\
    module_lpse.f90 \
    module_lpse_eqn.o
module_lst_eqn_IR.o: \
    module_lst_eqn_IR.f90 \
    module_baseflow.o \
    module_baseflow_normal.o \
    module_coordinate.o \
    module_difference.o \
    module_dis.o \
    module_dis_flux.o \
    module_dis_shape.o \
    module_dis_wavenum.o \
    module_grid.o \
    module_inverse_rayleigh.o \
    module_lns_op_normal.o \
    module_lpse_dis_normal.o \
    module_lst_dis_op_normal.o \
    module_lst_dis_op_point.o \
    module_sparsematrix.o \
    module_vector_cmplx.o \
    penf.o
module_solver.o: \
    module_solver.f90 \
    module_grid.o \
    module_baseflow.o \
    module_coordinate.o \
    module_dis.o
module_gmres_adapter.o: \
    module_gmres_adapter.f90 \
    module_mkl_gmres.o \
    module_sparsematrix.o \
    module_pardiso_adapter.o
module_mkl_gmres.o: \
    module_mkl_gmres.f90 \
    module_solver.o
module_lpse_eqn.o: \
    module_lpse_eqn.f90 \
    module_gmres_adapter.o \
    module_solver.o \
    module_baseflow.o \
    module_baseflow_org.o \
    module_baseflow_normal.o \
    module_coordinate.o \
    module_difference.o \
    module_dis.o \
    module_dis_flux.o \
    module_dis_wavenum.o \
    module_grid.o \
    module_lns_op_normal.o \
    module_lpse_bc.o \
    module_lpse_dis_normal.o \
    module_lpse_dis_op_normal.o \
    module_lpse_dis_op_point.o \
    module_lst_eqn_IR.o \
    module_parameter.o \
    module_pardiso_adapter.o \
    module_sparsematrix.o \
    module_vector_cmplx.o \
    penf.o \
    stringifor.o
module_lame.o: \
    module_lame.f90 \
    module_vector.o \
    penf.o
module_baseflow_normal.o: \
    module_baseflow_normal.f90 \
    module_baseflow.o \
    module_baseflow_org.o \
    module_difference.o \
    module_grid.o \
    penf.o
module_coordinate.o: \
    module_coordinate.f90 \
    module_basis.o \
    module_grid.o \
    module_lame.o \
    penf.o \
    stringifor.o
module_lns_op_normal.o: \
    module_lns_op_normal.f90 \
    module_baseflow_normal.o \
    module_grid.o \
    module_lns_op_point.o \
    module_coordinate.o
module_baseflow_org.o: \
    module_baseflow_org.f90 \
    module_vector.o \
    penf.o
module_parameter.o: \
    module_parameter.f90 \
    penf.o
module_sparsematrix.o: \
    module_sparsematrix.f90 \
    module_sparsekit.o \
    penf.o
module_lst_dis_op_point.o: \
    module_lst_dis_op_point.f90 \
    module_dis_wavenum.o \
    module_lns_op_point.o \
    module_coordinate.o \
    module_op_mat.o \
    module_op_mat_cmplex.o \
    penf.o
befor64_pack_data_m.o: \
    befor64_pack_data_m.F90 \
    penf.o
befor64.o: \
    befor64.F90 \
    befor64_pack_data_m.o \
    penf.o
module_baseflow.o: \
    module_baseflow.f90 \
    module_baseflow_org.o \
    module_difference.o \
    module_grid.o \
    module_vector.o \
    penf.o \
    stringifor.o
module_lpse_bc.o: \
    module_lpse_bc.f90 \
    module_dis_wavenum.o \
    module_generic_bc.o \
    module_lns_op_point.o \
    module_coordinate.o \
    penf.o
module_context.o: \
    module_context.f90 \
    module_apse_lpse.o \
    module_lpse_apse.o \
    module_apse_eqn.o \
    module_baseflow.o \
    module_cfgio_adapter.o \
    module_coordinate.o \
    module_difference.o \
    module_dis.o \
    module_dis_wavenum.o \
    module_grid.o \
    module_lpse.o \
    module_lst_eqn_IR.o \
    module_lst.o
string_conv_mod.o: \
    string_conv_mod.f90
module_vector.o: \
    module_vector.f90 \
    penf.o
module_lst_dis_op_normal.o: \
    module_lst_dis_op_normal.f90 \
    module_difference.o \
    module_dis_wavenum.o \
    module_lns_op_normal.o \
    module_coordinate.o \
    module_lst_dis_op_point.o
module_gas.o: \
    module_gas.f90 \
    module_parameter.o \
    penf.o
module_dis_flux.o: \
    module_dis_flux.f90 \
    module_vector_cmplx.o \
    module_fft.o \
    penf.o
module_op_mat_cmplex.o: \
    module_op_mat_cmplex.f90 \
    module_op_mat.o \
    module_parameter.o \
    penf.o
main.o: \
    main.f90 \
    module_cfgio_adapter.o \
    module_context.o
module_op_mat.o: \
    module_op_mat.f90 \
    module_parameter.o \
    penf.o
module_apse_eqn.o: \
    module_apse_eqn.f90 \
    module_baseflow.o \
    module_baseflow_org.o \
    module_baseflow_normal.o \
    module_coordinate.o \
    module_difference.o \
    module_dis.o \
    module_dis_flux.o \
    module_dis_wavenum.o \
    module_grid.o \
    module_lns_op_normal.o \
    module_lpse_bc.o \
    module_lpse_dis_normal.o \
    module_lpse_dis_op_normal.o \
    module_lpse_dis_op_point.o \
    module_lst_eqn_IR.o \
    module_parameter.o \
    module_pardiso_adapter.o \
    module_sparsematrix.o \
    module_vector_cmplx.o \
    penf.o \
    stringifor.o
module_dis_normal.o: \
    module_dis_normal.f90 \
    module_difference.o \
    module_dis_flux.o \
    module_dis_shape.o \
    module_dis_wavenum.o
module_dis_shape.o: \
    module_dis_shape.f90 \
    module_difference.o \
    module_dis_flux.o \
    module_grid.o \
    penf.o
module_tensor.o: \
    module_tensor.f90 \
    module_vector.o \
    penf.o
module_vector_cmplx.o: \
    module_vector_cmplx.f90 \
    penf.o
cfgio_mod.o: \
    cfgio_mod.f90 \
    string_conv_mod.o
penf_global_parameters_variables.o: \
    penf_global_parameters_variables.F90
module_lns_op_point.o: \
    module_lns_op_point.f90 \
    module_baseflow_org.o \
    module_baseflow_normal.o \
    module_gas.o \
    module_coordinate.o \
    module_op_mat.o \
    module_op_mat_cmplex.o \
    module_parameter.o \
    module_vector.o
module_pardiso_adapter.o: \
    module_pardiso_adapter.f90 \
    module_pardiso.o \
    module_sparsematrix.o \
    penf.o
module_sparsekit.o: \
    module_sparsekit.f90 \
    penf.o
module_pardiso.o: \
    module_pardiso.f90 \
    module_sparsematrix.o \
    penf.o
penf.o: \
    penf.F90 \
    penf_b_size.o \
    penf_global_parameters_variables.o \
    penf_stringify.o
module_grid.o: \
    module_grid.f90 \
    module_point.o \
    penf.o \
    stringifor.o
module_print.o: \
    module_print.f90 \
    penf.o \
    stringifor.o
module_basis.o: \
    module_basis.f90 \
    module_vector.o \
    penf.o
module_cfgio_adapter.o: \
    module_cfgio_adapter.f90 \
    cfgio_mod.o \
    module_dis_wavenum.o \
    module_parameter.o \
    penf.o \
    stringifor.o
module_generic_bc.o: \
    module_generic_bc.f90 \
    module_dis_wavenum.o \
    module_lns_op_point.o \
    module_coordinate.o \
    penf.o
stringifor_string_t.o: \
    stringifor_string_t.F90 \
    befor64.o \
    penf.o
stringifor.o: \
    stringifor.F90 \
    penf.o \
    stringifor_string_t.o
module_lpse_dis_normal.o: \
    module_lpse_dis_normal.f90 \
    module_difference.o \
    module_dis_flux.o \
    module_dis_normal.o \
    module_dis_shape.o \
    module_dis_wavenum.o \
    module_lpse_dis_op_normal.o \
    module_lpse_dis_op_point.o \
    penf.o
module_lpse_dis_op_point.o: \
    module_lpse_dis_op_point.f90 \
    module_dis_wavenum.o \
    module_lns_op_point.o \
    module_coordinate.o \
    module_lpse_bc.o \
    module_op_mat.o \
    module_op_mat_cmplex.o \
    penf.o
module_dis_wavenum.o: \
    module_dis_wavenum.f90 \
    penf.o
penf_b_size.o: \
    penf_b_size.F90 \
    penf_global_parameters_variables.o
module_difference.o: \
    module_difference.f90 \
    module_baseflow_org.o \
    module_dis_flux.o \
    module_grid.o \
    module_vector.o \
    penf.o \
    module_parameter.o
module_dis.o: \
    module_dis.f90 \
    module_dis_flux.o \
    module_dis_wavenum.o \
    module_grid.o \
    module_lpse_dis_normal.o \
    module_print.o \
    module_vector_cmplx.o \
    penf.o \
    stringifor.o
penf_stringify.o: \
    penf_stringify.F90 \
    penf_b_size.o \
    penf_global_parameters_variables.o
module_point.o: \
    module_point.f90 \
    penf.o
module_inverse_rayleigh.o: \
    module_inverse_rayleigh.f90 \
    module_pardiso.o \
    module_sparsematrix.o \
    penf.o
module_lpse_dis_op_normal.o: \
    module_lpse_dis_op_normal.f90 \
    module_difference.o \
    module_dis_wavenum.o \
    module_lns_op_normal.o \
    module_coordinate.o \
    module_lpse_bc.o \
    module_lpse_dis_op_point.o
module_fft.o: \
    module_fft.f90 \
    penf.o
module_apse.o: \
    module_apse.f90 \
    module_apse_eqn.o \
    module_solver.o \
    module_lst_eqn_IR.o
module_lst.o: \
    module_lst.f90
module_lpse_apse.o: \
    module_lpse_apse.f90 \
    module_solver.o \
    module_lpse_eqn.o \
    module_apse_eqn.o
module_apse_lpse.o: \
    module_apse_lpse.f90 \
    module_solver.o
module_npse.o: \
    module_npse.f90 \
    module_solver.o \
    module_npse_eqn.o
module_nonlinear_term.o:\
    module_nonlinear_term.f90 \
    module_fft.o \
    module_baseflow_normal.o \
    module_lpse_dis_normal.o \
    module_spectral_method.o \
    module_cfgio_adapter.o
module_npse_eqn.o: \
    module_npse_eqn.f90 \
    module_lpse_eqn.o \
    module_cfgio_adapter.o \
    module_nonlinear_term.o
module_spectral_method.o: \
    module_spectral_method.f90 \
    penf.o \
    module_fft.o \
    module_parameter.o
module_debug.o: \
    module_debug.f90
