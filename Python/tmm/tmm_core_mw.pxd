#import numpy as np
cimport numpy as np

cdef class Tmm_data_internalsource:
    cdef readonly np.ndarray vw_list_l1, vw_list_l2
    cdef readonly np.ndarray vw_list_r1, vw_list_r2
    cdef readonly np.ndarray kz_list
    cdef readonly np.ndarray th_list
    cdef readonly np.ndarray cos_th_list
    cdef readonly int pol
    cdef readonly np.ndarray  n_list
    cdef readonly np.ndarray  d_list
    cdef readonly complex th_i
    cdef readonly double lam_vac
    cdef readonly int source_layer
    cdef readonly double yA
    cdef readonly power_output1, power_output2
    cdef readonly r_list, t_list
    cdef readonly M_list_l, M_list_r
    cdef readonly double phi_i
    cdef readonly double front_loss1, front_loss2
    cdef readonly double back_loss1, back_loss2

cpdef Tmm_data_internalsource coh_tmm_internalsource(int pol, np.ndarray[complex, ndim=1, negative_indices=False] n_list,
            np.ndarray[double, ndim=1, negative_indices=False] d_list, int li, double yA, double phi_i, double lam_vac)


cdef double position_resolved_a(int layer, double dist, Tmm_data_internalsource data)

cdef get_k_vec(complex nk, double lambda0, double phi)

cdef list position_resolved_a_test_h5(int layer, double dist, Tmm_data_internalsource data)

cdef float calc_a_def_int(int layer,  double yB1, double yB2, Tmm_data_internalsource data)
