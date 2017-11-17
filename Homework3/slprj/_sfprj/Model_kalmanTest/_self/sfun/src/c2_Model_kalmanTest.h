#ifndef __c2_Model_kalmanTest_h__
#define __c2_Model_kalmanTest_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc2_Model_kalmanTestInstanceStruct
#define typedef_SFc2_Model_kalmanTestInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c2_sfEvent;
  boolean_T c2_isStable;
  boolean_T c2_doneDoubleBufferReInit;
  uint8_T c2_is_active_c2_Model_kalmanTest;
  real_T (*c2_H)[12];
  real_T (*c2_P_0)[36];
  real_T (*c2_Q)[4];
  real_T (*c2_R)[4];
  real_T (*c2_PHI)[36];
  real_T (*c2_DELTA)[6];
  real_T (*c2_GAMMA)[12];
  real_T (*c2_y)[2];
  real_T *c2_u;
  real_T (*c2_x_barprev)[6];
  real_T *c2_i_prev;
  real_T (*c2_P_prev)[36];
  real_T (*c2_x_hat)[6];
  real_T (*c2_x_barnext)[6];
  real_T *c2_i_next;
  real_T (*c2_P_next)[36];
} SFc2_Model_kalmanTestInstanceStruct;

#endif                                 /*typedef_SFc2_Model_kalmanTestInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c2_Model_kalmanTest_get_eml_resolved_functions_info
  (void);

/* Function Definitions */
extern void sf_c2_Model_kalmanTest_get_check_sum(mxArray *plhs[]);
extern void c2_Model_kalmanTest_method_dispatcher(SimStruct *S, int_T method,
  void *data);

#endif
