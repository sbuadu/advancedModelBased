#ifndef __c5_Model_kalmanTest_h__
#define __c5_Model_kalmanTest_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc5_Model_kalmanTestInstanceStruct
#define typedef_SFc5_Model_kalmanTestInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c5_sfEvent;
  boolean_T c5_isStable;
  boolean_T c5_doneDoubleBufferReInit;
  uint8_T c5_is_active_c5_Model_kalmanTest;
  real_T (*c5_H)[12];
  real_T (*c5_P_0)[36];
  real_T (*c5_Q)[4];
  real_T (*c5_R)[4];
  real_T (*c5_PHI)[36];
  real_T (*c5_DELTA)[6];
  real_T (*c5_GAMMA)[12];
  real_T (*c5_y)[2];
  real_T *c5_u;
  real_T (*c5_x_barprev)[6];
  real_T *c5_i_prev;
  real_T (*c5_P_prev)[36];
  real_T (*c5_x_hat)[6];
  real_T (*c5_x_barnext)[6];
  real_T *c5_i_next;
  real_T (*c5_P_next)[36];
} SFc5_Model_kalmanTestInstanceStruct;

#endif                                 /*typedef_SFc5_Model_kalmanTestInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c5_Model_kalmanTest_get_eml_resolved_functions_info
  (void);

/* Function Definitions */
extern void sf_c5_Model_kalmanTest_get_check_sum(mxArray *plhs[]);
extern void c5_Model_kalmanTest_method_dispatcher(SimStruct *S, int_T method,
  void *data);

#endif
