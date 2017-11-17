#ifndef __c4_Model_kalmanTest_h__
#define __c4_Model_kalmanTest_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc4_Model_kalmanTestInstanceStruct
#define typedef_SFc4_Model_kalmanTestInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c4_sfEvent;
  boolean_T c4_isStable;
  boolean_T c4_doneDoubleBufferReInit;
  uint8_T c4_is_active_c4_Model_kalmanTest;
  real_T (*c4_C)[6];
  real_T *c4_u;
  real_T *c4_x_dotdot;
  real_T *c4_w;
  real_T *c4_f;
  real_T *c4_x;
  real_T *c4_x_dot;
  real_T *c4_theta;
  real_T *c4_theta_dot;
  real_T *c4_theta_dotdot;
} SFc4_Model_kalmanTestInstanceStruct;

#endif                                 /*typedef_SFc4_Model_kalmanTestInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c4_Model_kalmanTest_get_eml_resolved_functions_info
  (void);

/* Function Definitions */
extern void sf_c4_Model_kalmanTest_get_check_sum(mxArray *plhs[]);
extern void c4_Model_kalmanTest_method_dispatcher(SimStruct *S, int_T method,
  void *data);

#endif
