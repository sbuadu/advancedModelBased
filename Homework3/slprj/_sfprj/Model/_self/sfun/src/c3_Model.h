#ifndef __c3_Model_h__
#define __c3_Model_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc3_ModelInstanceStruct
#define typedef_SFc3_ModelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c3_sfEvent;
  boolean_T c3_isStable;
  boolean_T c3_doneDoubleBufferReInit;
  uint8_T c3_is_active_c3_Model;
  real_T *c3_u;
  real_T *c3_w;
  real_T *c3_f;
  real_T *c3_x;
  real_T *c3_x_dot;
  real_T (*c3_data)[6];
  real_T *c3_theta;
  real_T *c3_theta_dotdot;
} SFc3_ModelInstanceStruct;

#endif                                 /*typedef_SFc3_ModelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c3_Model_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c3_Model_get_check_sum(mxArray *plhs[]);
extern void c3_Model_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
