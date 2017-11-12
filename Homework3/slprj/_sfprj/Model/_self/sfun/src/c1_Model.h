#ifndef __c1_Model_h__
#define __c1_Model_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc1_ModelInstanceStruct
#define typedef_SFc1_ModelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c1_sfEvent;
  boolean_T c1_isStable;
  boolean_T c1_doneDoubleBufferReInit;
  uint8_T c1_is_active_c1_Model;
  real_T *c1_u;
  real_T *c1_w;
  real_T *c1_theta;
  real_T (*c1_data)[6];
  real_T *c1_x_dotdot;
  real_T *c1_x;
  real_T *c1_x_dot;
} SFc1_ModelInstanceStruct;

#endif                                 /*typedef_SFc1_ModelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */
extern struct SfDebugInstanceStruct *sfGlobalDebugInstanceStruct;

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c1_Model_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c1_Model_get_check_sum(mxArray *plhs[]);
extern void c1_Model_method_dispatcher(SimStruct *S, int_T method, void *data);

#endif
