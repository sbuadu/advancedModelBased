/* Include files */

#include "Model_sfun.h"
#include "Model_sfun_debug_macros.h"
#include "c1_Model.h"
#include "c2_Model.h"
#include "c3_Model.h"
#include "c4_Model.h"
#include "c5_Model.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _ModelMachineNumber_;

/* Function Declarations */

/* Function Definitions */
void Model_initializer(void)
{
}

void Model_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_Model_method_dispatcher(SimStruct *simstructPtr, unsigned int
  chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_Model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_Model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==3) {
    c3_Model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==4) {
    c4_Model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==5) {
    c5_Model_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

extern void sf_Model_uses_exported_functions(int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[])
{
  plhs[0] = mxCreateLogicalScalar(0);
}

unsigned int sf_Model_process_check_sum_call( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(531020900U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2342136637U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(878341417U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4187411948U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3689752237U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(387893888U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1180458060U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3696567601U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_Model_get_check_sum(mxArray *plhs[]);
          sf_c1_Model_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_Model_get_check_sum(mxArray *plhs[]);
          sf_c2_Model_get_check_sum(plhs);
          break;
        }

       case 3:
        {
          extern void sf_c3_Model_get_check_sum(mxArray *plhs[]);
          sf_c3_Model_get_check_sum(plhs);
          break;
        }

       case 4:
        {
          extern void sf_c4_Model_get_check_sum(mxArray *plhs[]);
          sf_c4_Model_get_check_sum(plhs);
          break;
        }

       case 5:
        {
          extern void sf_c5_Model_get_check_sum(mxArray *plhs[]);
          sf_c5_Model_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1845332328U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(556278208U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3892017426U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2852207787U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1340808022U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(847957673U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3924104721U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3645928761U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_autoinheritance_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "TH7xSaNBwUOFLbrdMXJ5rC") == 0) {
          extern mxArray *sf_c1_Model_get_autoinheritance_info(void);
          plhs[0] = sf_c1_Model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "ZRNDPkmBH3AOPU8uGHPH1E") == 0) {
          extern mxArray *sf_c2_Model_get_autoinheritance_info(void);
          plhs[0] = sf_c2_Model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 3:
      {
        if (strcmp(aiChksum, "pE1JcGkWWBj37Qfq3ihjTD") == 0) {
          extern mxArray *sf_c3_Model_get_autoinheritance_info(void);
          plhs[0] = sf_c3_Model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 4:
      {
        if (strcmp(aiChksum, "KJ7WYei7pWoMZafeqGWcpD") == 0) {
          extern mxArray *sf_c4_Model_get_autoinheritance_info(void);
          plhs[0] = sf_c4_Model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 5:
      {
        if (strcmp(aiChksum, "xBWu5EfO8nD6kiLo3U4SCH") == 0) {
          extern mxArray *sf_c5_Model_get_autoinheritance_info(void);
          plhs[0] = sf_c5_Model_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_get_eml_resolved_functions_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray *sf_c1_Model_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_Model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray *sf_c2_Model_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_Model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 3:
      {
        extern const mxArray *sf_c3_Model_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c3_Model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 4:
      {
        extern const mxArray *sf_c4_Model_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c4_Model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 5:
      {
        extern const mxArray *sf_c5_Model_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c5_Model_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Model_third_party_uses_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the third_party_uses_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_third_party_uses_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "kn7CzH0P5qb04BCM9jYYKD") == 0) {
          extern mxArray *sf_c1_Model_third_party_uses_info(void);
          plhs[0] = sf_c1_Model_third_party_uses_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "mZNePuh5emkJBUqGEeJNfD") == 0) {
          extern mxArray *sf_c2_Model_third_party_uses_info(void);
          plhs[0] = sf_c2_Model_third_party_uses_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "Yalj2gptGXLLFJU5KxAhTF") == 0) {
          extern mxArray *sf_c3_Model_third_party_uses_info(void);
          plhs[0] = sf_c3_Model_third_party_uses_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "GKsPEmantSoy9iHiVPRmeF") == 0) {
          extern mxArray *sf_c4_Model_third_party_uses_info(void);
          plhs[0] = sf_c4_Model_third_party_uses_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "wsp9FZlSBj8oWynX1mJwOE") == 0) {
          extern mxArray *sf_c5_Model_third_party_uses_info(void);
          plhs[0] = sf_c5_Model_third_party_uses_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Model_jit_fallback_info( int nlhs, mxArray * plhs[], int nrhs,
  const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the jit_fallback_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_jit_fallback_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "kn7CzH0P5qb04BCM9jYYKD") == 0) {
          extern mxArray *sf_c1_Model_jit_fallback_info(void);
          plhs[0] = sf_c1_Model_jit_fallback_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "mZNePuh5emkJBUqGEeJNfD") == 0) {
          extern mxArray *sf_c2_Model_jit_fallback_info(void);
          plhs[0] = sf_c2_Model_jit_fallback_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "Yalj2gptGXLLFJU5KxAhTF") == 0) {
          extern mxArray *sf_c3_Model_jit_fallback_info(void);
          plhs[0] = sf_c3_Model_jit_fallback_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "GKsPEmantSoy9iHiVPRmeF") == 0) {
          extern mxArray *sf_c4_Model_jit_fallback_info(void);
          plhs[0] = sf_c4_Model_jit_fallback_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "wsp9FZlSBj8oWynX1mJwOE") == 0) {
          extern mxArray *sf_c5_Model_jit_fallback_info(void);
          plhs[0] = sf_c5_Model_jit_fallback_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

unsigned int sf_Model_updateBuildInfo_args_info( int nlhs, mxArray * plhs[], int
  nrhs, const mxArray * prhs[] )
{
  char commandName[64];
  char tpChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the updateBuildInfo_args_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  mxGetString(prhs[2], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_updateBuildInfo_args_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(tpChksum, "kn7CzH0P5qb04BCM9jYYKD") == 0) {
          extern mxArray *sf_c1_Model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c1_Model_updateBuildInfo_args_info();
          break;
        }
      }

     case 2:
      {
        if (strcmp(tpChksum, "mZNePuh5emkJBUqGEeJNfD") == 0) {
          extern mxArray *sf_c2_Model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c2_Model_updateBuildInfo_args_info();
          break;
        }
      }

     case 3:
      {
        if (strcmp(tpChksum, "Yalj2gptGXLLFJU5KxAhTF") == 0) {
          extern mxArray *sf_c3_Model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c3_Model_updateBuildInfo_args_info();
          break;
        }
      }

     case 4:
      {
        if (strcmp(tpChksum, "GKsPEmantSoy9iHiVPRmeF") == 0) {
          extern mxArray *sf_c4_Model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c4_Model_updateBuildInfo_args_info();
          break;
        }
      }

     case 5:
      {
        if (strcmp(tpChksum, "wsp9FZlSBj8oWynX1mJwOE") == 0) {
          extern mxArray *sf_c5_Model_updateBuildInfo_args_info(void);
          plhs[0] = sf_c5_Model_updateBuildInfo_args_info();
          break;
        }
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;
}

void sf_Model_get_post_codegen_info( int nlhs, mxArray * plhs[], int nrhs, const
  mxArray * prhs[] )
{
  unsigned int chartFileNumber = (unsigned int) mxGetScalar(prhs[0]);
  char tpChksum[64];
  mxGetString(prhs[1], tpChksum,sizeof(tpChksum)/sizeof(char));
  tpChksum[(sizeof(tpChksum)/sizeof(char)-1)] = '\0';
  switch (chartFileNumber) {
   case 1:
    {
      if (strcmp(tpChksum, "kn7CzH0P5qb04BCM9jYYKD") == 0) {
        extern mxArray *sf_c1_Model_get_post_codegen_info(void);
        plhs[0] = sf_c1_Model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 2:
    {
      if (strcmp(tpChksum, "mZNePuh5emkJBUqGEeJNfD") == 0) {
        extern mxArray *sf_c2_Model_get_post_codegen_info(void);
        plhs[0] = sf_c2_Model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 3:
    {
      if (strcmp(tpChksum, "Yalj2gptGXLLFJU5KxAhTF") == 0) {
        extern mxArray *sf_c3_Model_get_post_codegen_info(void);
        plhs[0] = sf_c3_Model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 4:
    {
      if (strcmp(tpChksum, "GKsPEmantSoy9iHiVPRmeF") == 0) {
        extern mxArray *sf_c4_Model_get_post_codegen_info(void);
        plhs[0] = sf_c4_Model_get_post_codegen_info();
        return;
      }
    }
    break;

   case 5:
    {
      if (strcmp(tpChksum, "wsp9FZlSBj8oWynX1mJwOE") == 0) {
        extern mxArray *sf_c5_Model_get_post_codegen_info(void);
        plhs[0] = sf_c5_Model_get_post_codegen_info();
        return;
      }
    }
    break;

   default:
    break;
  }

  plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
}

void Model_debug_initialize(struct SfDebugInstanceStruct* debugInstance)
{
  _ModelMachineNumber_ = sf_debug_initialize_machine(debugInstance,"Model",
    "sfun",0,5,0,0,0);
  sf_debug_set_machine_event_thresholds(debugInstance,_ModelMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(debugInstance,_ModelMachineNumber_,0);
}

void Model_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_Model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("Model", "Model");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_Model_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
