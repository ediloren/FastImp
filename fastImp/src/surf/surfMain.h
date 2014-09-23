/*-*-c++-*-*/
/******************************************************************************
 **                                                                          **
 ** Copyright (c) 2000 Massachusetts Institute of Technology, Cambridge, MA. **
 ** All rights reserved.                                                     **
 **                                                                          **
 ******************************************************************************
  ==========================================================================
  ============================= File Info ==================================

  Author: Zhenhai Zhu
  
  Description:

  Resources:

  See also:

  ==========================================================================
*/

#ifndef __SURF_MAIN_H_
#define __SURF_MAIN_H_

#include "surf.h"
#include "surfConst.h"
#include <getopt.h>
#include <string> // for exit()

namespace surf {

  const std::string SURF_VERSION = "fastimp v1.0";
  const std::string SURF_DATE = "07/01/03";

  enum cmdLineOptions {
    //    ACCELARATION_SWITCH        = 'a',
    FREQUENCY_BEGIN            = 'b',
    CURRENT_COMP_MODE          = 'c',
    SOLVER_TYPE                = 'd',
    FREQUENCY_END              = 'e',
    FORMULATION_TYPE           = 'f',
    GLOBAL_COORD_SWITCH        = 'g',
    HELP                       = 'h',
    INPUT_FILE_NAME            = 'i',
    ONE_COL_ADMMITANCE         = 'l',
    FREQUENCY_SAM_METHOD       = 'm',
    FREQUENCY_NUM_STEP         = 'n',
    OUTPUT_MESH_FILE           = 'o',
    PRE_CONDITIONER_TYPE       = 'p',
    QUAD_PANEL_SWITCH          = 'q',
    READ_MESH_FILE             = 'r',
    SIMULATION_TYPE            = 's',
    TRANSMISSION_LINE_SWITCH   = 't',
    UNIFORM_MESH               = 'u',
    VERSION_NUM                = 'v',
    WARNING_SWITCH             = 'w'
  };

  const char* simuTypeName[3] = {"MQS", "EMQS", "FULL_WAVE"};
  const char* formulationTypeName[2] = {"OLD", "NEW"};
  const char* freqSamMethodName[2] = {"LOGRITHMIC", "LINEAR"};
  const char* preCondTypeName[3] = {"NONE", "LEFT_PRE_CONDITIONER", 
				    "RIGHT_PRE_CONDITIONER"}; 
  const char* solverName[3] = {"DIRECT", "ITERATIVE", "ITERATIVE_WITH_pFFT"}; 
  const char* currentCompModeName[5] = {"AUTO", "LOW_FREQUENCY_MODE", 
					"HIGH_FREQUENCY_MODE", "2D_MODE",
					"ENERGY_MODE" };

  void cmdLineParsor(int argc, char *argv[]);
  void checkCmdLineParameter (const std::vector<surf::CondInfo*>& condInfoPtrList);
  void printUsage (void);
  void printVersion (void);
}

#endif
