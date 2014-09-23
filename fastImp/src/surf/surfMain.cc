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

  const static char cvsid[] = "$Id: surfMain.cc,v 1.21 2003/07/18 20:51:57 zhzhu Exp $";

  ==========================================================================
*/

#include "surfMain.h"

using namespace surf;
using namespace std;

static char* inputStructFile;
static char* inputMeshFile = 0;
static char* outputMeshFile = 0;
static SimulationType simuType = MQS;
static FrequencySample freqSam(1e9, 1e9, 1);
static int oneColIndex = -1;
static SolverType solverType = ITERATIVE_WITH_PFFT;
static PreconditionerType preCondType = RIGHT_PRE_CONDITIONER;
static FormulationType formulationType = OLD;
static CurrentCompMode currentCompMode = AUTO_MODE;
static bool useNonUniformMesh = true;
static bool useGlobalCoord = false;
static bool isTransmissionLine = false;
static bool useQuadPanel = true;
static bool checkPanelSize = false;

/**********************************************************************
 * main --
 **********************************************************************/
int
main(
     int argc,
     char *argv[])
{  
  printVersion();
  cmdLineParsor(argc, argv);

  vector<surf::CondInfo*> condInfoPtrList;
  readCondInfo(inputStructFile, condInfoPtrList);

  checkCmdLineParameter(condInfoPtrList);

  Surf surf(inputStructFile, outputMeshFile, inputMeshFile, 
	    simuType, preCondType, useGlobalCoord, formulationType, 
	    currentCompMode, solverType, useNonUniformMesh, 
	    checkPanelSize, isTransmissionLine, useQuadPanel, freqSam, 
	    oneColIndex, condInfoPtrList);
  surf.extract();
}

/**********************************************************************
 * cmdLineParsor --
 **********************************************************************/
void
surf::cmdLineParsor (
		     int argc,
		     char *argv[])
{
  if (argc == 1) {
    printUsage();
    exit(0);
  }	

  double freqBegin = freqSam.firstPoint();
  double freqEnd = freqSam.lastPoint();
  int numSam = freqSam.numSam();
  SampleMethod samMethod = freqSam.samMethod();
  opterr = 0;  /* Disable stderr Error Message by getopt */
  int option;
  while ((option = getopt(argc, argv, "b:c:d:e:f:g:hi:l:m:n:o:p:q:r:s:t:u:vw:")) != EOF) {
    switch (option) {
    case INPUT_FILE_NAME:
      inputStructFile = optarg;
      break;
    case READ_MESH_FILE:
      inputMeshFile = optarg;
      break;
    case OUTPUT_MESH_FILE:
      outputMeshFile = optarg;
      break;
    case SIMULATION_TYPE:
      simuType = static_cast<SimulationType>(atoi(optarg));
      break;
    case FREQUENCY_BEGIN:
      sscanf(optarg, "%lf", &freqBegin);
      break;
    case FREQUENCY_END:
      sscanf(optarg, "%lf", &freqEnd);
      break;
    case FREQUENCY_NUM_STEP:
      numSam = atoi(optarg);
      break;
    case FREQUENCY_SAM_METHOD:
      samMethod = static_cast<SampleMethod>(atoi(optarg));
      break;      
    case SOLVER_TYPE:
      if (atoi(optarg) == 1) {
	solverType = ITERATIVE;
      } else if (atoi(optarg) == 2) {
	solverType = ITERATIVE_WITH_PFFT;
      } else {
	solverType = DIRECT;
      }
      break;
    case ONE_COL_ADMMITANCE:
      oneColIndex = atoi(optarg);
      break;
    case TRANSMISSION_LINE_SWITCH:
      isTransmissionLine = (atoi(optarg) == 1);
      break;
    case WARNING_SWITCH:
      checkPanelSize = (atoi(optarg) == 1); 
      break;
    case HELP:
      printUsage();
      exit(0);
      break;
    case VERSION_NUM:
      printVersion();
      exit(0);
      break;
#ifdef DEBUG_VERSION
    case CURRENT_COMP_MODE:
      currentCompMode = static_cast<CurrentCompMode>(atoi(optarg));
      break;
    case PRE_CONDITIONER_TYPE:
      preCondType = static_cast<PreconditionerType>(atoi(optarg));
      break;
    case UNIFORM_MESH:
      useNonUniformMesh = (atoi(optarg) == 0);
      break;
    case GLOBAL_COORD_SWITCH:
      useGlobalCoord = (atoi(optarg) == 1);
      break;
    case FORMULATION_TYPE:
      formulationType = static_cast<FormulationType>(atoi(optarg));
      break;
    case QUAD_PANEL_SWITCH:
      useQuadPanel = (atoi(optarg) == 1);
      break;
#endif
    case '?': /* unknown option */
      cout << endl << endl
	   << "\t Warning in surf:" << endl
	   << "\t Wrong command line options" << endl
	   << "\t please check." << endl;
      printUsage();
      exit(1);
      break;
    default: /* should never happen */
      printUsage();
      exit(2);
    }
  }
  freqSam = FrequencySample(freqBegin, freqEnd, numSam, samMethod);
}

/**********************************************************************
 * checkCmdLineParameter --
 **********************************************************************/
void
surf::checkCmdLineParameter (
			     const vector<surf::CondInfo*>& condInfoPtrList)
{
  if (simuType >= NUM_SIMULATION_TYPE) {
    surf::errorMessage("surfMain.c : checkCmdLineParameter",
		 "Illegal simulation type, Please check command line options!! ");
  }

  if (!inputStructFile) {
    surf::errorMessage("surfMain.c : checkCmdLineParameter",
		       "No structure input file is specified!!");
  }

  if (oneColIndex >= 0) {
    if (oneColIndex >= condInfoPtrList.size()) {
      surf::errorMessage("surfMain.cc : main()",
			 "column index of admmitance matrix is too large");
    }
    if ( condInfoPtrList[oneColIndex]->condType() == CondInfo::GROUND) {
      surf::errorMessage("surfMain.cc : main()",
			 "column index is that of a grounded conductor");
    }
  }

  if (solverType == DIRECT) {
    // since no iterative solver is used the preconditioner is irrelevant.
    preCondType = NO_PRE_CONDITIONER;

    // the global coordinate must be used. Otherwise wrong entry would be used by 
    // port currect calculation.
    useGlobalCoord = true; 
  }

  if (formulationType >= NUM_FORMULATION_TYPE) {
    surf::errorMessage("surfMain.c : checkCmdLineParameter",
		 "Illegal formulation type, Please check command line options!!");
  }

  if (preCondType >= NUM_PRE_CONDITIONER_TYPE) {
    surf::errorMessage("surfMain.c : checkCmdLineParameter",
		 "Illegal pre-conditioner type, Please check command line options!!");
  }

  if (currentCompMode >= NUM_CURRENT_COMP_MODE) {
    surf::errorMessage("surfMain.c : checkCmdLineParameter",
		 "Illegal current computation mode, Please check command line options!!");
  }
}

/**********************************************************************
 * printUsage --
 **********************************************************************/
void
surf::printUsage (
		  void)
{
  cout << endl 
       << "\t Usage:  fastimp <options>" << endl
       << "\t Where the options are:" << endl
       << "\t -i <input structure file name>" << endl
       << "\t -r <input mesh file name>" << endl
       << "\t -o <output mesh file name>" << endl
       << "\t    A patran-like format mesh file will be generated." << endl
       << "\t    If combined with -r, it generates a fastcap input file. " << endl
       << "\t    This option could be used to visualize mesh" << endl
       << "\t -s <simulation type>" << endl
       << "\t    0: MQS, 1: EMQS, 2: Full_Wave" << endl
       << "\t    default is " << simuTypeName[simuType] << endl
       << "\t -b <begin of the frequncy range>"  << endl
       << "\t    default is " << freqSam.firstPoint() << "Hz" << endl
       << "\t -e <end of the frequncy range>" << endl
       << "\t    default is " << freqSam.lastPoint() << "Hz"  << endl
       << "\t -n <number of frequency sampling points>" << endl
       << "\t    default is " << freqSam.numSam() << endl
       << "\t -m <frequency sampling method>" << endl
       << "\t    0:LOGARITHMIC, 1:LINEAR" << endl
       << "\t    default is " << freqSamMethodName[freqSam.samMethod()] << endl
       << "\t -d <solver type>" << endl
       << "\t    1: gmres, 2: gmres+pfft, otherwise: superLU" << endl
       << "\t    default is " << solverName[solverType] << endl
       << "\t -l <column index of the admittance matrix>" << endl
       << "\t    default is extracting the whole matrix." << endl
       << "\t -t <transmission line switch>" << endl
       << "\t    1: treat the structure as a T-line and short one end" << endl
       << "\t    otherwise: treat it as a regular 3D structure" << endl
       << "\t    default is " << isTransmissionLine << endl
       << "\t -w <check panel size>" << endl
       << "\t    0: do not check, otherwise: check" << endl
       << "\t    default is " << checkPanelSize << endl
       << "\t -h Print this message and Exit" << endl
       << "\t -v Print Version Number and Exit" << endl 
#ifdef DEBUG_VERSION
       << "\t -p <preconditioner type>" << endl
       << "\t    0: No preconditioner, 1: LEFT, 2: RIGHT" << endl
       << "\t    default is " << preCondTypeName[preCondType] << endl
       << "\t -u <mesh type>" << endl
       << "\t    0: non-uniform, otherwise: uniform" << endl
       << "\t    default is " << !useNonUniformMesh << endl
       << "\t -q <quadrilateral panel switch>" << endl
       << "\t    1: use quad, otherwise: use triangle panels" << endl
       << "\t    default is " << useQuadPanel << endl
       << "\t    If mesh is read from a file, this option has no effect." << endl
       << "\t -c <port current computing mode>" << endl
       << "\t    0:AUTO, 1:LOW_FREQ, 2:HIGH_FREQ" << endl
       << "\t    default is " << currentCompModeName[currentCompMode] << endl
       << "\t -f <formulation index>" << endl
       << "\t    0: OLD, 1: NEW" << endl
       << "\t    default is "  << formulationTypeName[formulationType] << endl
       << "\t -g <pre-conditioner in global/local coordinate switch>" << endl
       << "\t    1: global, otherwise: local" << endl
       << "\t    default is " << useGlobalCoord << endl
#endif
       << endl;
}

/**********************************************************************
 * printVersion --
 **********************************************************************/
void
surf::printVersion (
		    void)
{
  cout << endl << "\t " << SURF_VERSION << "  " << SURF_DATE << endl;  
  cout << "\t FAST IMPedance extraction program" << endl;
}

