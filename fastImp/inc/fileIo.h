/*-*-c-*-*/
/*
  ==========================================================================
  ============================= File Info ==================================

  =========================== File Contents ================================

  Author: Zhenhai Zhu
  
  Description: 

  Resources:

  See also:

  ==========================================================================
*/

#ifndef __FILE_IO_H_
#define __FILE_IO_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <string.h>

  typedef enum FILE_IO_BOOLEAN {
    FILE_IO_FALSE = 0,
    FILE_IO_TRUE = 1
  } FILE_IO_BOOLEAN;

  extern FILE* Rd_OpenFile(char *fname);
  extern FILE* Rd_OpenFile_Simple(char *fname);
  extern void Rd_CloseFile(void);
  extern FILE_IO_BOOLEAN Rd_Int(int *val);
  extern FILE_IO_BOOLEAN Rd_Double(double *val);
  extern FILE_IO_BOOLEAN Rd_BlockStart(char* keyword);
  extern FILE_IO_BOOLEAN Rd_BlockEnd(char* keyword);

  extern FILE* Wr_OpenFile_Simple(char *fname);
  extern void Wr_CloseFile(void);
  extern void Wr_Space(void);
  extern void Wr_Comment(char *comment);
  extern void Wr_Int(int val, char *comment);
  extern void Wr_Double(double val, char *comment);
  extern void Wr_BlockStart(char* keyword);
  extern void Wr_BlockEnd(char* keyword);

#ifdef __cplusplus
} /* End of extern "C" { */
#endif

#endif
