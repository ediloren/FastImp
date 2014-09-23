/*-*-c-*-*/
/*
  ==========================================================================
  ============================= File Info ==================================

  =========================== File Contents ================================

  Author: Zhenhai Zhu
  
  Description: 

  (1) Writing files with hierarchical indentation format. 
      This module could handle multiple files.
  (2) Reading files with hierarchical indentation format. 
      This module could handle multiple files.

  Resources:

  See also:

  const static char cvsid[] = "$Id: fileIo.c,v 1.3 2002/07/18 15:05:44 zhzhu Exp $";

  ==========================================================================
*/


#include "fileIo.h"

typedef struct _FileIoDescs {
    char* fname;
    int lineNum;
    FILE *fd;
    int level;			/* Indentation level */
} FileIoDescs;

static FileIoDescs  WrFileDesc;
static FileIoDescs  RdFileDesc;
static FILE_IO_BOOLEAN currentWrFileDescInUse = FILE_IO_FALSE;
static FILE_IO_BOOLEAN currentRdFileDescInUse = FILE_IO_FALSE;

#define MAX_STACK_SIZE 20
typedef FileIoDescs StackElement;
static StackElement  WrFileDescStack[MAX_STACK_SIZE];
static StackElement  RdFileDescStack[MAX_STACK_SIZE];
static int topOfWrFileDescsStack = -1;
static int topOfRdFileDescsStack = -1;

static void
addStack (
	  int *top,
	  StackElement *stack,
	  StackElement element);
static StackElement
deleteStack (
	     int *top,
	     StackElement *stack);
static void
stackFull (void);
static void
stackEmpty (void);
static int
IsStackEmpty (
	      int top);
static int
IsStackFull (
	     int top);

static char* Rd_NextLine(void);

/**********************************************************************
 * Rd_OpenFile --
 **********************************************************************/
FILE
*Rd_OpenFile(
	     char *fname)
{
    if (fname == NULL) return NULL;

    if (currentRdFileDescInUse == FILE_IO_TRUE) {
	if (IsStackFull(topOfRdFileDescsStack)) {
	    stackFull();
	} else {
	    addStack(&topOfRdFileDescsStack, RdFileDescStack, RdFileDesc);
	}
    }

    /* setup current file description */
    RdFileDesc.level   = 0;
    RdFileDesc.fname   = fname;
    RdFileDesc.fd = fopen(fname, "r");
    RdFileDesc.lineNum = 0;       
    currentRdFileDescInUse = FILE_IO_TRUE;

    if ( (RdFileDesc.fd == NULL) ||
	 (!Rd_BlockStart("File")) ) {
	printf("\n\tFail to open File [%s] for read\n", fname);
	return NULL;
    } else {
	printf("\n\tFile [%s] has been opened for read\n", fname);
    }

    return RdFileDesc.fd;
}

/**********************************************************************
 * Rd_CloseFile --
 **********************************************************************/
void 
Rd_CloseFile(void)
{
    if ( (WrFileDesc.level != 0) ||
	 (!Rd_BlockEnd("File")) ) {
	fprintf(stderr, "\n\t Error in fileIo"
		"\n\t Opened Block in the file!  Closing the file anyway.\n");
    }
    fclose(RdFileDesc.fd);
    printf("\n\tFile [%s] has been closed for read\n", RdFileDesc.fname);

    if (!IsStackEmpty(topOfRdFileDescsStack)) {
	RdFileDesc = deleteStack(&topOfRdFileDescsStack, RdFileDescStack);

    } else {
	/* else reset file description */
	RdFileDesc.fname    = NULL;
	RdFileDesc.fd       = NULL;
	RdFileDesc.lineNum  =    0;
	RdFileDesc.level    =    0;
	currentRdFileDescInUse = FILE_IO_FALSE;
    }          
}                

/**********************************************************************
 * Rd_Int --
 * If these is a comment behind the integer, it will be ignored
 **********************************************************************/
FILE_IO_BOOLEAN   
Rd_Int(
       int *val)
{     
    char *p;

    p = Rd_NextLine();
    if (p == NULL) return FILE_IO_FALSE;
    sscanf(p, "%d", val);
    return FILE_IO_TRUE;
}

/**********************************************************************
 * Rd_Double --
 * If these is a comment behind the double, it will be ignored
 **********************************************************************/
FILE_IO_BOOLEAN   
Rd_Double(
	  double *val)
{  
    char *p;
   
    p = Rd_NextLine();
    if (p == NULL) return FILE_IO_FALSE;
    sscanf(p, "%lf", val);
    return FILE_IO_TRUE;
}

/**********************************************************************
 * Rd_BlockStart --
 * keyword could be NULL, blank or any comment
 **********************************************************************/
FILE_IO_BOOLEAN   
Rd_BlockStart(
	      char* keyword)
{  
    char bfr[300];
    char       *p;
    char blank[3] = "";

    p = Rd_NextLine();
    if (p == NULL) return FILE_IO_FALSE;

    if ( (keyword != NULL) && (strcmp(keyword, blank) != 0) ){
	sprintf(bfr, "{%s", keyword);
	if (strcmp(bfr, p) != 0) {
	    fprintf(stderr, "\n\t Error in fileIo:"
		    "string Compare failed p=%s, buffer=%s in %s line %d\n",
		    p, bfr, __FILE__,__LINE__);
	    exit(1); 
	    return FILE_IO_FALSE; 
	}   
    }
    RdFileDesc.level++;   
    return FILE_IO_TRUE;
}
      
/**********************************************************************
 * Rd_BlockEnd --
 * keyword could be NULL, blank or any comment
 **********************************************************************/
FILE_IO_BOOLEAN   
Rd_BlockEnd(
	    char* keyword)
{  
    char bfr[300];
    char       *p;
    char blank[3] = "";

    p = Rd_NextLine();
    if (p == NULL) return FILE_IO_FALSE;

    if ( (keyword != NULL) && (strcmp(keyword, blank) != 0) ){
	sprintf(bfr, "}%s", keyword);
	if (strcmp(bfr, p) != 0) {
	    fprintf(stderr, "\n\t Error in fileIo:"
		    "string Compare failed p=%s, buffer=%s in %s line %d\n",
		    p, bfr, __FILE__, __LINE__);
	    exit(1); 
	    return FILE_IO_FALSE;
	}
    }

    RdFileDesc.level--;   
    return FILE_IO_TRUE;
} 

/**********************************************************************
 * Rd_NextLine --
 **********************************************************************/
static char *
Rd_NextLine(void) 
{

    static char  lineBfr[300];
    char                   *p;
    char                  *p1;

    do {
	/* Read line from file, return NULL if eof */
	p = fgets(lineBfr, sizeof(lineBfr)-1, RdFileDesc.fd);
	if (p == NULL) return NULL;
	assert(p == lineBfr);

	RdFileDesc.lineNum++;

	/* Remove EOL and comments from the end of line buffer */
	for(p1=p;
	    (*p1)!= '\0' && (*p1) != '#' && (*p1) != '*' && (*p1) != '\n'; 
	    p1++) {
	}
	*p1 = '\0';

	/* Skip leading spaces */
	while((*p) && isspace(*p)) {
	    p++;
	}

	/* Remove trailing spaces */
	p1 = p + strlen(p);
	while (p1 > p && isspace(p1[-1])) {
	    p1--;
	}
	*p1 = '\0';

    } while (strlen(p) == 0);

    return (p);
}
	         
/**********************************************************************
 * Wr_OpenFile --
 **********************************************************************/
FILE
*Wr_OpenFile(
	     char *fname)
{
    char charBuffer[200];
    time_t timer; 

    if (fname == NULL) return NULL;
   
    if (currentWrFileDescInUse == FILE_IO_TRUE) {
	if (IsStackFull(topOfWrFileDescsStack)) {
	    stackFull();
	} else {
	    addStack(&topOfWrFileDescsStack, WrFileDescStack, WrFileDesc);
	}
    }

    /* setup current file description */
    WrFileDesc.level   = 0;
    WrFileDesc.fname   = fname;
    WrFileDesc.fd      = fopen(fname, "w");
    WrFileDesc.lineNum = 0;       
    currentWrFileDescInUse = FILE_IO_TRUE;

    /* proceed to open the file */
    if (WrFileDesc.fd == NULL) return NULL;
    printf("\n\tFile [%s] has been opened for write\n", fname);

    sprintf(charBuffer, "%s %s", "This file was generated by ",getenv("LOGNAME"));
    Wr_Comment(charBuffer);

    (void)time(&timer);
    sprintf(charBuffer, "%s %s", "date: ", ctime(&timer));
    Wr_Comment(charBuffer);

    Wr_BlockStart("File");
       
    return WrFileDesc.fd;
}                


/**********************************************************************
 * Wr_CloseFile --
 **********************************************************************/
void 
Wr_CloseFile(
	     void)
{
    Wr_BlockEnd("File");
    Wr_Space();
       
    Wr_Comment("--- End of file");
       
    if (WrFileDesc.level != 0) {
	fprintf(stderr, "\n\t Error in fileIo"
		"\n\t Opened Block in the file!  Closing the file anyway.\n");
    }
    fclose(WrFileDesc.fd);
    printf("\n\tFile [%s] has been closed for write\n", WrFileDesc.fname);

    if (!IsStackEmpty(topOfWrFileDescsStack)) {
	WrFileDesc = deleteStack(&topOfWrFileDescsStack, WrFileDescStack);

    } else {
	/* else reset file description */
	WrFileDesc.fname    = NULL;
	WrFileDesc.fd       = NULL;
	WrFileDesc.lineNum  =    0;
	WrFileDesc.level    =    0;
	currentWrFileDescInUse = FILE_IO_FALSE;
    }          
}      

/**********************************************************************
 * Wr_RawLine --
 **********************************************************************/
void
Wr_RawLine(
	   char *line)
{
    fprintf(WrFileDesc.fd, "%s\n", line);
    WrFileDesc.lineNum ++;
}

/**********************************************************************
 * Wr_IndentedLine --
 **********************************************************************/
void
Wr_IndentedLine(
		char *line)
{
    char    bfr[300];
    char indent[100];

    memset(indent, 0x20, sizeof(indent));
    indent[3*WrFileDesc.level] = 0;
    sprintf(bfr, "%s%s", indent, line);
    Wr_RawLine(bfr);
}

/**********************************************************************
 * Wr_Comment --
 **********************************************************************/
void
Wr_Comment(
	   char * comment)
{
    char bfr[300];

    if (!comment) return;
    sprintf(bfr, "   #%s", comment);
    Wr_IndentedLine(bfr);
}

/**********************************************************************
 * Wr_Int --
 **********************************************************************/
void
Wr_Int(
       int val,
       char * comment)
{
    char bfr[300];
    sprintf(bfr, "%-9d         %s", val, comment);
    Wr_IndentedLine(bfr);
}

/**********************************************************************
 * Wr_Double --
 **********************************************************************/
void
Wr_Double(
	  double val,
	  char * comment)
{
    char bfr[300];
    sprintf(bfr, "%12.6f        %s", val, comment);
    Wr_IndentedLine(bfr);
}

/**********************************************************************
 * Wr_Space --
 **********************************************************************/
void
Wr_Space(
	 void)
{
    Wr_RawLine("");
}

/**********************************************************************
 * Wr_BlockStart --
 **********************************************************************/
void
Wr_BlockStart(
	      char* keyword)
{
    char bfr[300];

    sprintf(bfr, "{%s", keyword);
    Wr_IndentedLine(bfr);
    WrFileDesc.level++;
}

/**********************************************************************
 * Wr_BlockEnd --
 **********************************************************************/
void 
Wr_BlockEnd (
	     char* keyword)
{
    char bfr[300];
   
    sprintf(bfr, "}%s", keyword);
    WrFileDesc.level--;
    Wr_IndentedLine(bfr);
}

/*
 * Stack Abstrct Data Structure
 */

/**********************************************************************
 * addStack --
 **********************************************************************/
static void
addStack (
	  int *top,
	  StackElement *stack,
	  StackElement element)
{
    if (IsStackFull(*top)) {
	stackFull();
	return;
    }
    ++*top;
    stack[*top] = element;
}

/**********************************************************************
 * deleteStack --
 **********************************************************************/
static StackElement
deleteStack (
	     int *top,
	     StackElement *stack)
{
    if(IsStackEmpty(*top)) {
	stackEmpty();
    }
    return (stack[(*top)--]);
}

/**********************************************************************
 * IsStackEmpty --
 **********************************************************************/
static int
IsStackEmpty (
	      int top)
{
    if (top == -1) {
	return 1;
    }
    return 0;
}

/**********************************************************************
 * IsStackFull --
 **********************************************************************/
static int
IsStackFull (
	     int top)
{
    if (top >= MAX_STACK_SIZE) {
	return 1;
    }
    return 0;
}

/**********************************************************************
 * stackFull --
 **********************************************************************/
static void
stackFull (void)
{
    fprintf(stderr, "\n\t Error in fileIo:"
	    "\t\n number of files to be written or read exceeds the internal limitation\n");
    exit(1);    
}

/**********************************************************************
 * stackEmpty --
 **********************************************************************/
static void
stackEmpty (void)
{
    fprintf(stderr, "\n\t Error in fileIo:"
	    "\t\n No more elements in the stack\n");
    exit(1);        
}


