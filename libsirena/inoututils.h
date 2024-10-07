/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      INOUTUTILS
*
*  File:       inoututils.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef INOUTUTILS_
#define INOUTUTILS_

// Utils module

#include "genutils.h"
	
using namespace std;

struct IOData
{
	fitsfile *inObject;
	char *nameTable;
	char *nameCol;
	char *unit;
	//MC char *type;
	int type;
	int iniCol;
	int endCol;
	long iniRow;
	long endRow;
};

// Structure to define input parameters
typedef struct
{
	string name;
	string description;
	string type;
	string ValStr;
	string defValStr;
	int ValInt;
	int minValInt;
	int maxValInt;
	int defValInt;
	double ValReal;
	double minValReal;
	double maxValReal;
	double defValReal;
} inparam;

#ifdef __cplusplus
extern "C"
#endif
int readFitsSimple(IOData obj, gsl_vector **result);
int readFitsComplex(IOData obj, gsl_matrix **result);

int writeFitsSimple (IOData obj,gsl_vector *vector);
int writeFitsComplex(IOData obj, gsl_matrix *matrix);

int toGslMatrix(void **buffer, gsl_matrix **matrix, long numCol,int numRow,int type, int eventini);
int toGslVector(void **buffer, gsl_vector **array, long nevent, int eventini, int type);

int fromGslVector(void **buffer, gsl_vector **array, int type);
int fromGslMatrix(void **buffer, gsl_matrix **matrix, int type);

int interactivePars(inparam *taskPars, int np, string task);

#endif /*INOUTUTILS_*/
