/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      TESLIB
*
*  File:       teslib.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef TESLIB_H
#define TELIB_H 1

#include "initSIRENA.h"

#include "gti.h"

#define TOOLSUB teslib_main
#include "headas_main.c"

#include "versionSIRENA.h"

#include <time.h>

int getpar_teslib(struct Parameters* const par);

#endif /* TESLIB_H */
