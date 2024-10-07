/**********************************************************************
*   This software is part of the grant PID2021-122955OB-C41
*   and by 'ERDF A way of making Europe'.
*
***********************************************************************
*                      TESRECONS
*
*  File:       tesrecons.h
*  Developers: Beatriz Cobo
* 	           cobo@ifca.unican.es
*              IFCA
*              Maite Ceballos
*              ceballos@ifca.unican.es
*              IFCA
*                                                                     
***********************************************************************/

#ifndef TESRECONS_H
#define TELIB_H 1

#include "initSIRENA.h"

#include "gti.h"

#define TOOLSUB tesrecons_main
#include "headas_main.c"

#include "versionSIRENA.h"

#include <time.h>

int getpar_tesrecons(struct Parameters* const par);

#endif /* TESRECONS_H */
