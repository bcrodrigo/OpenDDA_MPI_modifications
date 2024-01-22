/* Written by James Mc Donald 2006
   [Antispam: email in reverse] ei.yawlagiun.ti@semaj

   Function prototypes for efficiencies

   Copyright (C) 2006 James Mc Donald,
   Computational Astrophysics Laboratory,
   National University of Ireland, Galway
   This code is covered by the GNU General Public License 

   
   MODIFICATION:
   Efficiencies are written into a csv file, instead of the standard output
*/
#include <stdio.h>
#include <float.h>
#include "definitions.h"

void efficiencies(void);
FILE *q_output;