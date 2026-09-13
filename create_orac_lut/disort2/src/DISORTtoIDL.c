/* File:    FORTRANtoIDL.c
 * Purpose: IDL Backend to some simple Fortran functions
 * Author:  Randall Skelton <skelton@atm.ox.ac.uk>
 * Date:    Spring-Summer 2001 */

/* ANSI */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* IDL */
#include "export.h"

/* Local */
#include "DISORTtoIDL.h"

/* Definitions for Windows/MacOS */
#ifdef WIN32
#include <windows.h>
    __declspec(dllexport)
#endif
#ifdef __MACOS__
    __declspec(export)
#endif

/* Define message codes and their corresponding printf(3) format strings.
 * Note that message codes start at zero and each one is one less than the
 * previous one Codes must be monotonic and contiguous.  Must match to
 * corrisponding entries in FORTRANtoIDL.h */

static IDL_MSG_DEF msg_arr[] =
	{
		{ "DISORTtoIDL_Error", "%NError: %s." },
		{ "DISORTtoIDL_NOSTRINGARRAY", "%Nstring arrays not allowed %s"},
	};

/* The load function fills this message block handle with the opaque handle to
 * the message block used for this model.  The other routines can then us it
 * to throw errors from this block. */

IDL_MSG_BLOCK msg_block;

int IDL_Load(void)
{
	/* Define the message block */
	if (!(msg_block = IDL_MessageDefineBlock("DISORTtoIDL", ARRLEN(msg_arr),
				msg_arr))) {
		return IDL_FALSE;
	}

	/* Call the startup function to add the routines to IDL. */

	/* Routines */
	if (!DISORTtoIDL_Startup()) {
		IDL_MessageFromBlock(msg_block, DISORTtoIDL_ERROR,
			IDL_MSG_RET, "Unable to initialize DISORT to IDL interface");
	}

	return IDL_TRUE;
}


