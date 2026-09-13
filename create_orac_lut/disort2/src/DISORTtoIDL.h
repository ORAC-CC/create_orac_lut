/* These must match the corresponding entries in DISORTwrapper.c 	*/

#ifndef IDL_C_H
#define IDL_C_H

/* Message Numbers */
#define DISORTtoIDL_ERROR		 0
#define DISORTtoIDL_NOSTRINGARRAY	-1

/* Useful macro */
#define ARRLEN(arr) (sizeof(arr)/sizeof(arr[0]))

extern IDL_MSG_BLOCK msg_block;

/* Define the startup function that adds C functions to IDL
 * along with the exit handler */

/* IDL-Postres interface */
extern void DISORTtoIDL_exit_handler(void);
extern int  DISORTtoIDL_Startup(void);

#endif
