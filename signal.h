/**
* @file signal_h
* The control structure used to pass data to the SPU from PPU side
* @author Alexandru Paler and Fabio Campos
*/
#ifndef _SIGNAL_H_
#define _SIGNAL_H_

#define NUM_THREADS 6

#define EXIT_MESSAGE 123

typedef union
{
	unsigned long long ull;
	unsigned int ui[2];
}addr64;

typedef struct _control_block
{
	unsigned int rank; /*4bytes*/
	unsigned long speid;/*8bytes*/
	void* ls;/*1byte*/
	void* sig1;/*1byte*/
	unsigned int mpm_num[4];/*16bytes*/
	unsigned int r;/*4bytes*/
	unsigned int from;/*4bytes*/
	unsigned int until;/*4bytes*/

	unsigned char  pad[84];
}control_block;/*128bytes*/

#endif  /*_SIGNAL_H_*/
