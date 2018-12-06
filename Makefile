# @(#)06	1.6  src/samples/tutorial/simple/Makefile, sw.samples, sdk_pub 10/11/05 15:26:25 
# -------------------------------------------------------------- 
# (C) Copyright 2001,2005,                                       
# International Business Machines Corporation,                   
# Sony Computer Entertainment Incorporated,                      
# Toshiba Corporation.                                           
#                                                                
# All Rights Reserved.                                           
# -------------------------------------------------------------- 
# PROLOG END TAG zYx                                              

########################################################################
#			Subdirectories
########################################################################

DIRS		:= 	spu

########################################################################
#                       Target
########################################################################

PROGRAM_ppu	:= 	aks


########################################################################
#                       Local Defines
########################################################################
SYS_LIBS        += -lm -lgmp

IMPORTS         := spu/lib_spu_polydiv.a -lspe2 -lpthread -lm

#INSTALL_DIR	= $(SDKBIN)/samples
#INSTALL_FILES	= $(PROGRAM_ppu)

########################################################################
#			make.footer
########################################################################

include $(CELL_TOP)/make.footer

