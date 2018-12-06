# --------------------------------------------------------------  
# (C)Copyright 2001,2006,                                         
# International Business Machines Corporation,                    
# Sony Computer Entertainment, Incorporated,                      
# Toshiba Corporation,                                            
#                                                                 
# All Rights Reserved.                                            
# --------------------------------------------------------------  
# PROLOG END TAG zYx                                              

########################################################################
#                       Target
########################################################################

PROGRAMS_spu    = spu_polydiv

OBJS_spu_interrupt 	:= polyUtil.o

LIBRARY_embed	:= lib_spu_polydiv.a

########################################################################
#			Local Defines
########################################################################

IMPORTS		= $(SDKLIB_spu)/libmpm.a -lm
INCLUDE		:= -I ..

########################################################################
#			make.footer
########################################################################

include $(CELL_TOP)/make.footer