# Makefile for compiling the RR code on Linux systems

CC = gcc
CF = -O3 -Wall -Wextra -Wno-unused-parameter -Wuninitialized -Winit-self -pedantic -fopenmp -ffast-math -std=gnu11

INCLUDES = -I/data/home/mbreton/gsl-2.6/include/ -I/dec/users/mabreton/Cuba-4.2/
LIBRARIES = -L/data/home/mbreton/gsl-2.6/lib/ -lgsl -lgslcblas -lm -L/dec/users/mabreton/Cuba-4.2/ -lcuba

##### WARNING !!!!!! MUST WRITE THIS LINE IN TERMINAL BEFORE EXECUTING THE CODE IF GSL LIBRARY NOT IN PATH #############
#export LD_LIBRARY_PATH=/data/home/mbreton/gsl-2.6/lib:$LD_LIBRARY_PATH


#FINALIZE COMPILE FLAGS
CF += $(OPTIONS) #-g


## FINALIZE 
RR_INC = $(INCLUDES)
RR_LIB =  $(LIBRARIES) 



RR: compute_RR.c
	$(CC) compute_RR.c -o compute_RR $(CF) $(RR_INC) $(RR_LIB)


clean:
	rm -f compute_RR */*~ *~
	
