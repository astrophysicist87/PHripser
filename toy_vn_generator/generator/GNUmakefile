# =============================================================
#  Makefile                             Christopher J. Plumberg
# =============================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install	make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS= -O3 -std=c++11 -g -lgsl -lgslcblas -fopenmp

RM				=	rm -f
O               =	.o
LDFLAGS         =	$(CFLAGS)
SYSTEMFILES     =	$(SRCGNU)
INCDIR			=	include
SRCDIR			=	src
OBJDIR			=	obj

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	toy_vn_generator.e
endif

SRC			=	$(SRCDIR)/toy_vn_generator.cpp

INC			= 	$(INCDIR)/toy_vn_generator.h \
				$(INCDIR)/sampled_distribution.h

# -------------------------------------------------

OBJECTS				=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
						$(notdir $(basename $(SRC)))))
TARGET				=	$(MAIN)
INSTPATH			=	..

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	-@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(SRCDIR)/%.cc:
	#if [ -f $@ ] ; then touch $@ ; else false ; fi

$(TARGET):	$(OBJECTS)	
	$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 

# -------------------------------------------------

.PHONY:		all mkdirs clean distclean install

all:		mkdirs $(TARGET)

help:
	@grep '^##' GNUmakefile

mkdirs:	
	-@mkdir -p $(OBJDIR)

clean:		
	-rm -f $(OBJECTS)

distclean:	
	-rm -f $(TARGET)
	-rm -f $(OBJECTS)

install:	$(TARGET)
	cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
$(SRCDIR)/pair_particle_vn.cc: $(INCDIR)/toy_vn_generator.h \
                               $(INCDIR)/sampled_distribution.h


