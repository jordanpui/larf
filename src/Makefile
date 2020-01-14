SHELL=/bin/bash

OPT := -O3
# OPT     := -O3 -DWRITE -DDRAW
# OPT     := -O3 -DNDEBUG
LIBS := -pthread
ifeq ($(mode),debug)
	OPT := -O0 -g -DWRITE -DDRAW
	LIBS += -lgd
else ifeq ($(mode),light_debug)
	OPT := -DWRITE -DDRAW
	LIBS += -lgd
else ifeq ($(mode),release)
	OPT += -static
	LIBS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif

target = bin

WFLAG 	= -Wall
DEPFLAG = -MMD -MP

GUROBIFLAGS=-L$(GUROBI_HOME)/lib -lgurobi_c++ -lgurobi81
GUROBIDIR=$(GUROBI_HOME)/include

# CC = g++ -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 $(OPT) $(WFLAG) $(DEPFLAG) $(INCLUDE)
CC = g++-4.8 -std=c++11 $(OPT) $(WFLAG) $(DEPFLAG) $(INCLUDE)

INCLUDE = -I. -I../../boost_1_62_0 -I$(GUROBIDIR)
LIBS += $(addprefix alg/, patoh/libpatoh.a) \
		-lboost_system -ldl $(CCLNFLAGS) $(GUROBIFLAGS)

CC_OBJS = main
UT_OBJS = $(addprefix utils/, log draw)
DB_OBJS = $(addprefix db/, db db_draw db_bookshelf site instance net group swbox clkrgn)
GP_OBJS = $(addprefix gp/, gp gp_data gp_main gp_qsolve gp_spread gp_region gp_setting)
TDM_OBJS = $(addprefix tdm/, timing_graph tdm_db tdm_part tdm_net tdm_solve_lp tdm_solve_lag tdm_solve_lag_init tdm_solve_lag_update tdm_solve_lag_data tdm_leg tdm_refine_lp tdm_refine_greedy)
ALG_OBJS = $(addprefix alg/, matching bipartite)

OBJS = 	$(addsuffix .o, $(CC_OBJS) $(UT_OBJS) $(DB_OBJS) $(ALG_OBJS) $(TDM_OBJS) $(GP_OBJS))  

BFILE = larf

TAR_FILES = $(BFILE) ../scripts/*

.PHONY: all clean tags
all: $(BFILE)
	@date +'%D %T'
	mkdir -p ../$(target)
	cp -u $(TAR_FILES) ../$(target)/

$(BFILE): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LIBS)

%.o : %.cpp
	$(CC) -o $@ -c $<

-include $(OBJS:.o=.d)

clean:
	rm -f {*,*/*,*/*/*}.{o,d} $(BFILE)

tags:
	cscope -Rbq
	ctags -R *.cpp *.h
