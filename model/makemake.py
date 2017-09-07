
import discretize
prob = discretize.Problem
nF = discretize.nF              # include myfunctions.cpp if functions are used

template = """
# Makefile for %(probl)s IPOPT problem

# This Makefile depends on pkg-config correctly locating the IPOPT headers and
# libraries; if your installation doesn't support this, alter COIN_CFLAGS and
# COIN_LDFLAGS to point to the correct directories
#
# Compilation of large problems may require a fairly recent (>= 4.3) version of
# g++

CFLAGS = -O3 -pipe -DNDEBUG -pedantic-errors

COIN_CFLAGS = $(shell pkg-config --cflags ipopt)
COIN_LDFLAGS = $(shell pkg-config --libs ipopt)

SRCS = %(probl)s_main.cpp %(probl)s_nlp.cpp %(eqns)s
OBJS = $(SRCS:.cpp=.o)

all: %(probl)s

%%.o : %%.cpp
        $(CXX) $(CFLAGS) $(COIN_CFLAGS) -c $<

%(probl)s: $(OBJS)
        $(CXX) $(OBJS) $(COIN_LDFLAGS) -o %(probl)s
"""

open("Makefile","wt").write(template % dict(probl = prob.lower(),
                                            eqns = "myfunctions.cpp" if nF else ""))
