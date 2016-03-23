#lib_LTLIBRARIES = libdfa-1.0.la libedsm-1.0.la liblstar-1.0utilities.la
#
#libdfa_1_0_la_SOURCES = dfa.cpp dfa.h utility.h utility.cpp
#libdfa_1_0_la_CPPFLAGS = -DDEBUG  -DDEBUG_1 -DDEBUG_2 -DENGLISH_MSG -ggdb3 -O0 -std=gnu++11
#
#libedsm_1_0_la_SOURCES = edsm.cpp edsm.h utility.h utility.cpp
#libedsm_1_0_la_CPPFLAGS = -DDEBUG  -DDEBUG_1 -DDEBUG_2 -DENGLISH_MSG -ggdb3 -O0 -std=gnu++11
#
#liblstar_1_0_la_SOURCES = lstar.cpp lstar.h utility.h utility.cpp
#liblstar_1_0_la_CPPFLAGS = -DDEBUG  -DDEBUG_1 -DDEBUG_2 -DENGLISH_MSG -ggdb3 -O0 -std=gnu++11
#
#include_HEADERS = dfa.h edsm.h lstar.h messages.h
#
#bin_PROGRAMS = giLearning
#
#giLearning_SOURCES = GI-learning.cpp 
#giLearning_LDADD = $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(lib_LTLIBRARIES) 
#
#
#debug
#giLearning_CPPFLAGS = $(BOOST_CPPFLAGS) -DDEBUG  -DDEBUG_1 -DDEBUG_2 -DENGLISH_MSG -ggdb3 -O0 -std=gnu++11
#giLearning_LDFLAGS = $(BOOST_SYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)

#############################################################

#include_HEADERS = dfa.h dfaEDSM.h bluefringe.h blueStar.h edsm.h lstar.h messages.h utilities.h
#
#bin_PROGRAMS = giLearning
#
#giLearning_SOURCES = GI-learning.cpp bluefringe.cpp blueStar.cpp edsm.cpp dfa.cpp dfaEDSM.cpp utilities.cpp lstar.cpp 
#giLearning_LDADD = $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(LIBXMLXX_LIBS)
#
#giLearning_CPPFLAGS = $(BOOST_CPPFLAGS) $(LIBXMLXX_CFLAGS) -DDEBUG  -DDEBUG_1 -DDEBUG_2 -DENGLISH_MSG -ggdb3 -O0 -std=gnu++11
#giLearning_LDFLAGS = $(BOOST_SYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)

#############################################################

lib_LTLIBRARIES = libutilities-1.0.la libdfa-1.0.la libedsm-1.0.la liblstar-1.0.la 

libutilities_1_0_la_SOURCES = utilities.h utilities.cpp
libutilities_1_0_la_CPPFLAGS = --std=gnu++11 -O2 -fopenmp

libdfa_1_0_la_SOURCES = dfa.cpp dfa.h
libdfa_1_0_la_CPPFLAGS = -std=gnu++11 -O2 -fopenmp
libdfa_1_0_la_LIBADD = libutilities-1.0.la


libedsm_1_0_la_SOURCES = dfaEDSM.h dfaEDSM.cpp bluefringe.h blueStar.h bluefringe.cpp blueStar.cpp edsm.cpp edsm.h
libedsm_1_0_la_CPPFLAGS = -std=gnu++11 -O2 -fopenmp
libedsm_1_0_la_LIBADD = libdfa-1.0.la

liblstar_1_0_la_SOURCES = lstar.cpp lstar.h
liblstar_1_0_la_CPPFLAGS = -std=gnu++11 -O2 -fopenmp


include_HEADERS = dfa.h dfaEDSM.h bluefringe.h blueStar.h edsm.h lstar.h messages.h utilities.h

bin_PROGRAMS = giLearning

giLearning_SOURCES = GI-learning.cpp 
giLearning_LDADD = $(BOOST_SYSTEM_LIBS) $(BOOST_FILESYSTEM_LIBS) $(LIBXMLXX_LIBS) $(lib_LTLIBRARIES) 

giLearning_CPPFLAGS = $(BOOST_CPPFLAGS) $(LIBXMLXX_CFLAGS) -std=gnu++11 -O2 -fopenmp
giLearning_LDFLAGS = $(BOOST_SYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)