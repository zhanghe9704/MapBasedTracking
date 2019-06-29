INCDIR = ../include
CFLAGS = -O3 -Wall -std=c++11

TPSA_LIBS = tpsaso
TRACKING_LIBS = map_tracking
LIBS = -L. -l$(TPSA_LIBS) -l$(TRACKING_LIBS) -lstdc++ -lm

SRC = $(wildcard *.cc)
TARGET = $(SRC:.cc=)
DEPS = $(INCDIR)/map.hpp $(INCDIR)/tracking.hpp

.PHONY: all
all: $(TARGET)

.cc:
	$(CC) -o $@ $< $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TARGET)