CC = g++
NAME = libmap_tracking
VERSION = .so
TARGET_LIB = $(NAME)$(VERSION)

INCDIR = include
TPSA_INCDIR = tpsa_lib/include
CFLAGS = -O3 -Wall -shared -std=c++11 -fPIC -I$(INCDIR) -I$(TPSA_INCDIR)

LIBS = -lm

SRC = $(wildcard src/*.cc)
SRC += $(wildcard src/*.cpp)
OBJ_1 = $(SRC:.cpp=.o)
OBJ = $(OBJ_1:.cc=.o)
DEPS = $(wildcard $(INCDIR)/*.hpp)
DEPS += $(wildcard $(TPSA_INCDIR)/*.h)

$(info $$SRC is [${SRC}])
$(info $$OBJ is [${OBJ}])
$(info $$CFLAGS is [${CFLAGS}])

.PHONE: all
all = $(TARGET_LIB)

$(TARGET_LIB): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
	
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ)