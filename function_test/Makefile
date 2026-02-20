CXX      ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP

PKG_CONFIG   := $(shell command -v pkg-config 2>/dev/null)

ifeq ($(strip $(PKG_CONFIG)),)
  ARMA_CFLAGS :=
  ARMA_LIBS   :=
else
  ARMA_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags armadillo 2>/dev/null)
  ARMA_LIBS   := $(shell "$(PKG_CONFIG)" --libs armadillo 2>/dev/null)
endif

INCLUDES := -I/usr/local/include \
            -I/opt/homebrew/opt/boost/include \
            -I/opt/homebrew/include

ifneq ($(strip $(ARMA_CFLAGS)),)
  INCLUDES += $(ARMA_CFLAGS)
endif

LIBS := -larmadillo \
        -L/opt/homebrew/opt/openblas/lib -lopenblas \
        -llapack \
        -L/opt/homebrew/opt/superlu/lib -lsuperlu

ifneq ($(strip $(ARMA_LIBS)),)
  LIBS += $(ARMA_LIBS)
endif

# Source files for the test (exclude any main.cpp if present)
SOURCES := test_scoretest.cpp saige_test.cpp spa.cpp spa_binary.cpp UTIL.cpp cct.cpp getMem.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  := test_scoretest

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run

-include $(DEPS)
