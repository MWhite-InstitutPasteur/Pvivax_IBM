CXX      := -c++
CXXFLAGS := -pedantic-errors -Wall -Wextra
LDFLAGS  := -L/usr/lib -lstdc++ -lm
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)
TARGET   := sim
INCLUDE  := -Iinclude/ -Icontrib/
HDR      := \
	$(wildcard contrib/*.h) \
	$(wildcard Pv_mod/*.h) \

SRC      :=                      \
	$(wildcard contrib/*.cpp) \
	$(wildcard Pv_mod/*.cpp) \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

# We make all sources depend on all headers. Not necessary, but ensures a
# correct build, and low-cost if using ccache.
# This is a hacky alternative to using a more complex build system.
$(OBJ_DIR)/%.o: %.cpp $(HDR)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
