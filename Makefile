# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3

# Source files
SRCS = main.cc routines.cc util.cc models.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Executable name
TARGET = bootstrap

# Compile rule
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: test.cc
	g++ test.cc -o test
	./test

cal: $(TARGET)
	./main

icat: plt.py
	python -OO plt.py
	kitten icat plot.png

debug:
	@echo "OBJS: $(OBJS)"

clean:
	rm -f bootstrap plot.png data.csv debug

.PHONY: compile cal plot icat clean
