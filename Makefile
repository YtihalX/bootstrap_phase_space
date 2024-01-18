# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3

FOO = foo

# Source files
SRCS = main.cc routines.cc util.cc models.cc

HEADERS = routines.h util.h models.h

# Object files
OBJS = $(SRCS:.cc=.o)

# Executable name
TARGET = bootstrap

# Compile rule
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cc foo
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: test.cc
	g++ test.cc -o test
	./test

cal: $(TARGET)
	./$(TARGET)

icat: plt.py
	python -OO plt.py
	kitten icat plot.png

$(FOO): $(HEADERS)
	touch $@

debug:
	@echo "OBJS: $(OBJS)"
	@echo "HEADERS: $(HEADERS)"

clean:
	rm -f $(TARGET) plot.png data.csv debug $(OBJS)

.PHONY: compile cal plot icat clean
