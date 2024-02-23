# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3 -lmpfr

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

plot0.png: plt0.py data0.csv
	python -OO plt0.py

icat0: plot0.png
	kitten icat plot0.png

plot1.png: plt1.py data1.csv
	python -OO plt1.py

icat1: plot1.png
	kitten icat plot1.png

plot_two.png: plt_two.py
	python -OO plt_two.py

icat_two: plot_two.png
	kitten icat plot_two.png

$(FOO): $(HEADERS)
	touch $@

debug:
	@echo "OBJS: $(OBJS)"
	@echo "HEADERS: $(HEADERS)"

clean:
	rm -f $(TARGET) plot.png data.csv debug $(OBJS)

.PHONY: cal icat0 icat1 icat_two clean 
