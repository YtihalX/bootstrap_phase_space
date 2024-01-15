compile: main.cc
	g++ main.cc -o main -O3
test: test.cc
	g++ test.cc -o test
	./test

cal: compile
	./main

icat: plt.py
	python -OO plt.py
	kitten icat plot.png

clean:
	rm -f main plot.png data.csv debug

.PHONY: compile cal plot icat clean
