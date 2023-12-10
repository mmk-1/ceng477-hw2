CC = g++
src = *.cpp

all: clean rasterizer_cpp

rasterizer_cpp:
	$(CC) -Wfatal-errors $(src) -std=c++11 -O3 -o rasterizer

clean:
		rm -f rasterizer

run: clean all
	./rasterizer $$(cat scene.txt)

