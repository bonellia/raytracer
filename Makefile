all:
	g++ *.cpp -std=c++17 -o raytracer -O3 -pthread -flto
