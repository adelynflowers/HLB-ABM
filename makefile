abmmake: src/*
	g++ -std=c++17 -Ofast -I headers -o serial.out src/*
abmdebug:
	g++ -std=c++17 -g -I headers -o serial.out src/*
