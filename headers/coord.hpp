#ifndef COORD_HPP
#define COORD_HPP
#include<vector>
using namespace std;
class Coord {
public:
    int x;
    int y;
    Coord(int x, int y);
    Coord();

    vector<Coord> getNeighbors(int length, int width);
};

#endif