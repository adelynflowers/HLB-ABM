#include "../headers/coord.hpp"


//Constructors

Coord::Coord() {
    this->x = -1;
    this->y = -1;
}

Coord::Coord(int x, int y) {
    this->x = x;
    this->y = y;
}

vector<Coord> Coord::getNeighbors(int length, int width) {
    Coord neighbors[4];
    neighbors[0] = Coord(this->x+1,this->y);
    neighbors[1] = Coord(this->x,this->y+1);
    neighbors[2] = Coord(this->x-1,this->y);
    neighbors[3] = Coord(this->x,this->y-1);
    vector<Coord> cleanNeighbors;
    for (int i = 0; i < 4; i++) {
        if (neighbors[i].x < length && neighbors[i].x >= 0 && neighbors[i].y < width && neighbors[i].y >= 0) {
            cleanNeighbors.push_back(neighbors[i]);
        }
    }
    return cleanNeighbors;
}