#include "../headers/grove.hpp"

/*
Grove::Grove() {
    //Please dont use this
    setAge(0);
}*/



Grove::Grove() {
}

Grove::Grove(Commodity crop, bool agency, int i_lb, int i_ub, int j_lb, int j_ub) {
    this->crop = crop;
    setAgency(agency);
    ibounds[0] = i_lb;
    ibounds[1] = i_ub;
    jbounds[0] = j_lb;
    jbounds[1] = j_ub;
}

int* Grove::getIBounds() {
    return this->ibounds;
}

int* Grove::getJBounds() {
    return this->jbounds;
}

//set agency
void Grove::setAgency(bool agency) {
    this->agency = agency;
}

