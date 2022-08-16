#include "../headers/commodity.hpp"
#include <iostream>
#include <assert.h>



double Commodity::getYieldAtAge(int i) {
    assert(!yieldProfile.empty());
    if (i >= yieldProfile.size()) {
        return yieldProfile.back();
    }
    else if (i > maxAge) {
        return 0;
    }
    else {
        return yieldProfile[i];
    }

}

Commodity::Commodity() {
   
}

// Constructor, calls setters
Commodity::Commodity(double price, int maxAge, std::vector<double> yieldProfile, double variableCost, std::vector<int> harvests) {
    setPrice(price);
    setMaxAge(maxAge);
    this->yieldProfile = yieldProfile;
    setVariableCost(variableCost);
    this->harvestPeriods = harvests;
}

Commodity::Commodity(double freshYield, double juiceYield, double freshPrice, double juicePrice, double costs, std::vector<int> harvests) {
    this->harvestPeriods = harvests;
    this->freshYield = freshYield;
    this->juiceYield = juiceYield;
    this->freshPrice = freshPrice;
    this->juicePrice = juicePrice;
    this->costs = costs;
}

double Commodity::getReturns() {
    return (freshYield * freshPrice) + (juiceYield * juicePrice);
}

//Sets MaxAge
void Commodity::setMaxAge(int age) {
    if (age > 0) {
        this->maxAge = age;
    }
    else {
        this->maxAge = 0;
    }
}

//Sets Price
void Commodity::setPrice(double price) {
    if (price > 0) {
        this->price = price;
    }
    else {
        this->price = 0;
    }
}

// Sets variable cost
void Commodity::setVariableCost(double cost) {
    if (cost > 0) {
        this->variableCost = cost;
    }
    else {
        this->variableCost = 0;
    }
}

