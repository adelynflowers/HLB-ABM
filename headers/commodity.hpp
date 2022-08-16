#ifndef COMMODITY_HPP
#define COMMODITY_HPP
#include <vector>
#include <algorithm>
class Commodity;



class Commodity {
    
private:
    double price;
    int maxAge;
    std::vector<double> yieldProfile;
    double variableCost; // per year
    std::vector<int> harvestPeriods;
   
    
    //Setters; Private because they should ONLY be accessed via constructor
    void setPrice(double);
    void setMaxAge(int);
    void setVariableCost(double);
public:
    double freshYield;
    double juiceYield;
    double freshPrice;
    double juicePrice;
    double costs; // per acre

    //Constructors
    Commodity();
    Commodity(double price, int maxAge, std::vector<double> yieldProfile, double variableCost, std::vector<int> harvests);
    Commodity(double freshYield, double juiceYield, double freshPrice, double juicePrice, double costs, std::vector<int> harvests);

    // Getters
    double getPrice() {return this->price; }
    int getMaxAge() { return this->maxAge; }
    double getYieldAtAge(int i);
    double getVariableCost() { return this->variableCost; }
    bool isHarvestPeriod(int t) { return std::find(harvestPeriods.begin(), harvestPeriods.end(), t) != harvestPeriods.end(); }
    double getReturns();

    

};
#endif