#ifndef GROVE_HPP
#define GROVE_HPP
#include "commodity.hpp"
#include "planningFunc.hpp"
#include <vector>
using namespace std;

class Behavior;

class Grove {
private:
    Commodity crop; // The crop that is being grown
    bool agency; // 1: This grove makes decisions, 0: The grove has a fixed behavior pattern
    int ibounds[2]; // Lower inclusive, Upper Exclusive
    int jbounds[2]; // Lower inclusive, Upper Exclusive
    double lambda; // Trust in extension agent
    double alpha; // Expectation of neighbors coordination
    double sprayEfficacy;

public:
    double fixedCosts; //Fixed cost per year associated with the grove
    double costs = 0;
    double returns = 0;
    double lastExtensionRisk = 0;
    double lastGrowerRisk = 0;
    double lastAdjustedRisk = 0;
    double lastRowPsyllids = 0;
    double lastColPsyllids = 0;
    double maxE_i = 0;
    double maxE_j = 0;
    double lastNAEV = 0;
    double lastISEV = 0;
    double lastGSEV = 0;
    bool foundHLB = false;
    int foundHLB_day = -1;
    double rogueThreshold = 0;
    int rogueFreq = 0;
    int rogueRadius = 0;
    vector<Behavior*> behaviorPatterns;

    Grove();
    Grove(Commodity crop, bool agency, int i_lb, int i_ub, int j_lb, int j_ub);
    //Grove();
    // Getters 
    bool hasAgency() { return this->agency; }
    Commodity* getCrop() { return &crop; }
    double getFixedCosts() { return this->fixedCosts; }
    int* getIBounds();
    int* getJBounds();
    //Setters
    void setAgency(bool);
    //Get a planned action
    plan_func getAction(int relativePeriod);
    double getLambda() { return this->lambda; }
    double getAlpha() { return this->alpha; }
    double getSprayEfficacy() { return this->sprayEfficacy; }
};

#endif