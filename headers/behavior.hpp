#ifndef BEHAVIOR_HPP
#define BEHAVIOR_HPP
#include "grove.hpp"
#include "math.h"
#include "planningFunc.hpp"
#include "bioABM.h"
#include<vector>

using namespace std;

void sprayGrove(Grove * g);
double getInfectedYield(double severity);
class Behavior {
public:
    //Fill up the planned actions vector for the next season
    virtual void PlanActions() = 0;
    //Simulate outcomes of a grove under this behavior pattern
    double getExpectedValue(Grove g, double risk, int simulationLength, 
                            int startingPeriod, int planningLength,
                            double sprayEfficacy, double alpha, double additionalCosts);

    double* getExpectedValueTester(Grove g, double risk, int simulationLength, 
                            int startingPeriod, int planningLength,
                            double sprayEfficacy, double alpha);    
    //Get the expected mean infection at time t
    virtual double hlbSpread(int t, Grove *g) = 0;

    //Returns the variable costs of behavior per planning period
    virtual double getVariableCosts() = 0;

    virtual void executeAction(Grove *g) = 0;
    virtual string getName() = 0;
    virtual string getParams() = 0;

    bool actionPlannedOnDay(int day) {
        assert(day >= 0 && day < 365);
        return q[day];
    }

    //Action queue
    bool q[365] = { 0 };

};

//NOTE: PERMANENTLY INCREASES YIELD AND COSTS
class DensePlanting: public Behavior {
private:
    bool notModified = true;
public:
    double annualCosts; // the addition to annual costs based on modified density
    double yieldMultiplier; // multiplies yield by (1 + yieldMultiplier)

    DensePlanting(double annualCosts, double ym) {
        this->annualCosts = annualCosts;
        this->yieldMultiplier = ym;
    }
    void PlanActions();
    double getVariableCosts() { return annualCosts; }
    double hlbSpread(int t, Grove *g) { return 0; }
    void executeAction(Grove *g);
    string getName() { return "DensePlanting"; }
    string getParams() {
        stringstream ss;
        ss << annualCosts << ";" << yieldMultiplier;
        return ss.str();
    }
    
};

class RogueTrees: public Behavior {
public:
    //Cost per tree remove
    double removalCost;
    //Cost per survey
    double surveyCost;
    //Number of surveys per year
    int frequency; 
    //Radius of surrounding trees removed
    int radius;

    RogueTrees(double removalCost, double surveyCost, int frequency, int radius) {
        this->removalCost = removalCost;
        this->radius = radius;
        this->frequency = frequency;
        this->surveyCost = surveyCost;
    }

    //Check trees to be rogued based on frequency parameter
    void PlanActions();

    void executeAction(Grove *g);

    double hlbSpread(int t, Grove *g) { return 0; };

    //Unused and incorrect
    double getVariableCosts() { return this->removalCost; }
    string getName() { return "RogueTrees"; }
    string getParams() { 
        stringstream ss;
        ss << frequency << ";" << radius << ";" << removalCost << ";" << surveyCost;
        return ss.str();
    }

};

class RectangularRogue: public Behavior {
public:
    //Width dimension of rogue
    int width;
    //Height dimension of rogue
    int height;
    //Cost per tree remove
    double removalCost;
    //Cost per survey
    double surveyCost;
    //Number of surveys per year
    int frequency; 

    RectangularRogue(double removalCost, double surveyCost, int frequency, int width, int height) {
        this->removalCost = removalCost;
        this->width = width;
        this->height = height;
        this->frequency = frequency;
        this->surveyCost = surveyCost;
    }

    //Check trees to be rogued based on frequency parameter
    void PlanActions();

    void executeAction(Grove *g);

    double hlbSpread(int t, Grove *g) { return 0; };

    //Unused and incorrect
    double getVariableCosts() { return this->removalCost; }
    string getName() { return "RectangularRogue"; }
    string getParams() { 
        stringstream ss;
        ss << frequency << ";" << width << ";" << height << ";" << removalCost << ";" << surveyCost;
        return ss.str();
    }


};

class SprayTrees: public Behavior {
 private:
    //window around target dates
    int windowSize;
    double efficacy;
    //harvest start days
    int start1;
    int start2;
    int start3;
    //Cost per spray
    double sprayCost;
 public:
        SprayTrees(double efficacy, double sprayCost, int start1, int start2, int start3) {  
            this->efficacy = efficacy;
            this->sprayCost = sprayCost;
            this->start1 = start1;
            this->start2 = start2;
            this->start3 = start3;
            this->windowSize = 21; //Hard coded for 2nd manuscript
        }
        // //Fills planning Q with sprays in a window
        void PlanActions();

        //Returns expected value of this behavior pattern if continued until the end of the simulation
        double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod, double efficacy);

        //Returns the expected mean infection
        double hlbSpread(int t, Grove *g) { return 0; }

        //Returns the variable costs per year
        double getVariableCosts() { return this->sprayCost; }

        string getName() { return "SprayTrees"; }
        string getParams() { 
            stringstream ss;
            ss << efficacy << ";" << sprayCost;
            return ss.str();
        }

        void executeAction(Grove *g);
};


/*
class NoAction : public Behavior {
public:
    //Fills planning Q with NULLs
    void PlanActions();

    //Returns expected value of this behavior pattern if continued until the end of the simulation
    double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod);

    //Returns the expected mean infection with no mitigation
    double hlbSpread(int t, Grove *g);

    //Does nothing
    void executeAction(Grove *g) { return; }

    //Returns the variable costs per planning period
    double getVariableCosts() { return 0; }

    string getName() { return "NoAction"; }
    string getParams() { return "NA"; }
};*/

/*
class IndividualAction: public Behavior {
    private:
        int windowSize;
        double vc;
    public:
        IndividualAction(int window, double sprayCost, int planningLength) { 
            windowSize = window; 
            this->vc = sprayCost *6 / 4; 
        }
        // //Fills planning Q with sprays in a window
        void PlanActions(Grove *g);

        //Returns expected value of this behavior pattern if continued until the end of the simulation
        double SimulateOutcome(Grove * g, double risk, int simulationLength,int startingPeriod, double efficacy);

        //Returns the expected mean infection
        double hlbSpread(int t, Grove *g);

        //Returns the variable costs per year
        double getVariableCosts() { return this->vc; }
};*/


#endif
