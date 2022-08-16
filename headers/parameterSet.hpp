#ifndef PARAMETERSET_HPP
#define PARAMETERSET_HPP
#include <string>

namespace ParameterSet {
    int planningLength = 13; //Length of a planned period
    double sprayingPopEff = 0.10; //% chance of killing
    double freshYield = 736; // fresh yield per acre / num cells
    double juiceYield = 636; // juice yield per acre / num cells
    double freshPrice = 14.0; // fresh price
    double juicePrice = 2.70; // juice price
    double costs = 200; //costs per acre at given yields
    bool biologicalRun = false; // run only the bio model
    int numIndividualSprays = 3; //Number of sprays for individual behavior pattern
    int numGroupSprays = 6; //Number of sprays for group behavior pattern
    int projectionLength = 1825; //Number of days growers plan into future when choosing a behavior, longer than model duration
    const int gridLength = 1; //Number of rows in grove grid
    const int gridWidth = 1; //Length of rows in grove grid
    const int numBehaviorPatterns = 3; //Number of possible behavior patterns, used to index list
    double sprayCost = 3; // cost per spray
    int groupWindow = 21;
    int individualWindow = 60;
    int rogueRadius = 1; //Radius of adjacent trees removed when infected one found
    int rogueFrequency = 20; // Frequency that rogueing is checked
    double groupPremium = 5;
    //Trust in extension agent
    double g00_lambda = 0.5;
    double g01_lambda = 0.5;
    double g02_lambda = 0.5;
    double g10_lambda = 0.5;
    double g11_lambda = 0.5;
    double g12_lambda = 0.5;
    double g20_lambda = 0.5;
    double g21_lambda = 0.5;
    double g22_lambda = 0.5;
    //Expectations of neighbors cooperation
    double g00_alpha = 0.5;
    double g01_alpha = 0.5;
    double g02_alpha = 0.5;
    double g10_alpha = 0.5;
    double g11_alpha = 0.5;
    double g12_alpha = 0.5;
    double g20_alpha = 0.5;
    double g21_alpha = 0.5;
    double g22_alpha = 0.5;
    //Initial behavior, indexes behavior list
    string g00_behavior = "0";
    string g01_behavior = "0";
    string g02_behavior = "0";
    string g10_behavior = "0";
    string g11_behavior = "0";
    string g12_behavior = "0";
    string g20_behavior = "0";
    string g21_behavior = "0";
    string g22_behavior = "0";
    //Agency
    bool g00_agency = true;
    bool g01_agency = true;
    bool g02_agency = true;
    bool g10_agency = true;
    bool g11_agency = true;
    bool g12_agency = true;
    bool g20_agency = true;
    bool g21_agency = true;
    bool g22_agency = true;
    //Group OC/premium
    double g00_premium = 0.5;
    double g01_premium = 0.5;
    double g02_premium = 0.5;
    double g10_premium = 0.5;
    double g11_premium = 0.5;
    double g12_premium = 0.5;
    double g20_premium = 0.5;
    double g21_premium = 0.5;
    double g22_premium = 0.5;
    //Rogue threshold
    double g00_threshold = 0.5;
    double g01_threshold = 0.5;
    double g02_threshold = 0.5;
    double g10_threshold = 0.5;
    double g11_threshold = 0.5;
    double g12_threshold = 0.5;
    double g20_threshold = 0.5;
    double g21_threshold = 0.5;
    double g22_threshold = 0.5;
    //Rogueing frequency
    int g00_removalFreq = 0;
    int g01_removalFreq = 0;
    int g02_removalFreq = 0;
    int g10_removalFreq = 0;
    int g11_removalFreq = 0;
    int g12_removalFreq = 0;
    int g20_removalFreq = 0;
    int g21_removalFreq = 0;
    int g22_removalFreq = 0;

    
};

    
#endif