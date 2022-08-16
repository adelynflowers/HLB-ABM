#include "../headers/behavior.hpp"




/***************************************************************
* Get Infected yield
* Relates HLB severity to yield. From Bassanezi et. al (2011)
******************************************************************/
double getInfectedYield(double hlbSeverity) {
    return exp(-1.85 * hlbSeverity);
}

/***************************************************************
* Spray Grove
***************************************************************/
void sprayGrove(int* ibounds, int* jbounds, double efficacy) {
    vector<boost::tuple<int, int>> coords;
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            coords.push_back(boost::tuple<int, int>(i, j));
        }
    }
    bioABM::sprayTrees(efficacy, coords);
}

/****************************************************************
 * Check and rogue
 * Checks all trees based on a frequency parameter, rogueing 
 * (w/o replacement) symptomatic trees. Returns number of trees
 * removed
 * **************************************************************/
int checkAndRogue(int* ibounds, int* jbounds, int width, int height) {
    int removalCount = 0;
    boost::uniform_int<> gen(0,1);
    boost::random::mt19937 econ_rng(std::time(0));
    for (int i = ibounds[0]; i < ibounds[1]; i++) {
        for (int j = jbounds[0]; j < jbounds[1]; j++) {
            if (bioABM::isSymptomatic(i,j) && gen(econ_rng) <= bioABM::getSeverityAt(i,j)) {
                //Rogue within a certain radius
                for (int k = -height; k <= height; k++) {
                    for (int l = -width; l <= width; l++) {
                        if (bioABM::isTreeAlive(i+k,j+l)) {
                            bioABM::rogueTreeAt(i+k,j+l);
                            removalCount++;
                        }
                    }
                }
            }
        }
    }
    return removalCount;
}

void DensePlanting::PlanActions() {

    for (int i = 0; i < 365; i++) {
        if (this->notModified) {
            q[i] = true;
        }
        else {
            q[i] = false;
        }
    }
}

void DensePlanting::executeAction(Grove *g) {
    Commodity* c = g->getCrop();
    c->freshYield = (1 + this->yieldMultiplier) * c->freshYield;
    c->juiceYield = (1 + this->yieldMultiplier) * c->juiceYield;
    g->fixedCosts += this->annualCosts;
    this->notModified = false;
    this->PlanActions();
}
void RogueTrees::executeAction(Grove *g) {
    int numRemoved = checkAndRogue(g->getIBounds(), g->getJBounds(), this->radius, this->radius);
    g->costs += this->surveyCost;
    g->costs += numRemoved * this->removalCost;
}
void RectangularRogue::executeAction(Grove *g) {
    int numRemoved = checkAndRogue(g->getIBounds(), g->getJBounds(), this->width, this->height);
    g->costs += this->surveyCost;
    g->costs += numRemoved * this->removalCost;
}
/*****************************************************************
 * 
 * ***************************************************************/
void RogueTrees::PlanActions() {
    //Fill with false flags
    for (int i = 0; i < 365; i++) {
        q[i] = false;
    }
    //Check once at start of year, and every frequency days after
    int daysPerSurvey = 365 / this->frequency;
    for (int i = 0; i < 365; i += daysPerSurvey) {
        q[i] = true;
    }
}

void RectangularRogue::PlanActions() {
    //Fill with false flags
    for (int i = 0; i < 365; i++) {
        q[i] = false;
    }
    //Check once at start of year, and every frequency days after
    int daysPerSurvey = 365 / this->frequency;
    for (int i = 0; i < 365; i += daysPerSurvey) {
        q[i] = true;
    }
}

void SprayTrees::executeAction(Grove *g) {
    sprayGrove(g->getIBounds(), g->getJBounds(), this->efficacy);
    g->costs += this->sprayCost;
}
/*******************************************************
* Group Action : Plan Actions
* CURRENTLY DEPRECATED PLANS TO REIMPLEMENT
********************************************************/

void SprayTrees::PlanActions() {
    //Fill with nulls
    for (int i = 0; i < 365; i++) {
        q[i] = false;
    }
    //Initialize RNG
    int target1 = 2;
    int target2 = 30;
    int gen_ub = floor(this->windowSize / 2);
    int gen_lb = -1 * gen_ub;
    boost::uniform_int<> gen(gen_lb, gen_ub);
    boost::random::mt19937 econ_rng(std::time(0));

    //Harvest 1
    int harvest1spray1day = this->start1 + target1;
    int harvest1spray2day = this->start1 + target2;
    q[harvest1spray1day] = true;
    q[harvest1spray2day] = true;

    //Harvest 2
    int harvest2spray1day = this->start2 + target1;
    int harvest2spray2day = this->start2 + target2;
    q[harvest2spray1day] = true;
    q[harvest2spray2day] = true;

    //Harvest 3
    int harvest3spray1day = this->start3 + target1;
    int harvest3spray2day = this->start3 + target2;
    q[harvest3spray1day] = true;
    q[harvest3spray2day] = true;

}

/*********************************************************
* Behavior : Get Expected Value
* Returns the expected value of a behavior pattern
**********************************************************/
double Behavior::getExpectedValue(Grove g, double risk, int projectionLength, int startingPeriod, 
                                  int planningLength, double sprayEfficacy, double alpha, double additionalCosts) {
    double EV = 0;
    double ui_outcome = 0;
    double i_outcome = 0;
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
   
    for (int t = startingPeriod; t < projectionLength; t++) {

        //Add variable costs from mitigation every planning period
        if (t % planningLength == 0) {
            ui_outcome -= this->getVariableCosts();
            i_outcome -= this->getVariableCosts();
            ui_outcome -= additionalCosts;
            i_outcome -= additionalCosts;
        }

        //Add crop income every harvest
        if (g.getCrop()->isHarvestPeriod(t % 365)) {
            int numCrops = (g.getIBounds()[1] - g.getIBounds()[0]) * (g.getJBounds()[1] - g.getJBounds()[0]);
            double severity;
            if (!g.foundHLB) {
                severity = this->hlbSpread(t - startingPeriod, &g);
            }
            else {
                severity = this->hlbSpread(t - g.foundHLB_day, &g);
            }
         
            double returns = g.getCrop()->getReturns();
            double infectedReturns = returns * getInfectedYield(severity);
            ui_outcome += returns * numCrops;
            i_outcome += infectedReturns * numCrops;
        }

        //Add variable costs from crops and fixed costs on land at the end of the year
        if (t % 365 == 0) {
            ui_outcome -= g.getCrop()->costs;
            i_outcome -= g.getCrop()->costs;
        }
    }
    EV = (risk * i_outcome) + ((1 - risk) * ui_outcome);
    return EV;
}

double* Behavior::getExpectedValueTester(Grove g, double risk, int projectionLength, int startingPeriod, 
                                  int planningLength, double sprayEfficacy, double alpha) {
    double EV = 0;
    double ui_outcome = 0;
    double i_outcome = 0;
    double costs = 0;
    int* ibounds = g.getIBounds();
    int* jbounds = g.getJBounds();
   
    for (int t = startingPeriod; t < projectionLength; t++) {

        //Add variable costs from mitigation every planning period
        if (t % planningLength == 0) {
            ui_outcome -= this->getVariableCosts();
            i_outcome -= this->getVariableCosts();
            costs += this->getVariableCosts();
        }

        //Add crop income every harvest
        if (g.getCrop()->isHarvestPeriod(t % 365)) {
            int numCrops = (g.getIBounds()[1] - g.getIBounds()[0]) * (g.getJBounds()[1] - g.getJBounds()[0]);
            double severity;
            if (!g.foundHLB) {
                severity = this->hlbSpread(t - startingPeriod, &g);
            }
            else {
                severity = this->hlbSpread(t - g.foundHLB_day, &g);
            }
         
            double returns = g.getCrop()->getReturns();
            double infectedReturns = returns * getInfectedYield(severity);
            ui_outcome += returns * numCrops;
            i_outcome += infectedReturns * numCrops;
        }

        //Add variable costs from crops and fixed costs on land at the end of the year
        if (t % 365 == 0) {
            ui_outcome -= g.getCrop()->costs;
            i_outcome -= g.getCrop()->costs;
            costs += g.getCrop()->costs;
        }
    }
    EV = (risk * i_outcome) + ((1 - risk) * ui_outcome);
    i_outcome += costs;
    ui_outcome += costs;
    double ER = (risk * i_outcome) + ((1 - risk) * ui_outcome);

    double* retval = new double[3];
    retval[0] = EV;
    retval[1] = costs;
    retval[2] = ER;
    return retval;
}

/********************************************************
 * Beta Spread
 * Measures HLB Spread using beta regression results
 * *****************************************************/
double betaSpread(int relT, double efficacy, string strategy, double alpha) {
    // Coefficients
    double intercept = -2.44978330;
    double e75_coef = -0.20743718;
    double e85_coef = -0.07299259;
    double alpha1_coef = 0.01494436;
    double indv_coef = -0.15463509;
    double group_coef = -0.10212705;
    double maxT_coef = 0.00228886;
    double e75alpha1_coef = -0.06266322;
    double e85alpha1_coef = -0.79121443;
    double e75indv_coef = -0.16353347;
    double e85indv_coef = -0.64362937;
    double e75group_coef = -0.13874484;
    double e85group_coef = -0.85018937;
    double e75maxT_coef = 0.00018147;
    double e85maxT_coef = 0.00020873;
    double indvmaxT_coef = -0.00009667;
    double groupmaxT_coef = -0.00013745;
    double e75indvmaxT_coef = -0.00003253;
    double e85indvmaxT_coef = -0.00017724;
    double e75groupmaxT_coef = -0.00005853;
    double e85groupmaxT_coef = -0.00014676;
    
    //covariates
    int e75 = (efficacy == 0.75);
    int e85 = (efficacy == 0.85);
    int indv = (strategy == "Individual Action");
    int group = (strategy == "Group Action");

    //sum em up
    double sum_alpha1 = 0;
    double sum_alpha0 = 0;
    //sum the alpha0 case
    sum_alpha0 += intercept + 
                  (e75_coef * e75) + 
                  (e85_coef * e85) + 
                  (indv_coef * indv) + 
                  (group_coef * group) +
                  (maxT_coef * relT) + 
                  (e75indv_coef * (e75 * indv)) +
                  (e85indv_coef * (e85 * indv)) + 
                  (e75group_coef * (e75 * group)) +
                  (e85group_coef * (e85 * group)) +
                  (e75maxT_coef * (e75 * relT)) +
                  (e85maxT_coef * (e85 * relT)) +
                  (indvmaxT_coef * (indv * relT)) +
                  (groupmaxT_coef * (group * relT)) +
                  (e75indvmaxT_coef * (e75 * indv * relT)) +
                  (e85indvmaxT_coef * (e85 * indv * relT)) +
                  (e75groupmaxT_coef * (e75 * group * relT)) +
                  (e85groupmaxT_coef * (e85 * group * relT));
    // add the additional alpha terms
    sum_alpha1 = sum_alpha0 + 
                 alpha1_coef + 
                 (e75alpha1_coef * e75) +
                 (e85alpha1_coef * e85);
    //Take the weighted sum based on alpha
    double sum = sum_alpha0 + alpha*(sum_alpha1 - sum_alpha0);
    //Apply the logistic transformation 
    double transformed = 1 / (1 + exp(-1*sum));
    return transformed;
}

/********************************************************
* No Action : HLB Spread
* Returns the mean hlb spread t days after infection
********************************************************/
/*
double NoAction::hlbSpread(int t, Grove *g) {
    double res = betaSpread(t, g->getSprayEfficacy(), "No Action", g->getAlpha());
    return res;
}*/

/*******************************************************
* No Action : Plan Actions
* Fills a planning queue with empty actions
********************************************************/
/*void NoAction::PlanActions() {
        for (int i = 0; i < 365; i++) {
            q[i] = false;
        }
}*/



/********************************************************
 * Individual Action: HLB Spread
 * Returns the mean hlb spread t days after infection
*********************************************************/
/*
double IndividualAction::hlbSpread(int t, Grove *g) {
    double res = betaSpread(t, g->getSprayEfficacy(), "Individual Action", g->getAlpha());
    return res;
}*/

/*******************************************************
* Individual Action : Plan Actions
* CURRENTLY DEPRECATED PLANS TO REIMPLEMENT
********************************************************/
/*
void IndividualAction::PlanActions() {
      return;
}*/

/********************************************************
 * Individual Action: HLB Spread
 * Returns the mean hlb spread t days after infection
*********************************************************/
/*
double GroupAction::hlbSpread(int t, Grove *g) {
    double res = betaSpread(t, g->getSprayEfficacy(), "Group Action", g->getAlpha());
    return res;
}*/

/*******************************************************
* Group Action : Plan Actions
* CURRENTLY DEPRECATED PLANS TO REIMPLEMENT
********************************************************/
/*
void GroupAction::PlanActions() {
    /*
        //Fill with nulls
        for (int i = 0; i < 365; i++) {
            q[i] = NULL;
        }
        //Initialize RNG
        plan_func spray = sprayGrove;
        int target1 = 16;
        int target2 = 30;
        int gen_ub = floor(this->windowSize / 2);
        int gen_lb = -1 * gen_ub;
        boost::uniform_int<> gen(gen_lb, gen_ub);
        boost::random::mt19937 econ_rng(std::time(0));

        //Harvest 1
        int harvest1spray1day = this->start1 + target1 + gen(econ_rng);
        int harvest1spray2day = this->start1 + target2 + gen(econ_rng);
        q[harvest1spray1day] = sprayGrove;
        q[harvest1spray2day] = sprayGrove;

        //Harvest 2
        int harvest2spray1day = this->start2 + target1 + gen(econ_rng);
        int harvest2spray2day = this->start2 + target2 + gen(econ_rng);
        q[harvest2spray1day] = sprayGrove;
        q[harvest2spray2day] = sprayGrove;

         //Harvest 3
        int harvest3spray1day = this->start3 + target1 + gen(econ_rng);
        int harvest3spray2day = this->start3 + target2 + gen(econ_rng);
        q[harvest3spray1day] = sprayGrove;
        q[harvest3spray2day] = sprayGrove;

}*/












