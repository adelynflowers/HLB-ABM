/*
#include <iostream>
#include "../headers/bioABM.h"
using namespace std;


int main(int argc, char ** argv) {
	auto start = std::chrono::high_resolution_clock::now();
	bioABM::parseParameterFile(string("configs/modified_default.json"));
	while (bioABM::getModelDay() < bioABM::getModelDuration()) {
		cout << "Period " << bioABM::getModelDay() << endl;
		bioABM::advanceBiologicalModel_parallel();
	}
	bioABM::finishRun();
	auto stop = chrono::high_resolution_clock::now();
	cout << "Total run time: " << chrono::duration_cast<chrono::seconds>(stop - start).count() << "s " << endl;
	return 0;
}
*/