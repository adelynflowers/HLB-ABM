#include <iostream>
#include "headers/bioABM.h"

using namespace std;


int main() {
	for (int i = 0; i < 365; i++) {
		bioABM::advanceBiologicalModel();
	}
	bioABM::finishRun();
	return 0;
}