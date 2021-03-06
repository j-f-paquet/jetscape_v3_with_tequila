#include <vector>
#include <string>
#include <iostream>
#include <memory>

#include "JetScapeTask.h"
#include "TequilaMutex.h"

using namespace std;
using std::shared_ptr;
using namespace Jetscape;

TequilaMutex::TequilaMutex() {}

TequilaMutex::~TequilaMutex() {}

bool TequilaMutex::CheckMutex(vector<shared_ptr<JetScapeTask>> modules) {
  bool isLbt = false;
  bool isAdscft = false;

  for (auto module : modules) {
    string name = module->GetId();
    if (!name.compare("LBT"))
      isLbt = true;
    if (!name.compare("AdSCFT"))
      isAdscft = true;
  }

  if (isLbt || isAdscft)
    return false;
  return true;
}
