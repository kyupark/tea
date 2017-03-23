using namespace std;

#include "ReadInfo.h"

ReadInfo::ReadInfo() : used(false) { }

bool ReadInfo::is_used() const { return used; }

void ReadInfo::set_used() { used = true; }
