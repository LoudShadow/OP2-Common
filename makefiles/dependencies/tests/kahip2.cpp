#include <parhip_interface.h>
#include <stdio.h>
int main() {

  // Dummy function call to make sure linker doesn't do anything strange -
  // this executable will never actually be run
ParHIPPartitionKWay(nullptr, nullptr, nullptr, nullptr, nullptr,
                         nullptr, nullptr, true, 0, 0, nullptr, nullptr,
                         nullptr);
  return 0;
}
