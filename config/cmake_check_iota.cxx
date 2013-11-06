#include <algorithm>
#include <numeric>
#include <vector>

using namespace std;

int main ( void ) 
{ 
  vector<int> ipVec( 5, 0 );

  iota( ipVec.begin(), ipVec.end(), 1 );

  return 0; 
}
