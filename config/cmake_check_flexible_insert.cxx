#include <vector>
#include <list>

using namespace std;

int main ( void ) 
{ 
  list<int> ipList(1,0); 
  vector<int> ipVec; 

  ipVec.insert(ipVec.begin(),ipList.begin(),ipList.end()); 

  return 0; 
}
