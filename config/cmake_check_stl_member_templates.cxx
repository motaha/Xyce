#include <list>
#include <vector>

using namespace std;

int main ( void ) 
{ 
  vector<int> v(2);

  v[0] = v[1] = 0;

  list<char> l( v.begin(), v.end() );

  return 0; 
}
