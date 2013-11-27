
#include <unistd.h>

#include <ctime>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct triple
{
  int x, y, z;
  
  triple() : x(0), y(0), z(0) {
  }
  
  bool operator<(const triple &other) const {
    return (x < other.x && y < other.y && z < other.z);
  }
};

typedef map<triple, double> tensor_t;

int bounded_rand(int bound)
{
  return rand() % bound;
}

template<class T>
void add_point(int i, T &data, int bound = 128)
{
  triple t;
  t.x = bounded_rand(bound);
  t.y = bounded_rand(bound);
  t.z = bounded_rand(bound);
  data.insert(make_pair(t, 1.0));
  printf("%5d %d %d %d\n", i, t.x, t.y, t.z);
}

static int maxdim = 512;
static int nnz    = (int)(0.01 * (maxdim * maxdim * maxdim));

int main(int argc, char **argv)
{
  tensor_t data;
  
  //cout << nnz << " / " << (maxdim * maxdim * maxdim) << endl;
  //return 0;
  
  srand(time(0));
  for (int i = 0; i < nnz; ++i) {
    add_point(i, data, maxdim);
  }
  while (data.size() < nnz) {
    add_point(data.size(), data, maxdim);
  }
  
  return 0;
}
