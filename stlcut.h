#include <iostream>
#include <deque>
#include <vector>
#include <set>
#include <math.h>
#include <map>
#include <algorithm>
#include <admesh/stl.h>
#include <poly2tri/poly2tri.h>
#include <string>
#include <limits>
using namespace std;
using namespace p2t;

enum stl_position { above, on, below };

// makes more sense in this program
typedef stl_vertex stl_vector;
struct stl_plane
{
  float x;
  float y;
  float z;
  float d;
  stl_plane(float x, float y, float z,float d);
};

struct setVertComp 
 {
  bool operator() (const stl_vertex& lhs, const stl_vertex& rhs) const;
 };

class Mesh
{
public:
   //Mesh();
  ~Mesh();
  void cut(stl_plane plane);
  stl_position vertex_position(stl_vertex vertex);
  stl_vertex intersection(stl_vertex a, stl_vertex b);
  void open(char* name);
  void open( stl_file file);
  void export_stl(deque<stl_facet> facets, const char* name);
  stl_file* export_stl2(deque<stl_facet> facets, const char* name);
  void save();
  std::array<stl_file*,2> save2();
  void close(){stl_close(&mesh_file);}
  bool createBorderPolylines();
  void tmpBorder();
  void findHoles();
  void triangulateCut();
  std::array<string,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes);
private:
  double calculatePolygonArea(vector<p2t::Point*> polygon);
  bool vertexInPolygon(int nvert, const vector<Point* >& vertex,  const double &testx, const double testy);
  void createFaces(vector<p2t::Triangle*> &triangles);
  void setMissingCoordinate(const p2t::Point* a,stl_vertex & b);
  void pushToPolylines(vector<p2t::Point*> &vec,stl_vertex vert);
  void volumeTest();
  stl_file mesh_file;
  stl_plane plane=stl_plane(0,0,0,0);
  deque<stl_facet>top_facets,bot_facets;
  vector<stl_vertex>border;
  set<stl_vertex,setVertComp> border2;
  deque<stl_vertex>remainingBorder;
  vector<vector <p2t::Point*> > polylines;
  vector< vector< pair <vector<p2t::Point*>,int> > > polygonsWithHoles; 
  // this vector contains vectors which contains pair<polygon, -1 for polygon and positive number representing in which polygon is this hole>
  int numOfPolylines=0;
  float zCoord;
  
};


// vertex position related to the plane


 std::array<stl_file*,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes);



