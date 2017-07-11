#include <iostream>
#include <deque>
#include <vector>
#include <set>
#include <math.h>
#include <map>
#include <algorithm>
#include <admesh/stl.h>
#include "poly2tri/poly2tri.h"
//#include <poly2tri/poly2tri.h> tohle tam ma byt kdyz pouzivam -lpoly2tri
#include <string>
#include <limits>
#include <tuple>
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
  bool operator() (const p2t::Point* lhs, const p2t::Point* rhs) const;
  bool operator() (const p2t::Point lhs, const p2t::Point rhs) const;
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
  stl_file* export_stl2(deque<stl_facet> facets);
  void save();
  std::array<stl_file*,2> save2();
  void close(){stl_close(&mesh_file);}
  bool createBorderPolylines(bool processOnFac = true);
  //void tmpBorder();
  void findHoles();
  void triangulateCut(int topOrBot=0);
  void poly2triTest();
  std::array<string,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes);
  bool runUnitTests();
  
private:
  double calculatePolygonArea(vector<p2t::Point*> polygon);
  bool vertexInPolygon(int nvert, const vector<Point* >& vertex,  const double &testx, const double testy);
  void createFacet(vector<p2t::Triangle*> &triangles, int side = 0);
  stl_facet createFacet(stl_facet facet, int s, int i, stl_vertex intersect);
  stl_facet createFacet(stl_facet facet, int s, int i, stl_vertex intersect1, stl_vertex intersect2);
  void setMissingCoordinate(const p2t::Point* a,stl_vertex & b);
  void pushToPolylines(vector<p2t::Point*> &vec,stl_vertex vert);
  void pushToPolylinesFront(vector<p2t::Point*> &vec,stl_vertex vert);
  void volumeTest();
  void fixNonsimplePolygon(vector<p2t::Point*>& npolygon);
  bool checkForNewPoints(vector<p2t::Triangle*> &triangles, vector<p2t::Point*>& npolygon);
  void repairIfNonsimplePolygon();
  void setRemovedAxis();
  void setPlane(stl_plane plane);
  bool ccw(p2t::Point* a, p2t::Point* b, p2t::Point* c);
  bool intersect (p2t::Point* a, p2t::Point* b, p2t::Point* c, p2t::Point* d);
  void removeNonsimplePolygonPoints(vector<p2t::Point*> & p);
  bool acquireSaveName(string& name);
  //void checkPoly2triResult();
  void checkPoly2triResult( vector<p2t::Triangle*>& triangles );
  void initializeStl(stl_file * stl,int numOfFacets);
  void setVertex(stl_vertex& a, stl_vertex& b,const int &s,const stl_facet & facet);
  void processOnFacets();
  bool haveEqualEdges(tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet1, tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet2);
  void insertTo(stl_vertex x, stl_vertex y, vector<stl_vertex>& a, set<stl_vertex,setVertComp> & b);
  void sortPolylines();
  //Unit test methods
  bool t_setRemovedAxis();

  // Integration tests
  bool t_minMaxPointsSame();

  stl_file mesh_file;
  stl_plane plane=stl_plane(0,0,0,0);
  deque<stl_facet>top_facets,bot_facets;
  vector<stl_vertex>border;
  set<stl_vertex,setVertComp> border2; //border as set, used to calculate missing coordinate during 2d->3d conversion
  deque<stl_vertex>remainingBorder;
  vector<vector <p2t::Point*> > polylines;
  vector< vector< pair <vector<p2t::Point*>,int> > > polygonsWithHoles;
  vector<tuple<stl_facet,stl_position,stl_vertex,stl_vertex>> facetsOnPlane;
  // this vector contains vectors which contains pair<polygon, -1 for polygon and positive number representing in which polygon is this hole>
  int numOfPolylines=0;
  //float zCoord;
  
};


// vertex position related to the plane


 std::array<stl_file*,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes);



