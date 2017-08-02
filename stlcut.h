#include <iostream>
#include <deque>
#include <vector>
#include <set>
#include <math.h>
#include <map>
#include <algorithm>
#include <admesh/stl.h>
#include "poly2tri/poly2tri.h"
#include <string>
#include <limits>
#include <tuple>
#include <signal.h>
#include <setjmp.h>
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


struct comparatorStruct 
 {
  bool operator() (const stl_vertex& lhs, const stl_vertex& rhs) const;
  bool operator() (const p2t::Point* lhs, const p2t::Point* rhs) const;
  bool operator() (const p2t::Point lhs, const p2t::Point rhs) const;
  bool operator()(const stl_facet & a, const stl_facet & b) const;

 };

class Mesh
{
public:

  bool cut(stl_plane plane);
  void openStl(char* name);
  void setStl( stl_file file); 
  void save(string name = "");
  std::array<stl_file*,2> getFinalStls(); 
  void close();
  std::array<string,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes);
  void setOptions(bool silent, bool error_recovery); 
  bool runUnitTests();
  bool runIntegrationTests();
  
private:
  void exportStl(deque<stl_facet> facets, const char* name);
  stl_file* getExportedStl(deque<stl_facet> facets);
  bool createBorderPolylines(bool firstCall = true); 
  void findHoles();
  stl_position vertexPosition(stl_vertex vertex); 
  stl_vertex intersection(stl_vertex a, stl_vertex b); 
  void divideFacets();
  void triangulateCut(int topOrBot=0);
  void checkDuplicity(); 
  double calculatePolygonArea(vector<p2t::Point*> polygon); 
  bool vertexInPolygon(const vector<Point* >& polygon,  const double &testx, const double &testy); 
  void createFacets(vector<p2t::Triangle*> &triangles, int side = 0);
  stl_facet createFacet(stl_facet facet, int s, int i, stl_vertex intersect); 
  stl_facet createFacet(stl_facet facet, int s, int i, stl_vertex intersect1, stl_vertex intersect2); 
  stl_vertex getMissingCoordinate(const p2t::Point* a); 
  void pushBackToPolylines(vector<p2t::Point*> &vec,stl_vertex vert); 
  void pushFrontToPolylines(vector<p2t::Point*> &vec,stl_vertex vert); 
  void repairIfNonsimplePolygon(); 
  void setRemovedAxis(); 
  void setPlane(stl_plane plane); 
  bool ccw(p2t::Point* a, p2t::Point* b, p2t::Point* c); 
  bool edgesIntersect (p2t::Point* a, p2t::Point* b, p2t::Point* c, p2t::Point* d); 
  void removeNonsimplePolygonPoints(vector<p2t::Point*> & p); 
  void checkPoly2triResult( vector<p2t::Triangle*>& triangles ); 
  void initializeStl(stl_file * stl,int numOfFacets); 
  void setVertexFromFacet(stl_vertex& a, stl_vertex& b,const int &s,const stl_facet & facet); 
  bool processOnFacets();
  bool processOnBorder();
  bool haveEqualEdges(tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet1, tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet2); 
  void insertTo(stl_vertex x, stl_vertex y, vector<stl_vertex>& a, set<stl_vertex,comparatorStruct> & b);
  void sortPolylines(); 
  bool isStringValid(const string &str); 
  void pushOns(const int ons, stl_vertex& a, stl_vertex& b,const stl_facet &facet, const stl_position* pos);
  void pushAboveBelow(const int aboves,stl_vertex& a,stl_vertex& b,const stl_facet &facet, const stl_position* pos);
  void popTo(stl_vertex& a, stl_vertex& b); 
  void writeFails();
  void cleanupVariables();
  void deletePolygonsWithHoles();
  stl_vertex getVertex(double x, double y, double z);

  //Unit test methods
  bool t_setRemovedAxis();
  bool t_intersection();
  bool t_sortPolylines();
  bool t_calculatePolygonArea();
  bool t_vertexPosition(); 
  bool t_checkPoly2triResult();
  bool t_vertexInPolygon();
  bool t_getMissingCoordinate();
  bool t_checkDuplicity();
  bool t_createFacet();
  bool t_pushBackToPolylines();
  bool t_pushFrontToPolylines();
  bool t_ccw();
  bool t_edgesIntersect ();
  bool t_removeNonsimplePolygonPoints();
  bool t_setVertexFromFacet();
  bool t_isStringValid();
  bool t_setOptions();
  bool t_haveEqualEdges();
  bool t_popTo();

  // Integration tests
  bool t_minMaxPointsSameAfterCut();
  bool t_sameVolume();
  bool t_allTrianglesOriginalOrPartOfCut();
  bool t_noVertexOnOpositeSide();

  stl_file meshFile;
  stl_plane plane=stl_plane(0,0,0,0);
  deque<stl_facet>topFacets,botFacets;
  vector<stl_vertex>border,botBorder,topBorder;
  set<stl_vertex,comparatorStruct> originalVertices; //border as set, used to calculate missing coordinate during 2d->3d conversion
  vector<vector <p2t::Point*> > polylines;
  vector< vector< pair <vector<p2t::Point*>,int> > > polygonsWithHoles;
  vector<tuple<stl_facet,stl_position,stl_vertex,stl_vertex>> facetsOnPlane;
  // this vector contains vectors which contains pair<polygon, -1 for polygon and positive number representing in which polygon is this hole>
  int numOfPolylines=0;
  vector<int> fails;
  bool silent = false;
  bool errorRecovery = true;
  
};



