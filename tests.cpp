
#include "stlcut.h"
#include <iostream>
using namespace std;

extern char removedAxis;

extern bool operator == (const stl_vertex& a, const stl_vertex& b);

bool meshIntegrationTests(Mesh& mesh, int pos, char* name, vector <pair<char* , int >> &testFails)
{
  bool result = mesh.runIntegrationTests();
  if(result)
  {
    cout<<endl<<"All tests completed without a problem."<<endl<<endl;
  	return true;
  }

  else
  {
    cout<<endl<<"Some tests failed."<<endl<<endl;
    testFails.push_back(make_pair(name,pos));
    return false;
  }
}


bool Mesh::t_intersection()
{
  fails.clear();
  int num=1;
  stl_vertex vert, a, b ;
  cout<<"Testing intersection"<<endl;
  setPlane(stl_plane(0, 0, 1, 0));
  setRemovedAxis();

  vert = intersection(getVertex(0, 0, 2), getVertex(0, 0, -2));
  if(!(vert == getVertex(0, 0, 0))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 0, 2), getVertex(1, 0, -2));
  if(!(vert == getVertex(1, 0, 0))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 1, 1), getVertex(-1, -1, -1));
  if(!(vert == getVertex(0, 0, 0))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 1, 1), getVertex(0, 0, 0));
  if(!(vert == getVertex(0, 0, 0))) 
    fails.push_back(num);
  num++;

  setPlane(stl_plane(1, 0, 0, 5));
  setRemovedAxis();

  vert = intersection(getVertex(6, 0, 0), getVertex(1, 0, 0));
  if(!(vert == getVertex(5, 0, 0))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(7, 1, 2), getVertex(-3, 5, -2));
  if(!(vert == getVertex(5, 1.8, 1.2))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(5, 0, 2), getVertex(10, 3, 5));
  if(!(vert == getVertex(5, 0, 2))) 
    fails.push_back(num);
  num++;

  setPlane(stl_plane(0, 0.5, 0, -5));
  setRemovedAxis();

  vert = intersection(getVertex(3, 0, -2.5), getVertex(11, 0.8, 22));
  if(!(vert == getVertex(-47, -5, -155.625))) 
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(7, 1, 2), getVertex(-3, 5, -2));
  if(!(vert == getVertex(22, -5, 8))) 
    fails.push_back(num);
  num++;


  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_edgesIntersect ()
{
  cout<<"Testing edgesIntersect"<<endl;
  fails.clear();
  int num=1;
  p2t::Point p1,p2,p3,p4;
  p1 = Point(0,0);
  p2 = Point(0,2);
  p3 = Point(-1,-1);
  p4 = Point(1,1);

  if(edgesIntersect(&p1, &p2, &p3, &p4) != true)
    fails.push_back(num);
  num++;

  if(edgesIntersect(&p3, &p2, &p1, &p4) != false)
    fails.push_back(num);
  num++;

  p3 = Point(-2,-2);
  p4 = Point(2,2);
  if(edgesIntersect(&p1, &p2, &p3, &p3) != false)
    fails.push_back(num);
  num++;

  p3 = Point(-2 ,-2 + (1e-25));
  p4 = Point(2,2 - (1e-25));
  if(edgesIntersect(&p1, &p2, &p3, &p4) != true)
    fails.push_back(num);
  num++;


  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_calculatePolygonArea()
{
  cout<<"Testing calculatePolygonArea"<<endl;
  fails.clear();
  int num=1;
  double area;
  vector<p2t::Point*> pol;
  pol.push_back(new p2t::Point(0,0));
  pol.push_back(new p2t::Point(0,3));
  pol.push_back(new p2t::Point(3,3));
  pol.push_back(new p2t::Point(3,0));
  area = calculatePolygonArea(pol);
    if(area != 9)
     fails.push_back(num);
  num++;
  pol.clear();

  pol.push_back(new p2t::Point(3,0));
  pol.push_back(new p2t::Point(6,6));
  pol.push_back(new p2t::Point(3,12));
  pol.push_back(new p2t::Point(-3,12));
  pol.push_back(new p2t::Point(-6,6));
  pol.push_back(new p2t::Point(-3,0));
  area = calculatePolygonArea(pol);
    if(!area == 108)
     fails.push_back(num);
  num++;
  pol.clear();

  pol.push_back(new p2t::Point(0,0));
  pol.push_back(new p2t::Point(0,3));
  pol.push_back(new p2t::Point(3,3));
  area = calculatePolygonArea(pol);
    if(area != 4.5)
     fails.push_back(num);
  num++;
  pol.clear();

  pol.push_back(new p2t::Point(0,0));
  pol.push_back(new p2t::Point(3,3));
  pol.push_back(new p2t::Point(6,6));
  pol.push_back(new p2t::Point(9,9));
  area = calculatePolygonArea(pol);
    if(area != 0)
     fails.push_back(num);
  num++;
  pol.clear();

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;


}

bool Mesh::t_sortPolylines()
{
  cout<<"Testing sortPolylines"<<endl;
  fails.clear();
  int num=1;
  vector < p2t::Point*> pol;
  pol.push_back(new p2t::Point(1,1));
  pol.push_back(new p2t::Point(3,3));
  pol.push_back(new p2t::Point(2,4));
  polylines.push_back(pol);
  pol.clear();
  pol.push_back(new p2t::Point(0,0));
  pol.push_back(new p2t::Point(0,5));
  pol.push_back(new p2t::Point(5,5));
  pol.push_back(new p2t::Point(10,10));
  pol.push_back(new p2t::Point(-10,8));
  polylines.push_back(pol);
  pol.clear();
  pol.push_back(new p2t::Point(2,2));
  pol.push_back(new p2t::Point(3,3));
  polylines.push_back(pol);
  sortPolylines();
  if(!(polylines.size() == 2  && (*polylines[0][0]) == p2t::Point(1,1) && (*polylines[1][0]) == p2t::Point(0,0) ))
    fails.push_back(num);  
  num++;
  pol.clear();

  pol.push_back(new p2t::Point(22,22));
  pol.push_back(new p2t::Point(25,20));
  pol.push_back(new p2t::Point(20,20));
  polylines.push_back(pol);
  pol.clear();
  pol.push_back(new p2t::Point(-1,-1));
  pol.push_back(new p2t::Point(80,80));
  pol.push_back(new p2t::Point(-80,0));
  polylines.push_back(pol);
  sortPolylines();
  if(!(polylines.size() == 4  && (*polylines[0][0]) == p2t::Point(1,1) && (*polylines[3][0]) == p2t::Point(-1,-1) ))
    fails.push_back(num); 
  num++;
  pol.clear();

  pol.push_back(new p2t::Point(-8,-8));
  pol.push_back(new p2t::Point(-9,-9));
  pol.push_back(new p2t::Point(-15,-15));
  polylines.push_back(pol);
  pol.clear();
  pol.push_back(new p2t::Point(10,0));
  pol.push_back(new p2t::Point(50,0));
  pol.push_back(new p2t::Point(100,0));
  polylines.push_back(pol);
  sortPolylines();
  if(!(polylines.size() == 4  && (*polylines[0][0]) == p2t::Point(1,1) && (*polylines[3][0]) == p2t::Point(-1,-1) ))
    fails.push_back(num); 

  for (int i = 0; i < polylines.size(); ++i)
  {
    for (int j = 0; j < polylines[i].size(); ++j)
    {
      delete polylines[i][j];
    }
  }
  polylines.clear();

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

void Mesh::writeFails()
{
  cout<<"   Test/s number: ";
  for (int i = 0; i < fails.size(); ++i)
  {
    cout<<fails[i]<<" ";
  }
  cout<<"failed."<<endl;
}
bool Mesh::t_vertexPosition()
{
  fails.clear();
  int num=1;
  cout<<"Testing vertexPosition"<<endl;
  setPlane(stl_plane(1, 1, -1, 5));
  setRemovedAxis();

  if(!(vertexPosition(getVertex(0,0,10)) == below))
    fails.push_back(num);
  num++;
  if(!(vertexPosition(getVertex(0,0,-10)) == above))
    fails.push_back(num);
  num++;

  setPlane(stl_plane(1, 0, 0, 0));
  setRemovedAxis();
  if(!(vertexPosition(getVertex(1,5,10)) == above))
    fails.push_back(num);
  num++;
  if(!(vertexPosition(getVertex(1,5,-10)) == above))
    fails.push_back(num);
  num++;
  if(!(vertexPosition(getVertex(0,0,0)) == on))
    fails.push_back(num);
  num++;
  if(!(vertexPosition(getVertex(0,0,5)) == on))
    fails.push_back(num);
  num++;

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
} 

bool Mesh::t_checkPoly2triResult()
{
  fails.clear();
  int num=1;
  cout<<"Testing checkPoly2triResult"<<endl;
  setPlane(stl_plane(0, 0, 1 ,0));
  setRemovedAxis();
  vector<p2t::Triangle*> triangles;
  p2t::Triangle* triangle;
  vector<p2t::Point> point;
  point.push_back(p2t::Point(0,0));
  point.push_back(p2t::Point(2,0));
  point.push_back(p2t::Point(2,2));
  point.push_back(p2t::Point(3,3));
  point.push_back(p2t::Point(0,4));
  point.push_back(p2t::Point(4,4));
  point.push_back(p2t::Point(5,0));
  triangle = new Triangle(point[0],point[1],point[2]);  triangles.push_back(triangle);
  triangle = new Triangle(point[2],point[5],point[6]);  triangles.push_back(triangle);
  triangle = new Triangle(point[0],point[2],point[4]);  triangles.push_back(triangle);
  triangle = new Triangle(point[3],point[2],point[5]);  triangles.push_back(triangle);
  originalVertices.insert(getVertex(0,0,0));
  originalVertices.insert(getVertex(2,0,0));
  originalVertices.insert(getVertex(2,2,0));
  originalVertices.insert(getVertex(0,4,0));
  originalVertices.insert(getVertex(5,0,0));
  originalVertices.insert(getVertex(4,4,0));

  checkPoly2triResult(triangles);
  if(!(triangles.size() == 3) )
    fails.push_back(num);
  num++;

  triangle = new Triangle(point[0],point[1],point[3]); triangles.push_back(triangle);
  triangle = new Triangle(point[3],point[1],point[5]); triangles.push_back(triangle);
  checkPoly2triResult(triangles);
  if(!(triangles.size() == 3) )
    fails.push_back(num);
  num++;

  originalVertices.erase(getVertex(4,4,0));
  checkPoly2triResult(triangles);
  if(!(triangles.size() == 2 && triangles[0]->GetPoint(0)->x == 0 && triangles[0]->GetPoint(0)->y == 0 && triangles[1]->GetPoint(2)->x == 0 && triangles[1]->GetPoint(2)->y == 4) )
    fails.push_back(num);
  num++;

  setPlane(stl_plane(1, 0, 1 ,0));
  originalVertices.insert(getVertex(4,4,0));
  triangle = new Triangle(point[0],point[2],point[4]);  triangles.push_back(triangle);
  checkPoly2triResult(triangles);
  if(!(triangles.size() == 3) )
    fails.push_back(num);
  num++;

  for (int i = 0; i < triangles.size(); ++i)
  {
    delete triangles[i];
  }
  originalVertices.clear();
  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_vertexInPolygon()
{
  fails.clear();
  int num=1;
  cout<<"Testing vertexInPolygon"<<endl;
  vector<Point*> polygon;
  polygon.push_back(new p2t::Point(0,0));
  polygon.push_back(new p2t::Point(0,3));
  polygon.push_back(new p2t::Point(3,3));
  polygon.push_back(new p2t::Point(3,0));
  if(!(vertexInPolygon(polygon, 1, 1) == true))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 2, 2) == true))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 2.99999999999, 2.999999999) == true))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 1e-25, 1e-25) == true))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, -5, -1) == false))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 3, 1) == false)) // points on the edge
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 1, 3) == false))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 3, 0) == false)) // points on the edge
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, 0, 3) == false))
    fails.push_back(num);
  num++;

  for (std::vector<p2t::Point*>::iterator i = polygon.begin(); i != polygon.end(); ++i)
  {
    delete (*i);
  }
  polygon.clear();

  polygon.push_back(new p2t::Point(5,5));
  polygon.push_back(new p2t::Point(0,4));
  polygon.push_back(new p2t::Point(-3,3));
  polygon.push_back(new p2t::Point(-5,5));
  polygon.push_back(new p2t::Point(-6,0));
  polygon.push_back(new p2t::Point(-2,-2));

  if(!(vertexInPolygon(polygon, -5, 2.5) == true))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, -6, 2.5) == false))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, -3, (3+1e-25) == false)))
    fails.push_back(num);
  num++;

  if(!(vertexInPolygon(polygon, (-6 +1e-25), 0) == true))
    fails.push_back(num);
  num++;

  for (std::vector<p2t::Point*>::iterator i = polygon.begin(); i != polygon.end(); ++i)
  {
    delete (*i);
  }

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true; 
}

bool Mesh::t_getMissingCoordinate()
{
  fails.clear();
  int num=1;
  stl_vertex x;
  p2t::Point p; 
  vector<stl_vertex> vert;
  cout<<"Testing getMissingCoordinate"<<endl;
  setPlane(stl_plane (0, 0, 1, 0));
  setRemovedAxis();
  vert.push_back(getVertex(0, 0, 0));
  vert.push_back(getVertex(2, 0, 0));
  vert.push_back(getVertex(2, 2, 0));
  vert.push_back(getVertex(10, 0, 0));
  vert.push_back(getVertex(0, 5, 0));
  vert.push_back(getVertex(-3, -8, 0));
  vert.push_back(getVertex(-3, -9, 0));
  for (std::vector<stl_vertex>::iterator i = vert.begin(); i != vert.end(); ++i)
  {
    originalVertices.insert((*i));
  }

  p = p2t::Point(0,0);
  x = getMissingCoordinate(&p);
  if(!(x == vert[0]))
    fails.push_back(num);
  num++;

  p = p2t::Point(2,0);
  x = getMissingCoordinate(&p);
  if(!(x == vert[1]))
    fails.push_back(num);
  num++;

  p = p2t::Point(10.5,0);
  x = getMissingCoordinate(&p);
  if(x == vert[3] || x.x != 10.5f || x.y != 0)
    fails.push_back(num);
  num++;

  p = p2t::Point(-3,-8);
  x = getMissingCoordinate(&p);
  if(!(x == vert[5]))
    fails.push_back(num);
  num++;

  p = p2t::Point(-3,-8.1);
  x = getMissingCoordinate(&p);
  if(x == vert[5])
    fails.push_back(num);
  num++;

  vert.clear();
  vert.push_back(getVertex(6, 0 ,5));
  vert.push_back(getVertex(6, 1 ,5));
  vert.push_back(getVertex(6, 2 ,5));
  vert.push_back(getVertex(6, 3 ,5));
  vert.push_back(getVertex(6, 5 ,5));
  setPlane(stl_plane (1, 0, 0, 6));
  setRemovedAxis();

  p = p2t::Point(-3,-8.1);
  x = getMissingCoordinate(&p);
  if(!(x.z == -3 && x.y == -8.1f))
    fails.push_back(num);
  num++;

  p = p2t::Point(3,5);
  x = getMissingCoordinate(&p);
  if(!(x.x == 6 && x.y == 5 &&  x.z == 3))
    fails.push_back(num);
  num++;

  p = p2t::Point(5,3);
  x = getMissingCoordinate(&p);
  if(!(x == vert[3]))
    fails.push_back(num);
  num++;


 if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true; 
}

bool Mesh::t_checkDuplicity()
{
  fails.clear();
  int num = 1;
  cout<<"Testing checkDuplicity"<<endl;
  border.push_back(getVertex(1, 0 ,0));
  border.push_back(getVertex(2, 0 ,0));
  border.push_back(getVertex(0, 5 ,5));
  border.push_back(getVertex(1, 1.11 ,0));
  border.push_back(getVertex(1, 1.11 ,0));
  border.push_back(getVertex(0, 5 ,5));

  checkDuplicity();
  if(! (border.size() == 4))
    fails.push_back(num);
  num++;
  border.clear();

  border.push_back(getVertex(1, 0 ,1));
  border.push_back(getVertex(5, 0 ,2));
  border.push_back(getVertex(5, 0 ,2));
  border.push_back(getVertex(1, 0 ,1));
  border.push_back(getVertex(1, 1.11 ,0));
  border.push_back(getVertex(2, 1.11 ,3));
  border.push_back(getVertex(3, 3, 3));
  border.push_back(getVertex(1, 1 ,1));
  checkDuplicity();
  if(! (border.size() == 6))
    fails.push_back(num);
  num++;
  border.clear();


  border.push_back(getVertex(1, 0 ,1));
  border.push_back(getVertex(5, 0 ,2));
  border.push_back(getVertex(5, 0 ,3));
  border.push_back(getVertex(1, 0 ,1));
  border.push_back(getVertex(1, 1.11 ,0));
  border.push_back(getVertex(2, 1.11 ,3));
  border.push_back(getVertex(3, 3, 3));
  border.push_back(getVertex(1, 1 ,1));
  checkDuplicity();
  if(! (border.size() == 8))
    fails.push_back(num);
  num++;
  border.clear();

  border.push_back(getVertex(1, 1 ,1));
  border.push_back(getVertex(2, 2 ,2));
  border.push_back(getVertex(1, 1 ,1));
  border.push_back(getVertex(3, 3 ,3));
  border.push_back(getVertex(0, 0 ,0));
  border.push_back(getVertex(0, 0 ,0));
  border.push_back(getVertex(0, 0, 0));
  border.push_back(getVertex(1, 1 ,1));
  border.push_back(getVertex(1, 1, 1));
  border.push_back(getVertex(2, 2 ,2));
  checkDuplicity();
  if(! (border.size() == 8))
    fails.push_back(num);
  num++;
  border.clear();

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_createFacet()
{
  fails.clear();
  int num=1;
  cout<<"Testing createFacet"<<endl;
  stl_facet fac,testFac;
  fac.vertex[0] = getVertex(0 ,0 , 0);
  fac.vertex[1] = getVertex(1 ,0, -1);
  fac.vertex[2] = getVertex(2.2, 0 ,10);
  testFac=createFacet(fac,0,1,getVertex(2.2, 0, 10));

  if(!(testFac.vertex[0] == fac.vertex[0] && testFac.vertex[1] == fac.vertex[1] && testFac.vertex[2] == fac.vertex[2]) )
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,0,2,getVertex(2.2, 0, 10));
  if(!(testFac.vertex[0] == fac.vertex[0] && testFac.vertex[1] == fac.vertex[2] && testFac.vertex[2] == getVertex(2.2, 0, 10)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,1,1,getVertex(2.2, 0, 10));
  if(!(testFac.vertex[0] == fac.vertex[1] && testFac.vertex[1] == fac.vertex[2] && testFac.vertex[2] == getVertex(2.2, 0, 10)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,1,2,getVertex(2.2, 0, 10));
  if(!(testFac.vertex[0] == fac.vertex[1] && testFac.vertex[1] == fac.vertex[0] && testFac.vertex[2] == getVertex(2.2, 0, 10)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,2,1,getVertex(2.2, 0, 10));
  if(!(testFac.vertex[0] == fac.vertex[2] && testFac.vertex[1] == fac.vertex[0] && testFac.vertex[2] == getVertex(2.2, 0, 10)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,2,2,getVertex(2.2, 0, 10));
  if(!(testFac.vertex[0] == fac.vertex[2] && testFac.vertex[1] == fac.vertex[1] && testFac.vertex[2] == getVertex(2.2, 0, 10)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,0,0,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == fac.vertex[0] && testFac.vertex[1] == getVertex(2.2,0,10) && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,1,0,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == fac.vertex[1] && testFac.vertex[1] == getVertex(2.2,0,10) && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,2,0,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == fac.vertex[2] && testFac.vertex[1] == getVertex(2.2,0,10) && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,0,1,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(2.2,0,10) && testFac.vertex[1] == fac.vertex[1] && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,1,1,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(2.2,0,10) && testFac.vertex[1] == fac.vertex[2] && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,2,1,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(2.2,0,10) && testFac.vertex[1] == fac.vertex[0] && testFac.vertex[2] == getVertex(1,2,3)))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,0,2,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(1,2,3) && testFac.vertex[1] == fac.vertex[1] && testFac.vertex[2] == fac.vertex[2]))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,1,2,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(1,2,3) && testFac.vertex[1] == fac.vertex[2] && testFac.vertex[2] == fac.vertex[0]))
    fails.push_back(num);
  num++;

  testFac=createFacet(fac,2,2,getVertex(2.2, 0, 10),getVertex(1, 2 ,3));
  if(!(testFac.vertex[0] == getVertex(1,2,3) && testFac.vertex[1] == fac.vertex[0] && testFac.vertex[2] == fac.vertex[1]))
    fails.push_back(num);
  num++;

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}


bool Mesh::t_setVertexFromFacet()
{
  fails.clear();
  int num=1;
  cout<<"Testing setVertexFromFacet"<<endl;
  stl_vertex a, b;
  stl_facet f;
  a = getVertex( 1,    2,    3);
  b = getVertex(-1.1, -2.2, -3.3);
  f.vertex[0] = getVertex(0  , -3.1,   99.9);
  f.vertex[1] = getVertex(1  ,  5.2,  -1.5);
  f.vertex[2] = getVertex(2.2, -11.1 , 10);

  setVertexFromFacet(a,b,0,f);
  if(!(a == f.vertex[1] && b == f.vertex[2]))
    fails.push_back(num);
  num++;

  setVertexFromFacet(a,b,1,f);
  if(!(a == f.vertex[2] && b == f.vertex[0]))
    fails.push_back(num);
  num++;

  setVertexFromFacet(a,b,2,f);
  if(!(a == f.vertex[0] && b == f.vertex[1]))
    fails.push_back(num);
  num++;


  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}
bool Mesh::t_haveEqualEdges()
{
  fails.clear();
  int num=1;
  cout<<"Testing haveEqualEdges"<<endl;
  stl_facet f1, f2;
  f1.vertex[0] = getVertex(0  , 10,   5);
  f1.vertex[1] = getVertex(0  ,  0,  5);
  f1.vertex[2] = getVertex(20, 1 , 5);
  f2.vertex[0] = getVertex(0  , 10,   5);
  f2.vertex[1] = getVertex(-5  ,  12,  6);
  f2.vertex[2] = getVertex(0, -0 , 5);
  tuple<stl_facet, stl_position, stl_vertex, stl_vertex> t1, t2;
  t1 = make_tuple(f1,on,f1.vertex[0],f1.vertex[2]);
  t2 = make_tuple(f2,above, f2.vertex[0],f2.vertex[2]);

  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,on, f1.vertex[1], f1.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,on, f1.vertex[1], f1.vertex[0]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t2 = make_tuple(f2,above, f2.vertex[0], f2.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  f1.vertex[2] = getVertex(20, 1 , 1);
  t1 = make_tuple(f1,below, f1.vertex[0], f1.vertex[1]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,below, f1.vertex[1], f1.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != false)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,below, f1.vertex[2], f1.vertex[0]);
  if( haveEqualEdges( t1,t2 ) != false)
    fails.push_back(num);
  num++;


 if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

//removeNonSimplePolygons is designed to handle very small and specific amount of errorrs in polygon, so this test is taking it into account and tests only what
//this method is designed to handle.
bool Mesh::t_removeNonsimplePolygonPoints()
{
  fails.clear();
  int num=1;
  cout<<"Testing removeNonsimplePolygonPoints"<<endl;
  vector<p2t::Point*> p;
  p.push_back(new p2t::Point(0,0));
  p.push_back(new p2t::Point(10.1,0));
  p.push_back(new p2t::Point(12.1,3.3));
  p.push_back(new p2t::Point(20.123,19.23));
  p.push_back(new p2t::Point(-5,10));
  p.push_back(new p2t::Point(-20.2,-15.5));

  removeNonsimplePolygonPoints(p);
  if(p.size() != 6 )
    fails.push_back(num);
  num++;

  p.push_back(new p2t::Point(-21.1,-12.89));
  removeNonsimplePolygonPoints(p);
  if(!(p.size() == 6 && p[p.size()-1]->x == -20.2 && p[p.size()-1]->y == -15.5))
    fails.push_back(num);
  num++;

  for (int i = 0; i < p.size(); ++i)
  {
    delete p[i];
  }
  p.clear();

  p.push_back(new p2t::Point( 0, 0));
  p.push_back(new p2t::Point( 2.2, 1.1 ));
  p.push_back(new p2t::Point( 2,1 ));
  p.push_back(new p2t::Point( 0, 5));
  p.push_back(new p2t::Point(-5, 5));
  p.push_back(new p2t::Point( -3,-3 ));
  p.push_back(new p2t::Point( -1,-1 ));

  removeNonsimplePolygonPoints(p);
  if(!(p.size() == 6 && p[2]->x == 0 && p[2]->y == 5))
    fails.push_back(num);
  num++;

  for (int i = 0; i < p.size(); ++i)
  {
    delete p[i];
  }

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_ccw()
{ 
  fails.clear();
  int num=1;
  cout<<"Testing ccw"<<endl;
  p2t::Point p1,p2,p3;
  p1 = Point(0,0);
  p2 = Point(0,1);
  p3 = Point(1,0);

  if(ccw(&p1, &p2, &p3) != false)
    fails.push_back(num);
  num++;

  if(ccw(&p1, &p3, &p2) != true)
    fails.push_back(num);
  num++;

  p1 = Point(-10,-10);
  p2 = Point(0,-10);
  p3 = Point(-7.1,-3);

  if(ccw(&p1, &p2, &p3) != true)
    fails.push_back(num);
  num++;

  if(ccw(&p1, &p3, &p2) != false)
    fails.push_back(num);
  num++;

  if(ccw(&p3, &p2, &p1) != false)
    fails.push_back(num);
  num++;

  if(ccw(&p2, &p1, &p3) != false)
    fails.push_back(num);
  num++;

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_popTo()
{
  fails.clear();
  int num=1;
  cout<<"Testing popTo"<<endl;

  border.push_back(getVertex(1,123.3 ,0));
  border.push_back(getVertex(213,56 ,-312));
  border.push_back(getVertex(-129, 9384 ,214.4217));
  border.push_back(getVertex(9571,-623.98 ,0.56));

  stl_vertex a;
  stl_vertex b;

  popTo(a,b);
  if(!(a == getVertex(9571,-623.98 ,0.56) && b == getVertex(-129, 9384 ,214.4217) && border.size() == 2))
    fails.push_back(num);
  num++;

  popTo(a,b);
  if(!(b == getVertex(1,123.3 ,0) && a == getVertex(213,56 ,-312) && border.size() == 0))
    fails.push_back(num);
  num++;

 if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_pushBackToPolylines()
{
  fails.clear();
  int num=1;
  cout<<"Testing pushBackToPolylines"<<endl;
  setPlane(stl_plane (0, 0, 1, 0));
  setRemovedAxis();
  polylines.resize(1);

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 1 && polylines[0][0]->x == 1.1f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 2 && polylines[0][1]->x == 1.1f && polylines[0][1]->y == -2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3 && polylines[0][2]->x == -1.1f && polylines[0][2]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (1, 0, 0, 0));
  setRemovedAxis();

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 4 && polylines[0][3]->x == 3.3f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 5 && polylines[0][4]->x == 3.3f && polylines[0][1]->y == -2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 6 && polylines[0][5]->x == 3.3f && polylines[0][2]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 6))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (0, 1, 0, 0));
  setRemovedAxis();

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 7 && polylines[0][6]->x == 1.1f && polylines[0][6]->y == 3.3f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,-3.3));
  if(!(polylines[0].size() == 8 && polylines[0][7]->x == 1.1f && polylines[0][7]->y == -3.3f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 9 && polylines[0][8]->x == -1.1f && polylines[0][8]->y == 3.3f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,-2.2,3.3));
  if(!(polylines[0].size() == 9))
    fails.push_back(num);
  num++;

  for (int i = 0; i < polylines[0].size(); ++i)
  {
    delete polylines[0][i];
  }
  polylines.clear();

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true; 
}

bool Mesh::t_pushFrontToPolylines()
{
  fails.clear();
  int num=1;
  cout<<"Testing pushFrontToPolylines"<<endl;
  setPlane(stl_plane (0, 0, 1, 0));
  setRemovedAxis();
  polylines.resize(1);

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 1 && polylines[0][0]->x == 1.1f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 2 && polylines[0][0]->x == 1.1f && polylines[0][0]->y == -2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3 && polylines[0][0]->x == -1.1f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (1, 0, 0, 0));
  setRemovedAxis();

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 4 && polylines[0][0]->x == 3.3f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 5 && polylines[0][0]->x == 3.3f && polylines[0][0]->y == -2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 6 && polylines[0][0]->x == 3.3f && polylines[0][0]->y == 2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 6))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (0, 1, 0, 0));
  setRemovedAxis();

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 7 && polylines[0][0]->x == 1.1f && polylines[0][0]->y == 3.3f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,-3.3));
  if(!(polylines[0].size() == 8 && polylines[0][0]->x == 1.1f && polylines[0][0]->y == -3.3f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 9 && polylines[0][0]->x == -1.1f && polylines[0][0]->y == 3.3f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,-2.2,3.3));
  if(!(polylines[0].size() == 9))
    fails.push_back(num);
  num++;

  for (int i = 0; i < polylines[0].size(); ++i)
  {
    delete polylines[0][i];
  }
  polylines.clear();

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true; 

}


bool Mesh::t_setOptions()
{
  fails.clear();
  int num=1;
  cout<<"Testing setOptions"<<endl;

  setOptions(false, false);
  if(silent != false && errorRecovery != false)
    fails.push_back(num);
  num++;

  setOptions(true, true);
  if(silent != true && errorRecovery != true)
    fails.push_back(num);
  num++;

  setOptions(true, false);
  if(silent != true && errorRecovery != false)
    fails.push_back(num);
  num++;

  setOptions(false, true);
  if(silent != false && errorRecovery != true)
    fails.push_back(num);
  num++;


if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true; 
}


bool Mesh::t_setRemovedAxis()
{ 
  fails.clear();
  int num=1;
  cout<<"Testing setRemovedAxis"<<endl;
  stl_plane plane = stl_plane(1 , 0 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
     fails.push_back(num);
  num++;
  plane = stl_plane(0 , 1 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'y')
    fails.push_back(num);
  num++;
  plane = stl_plane(0 , 0 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'z')
    fails.push_back(num);
  num++;
  plane = stl_plane(1 , 0.5 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
    fails.push_back(num);
  num++;
  plane = stl_plane(1 , 1 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'y')
    fails.push_back(num);
  num++;
  plane = stl_plane(1 , 0 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
    fails.push_back(num);

  if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_isStringValid()
{
  fails.clear();
  int num=1;
  cout<<"Testing isStringValid"<<endl;
  string s;

  s = "abcASDH124526_ 912830asdA";
  if(!(isStringValid(s)))
    fails.push_back(num);
  num++;

  s = "abcASDH124526_ 91@##$#@(*)&2830asdA";
  if(isStringValid(s))
    fails.push_back(num);
  num++;

  s = "abcASDH12  //asd4526_ 912830asdA";
  if(isStringValid(s))
    fails.push_back(num);
  num++;

if(fails.size()>0)
  {
    writeFails();
    return false;
  }
  else return true;
}

bool Mesh::t_sameVolume()
{
    stl_calculate_volume(&meshFile);
    double org_volume = meshFile.stats.volume;
    double volume = 0.0;
    stl_file * botMesh;
    stl_file * topMesh;
    auto meshes = getFinalStls();
    topMesh = meshes[0];
    botMesh = meshes[1];
    stl_calculate_volume(topMesh);
    stl_calculate_volume(botMesh);
    volume += topMesh->stats.volume;
    volume += botMesh->stats.volume;
    if(abs(volume) <= org_volume * 1.01 && abs(volume) >= org_volume * 0.99)
      return true;
    return false;
}

bool Mesh::t_noVertexOnOpositeSide()
{
  auto meshes = getFinalStls();
  auto topMesh = meshes[0];
  stl_get_size(topMesh);
  auto botMesh = meshes[1];
  stl_get_size(botMesh);

  for (int i = 0; i < topMesh->stats.number_of_facets; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      stl_vertex & f = topMesh->facet_start[i].vertex[j];
      auto pos = vertexPosition(f);
      if(pos == below)
          if (abs(plane.x*f.x + plane.y*f.y + plane.z*f.z + plane.d) > 0.001)
            return false;
    }
  }
  for (int i = 0; i < botMesh->stats.number_of_facets; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      stl_vertex & f = botMesh->facet_start[i].vertex[j];
      auto pos = vertexPosition(f);
      if(pos == above)
        if  (abs(plane.x*f.x + plane.y*f.y + plane.z*f.z + plane.d) > 0.001)
          return false;
    }
  }

  return true;
}


/*
* Tests if every facet is from the original model or was modified due to being cut through.
*/
bool Mesh::t_allTrianglesOriginalOrPartOfCut()
{
  vector<stl_facet> originalFacets,topMeshFacets,botMeshFacets;
  for (int i = 0; i < meshFile.stats.number_of_facets; ++i)
  {
    originalFacets.push_back(meshFile.facet_start[i]);
  }

  sort(originalFacets.begin(),originalFacets.end(),[](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});

  auto meshes = getFinalStls();
  auto topMesh = meshes[0];
  stl_get_size(topMesh);
  auto botMesh = meshes[1];
  stl_get_size(botMesh);

  for (int i = 0; i < topMesh->stats.number_of_facets; ++i)
  {
    topMeshFacets.push_back(topMesh->facet_start[i]);
  }
  sort(topMeshFacets.begin(), topMeshFacets.end(),[](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});
  for (int i = 0; i < botMesh->stats.number_of_facets; ++i)
  {
    botMeshFacets.push_back(botMesh->facet_start[i]);
  }
  sort(botMeshFacets.begin(), botMeshFacets.end(),[](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});

  for (int i = 0; i < botMeshFacets.size(); ++i)
  {
    bool found = false;
    auto it = lower_bound(originalFacets.begin(),originalFacets.end(),botMeshFacets[i], [](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});
    while(it->vertex[0].x == botMeshFacets[i].vertex[0].x) 
    {
      if(botMeshFacets[i].vertex[0] == it->vertex[0])
        if(botMeshFacets[i].vertex[1] == it->vertex[1])
          if(botMeshFacets[i].vertex[2] == it->vertex[2])
          {
            found = true;
            break;
          }
          it++;
    }
    if(found == false)
    {
      int ons = 0;
      for (int k = 0; k < 3; ++k)
      {
        //if we  couldnt find the facet, lets test if one of its points is in the original border
        if(find(originalVertices.begin(),originalVertices.end(), botMeshFacets[i].vertex[k]) != originalVertices.end())
          ons++;
      }
      if(ons >= 1)
        continue;
      else
        return false;
    }
  }

  for (int i = 0; i < topMeshFacets.size(); ++i)
  {
    bool found = false;
    auto it = lower_bound(originalFacets.begin(),originalFacets.end(),topMeshFacets[i] ,[](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});
    while(it->vertex[0].x == topMeshFacets[i].vertex[0].x) 
    {
      if(topMeshFacets[i].vertex[0] == it->vertex[0])
        if(topMeshFacets[i].vertex[1] == it->vertex[1])
          if(topMeshFacets[i].vertex[2] == it->vertex[2])
          {
            found = true;
            break;
          }
          it++;
    }
    if(found == false)
    {
      int ons = 0;
      for (int k = 0; k < 3; ++k)
      {
        //if we  couldnt find the facet, lets test if one of its points is in the original border
        if(find(originalVertices.begin(),originalVertices.end(), topMeshFacets[i].vertex[k]) != originalVertices.end())
          ons++;
      }
      if(ons >= 1)
        continue;
      else 
        return false;
    }
  }

return true;
}



bool Mesh::t_minMaxPointsSameAfterCut()
{
  stl_get_size(&meshFile);
  stl_vertex maxOriginal;
  maxOriginal.x = meshFile.stats.max.x;
  maxOriginal.y = meshFile.stats.max.y;
  maxOriginal.z = meshFile.stats.max.z;
  stl_vertex minOriginal;
  minOriginal.x = meshFile.stats.min.x;
  minOriginal.y = meshFile.stats.min.y;
  minOriginal.z = meshFile.stats.min.z;

  stl_file * botMesh;
  stl_file * topMesh;
  auto meshes = getFinalStls();
  topMesh = meshes[0];
  botMesh = meshes[1];
  stl_get_size(botMesh);
  stl_get_size(topMesh);
  stl_vertex maxNew,minNew;

  maxNew.x = botMesh->stats.max.x > topMesh->stats.max.x ? botMesh->stats.max.x: topMesh->stats.max.x;
  maxNew.y = botMesh->stats.max.y > topMesh->stats.max.y ? botMesh->stats.max.y: topMesh->stats.max.y;
  maxNew.z = botMesh->stats.max.z > topMesh->stats.max.z ? botMesh->stats.max.z: topMesh->stats.max.z;
  minNew.x = botMesh->stats.min.x < topMesh->stats.min.x ? botMesh->stats.min.x: topMesh->stats.min.x;
  minNew.y = botMesh->stats.min.y < topMesh->stats.min.y ? botMesh->stats.min.y: topMesh->stats.min.y;
  minNew.z = botMesh->stats.min.z < topMesh->stats.min.z ? botMesh->stats.min.z: topMesh->stats.min.z;

  if(!(maxOriginal == maxNew && minOriginal == minNew))
    return false;
  return true;
}

bool Mesh::runIntegrationTests()
{
  bool success = true;
  bool result;

  result = t_minMaxPointsSameAfterCut(); 
  if(!result) 
    success = false;
  cout << "MinMax test was " <<(string)(!result ? "un" :"" )<<"succesfull"<<endl;

  result = t_sameVolume(); 
  if(!result) 
    success = false;
  cout << "Volume test was " <<(string)(!result ? "un" :"" )<<"succesfull"<<endl;

  result = t_allTrianglesOriginalOrPartOfCut();
  if(!result) 
    success = false;
  cout << "Triangle test was " <<(string)(!result ? "un" :"")<<"succesfull"<<endl;

  result = t_noVertexOnOpositeSide(); 
  if(!result) 
    success = false;
  cout << "Vertex test was " <<(string)(!result ? "un":"" )<<"succesfull"<<endl;

  return success;
}
bool Mesh::runUnitTests()
{
  cleanupVariables();
  bool success = true;
  if( !(t_setRemovedAxis()) ) 
    success = false;
  if( !(t_intersection()) ) 
    success = false;
  if( !(t_sortPolylines()) ) 
    success = false;
  if( !(t_calculatePolygonArea()) ) 
    success = false;
  if( !(t_vertexPosition()) ) 
    success = false;
  if( !(t_checkPoly2triResult()) ) 
    success = false;
  if( !(t_vertexInPolygon()) ) 
    success = false;
  if( !(t_getMissingCoordinate()) ) 
    success = false;
  if( !(t_checkDuplicity()) ) 
    success = false;
  if( !(t_createFacet()) ) 
    success = false;
  if( !(t_pushBackToPolylines()) ) 
    success = false;
  if( !(t_pushFrontToPolylines()) ) 
    success = false;
  if( !(t_ccw()) ) 
    success = false;
  if( !(t_edgesIntersect()) ) 
    success = false;
  if( !(t_removeNonsimplePolygonPoints()) ) 
    success = false;
  if( !(t_setVertexFromFacet()) ) 
    success = false;
  if( !(t_isStringValid()) ) 
    success = false;
  if( !(t_setOptions()) ) 
    success = false;
  if( !(t_haveEqualEdges()) ) 
    success = false;
  if( !(t_popTo()) ) 
    success = false;
  

  return success; 
}


string planeToString(stl_plane plane)
{
  string text;
  text += to_string(plane.x)+" ";
  text += to_string(plane.y)+" ";
  text += to_string(plane.z)+" ";
  text += to_string(-1 * plane.d);
  return text;
}

int main(int argc, char ** argv)
{
  if(argc >1 && (string(argv[1]) == "--help" || string(argv[1]) == "-h"))
  {
  	cout << "\e[1mUse:\e[0m "<<endl;
  	cout <<argv[0]<<" for unit tests."<<endl;
  	cout <<argv[0]<<" <stl_files> for also runing integration tests."<<endl;
  	cout << "Common use would be "<<argv[0]<<" ./stl_files/*.stl"<<endl;
  	return 0;
  }

  Mesh mesh;
  vector <pair<char* , int >> testFails;
  vector<stl_plane> planes;
  planes.push_back(stl_plane(1, 0 ,0 ,0));
  planes.push_back(stl_plane(0, 1 ,0 ,0));
  planes.push_back(stl_plane(0, 0 ,1 ,0));
  planes.push_back(stl_plane(1, 0 ,0 ,1));
  planes.push_back(stl_plane(0, 1 ,0 ,1));
  planes.push_back(stl_plane(0, 0 ,1 ,1));
  planes.push_back(stl_plane(1, 0 ,0 ,-1));
  planes.push_back(stl_plane(0, 1 ,0 ,-1));
  planes.push_back(stl_plane(0, 0 ,1 ,-1));
  planes.push_back(stl_plane(1, 1 ,0 ,0));
  planes.push_back(stl_plane(1, 0 ,0.5 ,0));
  planes.push_back(stl_plane(0, 1 ,1 ,0));

  bool unitTestResult = mesh.runUnitTests();

  if(argc > 1)
  	for (int i = 1; i < argc; ++i)
  	{
  		cout << "Testing file: "<<argv[i]<<endl;
  		for (int j = 0; j < planes.size(); ++j)
  		{
  		  cout << planeToString(planes[j])<<endl;
  		  mesh.openStl(argv[i]);
  		  if (mesh.cut(planes[j]))
  	    {
  	      meshIntegrationTests(mesh, j, argv[i],testFails);
  	    }
  		}
  	}
  
  if(unitTestResult)
    cout<<endl<<"All unit tests \e[32mcompleted without a problem\e[0m."<<endl;
  else
    cout<<endl<<"Some unit tests \e[31mfailed\e[0m."<<endl;

  if(argc > 1)
  	if(testFails.size()!=0)
  	  for (int i = 0; i < testFails.size(); ++i)
  	    cout << "Tests of: "<< (string)testFails[i].first << " using plane: "<<planeToString( planes[testFails[i].second] )<<" \e[31mfailed\e[0m."<<endl;
  	else
  	  cout << "All intergration tests \e[32mcompleted without a problem\e[0m."<<endl;

  return 0;
}