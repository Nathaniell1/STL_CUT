#include "stlcut.h"
using namespace p2t;
using namespace std;

double eps = 1e-25;
double minEps = numeric_limits<double>::max();
char removedAxis = 'z';

#ifndef SEGV
#define SEGV
jmp_buf buf;
sigset_t signal_set;
int numOfSegv=0;
#endif


  stl_plane::stl_plane(float x, float y, float z,float d) 
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = -1*d;
  }
bool operator == (const stl_vertex& a, const stl_vertex& b) 
{
  if(a.x == b.x && a.y == b.y && a.z == b.z)
    return true;
  return false;
}


/*
*Compares two vertices for equality (with tolerance).
* Comparision is made in 2D, 2 out of 3 points are selected based on removed axis.
*@param [in] a First vertex
*@param [in] b Second vertex
*/
bool vertexEqual (stl_vertex a, stl_vertex b) 
{
  float tolerance = eps;
  double tol1 = ABS(a.x-b.x);
  double tol2 = ABS(a.y-b.y);
  double tol3 = ABS(a.z-b.z);
  if(removedAxis == 'x') { tol1 = tol3; a.x = a.z; b.x = b.z; } 
  if(removedAxis == 'y') { tol2 = tol3; a.y = a.z; b.y = b.z; }
  double tmp = tol1>tol2 ? tol1:tol2;
  if(tmp < minEps) minEps = tmp;
  return (ABS((float)(a.x-b.x))<tolerance && ABS((float)(a.y-b.y))<tolerance  );
}

/*
*Checks for duplicit vertices in border. Duplicits are erased.
*/
void Mesh::checkDuplicity()
{
   for (unsigned int i = 0; i < border.size()-2; i+=2)
  {
    for (unsigned  int j = i+2; j < border.size(); j+=2)
    {
      stl_vertex tmp1 = border[j];
      stl_vertex tmp2 = border[j+1];
      stl_vertex tmp3 = border[i];
      stl_vertex tmp4 = border[i+1];
      if( (tmp1 == tmp3 && tmp2 == tmp4 ) || (tmp2 == tmp3 && tmp1 == tmp4) )
      {
        border.erase(border.begin()+j,border.begin()+(j+2));
        j-=2;
      }
    }
  }
}
/*
* Function to help handle segmentation faults of poly2tri.
*/
void segv_handler(int s)
{
    switch(s)
    {
        case SIGSEGV:
        cerr<<"\nSegmentation fault signal caught! Attempting recovery.."<<endl<<endl;
        numOfSegv++;
        longjmp(buf, numOfSegv);
        break;
    }
}
/*
*Function to make segfault recovery possible.
*Sets function to call in case on segfault and unmasks segfault signal.
*/
void segvInit()
{
  sigemptyset(&signal_set);
  sigaddset(&signal_set, SIGSEGV); 
  sigprocmask(SIG_UNBLOCK, &signal_set, NULL); //clearnig segfault signal
  signal(SIGSEGV, segv_handler);
}

/*
* Comparator of two poly2tri points.
*/
bool comparatorStruct::operator() (const p2t::Point* lhs, const p2t::Point* rhs) const
{
  if(lhs->x < rhs->x) return true;
  else 
  {
    if(lhs->x != rhs->x) return false;
    if(lhs->y <  rhs->y) return true;
  } 
  return false;
}

/*
* Comparator of two poly2tri points.
*/
bool comparatorStruct::operator() (const p2t::Point lhs, const p2t::Point rhs) const
{
  if(lhs.x < rhs.x) return true;
  else 
  {
    if(lhs.x != rhs.x) return false;
    if(lhs.y <  rhs.y) return true;
  } 
  return false;
}

bool comparatorStruct::operator()(const stl_facet & a, const stl_facet & b) const
{
  vector<stl_vertex> v1,v2; 
  for (int i = 0; i < 3; ++i)
  {
    v1.push_back(a.vertex[i]);
    v2.push_back(b.vertex[i]);
  }

  sort(v1.begin(),v1.end(),comparatorStruct());
  sort(v2.begin(),v2.end(),comparatorStruct());
  if(v1 == v2)
    return true;
  return false;
}


/*
* Compares vertexes based on removed axis. No tolerance used.
* This comparator serves for sorting purpouses.
*/
  bool comparatorStruct::operator() (const stl_vertex& lhs, const stl_vertex& rhs) const
  {

    if(removedAxis == 'z')
    {
      if(lhs.x < rhs.x) return true;
      else 
        {
          if(lhs.x != rhs.x) return false;
          if(lhs.y < rhs.y) return true;
          else 
          {
            if(lhs.y != rhs.y) return false;
            if(lhs.z < rhs.z) return true;
          }
        } 
      return false;
    }

    if(removedAxis == 'x')
    {
      if(lhs.y < rhs.y) return true;
      else 
        {
          if(lhs.y != rhs.y) return false;
          if(lhs.z <  rhs.z) return true;
          else 
          {
            if(lhs.z != rhs.z) return false;
            if(lhs.x <  rhs.x) return true;
           }
        } 
      return false;
    }

    if(removedAxis == 'y')
    {
      if(lhs.x < rhs.x) return true;
      else 
        {
          if(lhs.x != rhs.x) return false;
          if(lhs.z <  rhs.z) return true;
          else 
          {
            if(lhs.z != rhs.z) return false;
            if(lhs.y <  rhs.y) return true;
          }
        } 
      return false;
    }

  }

/*
* Returns missing coordinate when reverting points from 2d back to 3d. Missing coordinate is first calculated, and then its real value is found in original vertices.
* Original 3D point is returned.
*/

stl_vertex Mesh::getMissingCoordinate(const p2t::Point* a)
{
  stl_vertex b;
  if(removedAxis == 'x')
  {
    b.y = a->y;
    b.z = a->x;
    b.x = (plane.y*b.y + plane.z*b.z + plane.d) / ((-1.0) * plane.x);
  }
  if(removedAxis == 'y')
  {
    b.x = a->x;
    b.z = a->y;
    b.y = (plane.x*b.x+plane.z*b.z+plane.d) / ((-1.0) * plane.y);
  }
  if(removedAxis == 'z')
  {
    b.x = a->x;
    b.y = a->y;
    b.z = (plane.x*b.x + plane.y*b.y + plane.d) / ((-1.0) * plane.z);
  }

  set<stl_vertex>::iterator itlow,itup;
  itlow = originalVertices.lower_bound (b);              
  itup  = originalVertices.upper_bound (b);  

  if(itlow != itup) 
      b = *itlow;
  else
  {
    if(removedAxis == 'z' && b.x == itlow->x && b.y == itlow->y || removedAxis == 'x' && b.z == itlow->z && b.y == itlow->y || removedAxis == 'y' && b.x == itlow->x && b.z == itlow->z   )
      b = *itlow;
    else 
    {
      if(itlow-- != originalVertices.begin())
        if(removedAxis == 'z' && b.x == itlow->x && b.y == itlow->y || removedAxis == 'x' && b.z == itlow->z && b.y == itlow->y || removedAxis == 'y' && b.x == itlow->x && b.z == itlow->z   )
          b = *itlow;
    }
  }

  return b;
}


/*
*Swap coordinates based on the one we ignore and push it to polylines
*/
void Mesh::pushBackToPolylines(vector<p2t::Point*> &vec,stl_vertex vert)
{

  if(removedAxis == 'x')
  {
    vert.x = vert.z;
  }
  if(removedAxis == 'y')
  {
    vert.y = vert.z;
  }
  if(vec.size() == 0 ||  !(vec.back()->x == vert.x && vec.back()->y == vert.y))
    vec.push_back(new p2t::Point(vert.x,vert.y));
}

/*
* Closes the mesh file.
*/
void Mesh::close()
{
  stl_close(&meshFile);
}

/*
* Self-explenatory
*/
double Mesh::calculatePolygonArea(vector<p2t::Point*> polygon)
{
    int n = polygon.size();
    int j = 0;
    double area = 0.0;
    //shoelace algorithm to calculate polygon area
    for (int i = 0; i < n; ++i)
    {
      j = (i+1)%n;
      area += polygon[i]->x * polygon[j]->y;
      area -= polygon[j]->x * polygon[i]->y;
    }
    area = abs(area) / 2.0;
    return area;
}


/*
*Calculates if 3 point are in counter clockwise order. Used in finding non-simple polygons.
*/
bool Mesh::ccw(p2t::Point* a, p2t::Point* b, p2t::Point* c)
{
    return ( (c->y - a->y) * (b->x - a->x) > (b->y - a->y) * (c->x - a->x) );
}

/*
* Calculates if two edges intersect.
*/
bool Mesh::edgesIntersect (p2t::Point* a, p2t::Point* b, p2t::Point* c, p2t::Point* d)
{
    return ( ccw(a, c, d) != ccw(b, c, d) && ccw(a, b, c) != ccw(a, b, d) );
}


/*
* Tests if given vertex is inside given polygon .
* Used to find holes in polygons.
*/
bool Mesh::vertexInPolygon( const vector<Point* >& polygon,  const double &testx, const double &testy)
{
  int i, j;
  bool in = false;
  //vector from point to the right(to infinity) vs edges
  for (i = 0, j = polygon.size()-1; i < polygon.size(); j = i++) 
  {
    if ( ((polygon.at(i)->y > testy) != (polygon.at(j)->y > testy)) &&
     (testx < (polygon.at(j)->x - polygon.at(i)->x) * (testy - polygon.at(i)->y) / (polygon.at(j)->y - polygon.at(i)->y) + polygon.at(i)->x) )
    {
       in = !in;
    }
  }
  return in;
}

/*
* Sorts polylines (polygons) in descending order based on their area and removes polygons with area of 0;
*/
void Mesh::sortPolylines()
{
  vector<double> polygonArea;
  polygonArea.resize(polylines.size());
  for (unsigned  int i = 0; i < polylines.size(); ++i)
  {
    polygonArea[i] = calculatePolygonArea(polylines[i]);    
  }

  vector<int> areaOrder;
  for (unsigned int i = 0; i < polylines.size(); ++i)
  {
    areaOrder.push_back(i);
  }
  sort(areaOrder.begin(), areaOrder.end(), [&polygonArea](const int &a, const int &b)->bool{return polygonArea[a] < polygonArea[b];});
  //AreaOrdernow contains sorted indexes of polygons based on their area
  
  vector<vector <p2t::Point*> > tmpPolylines;
  tmpPolylines.resize(polylines.size());
  for (unsigned  int i = 0; i < polylines.size(); ++i)
  {
    tmpPolylines[i] = polylines[areaOrder[i]];
  }
  vector<int> indexes;
  for (unsigned int i = 0; i < tmpPolylines.size(); ++i)
  {
    if(polygonArea[areaOrder[i]] == 0)
    {
      indexes.push_back(i);
    }
  }
  sort(indexes.begin(),indexes.end());
  for (int i = indexes.size()-1; i >= 0; --i)
  {
    for (unsigned  int n = 0; n < tmpPolylines[indexes[i]].size(); ++n)
    {
      delete tmpPolylines[indexes[i]][n];
    }
    //erasing wrong polygons with area = 0
    tmpPolylines.erase(tmpPolylines.begin()+indexes[i]);
  }
  polylines = tmpPolylines;
}

/*
* This method determines if polygons in polylines are polygons or holes and distributes them in the designated container.
* When this method returns polygonsWithHoles contains those polygons and holes.
*/
void Mesh::findHoles()
{
  vector<p2t::Point*> tmpPolygon = polylines.back();
  polylines.pop_back();
  pair<vector<p2t::Point*>,int > tmpPair = make_pair(tmpPolygon,-1);
  vector<pair<vector<p2t::Point*> , int>> tmpVecPair;
  tmpVecPair.push_back(tmpPair);
  polygonsWithHoles.push_back(tmpVecPair);
  int in;
  int out;
  int holeIn = -1; //-1 for polygon
  bool placeFound = false;
  int pos = 0;
  while(polylines.size() > 0)
  {
    tmpPolygon = polylines.back(); 
    polylines.pop_back();
    pos = 0;
    placeFound = false;
    holeIn = -1; 
    for(unsigned  int k = 0;k < polygonsWithHoles.size();k++)
    {
      in = 0;
      out = 0; 
      //test 3 points and based on majority decide if this polygons is inside biggerPolygon
      for(int j = 0;j < 3;j+= 1)
      {
        vector<p2t::Point*>  biggerPolygon = polygonsWithHoles[k][pos].first;
        double a = tmpPolygon[j]->x; 
        double b = tmpPolygon[j]->y;
        if(vertexInPolygon(biggerPolygon, a, b))
          in++;
        else
          out++;
      }
      if(in > out || placeFound)
      {
        if(in > out) 
        {
          if(holeIn != -1)
            holeIn   = -1;  //swaping between polygon
          else
            holeIn   = pos; //and hole (pos = in which polygon is this one as a hole)
        }  
        placeFound = true;
        //as long as we didnt get to the end of the vector, continue - basicaly
        if(pos < polygonsWithHoles[k].size()-1) // we dont want to test it vs itself
        {
          pos++;
          k--; // we need to stay at the same index, 
        }
        else //we get to the end, make a pair and push it to the end
        {
          tmpPair = make_pair(tmpPolygon, holeIn);
          polygonsWithHoles[k].push_back(tmpPair);
          break;
        }
      }
    }
    if(!placeFound)          // if we get here and point wasnt in
    {                        // we found a new polygon
      tmpPair = make_pair(tmpPolygon,-1);
      tmpVecPair.erase(tmpVecPair.begin(), tmpVecPair.end());
      tmpVecPair.push_back(tmpPair);
      polygonsWithHoles.push_back(tmpVecPair);
    }
  }
}

/*
*This method tryes to fix non-simple polygons, but in very basic way.
*Its designed to remove falsely made intersecting edges due to floating points calcuation errors.
*Thats the reason why it always removes first point of second edge.. because this points is usualy close to the second point of first edge which causes intersection.
*/
void Mesh::removeNonsimplePolygonPoints(vector<p2t::Point*> & p)
{
  bool pointRemoved;
  do
  {
    pointRemoved = false;
    for (int m = 0; m < p.size()-1; m+=1)
    {
      for (int j = m+2; j < p.size()-1; j+=1) //I assume that edges which share one point wont intersect - to prevent false intersection results
      {
        if(edgesIntersect(p[m], p[m+1], p[j], p[j+1] ))
        {
          p.erase(p.begin()+j,p.begin()+j+1);
          j--;
          pointRemoved = true;
          continue;
        }
      }
    }

    for (int k = 1; k < p.size()-2; k+=1)
    {
      if(edgesIntersect(p[k], p[k+1], p[0], p[ p.size()-1] ))
      {
        p.erase(p.begin()+p.size()-1,p.begin()+p.size());
        k--;
        pointRemoved=true;
        continue;
      }
    }
  }while(pointRemoved == true);


}


/*
*Tests every polygon if its non-simple, if not, some points are removed in order to (try to) make it a simple polygon.
*/
void Mesh::repairIfNonsimplePolygon()
{
  for (int i = 0; i < polygonsWithHoles.size(); ++i) 
  {
    for (int j = 0; j < polygonsWithHoles[i].size(); ++j)
    {
      removeNonsimplePolygonPoints(polygonsWithHoles[i][j].first );
    }
  }
}

/*
* Searches through newly triagulated polygon and tryes to find triangles with points that did not exist before triangulation 
* (poly2tri sometime make them when provided with nonsimple polygon). If those tringles are found, they are deleted.
* Its one of the many just-in-case-something-went-wrong methods. Can help to make succesful cut through model with small errors.
*/
void Mesh::checkPoly2triResult( vector<p2t::Triangle*>& triangles )
{
  set<p2t::Point,comparatorStruct> originalBorderPoints;
  for (std::set<stl_vertex>::iterator i = originalVertices.begin(); i != originalVertices.end(); ++i) 
  {
    double x = (*i).x;
    double y = (*i).y;
    if(removedAxis == 'x')
      x = (*i).z;
    if(removedAxis == 'y')
      y = (*i).z;
    p2t::Point p = Point(x,y);
    originalBorderPoints.insert(p);
  }

  for (int i = 0; i < triangles.size(); ++i)
  {
    for(int j = 0; j<3; j++)
    {
      p2t::Point* tmp = triangles[i]->GetPoint(j);
      p2t::Point p = Point(tmp->x,tmp->y);
      if( originalBorderPoints.count(p) == 0)// test if 2D point is in the original border created by cut
      {
        // if its not, erase the triangle
        triangles.erase(triangles.begin()+i,triangles.begin()+i+1);
        i--;
        break;
      }
    }
  }
}

/*
* Clean up method for polygonsWithHoles.
*/
void Mesh::deletePolygonsWithHoles()
{
  for (int i = 0; i < polygonsWithHoles.size(); ++i)
  {
    for (int j = 0; j < polygonsWithHoles[i].size(); ++j)
    {
      for (int k = 0; k < polygonsWithHoles[i][j].first.size(); ++k)
      {
        delete polygonsWithHoles[i][j].first[k];
      }
    }
    polygonsWithHoles[i].clear();
  }
  polygonsWithHoles.clear(); 
}

/*
* Triangulates all identified polygons with its holes.
*/
void Mesh::triangulateCut(int topOrBot)
{
  map<int, p2t::CDT*> polygons;
  repairIfNonsimplePolygon();
  for (unsigned int i = 0; i < polygonsWithHoles.size(); ++i)
  {
    for (unsigned  int j = 0; j < polygonsWithHoles[i].size(); ++j)
    {
      if(polygonsWithHoles[i][j].second == -1) // we found polygon, not a hole
      {   
        p2t::CDT* tmp = new  p2t::CDT(polygonsWithHoles[i][j].first);
        polygons.insert ( pair<int,p2t::CDT*>(j,tmp) );
      }
      else // we found a hole
      {
        int holeIn = polygonsWithHoles[i][j].second;
        map<int, p2t::CDT*>::iterator it;
        it = polygons.find(holeIn);
        it->second->AddHole(polygonsWithHoles[i][j].first);
      }
    }
    //all polygons processed and ready to triangulation
    if(polygons.size() > 0)
      for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
      {
        k->second->Triangulate();
        vector<p2t::Triangle*> triangles = k->second->GetTriangles();
        checkPoly2triResult(triangles);
        createFacets(triangles, topOrBot);   
      }
    for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
    {
      delete (*k).second;
    }
    polygons.erase(polygons.begin(),polygons.end()); 
    
  }
  deletePolygonsWithHoles();  
}

/*
* Makes facets in 3D from triangles in 2D provided by poly2tri triangulation.
*/
void Mesh::createFacets(vector<p2t::Triangle*> &triangles, int side)
{
  // for each triangle, create facet
  for (vector<p2t::Triangle*>::iterator i = triangles.begin(); i != triangles.end(); i++) 
  {
    stl_vertex vertex;
    stl_facet facet;
    for (size_t j = 0; j < 3; j++) 
    {
      p2t::Point* p = (*i)->GetPoint(j);
      facet.vertex[j] = getMissingCoordinate(p);
    }
    float test_norm[3];
    test_norm[0] = plane.x;
    test_norm[1] = plane.y;
    test_norm[2] = plane.z;
    stl_normalize_vector(test_norm);
    facet.normal.x = test_norm[0];
    facet.normal.y = test_norm[1];
    facet.normal.z = test_norm[2];
    if(side == 0 || side == -1)
      botFacets.push_back(facet);
    //reverse normals
    facet.normal.x *= -1.0;
    facet.normal.y *= -1.0;
    facet.normal.z *= -1.0;
    vertex = facet.vertex[1];
    facet.vertex[1] = facet.vertex[2];
    facet.vertex[2] = vertex;
    if(side == 0 || side == 1)
      topFacets.push_back(facet);
  }
   
}

/*
* Insert stl_vertex to the front of polylines. Used to clear up code. Used during find for continious edges after the cut.
*/
void Mesh::pushFrontToPolylines(vector<p2t::Point*> &vec,stl_vertex vert)
{

  if(removedAxis == 'x')
  {
    vert.x = vert.z;
  }
  if(removedAxis == 'y')
  {
    vert.y = vert.z;
  }
  if(vec.size() == 0 ||  !(vec.front()->x == vert.x && vec.front()->y == vert.y))
    vec.insert(vec.begin(),new p2t::Point(vert.x,vert.y));
}

/*
* Pops from border into the two provided stl_vertex variables.
*/
void Mesh::popTo(stl_vertex& a, stl_vertex& b)
{
  a = border.back();
  border.pop_back();
  b = border.back();
  border.pop_back();
}

/*
* Creates facet. Used to clean up code. Creates new triangles during cut.
*/
stl_facet Mesh::createFacet(stl_facet facet, int s, int i, stl_vertex intersect)
{
  stl_facet tmp_facet = facet;
  tmp_facet.vertex[0] = facet.vertex[s];
  tmp_facet.vertex[1] = facet.vertex[(s+i)%3];
  tmp_facet.vertex[2] = intersect;
  return tmp_facet;
}

/*
* Creates facet. Used to clean up code. Creates new triangles during cut.
*/
stl_facet Mesh::createFacet(stl_facet facet,int s, int i, stl_vertex intersect1, stl_vertex intersect2)
{
  stl_facet tmp_facet = facet;
  switch(i)
  {
    case 0:
    {
      tmp_facet.vertex[0] = facet.vertex[s];
      tmp_facet.vertex[1] = intersect1;
      tmp_facet.vertex[2] = intersect2;
      return tmp_facet;
    }
    case 1:
    {
      tmp_facet.vertex[0] = intersect1;
      tmp_facet.vertex[1] = facet.vertex[(s+1)%3];
      tmp_facet.vertex[2] = intersect2;
      return tmp_facet;
    }
    case 2:
    {
      tmp_facet.vertex[0] = intersect2;
      tmp_facet.vertex[1] = facet.vertex[(s+1)%3];
      tmp_facet.vertex[2] = facet.vertex[(s+2)%3];
      return tmp_facet;
    }
  }
}

/*
* Normalizes plane normal and sets up plane variable.
*/
void Mesh::setPlane(stl_plane plane)
{
  this->plane = plane;
  float norm[3];
  norm[0] = plane.x;
  norm[1] = plane.y;
  norm[2] = plane.z;
  stl_normalize_vector(norm);
  this->plane.x = norm[0];
  this->plane.y = norm[1];
  this->plane.z = norm[2];
}
/*
* Sets removed axis based on the plane normal vector.
*/
void Mesh::setRemovedAxis()
{
  removedAxis = 'z';
  if(ABS(plane.x) >= ABS(plane.y) && ABS(plane.x) >= ABS(plane.z) )
    removedAxis = 'x';
  if(ABS(plane.y) >= ABS(plane.x) && ABS(plane.y) >= ABS(plane.z) )
    removedAxis = 'y';
}


void Mesh::setVertexFromFacet(stl_vertex& a, stl_vertex& b,const int &s,const stl_facet & facet)
{
  a = facet.vertex[(s+1)%3];
  b = facet.vertex[(s+2)%3];
}

/* 
* Calculates if two facets have equal edge. Facets provided have both at least one edge "on" the cuting plane.
* Used during prorocessing facets that have whole edge on the plane. In case of comparing facets that have all 3 points on, the stl_vertexes provided in
* the touple are irelevant.
*@param [in] facet1 Its a tuple. stl_facet its facet. stl_position is "above/below" if 2 points are "on" and one is not, or "on" when three points are on the cuting plane. stl_vertices represent edge on the cuting plane.
*@param [in] facet2
*/
bool Mesh::haveEqualEdges(tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet1, tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet2)
{

  //ifs facet is  "on" , all 3 points were on and we have to try to match all 3 edges of facet1 to 3 edges of facet2 
  if( (get<1>(facet1) == on && get<1>(facet2)!= on))
  {
    pair<stl_vertex,stl_vertex> a,b;
    for (int i = 0; i < 3; ++i)
    {
      a = make_pair(get<0>(facet1).vertex[i] ,       get<0>(facet1).vertex[(i+1)%3]); 
      b = make_pair(get<0>(facet1).vertex[(i+1)%3] , get<0>(facet1).vertex[i]); 
      if(a == (make_pair(get<2>(facet2),get<3>(facet2))) || b == (make_pair(get<2>(facet2),get<3>(facet2))))
        return true;
    }
  }
  if( (get<1>(facet2) == on && get<1>(facet1)!= on))
  {
    pair<stl_vertex,stl_vertex> a,b;
    for (int i = 0; i < 3; ++i)
    {
      a = make_pair(get<0>(facet2).vertex[i] ,       get<0>(facet2).vertex[(i+1)%3]); 
      b = make_pair(get<0>(facet2).vertex[(i+1)%3] , get<0>(facet2).vertex[i]); 
      if(a == make_pair(get<2>(facet1),get<3>(facet1)) || b == make_pair(get<2>(facet1),get<3>(facet1)))
        return true;
    }
  }
  // tests if its a same edge
  if( ((get<2>(facet1) == get<2>(facet2)) && (get<3>(facet1) == get<3>(facet2))) || ((get<3>(facet1) == get<2>(facet2) && get<2>(facet1) == get<3>(facet2))))     
    return true;

  return false;
}
/*
* Inserts and puhes back vertex x and y. Used to clean up code from reocurring parts.
*/
void Mesh::insertTo(stl_vertex x, stl_vertex y, vector<stl_vertex>& a, set<stl_vertex,comparatorStruct> & b)
{
  a.push_back(x);
  a.push_back(y);
  b.insert(x);
  b.insert(y);
}

/*
* After the cut, this method identifies polylines/polygons made during cut.
* In case of edges from original mode  on cutting plane, this method recursively calls itself via processOnBorder.
*@param [in] firstCall Signalise if this method is called from cut method or processOnBorder method via recursion.
*/
bool Mesh::createBorderPolylines(bool firstCall)
{
  if(firstCall == true && processOnFacets() == true)
    if(processOnBorder() == true)
      return false; // this means everything worked well, but we dont want to triangulate again, because processOnBorder already did that.

  if(border.size() == 0 )
  {
    topFacets.clear();
    botFacets.clear();
    if(!silent) cerr<<"Nothing to cut"<<endl;
    return false;
  }

  checkDuplicity();
  stl_vertex cont,end,tmp1,tmp2;
  popTo(cont,end);
  polylines.resize(20);
  numOfPolylines  = 0;
  eps             = 1e-24;
  minEps          = numeric_limits<double>::max(); 
  bool found      = true;
  pushBackToPolylines(polylines[numOfPolylines],end);
  pushBackToPolylines(polylines[numOfPolylines],cont);

  while(border.size()!=0)
  { 
    if(!found) 
      eps = minEps*1.005;
    if(eps > 0.25) // ignoring edges, if it gets to this, there was probably a problem with mesh
    {
      eps = 1e-24;
      minEps = numeric_limits<double>::max();
      popTo(cont,end);
      if(!polylines[numOfPolylines].empty()) 
        numOfPolylines++;
      pushBackToPolylines(polylines[numOfPolylines],end);
      pushBackToPolylines(polylines[numOfPolylines],cont);// if we didnt found it with 0.1 tolerance.. we will triangulate what we found 
    }
    for (int i = border.size()-1;i >= 0; i-= 2)
    {
      found = false;
      tmp1  = border[i]; 
      tmp2  = border[i-1];
      if(vertexEqual(tmp1, cont) || vertexEqual(tmp2, cont) || vertexEqual(tmp1, end) || vertexEqual(tmp2, end))// ve found next vertex in polyline
      {
        found = true;
        if(vertexEqual(tmp1,cont)) 
        {            
          pushBackToPolylines(polylines[numOfPolylines], tmp2);
          cont = tmp2;
          border.erase(border.begin()+(i-1), border.begin()+i+1); 
        }  
        else if(vertexEqual(tmp2,cont))
        {
          pushBackToPolylines(polylines[numOfPolylines], tmp1);
          cont = tmp1;
          border.erase(border.begin()+(i-1), border.begin()+i+1); 
        }
        else if(vertexEqual(tmp1, end))
        {
          pushFrontToPolylines(polylines[numOfPolylines], tmp2);
          end = tmp2;
          border.erase(border.begin()+(i-1), border.begin()+i+1); 
        }
        else if(vertexEqual(tmp2, end))
        {
          pushFrontToPolylines(polylines[numOfPolylines], tmp1);
          end = tmp1;
          border.erase(border.begin()+(i-1), border.begin()+i+1); 
        }

        if(vertexEqual(cont,end)) //after we found next point, we have to check if another point is an end
        {
          if(!polylines[numOfPolylines].empty()) 
          {
            delete polylines[numOfPolylines][polylines[numOfPolylines].size()-1];
            polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
          }
          if(border.size() > 0)
          {
            numOfPolylines++;
            if(numOfPolylines+5 > polylines.size())
              polylines.resize(numOfPolylines*10);
            popTo(end, cont);
            pushBackToPolylines(polylines[numOfPolylines], end); //start of new polyline
            pushBackToPolylines(polylines[numOfPolylines], cont);
          }
          minEps = numeric_limits<double>::max();
          break;
        }
        else 
        { 
          minEps = numeric_limits<double>::max();
          eps = 1e-24;
          break;
        }       
      }
    }
  }
  polylines.resize(numOfPolylines+1);
  sortPolylines();
  if (polylines.size() == 0)
    return false;
  return true;
}

/*
* Called if processOnFacets returns true. This method is called in special case that there are edges that were exactly on the cuting plane.
* In situation like this, its necessary to triangulate twice with addition of top-related and bot-related edges to the border.
* Returns false if cut was made through side of a model (and every triangle in on one side of the plane or on the plane).
* Returns true otherwise.
*/
bool Mesh::processOnBorder()
{
  if (!(topFacets.size()!=0 && botFacets.size() != 0 ))
    return false;
  vector<stl_vertex> borderBackUp = border;
  border.insert(border.end(), botBorder.begin(), botBorder.end());
  if(createBorderPolylines(false))
  {
    findHoles();
    triangulateCut(-1);
  }

  border = borderBackUp;
  border.insert(border.end(), topBorder.begin(), topBorder.end());
  if(createBorderPolylines(false))
  {
    findHoles();
    triangulateCut(1);
  }
  return true;
}

/*
* Sorts facets with edge on the cutting plane to the corresponding containers. 
* Returns true if there are facet pairs with one facet completely on the cutting plane and second with two points on and one below or above.
*/
bool Mesh::processOnFacets()
{
  for (int i = 0; i < facetsOnPlane.size(); ++i)
  {
    for (int j = i+1; j < facetsOnPlane.size(); ++j)
    {
     if( haveEqualEdges(facetsOnPlane[i],facetsOnPlane[j]) )
      {
        auto k = get<1>(facetsOnPlane[i]); auto l = get<1>(facetsOnPlane[j]);
        if( (k == below && l == above) || (k == above && l == below) )
          insertTo(get<2> (facetsOnPlane[i]), get<3> (facetsOnPlane[i]), border, originalVertices);
        if((k == on && l == below) || (l == on && k == below))
        {
          int pos = (k == on) ? pos=j:pos=i;
          insertTo(get<2>(facetsOnPlane[pos]), get<3>(facetsOnPlane[pos]), botBorder, originalVertices);
        }
         if((k == on && l == above) || (l == on && k == above))
        {
          int pos = (k == on) ? pos=j:pos=i;
          insertTo(get<2>(facetsOnPlane[pos]), get<3>(facetsOnPlane[pos]), topBorder, originalVertices);
        }
      }
    }
  }

  if( (botBorder.size() != 0 || topBorder.size() != 0))
    return true;
  return false;
}

/*
* Pushes facet to bot or top container accoring to its position to cuting plane.
*/
void Mesh::pushAboveBelow(const int aboves,stl_vertex& a,stl_vertex& b,const stl_facet &facet, const stl_position* pos)
{
  if(aboves == 1)
    for (int s = 0; s < 3; ++s)
      {
        if(pos[s] == above) 
        {
          setVertexFromFacet(a,b,s,facet);
          stl_vertex intersect1 = intersection(facet.vertex[s],a);
          stl_vertex intersect2 = intersection(facet.vertex[s],b);
          topFacets.push_back( createFacet(facet,s,0,intersect1,intersect2) ); // facet with above vertex of triangle and intersecting vertices added
          botFacets.push_back( createFacet(facet,s,1,intersect1,intersect2) ); // facet with intersecting vertices and below vertex 1            
          botFacets.push_back( createFacet(facet,s,2,intersect1,intersect2) ); // facet with below vertices and intersecting vertex 2
          border.push_back(intersect1);
          border.push_back(intersect2);
          break;
        }  
      }
    if(aboves == 2)
      for (int s = 0; s < 3; ++s)
      {
        if(pos[s] == below) 
        {
          setVertexFromFacet(a,b,s,facet);
          stl_vertex intersect1 = intersection(facet.vertex[s],a);
          stl_vertex intersect2 = intersection(facet.vertex[s],b);
          botFacets.push_back( createFacet(facet,s,0,intersect1,intersect2) );  // facet with above vertex and intersecting vertices added
          topFacets.push_back( createFacet(facet,s,1,intersect1,intersect2) );  // facet with intersecting vertices and below vertex 1  
          topFacets.push_back( createFacet(facet,s,2,intersect1,intersect2) );  // facet with below vertices and intersecting vertex 2
          border.push_back(intersect1);
          border.push_back(intersect2);
          break;
        }
      }   
}

/*
* Pushes vertices to border nd facets to bot and top border.
* It processes facets with one points on one above and one below or two points on and one above/below.
*/
void Mesh::pushOns(const int ons,stl_vertex& a,stl_vertex& b,const stl_facet &facet, const stl_position* pos)
{
  if(ons == 1)
  {
    for (int s = 0; s < 3; ++s)
    {
      if(pos[s] == on)
      {
        setVertexFromFacet(a,b,s,facet);
        stl_vertex intersect = intersection(a,b);
        if(pos[(s+1)%3] == below) 
        { 
          botFacets.push_back(createFacet(facet,s,1,intersect));
          topFacets.push_back(createFacet(facet,s,2,intersect));
          border.push_back(facet.vertex[s]);
          border.push_back(intersect);
        }
        else // pos[(s+1)%3] == above
        {
          botFacets.push_back(createFacet(facet,s,2,intersect));
          topFacets.push_back(createFacet(facet,s,1,intersect));
          border.push_back(facet.vertex[s]);
          border.push_back(intersect);
        } 
        break;
      }
    }
  }
  else 
  if(ons == 2)
  {
    for (int n = 0; n < 3; ++n)
    { 
      if(pos[n] == below)  
      { 
        facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));
        botFacets.push_back(facet); break;   
      }
      if(pos[n] == above) 
      { 
        facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));      
        topFacets.push_back(facet); break;
      }
    }
  }
}



/*
* Cleans up private variables to make sure its possible to call cut again on same instance of Mesh.
*/
void Mesh::cleanupVariables()
{
  botFacets.clear();
  topFacets.clear();
  originalVertices.clear();
  border.clear();
  polylines.clear();
  polygonsWithHoles.clear();
  facetsOnPlane.clear();
  numOfPolylines = 0;
  botBorder.clear();
  topBorder.clear();
}

/*
* Sets up options to run in silent mode and to disable recovery from error.
*/
void Mesh::setOptions(bool sil, bool err)
{
  errorRecovery = err; 
  silent = sil;
}

/*
* Makes cut through stl_file. Throws exception if input file wasnt provided before calling this method. Returns true if succesfull otherwise false.
*@param [in] plane The stl_plane used to cut the mesh.
*/
bool Mesh::cut(stl_plane plane)
{
  if(stl_get_error(&meshFile) != 0)
    throw std::runtime_error("Mesh to cut wasn't provided.");

  double error_correction = 0.00015;
  numOfSegv = 0;
  if(errorRecovery != false)
  { 
    setjmp(buf);
    segvInit();
  }
  if(numOfSegv > 6)
  {
    if(!silent) cerr<<"STLCUT wasnt able to made this cut. Try changing the plane position slightly and make sure that your model is 2-manifold."<<endl;
    return false;
  }
  else
  {
    if(numOfSegv > 0)
    {
      plane.d = plane.d + error_correction;
      error_correction = error_correction > 0?(-1)*error_correction:error_correction*10;
      if(!silent) cerr<<endl<<"Recovered from segmentation fault."<<endl;
    }

    cleanupVariables(); // makes using cut multiple times in row possible
    setPlane(plane);
    divideFacets();
    if(createBorderPolylines())
    {
      findHoles();
      triangulateCut();
    }
    
    if(topFacets.size() != 0 && botFacets.size() != 0)
      return true;
    return false;
  }
}

/*
* Divides facets below above and on cutting plane into respective containers.
*/
void Mesh::divideFacets()
{
  setRemovedAxis();
  size_t aboves = 0;
  size_t belows = 0;
  size_t ons = 0;
  stl_position pos[3];
  stl_vertex a,b;
  for (size_t k = 0; k < meshFile.stats.number_of_facets; k++)
  {
    stl_facet facet = meshFile.facet_start[k];
    aboves = belows = ons = 0;
    for (int i = 0; i < 3; ++i)
    {
      pos[i] = vertexPosition(facet.vertex[i]);
      if(pos[i] == above) aboves++;
      if(pos[i] == on)    ons++;
      if(pos[i] == below) belows++;
    }

    if(aboves == 3 || (aboves == 2 && ons == 1))
    {
      topFacets.push_back(facet);
      continue;
    }
    if(belows == 3 || (belows == 2 && ons == 1))
    {
      botFacets.push_back(facet);
      continue;
    }
    if(ons == 3)
    { // last 2 vertices in this tuple are not important in case of ons == 3
      facetsOnPlane.push_back(make_tuple(facet,pos[0],facet.vertex[0],facet.vertex[1]));
      continue;
    }
    if(ons == 1 && belows == 1 && aboves == 1)
    {
      pushOns(1,a,b,facet,pos);
      continue;
    }
    if(ons == 2)
    {
      pushOns(2,a,b,facet,pos);
      continue;
    }

    if(aboves == 1 && belows == 2) // last possibility... the plane cuts the triangle and doesnt intersect with any (already given) vertex,
    {
      pushAboveBelow(1,a,b,facet,pos);  
    }
    else // below == 1, above == 2
    {
      pushAboveBelow(2,a,b,facet,pos);
    }
  }
  for (unsigned int i = 0; i < border.size(); ++i)
  {
    originalVertices.insert(border[i]);
  }
}

/*
* Returns the position of the vertex relative to the plane (on, below, above)
@param [in] vertex 
*/
stl_position Mesh::vertexPosition(stl_vertex vertex) 
{
  double result = plane.x*vertex.x + plane.y*vertex.y + plane.z*vertex.z + plane.d;
  if (result > 0) return above;
  if (result < 0) return below;
  return on;
}

/*
* Calculates intersection between edge and plane. Returns vertex where intersection occurs.
@param [in] a First vertex of the edge
@param [in] b Second vertex of the edge
*/
stl_vertex Mesh::intersection(stl_vertex a, stl_vertex b) 
{
  stl_vertex tmp = b;
  if(a.x < b.x)
    {b = a; a = tmp;}
  else
    if(a.x == b.x && a.y < b.y)
      {b = a; a = tmp;}
    else
      if(a.x == b.x && a.y == b.y && a.z < b.z)
        {
          b = a; a = tmp;
        }

  stl_vector ab; // vector from A to B
  ab.x = b.x - a.x;
  ab.y = b.y - a.y;
  ab.z = b.z - a.z;
  double t = - (a.x*plane.x + a.y*plane.y + a.z*plane.z + plane.d) / (ab.x*plane.x + ab.y*plane.y + ab.z*plane.z);
  stl_vertex result;
  result.x = a.x + ab.x*t;
  result.y = a.y + ab.y*t;
  result.z = a.z + ab.z*t;
  return result;
}

/*
* Sets stl file which will be cut. Throws runtime error if file is in error state when provided.
*@param [in] file Takes stl_file type.
*/
void Mesh::setStl(stl_file file)
{
  meshFile = file;
  if(stl_get_error(&meshFile) != 0)
    throw std::runtime_error("Provided Stl file is in error state.");
}

/*
* Sets stl file which will be cut. Throws runtime error if file cannot be opened or there is other problem with it.
*@param [in] file Takes name of stl_file.
*/
void Mesh::openStl(char * name)
{
  stl_open(&meshFile, name);
  if(stl_get_error(&meshFile) != 0)
    throw std::runtime_error("Can't open file.");
}

/*
* Returns new meshes (top or bottom part) as pointer to stl_file.
*/
stl_file* Mesh::getExportedStl(deque<stl_facet> facets) 
{
  stl_file* stl_out = new stl_file;
  initializeStl(stl_out,facets.size());
  int first = 1;
  for (deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++) 
  {
    stl_out->facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(stl_out, *facet, first);
    first = 0;
  }
  stl_repair(stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//reverse normal directions
  return stl_out;
}

/*
*Tests if provided string is valid - Can contains alhanumeric characters, underscore and space.
*/
bool Mesh::isStringValid(const std::string &str)
{
    return find_if(str.begin(), str.end(), 
        [](char c) { return !(isalnum(c) || (c == ' ') || (c == '_' )); }) == str.end();
}

/*
*Saves 2 files. with name_1.stl and name_2.stl. If provided name is invalid (or is not provded) it uses Cut_Mesh_1.stl and Cut_Mesh_2.stl
*Throws an runtime_error exception if file cannot be created or saved.
*@param [in] name The name of new saved files (optional). Name can only contain alhanumeric characters, underscore and space.
*/
void Mesh::save(string name)
{
  if(! (name != "" && isStringValid(name) ) )
  {
    name = "Cut_Mesh";
  }
  exportStl(topFacets,(name + "_1.stl").c_str());
  exportStl(botFacets,(name + "_2.stl").c_str());
  if(!silent)cout<<"Files saved to "<<name<<"_1.stl and "<<name<<"_2.stl"<<endl;
}

/*
* Returns array of 2 stl_file pointers to new cut meshes.
*/
std::array<stl_file*,2> Mesh::getFinalStls()
{
  return{getExportedStl(topFacets),getExportedStl(botFacets)};
}

void Mesh::initializeStl(stl_file * stl,int numOfFacets)
{
  stl->stats.type = inmemory;
  stl->stats.number_of_facets = numOfFacets;
  stl->stats.original_num_facets = stl->stats.number_of_facets;
  stl->v_indices = NULL;
  stl->v_shared = NULL;
  stl->neighbors_start = NULL;
  stl_clear_error(stl);
  stl_allocate(stl);
  stl->stats.degenerate_facets = 0;
  stl->stats.edges_fixed = 0;
  stl->stats.facets_removed = 0;
  stl->stats.facets_added = 0;
  stl->stats.facets_reversed = 0;
  stl->stats.backwards_edges = 0;
  stl->stats.normals_fixed = 0;
}

/*
*Creates new file and exports STL to it.
*Throws runtime_error if fail can't be created.
*/
void Mesh::exportStl(deque<stl_facet> facets, const char* name) 
{
  stl_file stl_out;
  initializeStl(&stl_out,facets.size());
  
  int first = 1;
  for (deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++)
  {
    stl_out.facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(&stl_out, *facet, first);
    first = 0;
  }
  
  stl_write_ascii(&stl_out, name, "stlcut");
  if(stl_get_error(&stl_out) != 0)
    throw std::runtime_error("Can't create new file.");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

/*
* Return stl_vertex from x,y,z input.
*/
stl_vertex Mesh::getVertex(double x, double y, double z)
{
  stl_vertex a;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}


