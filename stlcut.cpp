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
  float tolerance = eps;//1e-25;//1e-15;
  double tol1 = ABS(a.x-b.x);
  double tol2 = ABS(a.y-b.y);
  double tol3 = ABS(a.z-b.z);
  if(removedAxis == 'x') { tol1 = tol3; a.x = a.z; b.x = b.z; } 
  if(removedAxis == 'y') { tol2 = tol3; a.y = a.z; b.y = b.z; }
  //TODO NEMELO By to resit i treti bod? pouzivam to na hledani navazujicich hran, co kdyz mam krychli s 2 hranama co jsou stejne az na zanedbanou hranu
  //prooomyslet.. asi je to ale blbost
  double tmp = tol1>tol2 ? tol1:tol2;
  if(tmp < minEps) minEps = tmp;
  return (ABS((float)(a.x-b.x))<tolerance && ABS((float)(a.y-b.y))<tolerance  );//&& ABS((float)(a.z-b.z))<tolerance);
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
      //TODO FIX with new == operator in mind
      if( (tmp1 == tmp3 && tmp2 == tmp4 ) || (tmp2 == tmp3 && tmp1 == tmp4) )
      //if ((tmp1.x == tmp3.x && tmp1.y == tmp3.y &&tmp1.z == tmp3.z && tmp2.x == tmp4.x && tmp2.y == tmp4.y &&  tmp2.z == tmp4.z) || ( tmp2.x == tmp3.x && tmp2.y == tmp3.y &&tmp2.z == tmp3.z && tmp1.x == tmp4.x && tmp1.y == tmp4.y && tmp1.z == tmp4.z ) )
      {
        border.erase(border.begin()+j,border.begin()+(j+2));
        j-=2;
        //cout<<"Mazu duplicitu"<<endl;
       // break;
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
  
/*
    if(a.vertex[0] == b.vertex[j])
      if(a.vertex[1] == b.vertex[(j+1)%3]);
        if(a.vertex[2] == b.vertex[(j+2)%3])
          return true;
  
  return false;*/
  
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
  /*
  bool x,y,z;
  x=y=z=false;
  for (int j = 0; j < 3; ++j)
  {
    if(a.vertex[0] == b.vertex[j])
      x = true;
    if(a.vertex[1] == b.vertex[j])
      y = true;
    if(a.vertex[2] == b.vertex[j])
      z = true;
  }
  return(x && y && z);
*/
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
stl_vertex Mesh::getMissingCoordinate(const p2t::Point* a)
{
  stl_vertex b;
 // cout<<"x: "<<a->x<<" z: "<<a->y<<"  calculated to:"<<endl;
  if(removedAxis == 'x')
  {
    b.y = a->y;
    b.z = a->x;
    b.x= (plane.y*b.y+plane.z*b.z+plane.d) / ((-1.0)* plane.x);
  }
  if(removedAxis == 'y')
  {
    b.x = a->x;
    b.z = a->y;
    b.y= (plane.x*b.x+plane.z*b.z+plane.d) / ((-1.0)* plane.y);
  }
  if(removedAxis == 'z')
  {
    b.x = a->x;
    b.y = a->y;
    //b.z = 0;
    b.z= (plane.x*b.x+plane.y*b.y+plane.d) / ((-1.0)* plane.z);
  }
  cout<<"zac"<<endl;

  set<stl_vertex>::iterator itlow,itup;
  itlow = originalVertices.lower_bound (b);              
  itup = originalVertices.upper_bound (b);   
  itup--;
  if (itlow == itup||itlow == originalVertices.begin()) // we found it
  {
    cout << "z borderu"<<endl;
    b = (*itlow);
  }
  else
  {
    itup++;
    itlow--;
    double dif1 = ABS(max(  (max(ABS((*itlow).x-b.x),ABS((*itlow).y-b.y))) , ABS(((*itlow).z-b.z) )));
    double dif2 = ABS(max(  (max(ABS((*itup).x-b.x) ,ABS((*itup).y-b.y ))) , ABS(((*itup).z-b.z )) ));
    if(dif1<dif2 && itlow!=originalVertices.end() && dif1<0.00001) // find the closest
    {
      cout << "z borderu"<<endl;
      b = (*itlow); 
     
    }
    else 
    {
      if(dif2 <0.00001)     
        {
          cout << "z borderu"<<endl;
          b = (*itup);

        }
    }
  }

  //cout<<"x: "<<b.x<<" z: "<<b.z<<endl;
  return b;
}
*/

/*
* Returns missing coordinate when reverting points from 2d back to 3d. Missing coordinate is first calculated, and then its real value is found in original vertices.
* Original 3D point is returned.
*/

stl_vertex Mesh::getMissingCoordinate(const p2t::Point* a)
{
  stl_vertex b;
 // cout<<"x: "<<a->x<<" z: "<<a->y<<"  calculated to:"<<endl;
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
      b = *(--itlow);
  }


  return b;
}


void Mesh::volumeTest()
{
    stl_calculate_volume(&meshFile);
    double org_volume = meshFile.stats.volume;
    double volume = 0.0;
    stl_file cut_mesh;
    stl_open(&cut_mesh, (char*)"Cut_Mesh_1.stl");
    stl_calculate_volume(&cut_mesh);
    volume += meshFile.stats.volume;
    stl_close(&cut_mesh);
    stl_open(&cut_mesh, (char*)"Cut_Mesh_2.stl");
    stl_calculate_volume(&cut_mesh);
    volume += cut_mesh.stats.volume;
    stl_close(&cut_mesh);
    if(abs(org_volume-volume) > abs(org_volume/1000.0)) 
      cerr<<"Volume test failed, cut might be wrong! "<<endl;
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
bool Mesh::vertexInPolygon( const vector<Point* >& polygon,  const double &testx, const double testy)
{
  int i, j;
  bool in = false;
  //vector from point to the right(to infinity) vs edges
  for (i = 0, j = polygon.size()-1; i < polygon.size(); j = i++) {
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
  //cout<<"Polyline size:"<<polylines.size()<<endl;
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
    //if(polygonArea[areaOrder[i]]!=0)
    {
      tmpPolylines[i] = polylines[areaOrder[i]];
    }  
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
  polygonsWithHoles.push_back(tmpVecPair);//tmplist);
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
      for(int j = 0;j < 3;j+= 1)//ceil(tmpPolygon.size()/3.0))
      {
        vector<p2t::Point*>  biggerPolygon = polygonsWithHoles[k][pos].first;//.front();
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
  }
  polygonsWithHoles.clear(); 
}

/*
* Triangulates all identified polygons with its holes.
*/
void Mesh::triangulateCut(int topOrBot)
{
  map<int, p2t::CDT*> polygons;//vector<p2t::CDT>> polygons;
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
      bot_facets.push_back(facet);
    //reverse normals
    facet.normal.x *= -1.0;
    facet.normal.y *= -1.0;
    facet.normal.z *= -1.0;
    vertex = facet.vertex[1];
    facet.vertex[1] = facet.vertex[2];
    facet.vertex[2] = vertex;
    if(side == 0 || side == 1)
      top_facets.push_back(facet);
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
bool Mesh::createBorderPolylines(bool firstCall)//bool processOnFac)
{
  if(firstCall == true && processOnFacets() == true)
    if(processOnBorder() == true)
      return false; // this means everything worked well, but we dont want to triangulate again, because processOnBorder already did that.

  if(border.size() == 0 )
  {
    if(!silent)cerr<<"Nothing to cut"<<endl;
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
    if(eps > 0.1) // ignoring edges, if it gets to this, there was probably a problem with mesh
    {
      eps = 1e-24;
      minEps = numeric_limits<double>::max();
      //if(processOnFac == true) 
        //std::cerr<<"Unable to find a connected edge, mesh might be invalid"<<endl;
      popTo(cont,end);
      if(!polylines[numOfPolylines].empty()) 
        numOfPolylines++;
      pushBackToPolylines(polylines[numOfPolylines],end);
      pushBackToPolylines(polylines[numOfPolylines],cont);// if we didnt found it with 0.1 tolerance.. we will triangulate what we found 
    }
    for (int i = border.size()-1;i >= 0; i-=2)
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
  if (!(top_facets.size()!=0 && bot_facets.size() != 0 ))
    return false;

  vector<stl_vertex> borderBackUp = border;
  border.insert(border.end(),botBorder.begin(),botBorder.end());
  if(createBorderPolylines(false))
  {
    findHoles();
    triangulateCut(-1);
  }

  border = borderBackUp;
  border.insert(border.end(),topBorder.begin(),topBorder.end());
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
   // udelat privatni
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
          /*Fail safe.. if both intersection were same points... floating point calculation failsafe*/

          top_facets.push_back( createFacet(facet,s,0,intersect1,intersect2) ); // facet with above vertex of triangle and intersecting vertices added
          bot_facets.push_back( createFacet(facet,s,1,intersect1,intersect2) ); // facet with intersecting vertices and below vertex 1            
          bot_facets.push_back( createFacet(facet,s,2,intersect1,intersect2) ); // facet with below vertices and intersecting vertex 2
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
          /*stl_vertex a,b;
          a = facet.vertex[(s+1)%3];
          b = facet.vertex[(s+2)%3];*/
          setVertexFromFacet(a,b,s,facet);
          stl_vertex intersect1 = intersection(facet.vertex[s],a);
          stl_vertex intersect2 = intersection(facet.vertex[s],b);
          bot_facets.push_back( createFacet(facet,s,0,intersect1,intersect2) );  // facet with above vertex and intersecting vertices added
          top_facets.push_back( createFacet(facet,s,1,intersect1,intersect2) );  // facet with intersecting vertices and below vertex 1  
          top_facets.push_back( createFacet(facet,s,2,intersect1,intersect2) );  // facet with below vertices and intersecting vertex 2
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
          bot_facets.push_back(createFacet(facet,s,1,intersect));
          top_facets.push_back(createFacet(facet,s,2,intersect));
          border.push_back(facet.vertex[s]);
          border.push_back(intersect);
        }
        else // pos[(s+1)%3] == above
        {
          bot_facets.push_back(createFacet(facet,s,2,intersect));
          top_facets.push_back(createFacet(facet,s,1,intersect));
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
        bot_facets.push_back(facet); break;   
      }
      if(pos[n] == above) 
      { 
        facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));      
        top_facets.push_back(facet); break;
      }
    }
  }
}



/*
* Cleans up private variables to make sure its possible to call cut again on same instance of Mesh.
*/
void Mesh::cleanupVariables()
{
  bot_facets.clear();
  top_facets.clear();
  originalVertices.clear();
  border.clear();
  polylines.clear();
  polygonsWithHoles.clear();
  facetsOnPlane.clear();
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
* Makes cut through stl_file. Returns true if succesfull otherwise false.
*@param [in] plane The stl_plane used to cut the mesh.
*/
bool Mesh::cut(stl_plane plane)
{
  double error_correction = 0.00015;
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

    cleanupVariables(); // with this its possible to use cut multiple times in row
    setPlane(plane);
    divideFacets();
    if(createBorderPolylines())
    {
      findHoles();
      triangulateCut();
      //return true;
    }
    
    cout<< (string)(t_minMaxPointsSameAfterCut() == true? "true": "false")<<endl;
    cout << (string)(t_sameVolume() == true? "true":"false")<<endl;
    cout << (string)(t_allTrianglesOriginalOrPartOfCut() == true? "true": "false")<<endl; 
    cout << (string)(t_noVertexOnOpositeSide() == true?"true": "false")<<endl; 
    if(top_facets.size() != 0 && bot_facets.size() != 0)
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
      top_facets.push_back(facet);
      continue;
    }
    if(belows == 3 || (belows == 2 && ons == 1))
    {
      bot_facets.push_back(facet);
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
  //stl_exit_on_error(&meshFile);
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
   //todo pridat repair co opravuje diry
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
  exportStl(top_facets,(name + "_1.stl").c_str());//"Cut_Mesh_1.stl");//"pokus1.stl");
  exportStl(bot_facets,(name + "_2.stl").c_str());//"Cut_Mesh_2.stl");//"pokus2.stl");
  if(!silent)cout<<"Files saved to "<<name<<"_1.stl and "<<name<<"_2.stl"<<endl;
  /*
  string name = "";
  if(!acquireSaveName(name))
    name="Cut_Mesh";
  exportStl(top_facets,(name+"_1.stl").c_str());//"Cut_Mesh_1.stl");//"pokus1.stl");
  exportStl(bot_facets,(name+"_2.stl").c_str());//"Cut_Mesh_2.stl");//"pokus2.stl");
  cout<<"Files saved to "<<name<<"_1.stl and "<<name<<"_2.stl"<<endl;*/
}

/*
* Returns array of 2 stl_file pointers to new cut meshes.
*/
std::array<stl_file*,2> Mesh::getFinalStls()
{
  return{getExportedStl(top_facets),getExportedStl(bot_facets)};

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
  
  //stl_repair(&stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//fix normal directions
  //stl_repair(&stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,0 ,0,0);//fix normal directions
  stl_write_ascii(&stl_out, name, "stlcut");
  if(stl_get_error(&stl_out) != 0)
    throw std::runtime_error("Can't create new file.");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

//used by ADMeshGUI
std::array<stl_file*,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & success)
{
  std::array<stl_file*,2> cutMesh; 
	stl_plane plane = stl_plane(a,b,c,d);
	Mesh mesh;
  mesh.setStl(*stlMesh);
  if(mesh.cut(plane))
  {
    cutMesh = mesh.getFinalStls();
    success = true;
  }
  else
    success = false;
  /*
                               
  if(mesh.createBorderPolylines())
    {
      mesh.findHoles();
      mesh.triangulateCut();
      cutMesh = mesh.getFinalStls();
      succes = true;
    }
    else succes = false;*/
  
  return{cutMesh[0],cutMesh[1]};
}
stl_vertex Mesh::getVertex(double x, double y, double z)
{
  stl_vertex a;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
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
  if(!(vert == getVertex(0, 0, 0))) //Test 1
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 0, 2), getVertex(1, 0, -2));
  if(!(vert == getVertex(1, 0, 0))) //Test 2
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 1, 1), getVertex(-1, -1, -1));
  if(!(vert == getVertex(0, 0, 0))) //Test 3
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(1, 1, 1), getVertex(0, 0, 0));
  if(!(vert == getVertex(0, 0, 0))) //Test 4
    fails.push_back(num);
  num++;

  setPlane(stl_plane(1, 0, 0, 5));
  setRemovedAxis();

  vert = intersection(getVertex(6, 0, 0), getVertex(1, 0, 0));
  if(!(vert == getVertex(5, 0, 0))) //Test 5
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(7, 1, 2), getVertex(-3, 5, -2));
  if(!(vert == getVertex(5, 1.8, 1.2))) //Test 6
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(5, 0, 2), getVertex(10, 3, 5));
  if(!(vert == getVertex(5, 0, 2))) //Test 7
    fails.push_back(num);
  num++;

  setPlane(stl_plane(0, 0.5, 0, -5));
  setRemovedAxis();

  vert = intersection(getVertex(3, 0, -2.5), getVertex(11, 0.8, 22));
   //cout << vert.x<<" "<<vert.y<<" "<<vert.z<<endl;
  if(!(vert == getVertex(-47, -5, -155.625))) //Test 8
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(7, 1, 2), getVertex(-3, 5, -2));
  if(!(vert == getVertex(22, -5, 8))) //Test 9
    fails.push_back(num);
  num++;

  setPlane(stl_plane(1, 1, 0.5, -5));
  setRemovedAxis();

  vert = intersection(getVertex(6, 0, 0), getVertex(-1, 0, 5));
  if(!(vert == getVertex(-15, 0, 15))) //Test 10
    fails.push_back(num);
  num++;
  vert = intersection(getVertex(7, 7, 7), getVertex(-7, -7, -7));
  if(!(vert == getVertex(-3, -3, -3))) //Test 11
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
    fails.push_back(num);  //Test 1
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
    fails.push_back(num); //Test 2
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
    fails.push_back(num); //Test 3

  //cleanup
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
  //stl_position pos;
  cout<<"Testing vertexPosition"<<endl;
  setPlane(stl_plane(1, 1, -1, 5));
  setRemovedAxis();

  if(!(vertexPosition(getVertex(0,0,10)) == below))
    fails.push_back(num);
  num++;
  if(!(vertexPosition(getVertex(0,0,-10)) == above))
    fails.push_back(num);
  num++;
  //cout << "V:" << vertexPosition(getVertex(8.66,8.66,-8.66))<<endl; //above 0, below 2, on 1

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
  //stl_plane plane = stl_plane(0 , 0 , 1 , 0);
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
  //stl_vertex getMissingCoordinate(const p2t::Point* a);
  //pushbacktopoly
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
  //cout<<"vys: "<<x.x<<" "<<x.y<<" "<<x.z<<endl;
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
  int num=1;
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

//stl_facet Mesh::createFacet(stl_facet facet, int s, int i, stl_vertex intersect)
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
  //tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet1, tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet2)
  stl_facet f1, f2;
  f1.vertex[0] = getVertex(0  , 10,   5);
  f1.vertex[1] = getVertex(0  ,  0,  5);
  f1.vertex[2] = getVertex(20, 1 , 5);
  f2.vertex[0] = getVertex(0  , 10,   5);
  f2.vertex[1] = getVertex(-5  ,  12,  6);
  f2.vertex[2] = getVertex(0, -0 , 5);
  tuple<stl_facet, stl_position, stl_vertex, stl_vertex> t1, t2;
  t1 = make_tuple(f1,on,f1.vertex[0],f1.vertex[2]);
  t2 = make_tuple(f2,above,f2.vertex[0],f2.vertex[2]);

  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,on,f1.vertex[1],f1.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,on,f1.vertex[1],f1.vertex[0]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t2 = make_tuple(f2,above,f2.vertex[0],f2.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  f1.vertex[2] = getVertex(20, 1 , 1);
  t1 = make_tuple(f1,below,f1.vertex[0],f1.vertex[1]);
  if( haveEqualEdges( t1,t2 ) != true)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,below,f1.vertex[1],f1.vertex[2]);
  if( haveEqualEdges( t1,t2 ) != false)
    fails.push_back(num);
  num++;

  t1 = make_tuple(f1,below,f1.vertex[2],f1.vertex[0]);
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
//this method is designed to handle
bool Mesh::t_removeNonsimplePolygonPoints()
{
  fails.clear();
  int num=1;
  cout<<"Testing removeNonsimplePolygonPoints"<<endl;
  //void Mesh::removeNonsimplePolygonPoints(vector<p2t::Point*> & p)
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
   //return ( (c->y - a->y) * (b->x - a->x) > (b->y - a->y) * (c->x - a->x) );
  // 0 - 0 * 0 - 0 > 1 - 0 * 1-0

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
  if(!(polylines[0].size() == 1 && polylines[0][0]->x == 1.1f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 2 && polylines[0][1]->x == 1.1f && polylines[0][1]->y ==-2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3 && polylines[0][2]->x == -1.1f && polylines[0][2]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (1, 0, 0, 0));
  setRemovedAxis();

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 4 && polylines[0][3]->x == 3.3f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 5 && polylines[0][4]->x == 3.3f && polylines[0][1]->y ==-2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 6 && polylines[0][5]->x == 3.3f && polylines[0][2]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 6))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (0, 1, 0, 0));
  setRemovedAxis();

  pushBackToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 7 && polylines[0][6]->x == 1.1f && polylines[0][6]->y ==3.3f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(1.1,-2.2,-3.3));
  if(!(polylines[0].size() == 8 && polylines[0][7]->x == 1.1f && polylines[0][7]->y ==-3.3f))
    fails.push_back(num);
  num++;

  pushBackToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 9 && polylines[0][8]->x == -1.1f && polylines[0][8]->y ==3.3f))
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
  if(!(polylines[0].size() == 1 && polylines[0][0]->x == 1.1f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 2 && polylines[0][0]->x == 1.1f && polylines[0][0]->y ==-2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3 && polylines[0][0]->x == -1.1f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 3))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (1, 0, 0, 0));
  setRemovedAxis();

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 4 && polylines[0][0]->x == 3.3f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,3.3));
  if(!(polylines[0].size() == 5 && polylines[0][0]->x == 3.3f && polylines[0][0]->y ==-2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 6 && polylines[0][0]->x == 3.3f && polylines[0][0]->y ==2.2f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 6))
    fails.push_back(num);
  num++;

  setPlane(stl_plane (0, 1, 0, 0));
  setRemovedAxis();

  pushFrontToPolylines(polylines[0], getVertex(1.1,2.2,3.3));
  if(!(polylines[0].size() == 7 && polylines[0][0]->x == 1.1f && polylines[0][0]->y ==3.3f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(1.1,-2.2,-3.3));
  if(!(polylines[0].size() == 8 && polylines[0][0]->x == 1.1f && polylines[0][0]->y ==-3.3f))
    fails.push_back(num);
  num++;

  pushFrontToPolylines(polylines[0], getVertex(-1.1,2.2,3.3));
  if(!(polylines[0].size() == 9 && polylines[0][0]->x == -1.1f && polylines[0][0]->y ==3.3f))
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
//return succes;
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
    cout << "vol1: "<<org_volume<<endl;
    cout << "vol2: "<<volume<<endl;
    if(abs(volume) <= org_volume * 1.001 && abs(volume) >= org_volume * 0.999)
      return true;
    return false;

}
/*
bool triangleArea(stl_facet &a)
}
    int n = 3;
    int j = 0;
    double area = 0.0;
    //shoelace algorithm to calculate polygon area
    for (int i = 0; i < n; ++i)
    {
      j = (i+1)%n;
      area += a.vertex[i]->x * a.vertex[j]->y;
      area -= a.vertex[j]->x * a.vertexi]->y;
    }
    area = abs(area) / 2.0;
    return area;
}
*/

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
      auto pos = vertexPosition(f);//(topMesh->facet_start[i].vertex[j]); 
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
      auto pos = vertexPosition(f);//(topMesh->facet_start[i].vertex[j]); 
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
    //cout << botMeshFacets[i].vertex[0].x<< " "<<botMeshFacets[i].vertex[0].y<< " "<<botMeshFacets[i].vertex[0].z<< " "<<endl;
    //cout << it->vertex[0].x<< " "<<it->vertex[0].y<< " "<<it->vertex[0].z<< " "<<endl;
    while(it->vertex[0].x == botMeshFacets[i].vertex[0].x) //dokud jsme v ramci tech serazenych
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
        //if we  couldnt find the facet, lets test if one of its points is in the
        if(find(originalVertices.begin(),originalVertices.end(), botMeshFacets[i].vertex[k]) != originalVertices.end())
          ons++;
      }
      if(ons >= 1)
        continue;
      else return false;
    }
  }

  for (int i = 0; i < topMeshFacets.size(); ++i)
  {
    bool found = false;
    auto it = lower_bound(originalFacets.begin(),originalFacets.end(),topMeshFacets[i] ,[](const stl_facet &a, const stl_facet &b)->bool{return a.vertex[0].x < b.vertex[0].x;});
    while(it->vertex[0].x == topMeshFacets[i].vertex[0].x) //dokud jsme v ramci tech serazenych
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
        //if we  couldnt find the facet, lets test if one of its points is in the
        if(find(originalVertices.begin(),originalVertices.end(), topMeshFacets[i].vertex[k]) != originalVertices.end())
          ons++;
      }
      if(ons >= 1)
        continue;
      else return false;
    }
  }

}





/*
bool Mesh::t_allTrianglesOriginalOrPartOfCut()
{
  //setPlane(stl_plane(0, 0 ,0 ,1));
  //setRemovedAxis();

  vector<stl_facet> originalFacets;
  //originalFacets.push_back(meshFile.facet_start[0]);
  comparatorStruct comp;
  bool flag;
  cout << "jdu na prvni for"<<endl;
  for (int i = 0; i < meshFile.stats.number_of_facets; ++i)
  {
    originalFacets.push_back(meshFile.facet_start[i]);
  }

  cout << "prvni for prosel"<<endl;

  auto meshes = getFinalStls();
  auto topMesh = meshes[0];
  stl_get_size(topMesh);
  auto botMesh = meshes[1];
  stl_get_size(botMesh);

  cout <<topMesh->stats.number_of_facets<<endl;
  for (int i = 0; i < topMesh->stats.number_of_facets; ++i)
  {
    cout << i<<endl;
    flag = false;
    for (int j = 0; j < originalFacets.size(); ++j)
    {
      if(comp( topMesh->facet_start[i], originalFacets[j] ) == true)
      {
        flag = true;
        break;
      }
    }
    if(flag == false)
    {   
       int ons = 0;
       for (int k = 0; k < 3; ++k)
       {
         if(find(originalVertices.begin(),originalVertices.end(), topMesh->facet_start[i].vertex[k]) != originalVertices.end())
          ons++;
       }
       if(ons >= 1)
        continue;
       else return false;
    }
  }
  cout << "druhy for prosel"<<endl;
  for (int i = 0; i < botMesh->stats.number_of_facets; ++i)
  {
    flag = false;
    for (int j = 0; j < originalFacets.size(); ++j)
    {
      if(comp( botMesh->facet_start[i], originalFacets[j] ) == false)
      {
        flag = true;
        break;
      }
    }
    if(flag == false)
    {   
       int ons = 0;;
       for (int k = 0; k < 3; ++k)
       {
         if(find(originalVertices.begin(),originalVertices.end(), botMesh->facet_start[i].vertex[k]) != originalVertices.end())
          ons++;
       }
       if(ons >= 1)
        continue;
       else return false;
    }

  }
  cout << "tret for prosel"<<endl;

  return true;

}*/

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
  cout << maxOriginal.x<<" "<<maxOriginal.y<<" "<<maxOriginal.z<<endl;
  cout << minOriginal.x<<" "<<minOriginal.y<<" "<<minOriginal.z<<endl;
  cout << maxNew.x<<" "<<maxNew.y<<" "<<maxNew.z<<endl;
  cout << minNew.x<<" "<<minNew.y<<" "<<minNew.z<<endl;

  if(!(maxOriginal == maxNew && minOriginal == minNew))
    return false;
  return true;
}


bool Mesh::runUnitTests()
{
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
  
  
  //if(!(t_minMaxPointsSame())) neni to unit test
    //success = false;

  return success; 
}
  
