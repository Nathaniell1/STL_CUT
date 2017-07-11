#include "stlcut.h"
using namespace p2t;
using namespace std;

double eps = 1e-25;
double minEps = numeric_limits<double>::max();
char removedAxis = 'z';


  stl_plane::stl_plane(float x, float y, float z,float d) 
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = -1*d;
  }
bool operator == (const stl_vertex& a, const stl_vertex& b) 
{
  //cout<<"use"<<endl;
  if(a.x == b.x && a.y == b.y && a.z == b.z)
    return true;
  return false;
}

// " == "
  /*
bool operator == (stl_vertex a, stl_vertex b) {
  float tolerance = eps;//1e-25;//1e-15;
  double tol1 = ABS(a.x-b.x);
  double tol2 = ABS(a.y-b.y);
  double tol3 = ABS(a.z-b.z);
  if(removedAxis == 'x') { tol1 = tol3; a.x = a.z; b.x = b.z; }
  if(removedAxis == 'y') { tol2 = tol3; a.y = a.z; b.y = b.z; }
  double tmp = tol1>tol2 ? tol1:tol2;
  if(tmp < minEps) minEps = tmp;
  return (ABS((float)(a.x-b.x))<tolerance && ABS((float)(a.y-b.y))<tolerance  );//&& ABS((float)(a.z-b.z))<tolerance);
}*/

/*Comparing function which also manipulates with eps nad minEps - this is used to increase tolerance during search for connected edges*/
  /*
bool vertexEqual (stl_vertex a, stl_vertex b) 
{
  float tolerance = eps;//1e-25;//1e-15;
  double tol1 = ABS(a.x-b.x);
  double tol2 = ABS(a.y-b.y);
  double tol3 = ABS(a.z-b.z);

  double tmp = tol1>tol2 ? tol1:tol2;
  if(tmp < minEps) minEps = tmp;
  return (ABS((float)(a.x-b.x))<tolerance && ABS((float)(a.y-b.y))<tolerance  );//&& ABS((float)(a.z-b.z))<tolerance);
}
*/
/*we need to check all 3 coordinates for comparing border, if we didnt, cube with same coordinates except for the one we erased would have 2 same edges*/
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
*Checks duplicit vertexes in border
*/
void checkDuplicity(vector<stl_vertex> &border)
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
      if ((tmp1.x == tmp3.x && tmp1.y == tmp3.y &&tmp1.z == tmp3.z && tmp2.x == tmp4.x && tmp2.y == tmp4.y &&  tmp2.z == tmp4.z) || ( tmp2.x == tmp3.x && tmp2.y == tmp3.y &&tmp2.z == tmp3.z && tmp1.x == tmp4.x && tmp1.y == tmp4.y && tmp1.z == tmp4.z ) )
      {
        border.erase(border.begin()+j,border.begin()+(j+2));
        j-=2;
        cout<<"Mazu duplicitu"<<endl;
       // break;
      }
    }
  }
}

bool setVertComp::operator() (const p2t::Point* lhs, const p2t::Point* rhs) const
{
  if(lhs->x < rhs->x) return true;
  else 
  {
    if(lhs->x != rhs->x) return false;
    if(lhs->y <  rhs->y) return true;
  } 
  return false;
}

bool setVertComp::operator() (const p2t::Point lhs, const p2t::Point rhs) const
{
  if(lhs.x < rhs.x) return true;
  else 
  {
    if(lhs.x != rhs.x) return false;
    if(lhs.y <  rhs.y) return true;
  } 
  return false;
}

/*
* Compares vertexes based of removed axis
*/
  bool setVertComp::operator() (const stl_vertex& lhs, const stl_vertex& rhs) const
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



Mesh::~Mesh()
{
  /*
  for (int i = 0; i < polylines.size(); ++i)
  {
    for (int j = 0; j < polylines[i].size(); ++j)
    {
      delete polylines[i][j];
    }
  }
*/

}

/*
*Calculates missing coordinate and then tryes to find it in the border to remove possible inaccuracy
*/
//TODO ODKOMENTOVAT
void Mesh::setMissingCoordinate(const p2t::Point* a,stl_vertex &b)
{
  
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


  set<stl_vertex>::iterator itlow,itup;
  itlow = border2.lower_bound (b);              
  itup = border2.upper_bound (b);   
  itup--;
  if (itlow == itup||itlow == border2.begin()) // we found it
  {
    b = (*itlow);
  }
  else
  {
    itup++;
    itlow--;
    double dif1 = ABS(max(  (max(ABS((*itlow).x-b.x),ABS((*itlow).y-b.y))) , ABS(((*itlow).z-b.z) )));
    double dif2 = ABS(max(  (max(ABS((*itup).x-b.x) ,ABS((*itup).y-b.y ))) , ABS(((*itup).z-b.z )) ));
    if(dif1<dif2 && itlow!=border2.end() && dif1<0.00001) // find the closest
    {
      b = (*itlow); 
     
    }
    else 
    {
      if(dif2 <0.00001)     
        {
          b = (*itup);

        }
    }
  }

  //cout<<"x: "<<b.x<<" z: "<<b.z<<endl;
  
}


void Mesh::volumeTest()
{
    
    stl_calculate_volume(&mesh_file);
    double org_volume = mesh_file.stats.volume;
    double volume = 0.0;
    stl_file cut_mesh;
    stl_open(&cut_mesh, (char*)"Cut_Mesh_1.stl");
    stl_calculate_volume(&cut_mesh);
    volume+=mesh_file.stats.volume;
    stl_close(&cut_mesh);
    stl_open(&cut_mesh, (char*)"Cut_Mesh_2.stl");
    stl_calculate_volume(&cut_mesh);
    volume+=cut_mesh.stats.volume;
    stl_close(&cut_mesh);
    if(abs(org_volume-volume) > abs(org_volume/1000.0)) 
      cerr<<"Volume test failed, cut might be wrong! "<<endl;
}
/*
*Swap coordinates based on the one we ignore and push it to polylines
*/
void Mesh::pushToPolylines(vector<p2t::Point*> &vec,stl_vertex vert)
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
*Self-explenatory
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
*Tests if given vertex is inside given polygon 
* used to find holes
*/

/*
*Calculates if 3 point are in counter clockwise order
*/
bool Mesh::ccw(p2t::Point* a, p2t::Point* b, p2t::Point* c)
{
    return ( (c->y - a->y) * (b->x - a->x) > (b->y - a->y) * (c->x - a->x) );
}

/*
* Calculates if two edges intersect
*/
bool Mesh::intersect (p2t::Point* a, p2t::Point* b, p2t::Point* c, p2t::Point* d)
{
    return ( ccw(a,c,d) != ccw(b,c,d) && ccw(a,b,c) != ccw(a,b,d) );
}



bool Mesh::vertexInPolygon(int nvert, const vector<Point* >& vertex,  const double &testx, const double testy)
{
  int i, j;
  bool in = false;
  //vector from point to the right(to infinity) vs edges
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((vertex.at(i)->y > testy) != (vertex.at(j)->y > testy)) &&
     (testx < (vertex.at(j)->x - vertex.at(i)->x) * (testy-vertex.at(i)->y) / (vertex.at(j)->y - vertex.at(i)->y) + vertex.at(i)->x) )
    {
       in = !in;
    }
  }
  return in;
}

void Mesh::sortPolylines()
{
  cout<<"Polyline size:"<<polylines.size()<<endl;
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
  sort(areaOrder.begin(),areaOrder.end(),[&polygonArea](const int &a, const int &b)->bool{return polygonArea[a]<polygonArea[b];});
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

void Mesh::findHoles()
{

  //////////////////////
  ///////////////////////
  //////////////--------------------------------------------------VYMAZAT--------------------------//////////////
  //////////////////////////
  /////////////////////////////
  //////////////////////////
  //TODO
  /*
  if(polylines.size() == 0) 
  {
    vector<p2t::Point*> xx;
    xx.push_back(new p2t::Point(0,0));
    xx.push_back(new p2t::Point(0,30));
    xx.push_back(new p2t::Point(30,30));
    //xx.push_back(new p2t::Point(20,30));
    polylines.push_back(xx);
    cout<<"Cheat"<<endl;
  }*/
  //--------------------------------------------------------------//
  //--------------------------------------------------------------//
  //--------------------------------------------------------------//
  cout<<"Polyline size po sorteni:"<<polylines.size()<<endl;
  vector<p2t::Point*> tmpPolygon = polylines.back();
  polylines.pop_back();
  pair<vector<p2t::Point*>,int >tmpPair = make_pair(tmpPolygon,-1);
  vector<pair<vector<p2t::Point*> , int> > tmpVecPair;
  tmpVecPair.push_back(tmpPair);
  polygonsWithHoles.push_back(tmpVecPair);//tmplist);
  int in;
  int out;
  int holeIn = -1; //-1 for polygon
  bool placeFound = false;
  int pos = 0;
  //
  for(;polylines.size()>0;)
  {
    tmpPolygon = polylines.back(); 
    polylines.pop_back();
    pos = 0;
    placeFound = false;
    holeIn = -1; 
    for(unsigned  int k = 0;k<polygonsWithHoles.size();k++)
    {
      in = 0;
      out = 0;
       
      //test 3 points and based on majority decide if this polygons is inside biggerPolygon
      for(int j = 0;j<3;j++)
      {
        vector<p2t::Point*>  biggerPolygon = polygonsWithHoles[k][pos].first;//.front();
        double a = tmpPolygon[j]->x; 
        double b = tmpPolygon[j]->y;
        if(vertexInPolygon((int)biggerPolygon.size(),biggerPolygon,a,b))
          in++;
        else
          out++;
      }
      if(in>out || placeFound)
      {
        if(in > out) 
        {
          if(holeIn!=-1)
            holeIn = -1;  //swaping between hole and not a hole
          else
            holeIn = pos;//pos; 
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
           tmpPair = make_pair(tmpPolygon,holeIn);
           polygonsWithHoles[k].push_back(tmpPair);
           break;
        }

      }
    }
    if(!placeFound)//out>in) // if we get here and point wasnt it, condition should be always true
    {                        // and we found a new polygon

        tmpPair = make_pair(tmpPolygon,-1);
        tmpVecPair.erase(tmpVecPair.begin(),tmpVecPair.end());
        tmpVecPair.push_back(tmpPair);
        polygonsWithHoles.push_back(tmpVecPair);//list2);
    }
  }
//musim projt kazdy polygonswithholes a asi testnout, ze tam nejsou "stejny" body
//polygonsWithHoles[0][0].first.erase(polygonsWithHoles[0][0].first.begin()+15,polygonsWithHoles[0][0].first.begin()+16);
/*cout<<"////////////////////"<<endl;
cout << "PolygonsWithholes"<<endl;
  for (int i = 0; i < polygonsWithHoles[0][0].first.size(); ++i)
  {
    //cout<<"translate(["<<polygonsWithHoles[0][0].first[i]->x<<", "<<polygonsWithHoles[0][0].first[i]->y<<",0]) sphere(r = 3);"<<endl;
    //polg.push_back(new p2t::Point(,));
    cout.precision(numeric_limits<double>::max_digits10);
    cout<<"polg.push_back(new p2t::Point("<<polygonsWithHoles[0][0].first[i]->x<<","<<polygonsWithHoles[0][0].first[i]->y<<"));"<<endl;
  }
  cout<<"Mame tu "<<polygonsWithHoles.size()<<" polygonu"<<endl;
  */
  /*
  polygonsWithHoles[0][0].first.resize(5);
  polygonsWithHoles[0][0].first[0]->x = 0; 
  polygonsWithHoles[0][0].first[0]->y = 0;
  polygonsWithHoles[0][0].first[1]->x = 5; 
  polygonsWithHoles[0][0].first[1]->y = 0;
  polygonsWithHoles[0][0].first[2]->x = 5; 
  polygonsWithHoles[0][0].first[2]->y = 5;
  polygonsWithHoles[0][0].first[3]->x = 5; 
  polygonsWithHoles[0][0].first[3]->y = 5+10e-10;
  polygonsWithHoles[0][0].first[4]->x = 0; 
  polygonsWithHoles[0][0].first[4]->y = 5;
*/
  }

void Mesh::poly2triTest()
{
  polygonsWithHoles.resize(1);
  polygonsWithHoles[0].resize(1);
  //polygonsWithHoles[0][0].first.resize(26);
  vector<p2t::Point*> polg;
  
  
  polg.push_back(new p2t::Point(0,00));
  polg.push_back(new p2t::Point(50,0));
  polg.push_back(new p2t::Point(60,0));
  polg.push_back(new p2t::Point(0,50));
  polg.push_back(new p2t::Point(50,50));
  polg.push_back(new p2t::Point(55,-50));
  polg.push_back(new p2t::Point(25,35));

  
  /*
  polg.push_back(new p2t::Point(109.271,3.61372e-08));
  polg.push_back(new p2t::Point(117.896,20.0236));
  polg.push_back(new p2t::Point(61.4449,101.42));
  polg.push_back(new p2t::Point(141.341,41.2166));
  polg.push_back(new p2t::Point(170.711,3.61372e-08));
  polg.push_back(new p2t::Point(190.61,61.1191));
  polg.push_back(new p2t::Point(178.605,87.4552));
  polg.push_back(new p2t::Point(169.867,141.42));
  polg.push_back(new p2t::Point(118.75,184.204));
  polg.push_back(new p2t::Point(82.2003,194.364));
  polg.push_back(new p2t::Point(-5.1589,207.326));
  polg.push_back(new p2t::Point(-68.8395,194.364));
  polg.push_back(new p2t::Point(-108.907,188.633));
  polg.push_back(new p2t::Point(-130.363,168.718));
  polg.push_back(new p2t::Point(-169.867,141.42));
  polg.push_back(new p2t::Point(-170.711,0));
  polg.push_back(new p2t::Point(-170.711,3.61372e-08));
  polg.push_back(new p2t::Point(-83.2822,92.4942));
  polg.push_back(new p2t::Point(-76.0951,94.5986));
  polg.push_back(new p2t::Point(36.1892,101.42));
  polg.push_back(new p2t::Point(-80.2957,89.1774));
  polg.push_back(new p2t::Point(-98.8727,47.4528));
  polg.push_back(new p2t::Point(-109.271,3.61372e-08));
  polg.push_back(new p2t::Point(-70.9228,3.61372e-08));
  polg.push_back(new p2t::Point(-3.25381e-08,3.61372e-08));
  polg.push_back(new p2t::Point(70.9228,3.61372e-08));
  */

  /*
  polg.push_back(new p2t::Point(-3.7082,-11.4127));
polg.push_back(new p2t::Point(-2.96459,-11.5708));
polg.push_back(new p2t::Point(-1.70277,-11.839));
polg.push_back(new p2t::Point(-1.25434,-11.9343));
polg.push_back(new p2t::Point(0.162123,-11.9343));
polg.push_back(new p2t::Point(1.25434,-11.9343));
polg.push_back(new p2t::Point(1.93465,-11.7897));
polg.push_back(new p2t::Point(3.19789,-11.5211));
polg.push_back(new p2t::Point(3.7082,-11.4127));
polg.push_back(new p2t::Point(4.54616,-11.0396));
polg.push_back(new p2t::Point(5.36919,-10.6732));
polg.push_back(new p2t::Point(6,-10.3923));
polg.push_back(new p2t::Point(7.08555,-9.6036));
polg.push_back(new p2t::Point(7.27636,-9.46497));
polg.push_back(new p2t::Point(7.81853,-9.07107));
polg.push_back(new p2t::Point(8.18999,-10.248));
polg.push_back(new p2t::Point(9.41904,-14.1421));
polg.push_back(new p2t::Point(11.4845,-13.469));
polg.push_back(new p2t::Point(11.9816,-13.3069));
polg.push_back(new p2t::Point(14.4535,-12.2027));
polg.push_back(new p2t::Point(15.9039,-11.5549));
polg.push_back(new p2t::Point(17.5742,-10.587));
polg.push_back(new p2t::Point(20.2701,-9.02482));
polg.push_back(new p2t::Point(20.9513,-8.52777));
polg.push_back(new p2t::Point(22.948,-7.07107));
polg.push_back(new p2t::Point(21.7855,-5.45827));
polg.push_back(new p2t::Point(21.081,-4.4809));
polg.push_back(new p2t::Point(19.1858,-2.36312));
polg.push_back(new p2t::Point(17.0711,-1.80688e-09));
polg.push_back(new p2t::Point(14.2811,-2.5246));
polg.push_back(new p2t::Point(13.8239,-2.93835));
polg.push_back(new p2t::Point(12.6019,-3.82989));
polg.push_back(new p2t::Point(11.8885,-4.35039));
polg.push_back(new p2t::Point(11.086,-4.93584));
polg.push_back(new p2t::Point(10.8527,-5.07107));
polg.push_back(new p2t::Point(10.9625,-4.88084));
polg.push_back(new p2t::Point(11.215,-4.10392));
polg.push_back(new p2t::Point(11.4207,-3.47089));
polg.push_back(new p2t::Point(11.7378,-2.49494));
polg.push_back(new p2t::Point(11.7811,-2.08265));
polg.push_back(new p2t::Point(12,-1.80688e-09));
polg.push_back(new p2t::Point(11.7378,-1.80688e-09));
polg.push_back(new p2t::Point(11.4755,-1.80688e-09));
polg.push_back(new p2t::Point(10.9625,-1.80688e-09));
polg.push_back(new p2t::Point(10.4495,-1.80688e-09));
polg.push_back(new p2t::Point(9.7082,-1.80688e-09));
polg.push_back(new p2t::Point(8.96686,-1.80688e-09));
polg.push_back(new p2t::Point(8.02957,-1.80688e-09));
polg.push_back(new p2t::Point(7.09228,-1.80688e-09));
polg.push_back(new p2t::Point(6,-1.80688e-09));
polg.push_back(new p2t::Point(4.90772,-1.80688e-09));
polg.push_back(new p2t::Point(3.7082,-1.80688e-09));
polg.push_back(new p2t::Point(2.50868,-1.80688e-09));
polg.push_back(new p2t::Point(1.25434,-1.80688e-09));
polg.push_back(new p2t::Point(-1.8991e-10,-1.80688e-09));
polg.push_back(new p2t::Point(-1.25434,-1.80688e-09));
polg.push_back(new p2t::Point(-2.50868,-1.80688e-09));
polg.push_back(new p2t::Point(-3.7082,-1.80688e-09));
polg.push_back(new p2t::Point(-4.90772,-1.80688e-09));
polg.push_back(new p2t::Point(-6,-1.80688e-09));
polg.push_back(new p2t::Point(-7.09228,-1.80688e-09));
polg.push_back(new p2t::Point(-8.02957,-1.80688e-09));
polg.push_back(new p2t::Point(-8.96686,-1.80688e-09));
polg.push_back(new p2t::Point(-9.7082,-1.80688e-09));
polg.push_back(new p2t::Point(-10.4495,-1.80688e-09));
polg.push_back(new p2t::Point(-10.9625,-1.80688e-09));
polg.push_back(new p2t::Point(-11.4755,-1.80688e-09));
polg.push_back(new p2t::Point(-11.7378,-1.80688e-09));
polg.push_back(new p2t::Point(-12,-1.80688e-09));
polg.push_back(new p2t::Point(-11.7974,-1.92804));
polg.push_back(new p2t::Point(-11.7378,-2.49494));
polg.push_back(new p2t::Point(-11.5485,-3.07756));
polg.push_back(new p2t::Point(-11.3106,-3.80964));
polg.push_back(new p2t::Point(-10.9625,-4.88084));
polg.push_back(new p2t::Point(-10.8527,-5.07107));
polg.push_back(new p2t::Point(-11.086,-4.93584));
polg.push_back(new p2t::Point(-12.213,-4.11361));
polg.push_back(new p2t::Point(-13.0753,-3.48446));
polg.push_back(new p2t::Point(-13.8239,-2.93835));
polg.push_back(new p2t::Point(-14.4596,-2.36312));
polg.push_back(new p2t::Point(-17.0711,-1.80688e-09));
polg.push_back(new p2t::Point(-21.081,-4.4809));
polg.push_back(new p2t::Point(-22.948,-7.07107));
polg.push_back(new p2t::Point(-20.2701,-9.02482));
polg.push_back(new p2t::Point(-18.7502,-9.90552));
polg.push_back(new p2t::Point(-15.9039,-11.5549));
polg.push_back(new p2t::Point(-12.7994,-12.9417));
polg.push_back(new p2t::Point(-11.9816,-13.3069));
polg.push_back(new p2t::Point(-9.41904,-14.1421));
polg.push_back(new p2t::Point(-8.35203,-10.7614));
//polg.push_back(new p2t::Point(-7.818531,-9.07107));
polg.push_back(new p2t::Point(-7.81853,-9.07107));
polg.push_back(new p2t::Point(-6,-10.3923));
polg.push_back(new p2t::Point(-5.09592,-10.7948));
polg.push_back(new p2t::Point(-4.24826,-11.1722));
*/
/*

polg.push_back(new p2t::Point(-3.7081999778747559,11.412699699401855));
polg.push_back(new p2t::Point(-2.9330453872680664,11.577456474304199));
polg.push_back(new p2t::Point(-1.7268095016479492,11.833859443664551));
polg.push_back(new p2t::Point(-1.2543400526046753,11.934300422668457));
polg.push_back(new p2t::Point(0.10308375209569931,11.934300422668457));
polg.push_back(new p2t::Point(1.2543400526046753,11.934300422668457));
polg.push_back(new p2t::Point(1.970956563949585,11.781937599182129));
polg.push_back(new p2t::Point(3.1763381958007812,11.525713920593262));
polg.push_back(new p2t::Point(3.7081999778747559,11.412699699401855));
polg.push_back(new p2t::Point(4.5800042152404785,11.024552345275879));
polg.push_back(new p2t::Point(5.351219654083252,10.681178092956543));
polg.push_back(new p2t::Point(6,10.392299652099609));
polg.push_back(new p2t::Point(7.1132984161376953,9.5834465026855469));
polg.push_back(new p2t::Point(7.2627043724060059,9.4748973846435547));
polg.push_back(new p2t::Point(7.6924142837524414,9.1626958847045898));
polg.push_back(new p2t::Point(7.9661688804626465,10.079180717468262));
polg.push_back(new p2t::Point(9.2224225997924805,14.28494930267334));
polg.push_back(new p2t::Point(11.490257263183594,13.54175853729248));
polg.push_back(new p2t::Point(12.032922744750977,13.363916397094727));
polg.push_back(new p2t::Point(14.472713470458984,12.267500877380371));
polg.push_back(new p2t::Point(15.963020324707031,11.597794532775879));
polg.push_back(new p2t::Point(17.603458404541016,10.640851974487305));
polg.push_back(new p2t::Point(20.328836441040039,9.0509929656982422));
polg.push_back(new p2t::Point(20.98573112487793,8.5680360794067383));
polg.push_back(new p2t::Point(22.92474365234375,7.1424946784973145));
polg.push_back(new p2t::Point(21.736728668212891,5.4699654579162598));
polg.push_back(new p2t::Point(21.025825500488281,4.469172477722168));
polg.push_back(new p2t::Point(19.169164657592773,2.3710436820983887));
polg.push_back(new p2t::Point(17.071100234985352,1.8251312683403853e-09));
polg.push_back(new p2t::Point(14.298360824584961,2.5319082736968994));
polg.push_back(new p2t::Point(13.847722053527832,2.9434208869934082));
polg.push_back(new p2t::Point(12.646105766296387,3.8268623352050781));
polg.push_back(new p2t::Point(11.91185188293457,4.3666973114013672));
polg.push_back(new p2t::Point(11.118185043334961,4.9501547813415527));
polg.push_back(new p2t::Point(10.823102951049805,5.1222929954528809));
polg.push_back(new p2t::Point(10.962499618530273,4.8808398246765137));
polg.push_back(new p2t::Point(11.212868690490723,4.1104645729064941));
polg.push_back(new p2t::Point(11.425088882446289,3.4573647975921631));
polg.push_back(new p2t::Point(11.737799644470215,2.4949400424957275));
polg.push_back(new p2t::Point(11.780765533447266,2.086097240447998));
polg.push_back(new p2t::Point(12,1.8251312683403853e-09));
polg.push_back(new p2t::Point(11.737799644470215,1.8251311573180828e-09));
polg.push_back(new p2t::Point(11.475545883178711,1.8251313793626878e-09));
polg.push_back(new p2t::Point(10.962499618530273,1.8251314903849902e-09));
polg.push_back(new p2t::Point(10.449520111083984,1.8251315793626878e-09)); // upravno 137 na 157
polg.push_back(new p2t::Point(9.7082004547119141,1.8251315793626879e-09));
//polg.push_back(new p2t::Point(8.9668588638305664,1.8251316573180828e-09));
/*polg.push_back(new p2t::Point(8.0295696258544922,1.8251309352734779e-09));
polg.push_back(new p2t::Point(7.0922760963439941,1.8251310462957804e-09));
polg.push_back(new p2t::Point(6,1.8251310462957804e-09));
polg.push_back(new p2t::Point(4.9077243804931641,1.8251314903849902e-09));
polg.push_back(new p2t::Point(3.7081999778747559,1.8251312683403853e-09));
polg.push_back(new p2t::Point(2.5086812973022461,1.8251318234518976e-09));
polg.push_back(new p2t::Point(1.2543400526046753,1.8251311573180828e-09));
polg.push_back(new p2t::Point(-1.9182831711983539e-10,1.8251324895857124e-09));
polg.push_back(new p2t::Point(-1.2543400526046753,1.8251311573180828e-09));
 polg.push_back(new p2t::Point(-2.508681058883667,1.8251324895857124e-09));
polg.push_back(new p2t::Point(-3.7081999778747559,1.8251312683403853e-09));
polg.push_back(new p2t::Point(-4.9077243804931641,1.8251314903849902e-09));
polg.push_back(new p2t::Point(-6,1.8251310462957804e-09));
polg.push_back(new p2t::Point(-7.0922760963439941,1.8251306022065705e-09));
polg.push_back(new p2t::Point(-8.0295696258544922,1.8251309352734779e-09));
polg.push_back(new p2t::Point(-8.9668588638305664,1.8251308242511755e-09));
polg.push_back(new p2t::Point(-9.7082004547119141,1.8251313793626878e-09));
polg.push_back(new p2t::Point(-10.449520111083984,1.8251318234518976e-09));
polg.push_back(new p2t::Point(-10.962499618530273,1.8251314903849902e-09));
polg.push_back(new p2t::Point(-11.475545883178711,1.8251313793626878e-09));
polg.push_back(new p2t::Point(-11.737799644470215,1.8251311573180828e-09));
polg.push_back(new p2t::Point(-12,1.8251311573180828e-09));
polg.push_back(new p2t::Point(-11.796916007995605,1.9324287176132202));
polg.push_back(new p2t::Point(-11.737799644470215,2.4949400424957275));
polg.push_back(new p2t::Point(-11.550782203674316,3.0703856945037842));
polg.push_back(new p2t::Point(-11.307831764221191,3.8180198669433594));
polg.push_back(new p2t::Point(-10.962499618530273,4.8808398246765137));
polg.push_back(new p2t::Point(-10.823099136352539,5.1222929954528809));
polg.push_back(new p2t::Point(-11.118185043334961,4.9501547813415527));
polg.push_back(new p2t::Point(-12.233139991760254,4.1304402351379395));
polg.push_back(new p2t::Point(-13.110909461975098,3.4850966930389404));
polg.push_back(new p2t::Point(-13.847722053527832,2.9434208869934082));
polg.push_back(new p2t::Point(-14.474519729614258,2.3710453510284424));
polg.push_back(new p2t::Point(-17.071100234985352,1.8251311573180828e-09));
polg.push_back(new p2t::Point(-21.025825500488281,4.469172477722168));
polg.push_back(new p2t::Point(-22.92474365234375,7.1424946784973145));
polg.push_back(new p2t::Point(-20.328836441040039,9.0509929656982422));
polg.push_back(new p2t::Point(-18.869289398193359,9.9024181365966797));
polg.push_back(new p2t::Point(-15.963020324707031,11.597794532775879));
polg.push_back(new p2t::Point(-12.921987533569336,12.964390754699707));
polg.push_back(new p2t::Point(-12.032922744750977,13.363916397094727));
polg.push_back(new p2t::Point(-9.2224225997924805,14.28494930267334));
polg.push_back(new p2t::Point(-8.1280031204223633,10.620997428894043));
polg.push_back(new p2t::Point(-7.6924114227294922,9.1626958847045898));
polg.push_back(new p2t::Point(-7.6924138069152832,9.1626958847045898));
polg.push_back(new p2t::Point(-6,10.392299652099609));
polg.push_back(new p2t::Point(-5.0698676109313965,10.806416511535645));
polg.push_back(new p2t::Point(-4.2703475952148438,11.162395477294922));
*/
/*
  polg.push_back(new p2t::Point(0.319643, 0.839787));
polg.push_back(new p2t::Point(0.319648, 0.8367));
polg.push_back(new p2t::Point(0.319437, 0.832799));
polg.push_back(new p2t::Point(0.319154, 0.830715));
polg.push_back(new p2t::Point(0.319073, 0.830124));
polg.push_back(new p2t::Point(0.318665, 0.827986));
polg.push_back(new p2t::Point(0.318566, 0.827472));
polg.push_back(new p2t::Point(0.318105, 0.825581));
polg.push_back(new p2t::Point(0.318011, 0.825198));
polg.push_back(new p2t::Point(0.317445, 0.823274));
polg.push_back(new p2t::Point(0.31735, 0.822953));
polg.push_back(new p2t::Point(0.316675, 0.821));
polg.push_back(new p2t::Point(0.316585, 0.82074));
polg.push_back(new p2t::Point(0.315796, 0.818764));
polg.push_back(new p2t::Point(0.315717, 0.818566));
polg.push_back(new p2t::Point(0.314802, 0.816555));
polg.push_back(new p2t::Point(0.31474, 0.816419));
polg.push_back(new p2t::Point(0.313702, 0.814395));
polg.push_back(new p2t::Point(0.313663, 0.814319));
polg.push_back(new p2t::Point(0.312497, 0.812289));
polg.push_back(new p2t::Point(0.312488, 0.812273));
polg.push_back(new p2t::Point(0.312145, 0.811736));
polg.push_back(new p2t::Point(0.311218, 0.810285));
polg.push_back(new p2t::Point(0.30989, 0.808407));
polg.push_back(new p2t::Point(0.309822, 0.808316));
polg.push_back(new p2t::Point(0.308408, 0.806504));
polg.push_back(new p2t::Point(0.308284, 0.806366));
polg.push_back(new p2t::Point(0.305254, 0.803004));
polg.push_back(new p2t::Point(0.30372, 0.801525));
polg.push_back(new p2t::Point(0.303446, 0.801273));
polg.push_back(new p2t::Point(0.301848, 0.79987));
polg.push_back(new p2t::Point(0.300233, 0.798578));
polg.push_back(new p2t::Point(0.299835, 0.798275));
polg.push_back(new p2t::Point(0.298168, 0.797064));
polg.push_back(new p2t::Point(0.296471, 0.795948));
polg.push_back(new p2t::Point(0.295931, 0.795611));
polg.push_back(new p2t::Point(0.294192, 0.794581));
polg.push_back(new p2t::Point(0.292453, 0.793659));
polg.push_back(new p2t::Point(0.291771, 0.793317));
polg.push_back(new p2t::Point(0.290001, 0.792482));
polg.push_back(new p2t::Point(0.288072, 0.791685));
polg.push_back(new p2t::Point(0.287158, 0.791335));
polg.push_back(new p2t::Point(0.285204, 0.790643));
polg.push_back(new p2t::Point(0.28327, 0.790068));
polg.push_back(new p2t::Point(0.28219, 0.789777));
polg.push_back(new p2t::Point(0.280244, 0.789306));
polg.push_back(new p2t::Point(0.278268, 0.788936));
polg.push_back(new p2t::Point(0.276991, 0.788731));
polg.push_back(new p2t::Point(0.275016, 0.788468));
polg.push_back(new p2t::Point(0.273624, 0.788394));
polg.push_back(new p2t::Point(0.269728, 0.788188));
polg.push_back(new p2t::Point(0.268183, 0.788271));
polg.push_back(new p2t::Point(0.264339, 0.788479));
polg.push_back(new p2t::Point(0.262491, 0.78873));
polg.push_back(new p2t::Point(0.260829, 0.789002));
polg.push_back(new p2t::Point(0.259013, 0.789349));
polg.push_back(new p2t::Point(0.25748, 0.789724));
polg.push_back(new p2t::Point(0.255996, 0.790123));
polg.push_back(new p2t::Point(0.254493, 0.790565));
polg.push_back(new p2t::Point(0.253031, 0.791071));
polg.push_back(new p2t::Point(0.251534, 0.791629));
polg.push_back(new p2t::Point(0.250106, 0.792199));
polg.push_back(new p2t::Point(0.248711, 0.792833));
polg.push_back(new p2t::Point(0.247217, 0.793556));
polg.push_back(new p2t::Point(0.24586, 0.794253));
polg.push_back(new p2t::Point(0.244547, 0.795006));
polg.push_back(new p2t::Point(0.243096, 0.795886));
polg.push_back(new p2t::Point(0.241825, 0.796698));
polg.push_back(new p2t::Point(0.240603, 0.797562));
polg.push_back(new p2t::Point(0.239222, 0.798589));
polg.push_back(new p2t::Point(0.238044, 0.799508));
polg.push_back(new p2t::Point(0.236757, 0.800668));
polg.push_back(new p2t::Point(0.234544, 0.802662));
polg.push_back(new p2t::Point(0.23353, 0.803714));
polg.push_back(new p2t::Point(0.232375, 0.804968));
polg.push_back(new p2t::Point(0.23141, 0.806067));
polg.push_back(new p2t::Point(0.230494, 0.807213));
polg.push_back(new p2t::Point(0.229468, 0.808559));
polg.push_back(new p2t::Point(0.228604, 0.809748));
polg.push_back(new p2t::Point(0.227782, 0.810998));
polg.push_back(new p2t::Point(0.226886, 0.812433));
polg.push_back(new p2t::Point(0.226121, 0.813724));
polg.push_back(new p2t::Point(0.225412, 0.81506));
polg.push_back(new p2t::Point(0.224671, 0.81654));
polg.push_back(new p2t::Point(0.224022, 0.817915));
polg.push_back(new p2t::Point(0.223379, 0.819471));
polg.push_back(new p2t::Point(0.222749, 0.821114));
polg.push_back(new p2t::Point(0.222183, 0.822712));
polg.push_back(new p2t::Point(0.221694, 0.82436));
polg.push_back(new p2t::Point(0.221254, 0.825989));
polg.push_back(new p2t::Point(0.220846, 0.827672));
polg.push_back(new p2t::Point(0.220512, 0.829457));
polg.push_back(new p2t::Point(0.22025, 0.831085));
polg.push_back(new p2t::Point(0.220008, 0.8329));
polg.push_back(new p2t::Point(0.219926, 0.834463));
polg.push_back(new p2t::Point(0.219728, 0.838188));
polg.push_back(new p2t::Point(0.219808, 0.839664));
polg.push_back(new p2t::Point(0.220019, 0.843577));
polg.push_back(new p2t::Point(0.220293, 0.845592));
polg.push_back(new p2t::Point(0.220504, 0.846889));
polg.push_back(new p2t::Point(0.22089, 0.848903));
polg.push_back(new p2t::Point(0.221322, 0.850676));
polg.push_back(new p2t::Point(0.221587, 0.851661));
polg.push_back(new p2t::Point(0.222106, 0.853423));
polg.push_back(new p2t::Point(0.222721, 0.855203));
polg.push_back(new p2t::Point(0.223036, 0.85605));
polg.push_back(new p2t::Point(0.223739, 0.857809));
polg.push_back(new p2t::Point(0.224551, 0.859594));
polg.push_back(new p2t::Point(0.224893, 0.860302));
polg.push_back(new p2t::Point(0.225793, 0.862056));
polg.push_back(new p2t::Point(0.226802, 0.863813));
polg.push_back(new p2t::Point(0.227143, 0.864376));
polg.push_back(new p2t::Point(0.228238, 0.86609));
polg.push_back(new p2t::Point(0.229444, 0.867796));
polg.push_back(new p2t::Point(0.229758, 0.868219));
polg.push_back(new p2t::Point(0.231048, 0.869872));
polg.push_back(new p2t::Point(0.231348, 0.870204));
polg.push_back(new p2t::Point(0.234202, 0.873372));
polg.push_back(new p2t::Point(0.235775, 0.874888));
polg.push_back(new p2t::Point(0.235964, 0.875062));
polg.push_back(new p2t::Point(0.237608, 0.876506));
polg.push_back(new p2t::Point(0.239365, 0.877911));
polg.push_back(new p2t::Point(0.239467, 0.877989));
polg.push_back(new p2t::Point(0.241289, 0.879312));
polg.push_back(new p2t::Point(0.243247, 0.8806));
polg.push_back(new p2t::Point(0.243248, 0.8806));
polg.push_back(new p2t::Point(0.243269, 0.880613));
polg.push_back(new p2t::Point(0.245211, 0.881763));
polg.push_back(new p2t::Point(0.245264, 0.881794));
polg.push_back(new p2t::Point(0.247226, 0.882835));
polg.push_back(new p2t::Point(0.247335, 0.882893));
polg.push_back(new p2t::Point(0.249289, 0.883815));
polg.push_back(new p2t::Point(0.249455, 0.883894));
polg.push_back(new p2t::Point(0.251578, 0.884771));
polg.push_back(new p2t::Point(0.25183, 0.884875));
polg.push_back(new p2t::Point(0.253929, 0.885618));
polg.push_back(new p2t::Point(0.254252, 0.885732));
polg.push_back(new p2t::Point(0.256321, 0.886348));
polg.push_back(new p2t::Point(0.256715, 0.886465));
polg.push_back(new p2t::Point(0.258749, 0.886957));
polg.push_back(new p2t::Point(0.259212, 0.887069));
polg.push_back(new p2t::Point(0.261265, 0.887455));
polg.push_back(new p2t::Point(0.261815, 0.887558));
polg.push_back(new p2t::Point(0.263821, 0.887825));
polg.push_back(new p2t::Point(0.264441, 0.887907));
polg.push_back(new p2t::Point(0.268267, 0.88811));
polg.push_back(new p2t::Point(0.26973, 0.88819));
polg.push_back(new p2t::Point(0.271359, 0.8881));
polg.push_back(new p2t::Point(0.275117, 0.887896));
polg.push_back(new p2t::Point(0.276867, 0.887659));
polg.push_back(new p2t::Point(0.277792, 0.887533));
polg.push_back(new p2t::Point(0.279478, 0.887211));
polg.push_back(new p2t::Point(0.280444, 0.887026));
polg.push_back(new p2t::Point(0.281854, 0.886682));
polg.push_back(new p2t::Point(0.282718, 0.886471));
polg.push_back(new p2t::Point(0.28408, 0.88607));
polg.push_back(new p2t::Point(0.284963, 0.88581));
polg.push_back(new p2t::Point(0.286277, 0.885356));
polg.push_back(new p2t::Point(0.287176, 0.885045));
polg.push_back(new p2t::Point(0.288442, 0.884539));
polg.push_back(new p2t::Point(0.28935, 0.884177));
polg.push_back(new p2t::Point(0.290578, 0.883618));
polg.push_back(new p2t::Point(0.291497, 0.8832));
polg.push_back(new p2t::Point(0.292679, 0.882594));
polg.push_back(new p2t::Point(0.293596, 0.882123));
polg.push_back(new p2t::Point(0.294732, 0.881471));
polg.push_back(new p2t::Point(0.295642, 0.880948));
polg.push_back(new p2t::Point(0.296733, 0.880252));
polg.push_back(new p2t::Point(0.297631, 0.879678));
polg.push_back(new p2t::Point(0.298675, 0.878939));
polg.push_back(new p2t::Point(0.299555, 0.878317));
polg.push_back(new p2t::Point(0.300555, 0.877536));
polg.push_back(new p2t::Point(0.301412, 0.876868));
polg.push_back(new p2t::Point(0.303287, 0.875178));
polg.push_back(new p2t::Point(0.304912, 0.873714));
polg.push_back(new p2t::Point(0.305771, 0.872823));
polg.push_back(new p2t::Point(0.306519, 0.872047));
polg.push_back(new p2t::Point(0.307337, 0.871115));
polg.push_back(new p2t::Point(0.308046, 0.870308));
polg.push_back(new p2t::Point(0.308824, 0.869335));
polg.push_back(new p2t::Point(0.309492, 0.8685));
polg.push_back(new p2t::Point(0.310229, 0.867485));
polg.push_back(new p2t::Point(0.310852, 0.866627));
polg.push_back(new p2t::Point(0.311557, 0.865555));
polg.push_back(new p2t::Point(0.31214, 0.864668));
polg.push_back(new p2t::Point(0.312801, 0.863552));
polg.push_back(new p2t::Point(0.313335, 0.862651));
polg.push_back(new p2t::Point(0.313952, 0.861489));
polg.push_back(new p2t::Point(0.314433, 0.86058));
polg.push_back(new p2t::Point(0.315004, 0.859371));
polg.push_back(new p2t::Point(0.315434, 0.85846));
polg.push_back(new p2t::Point(0.316004, 0.85708));
polg.push_back(new p2t::Point(0.316415, 0.856086));
polg.push_back(new p2t::Point(0.316925, 0.854646));
polg.push_back(new p2t::Point(0.317273, 0.853664));
polg.push_back(new p2t::Point(0.317718, 0.852164));
polg.push_back(new p2t::Point(0.318005, 0.851201));
polg.push_back(new p2t::Point(0.318383, 0.849641));
polg.push_back(new p2t::Point(0.31861, 0.848704));
polg.push_back(new p2t::Point(0.318923, 0.847032));
polg.push_back(new p2t::Point(0.319098, 0.846101));
polg.push_back(new p2t::Point(0.319329, 0.844364));
polg.push_back(new p2t::Point(0.319448, 0.843475));
*/
p2t::CDT* tmp = new  p2t::CDT(polg);
tmp->Triangulate();
vector<p2t::Triangle*> triangles = tmp->GetTriangles();
createFacet(triangles); 
  
}

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
        if(intersect(p[m], p[m+1], p[j], p[j+1] ))
        {

          cerr<<"CHYBA! POLYGON NENI JEDNODUCHY"<<endl;
          cerr<<"m= "<<m<<" j= "<<j<<endl;
          cerr.precision(numeric_limits<double>::max_digits10);
          cerr<<p[m]->x<<" "<<p[m]->y<<", "<<p[m+1]->x<<" "<<p[m+1]->y<<" vs "<<p[j]->x<<" "<<p[j]->y<<" , "<<p[j+1]->x<<" "<<p[j+1]->y<<endl;
          p.erase(p.begin()+j,p.begin()+j+1);
          pointRemoved = true;
          continue;
          //return;
        }
      }
    }
  
  
    /*in p is a polyline (points of polygon), we have to test last edge (from last to first point) as well*/
    /*TOHLE BY NEMELO NIKDY NASTAT, KDYZ JDE JEN O CHYBU V NEPRESNOSTI POCITANI REZU
    for (int k = 1; k < p.size()-2; k+=1)
    {
      if(intersect(p[k], p[k+1], p[0], p[ p.size()-1] ))
      {
        cerr<<"CHYBA! POLYGON NENI JEDNODUCHY2"<<endl;
        p.erase(p.begin()+j,p.begin()+j+1);
        pointRemoved = true;
        continue;
      }
    }
   */
  }while(pointRemoved == true);
}

/*
s novou funkci na removnuti spatnych pointu
bych ji mohl volat pouze v pripade, ze poly2tri spadne - nebo vzdy? kdyz je vse ok, nic by se nemelo vymaza
musim projit polygony.. a teoreticky mu predhodit kazdy polygon co mam.. i diry, a testnout ze nejsou non-simple
a kdyz jsou tak ty body vyhodit
a pak si pamatovat ze jsem neco vyhodil ( a jaka byla tolerance)
a zavolat spojovani pointu

*/
/*
*Tests every polygon if its non-simple, if not, some points are removed in order to make it simple polygon.
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

//TODO SMAZAT?, nefunguje uplne... nevymaze trojuhelniky navic, co vyuzivaj existujici body.. (bok u vlka)
void Mesh::checkPoly2triResult( vector<p2t::Triangle*>& triangles )
{
  /*prohledej vsechno co vytvoril
  porovnej to s kazdym bodem co mame
  uloz si vsechny novy body
  smaz vsechny trojuhelniky obsahujici tenhle novy bod*/
  set<p2t::Point,setVertComp> originalBorderPoints;
  //set<p2t::Point,setVertComp> extraPoints;
  // first we need to get set of all points, we have too keep in mind that coordinates are mixed
  for (std::set<stl_vertex>::iterator i = border2.begin(); i != border2.end(); ++i) 
  {
    /*
    kdyz je removle x, je y=y a v x je z... 
    kdzy je removedle y, x=x a v y je z
    kdyz je removedle z x=x a y=y
    */
    double x = (*i).x;
    double y = (*i).y;
    if(removedAxis == 'x')
      x = (*i).z;
    if(removedAxis == 'y')
      y=(*i).z;

    p2t::Point p = Point(x,y);
    originalBorderPoints.insert(p);
    //cout<<p.x<<" "<<p.y<<endl;
  }
  //cout<<"////////"<<endl;
  for (int i = 0; i < triangles.size(); ++i)
  {
    for(int j = 0; j<3; j++)
    {
      p2t::Point* tmp = triangles[i]->GetPoint(j);
      p2t::Point p = Point(tmp->x,tmp->y);//Point(triangles[i]->GetPoint(j) , triangles[i]->GetPoint(j)->y);
      cout << tmp->x<<" "<<tmp->y<<endl;
      //cout<<"OBPC: "<<originalBorderPoints.count(p)<<endl;
      if( originalBorderPoints.count(p) == 0)//(*p2t::Point)(triangles[i]->GetPoint(j))) == 0) // test if 2D point is in the original border created by cut
      {
        //its not in it
        //extraPoints.insert(p);
        //we found point which didnt exist, so we will erase this triangle
        //cout<<"MAZUUU"<<endl;
        triangles.erase(triangles.begin()+i,triangles.begin()+i+1);
        i--;
        break;
      }
    }
  }
}
/*
void Mesh::checkPoly2triResult()
{
    cout<<"Check zac"<<endl;
  vector<vector<p2t::Point*>> polygonsForTest;
  for (int i = 0; i < polygonsWithHoles.size(); ++i) 
  {
    for (int j = 0; j < polygonsWithHoles[i].size(); ++j)
    {
      if(polygonsWithHoles[i][j].second == -1)
      {
        vector<p2t::Point*> tmp=polygonsWithHoles[i][j].first;
        for (int k = 0; k < polygonsWithHoles[i][j].first.size(); ++k)
        {
          tmp[k]=new p2t::Point( tmp[k]->x, tmp[k]->y); // have to alocate those points again, because poly2tri will delete them via its destructor
        }
        polygonsForTest.push_back(tmp);
      }
    }
  }

  for (unsigned int i = 0; i < polygonsForTest.size(); ++i)
  {
    //isPolygonSimple(polygonsForTest[i]);
    p2t::CDT* tmp = new  p2t::CDT(polygonsForTest[i]);
    tmp->Triangulate();
    vector<p2t::Triangle*> triangles = tmp->GetTriangles();
    for(int k=0; k<3; k++)
    {
      if(checkForNewPoints(triangles,polygonsForTest[i]))
      {
        /////////////

        ////////////////
        //TODO UPravit, abych odstranoval body co najdu v ispolygonsimple
        fixNonsimplePolygon(polygonsForTest[i]);
        //TODO zapamatovat si toleranci v ktere se ubral bod a pridat funkci co zavola admesh opravu
        cerr<<"Mesh is non-manifold or plane cut made a non-simple polygon due to floating point calculations, trying to fix the problem, but result might be invalid or none."<<endl;
      }
      else break;
    }
    //else
      //cout<<"Poly2tri test completed without a problem, cut made a simple polygon."<<endl;
    delete tmp; // poly2tri destructor will clear points stored in polygonsForTest   
  }
  cout<<"check end"<<endl;
}
*/

void Mesh::triangulateCut(int topOrBot)
{
  map<int, p2t::CDT*> polygons;//vector<p2t::CDT>> polygons;
  repairIfNonsimplePolygon();
  //checkPoly2triResult();
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
    //all polygon processed and ready to triangulation
    cout<<"Polygons.size= "<<polygons.size()<<endl;
    if(polygons.size()>0)
      for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
      {
          cout<<"triangulujem"<<endl;
          k->second->Triangulate();
          cout<<"done1"<<endl;
          vector<p2t::Triangle*> triangles = k->second->GetTriangles();
          cout<<"done2"<<endl;
          checkPoly2triResult(triangles);
          createFacet(triangles, topOrBot);   
      }
    for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
    {
      delete (*k).second;
    }
    polygons.erase(polygons.begin(),polygons.end());   
  }
  polygonsWithHoles.clear(); 
}
bool Mesh::checkForNewPoints(vector<p2t::Triangle*> &triangles, vector<p2t::Point*>& npolygon)
{
  set<p2t::Point*,setVertComp> vertices;
  set<p2t::Point*,setVertComp> vertices2;
  int size = 0;
  int size2 = 0;
  for (int i = 0; i < npolygon.size(); ++i)
  {
    vertices.insert(npolygon[i]);
  }
  size = vertices.size();
  //vertices.clear();
  for (int i = 0; i < triangles.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      vertices2.insert(triangles[i]->GetPoint(j));
    }
  }
  size2 = vertices2.size();
  if(size!=size2)//if(size!=vertices.size())
    return true; // new points found
  return false;

}
void Mesh::fixNonsimplePolygon(vector<p2t::Point*>& npolygon)
{
  double minDif = numeric_limits<double>::max(); 
  int index = -1;
  for (int i = 1; i < npolygon.size(); ++i)
  {
    double dif1 = ABS((npolygon[i-1]->x - npolygon[i]->x));
    double dif2 = ABS((npolygon[i-1]->y - npolygon[i]->y));
    if(minDif > dif1 && minDif > dif2)
    {
      minDif = dif1 > dif2? dif1:dif2;
      index = i;
    }
  }
  
  if(minDif<0.01)
  {
    cout << "Mazu bod: x:" << npolygon[index]->x << " y: "<<npolygon[index]->y<<endl;
    cout << "Pred bod: x:" << npolygon[index-1]->x << " y: "<<npolygon[index-1]->y<<endl;
    npolygon.erase(npolygon.begin()+index-1,npolygon.begin()+index);
  }
}
void Mesh::createFacet(vector<p2t::Triangle*> &triangles, int side)
{
  // for each triangle, create facet
  for (vector<p2t::Triangle*>::iterator i = triangles.begin(); i != triangles.end(); i++) 
  {
    stl_vertex vertex;
    stl_facet facet;
    for (size_t j = 0; j < 3; j++) 
    {
      p2t::Point* p = (*i)->GetPoint(j);
      setMissingCoordinate(p,vertex);
      facet.vertex[j] = vertex;
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

void Mesh::pushToPolylinesFront(vector<p2t::Point*> &vec,stl_vertex vert)
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
bool Mesh::createBorderPolylines(bool processOnFac)
{
  if(processOnFac == true) 
    processOnFacets();
  if(border.size() == 0) 
  {
    if(!(bot_facets.size() != 0 && top_facets.size() != 0))//(bot_facets.size()!=0 || top_facets.size()!=0 ))
    {
      cerr<<"Nothing to cut"<<endl;
      return false;
    }
    else
    {
      this->save();
      return false;
    }
  }

  {
    //TODO rozmyslet, jestli to po zmene zachazeni s hranama NA potrebuju
  checkDuplicity(border);
  

  cout<<"Jsem v createborderpolylines"<<endl;
  cout<<"V border je: "<<endl;
  for (int i = 0; i < border.size(); ++i)
  {
    cout<<border[i].x << " "<<border[i].y<<" "<<border[i].z<<endl;
  }
  numOfPolylines = 0;
  stl_vertex cont = border.back();
  border.pop_back();
  stl_vertex end = border.back();
  border.pop_back();
  polylines.resize(20);
  eps = 1e-24;
  minEps = numeric_limits<double>::max(); 
  bool found = true;
  stl_vertex tmp1,tmp2;
  pushToPolylines(polylines[numOfPolylines],end);
  pushToPolylines(polylines[numOfPolylines],cont);

  while(border.size()!=0)
  { 
    if(!found) 
    {
      eps = minEps*1.005;
      /*
      if( eps <=0.5 && vertexEqual(cont,end))//(cont == end) )//&& border.size()>1)
      {
        if(!polylines[numOfPolylines].empty()) polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
        
        if(border.size()>0)
        {
          numOfPolylines++;
          if(numOfPolylines+5 > polylines.size()) 
            polylines.resize(numOfPolylines*10);
          end = border.back();
          border.pop_back();
          cont = border.back();
          border.pop_back();
          pushToPolylines(polylines[numOfPolylines],end);
          pushToPolylines(polylines[numOfPolylines],cont);
        }
        //nevim jestli jsou ty posledni 4 radky spravne.. byl tu break a nic, ale neungovalo to na duplu 0.1 0.1 1 0.9
        eps = 1e-24;
        minEps = numeric_limits<double>::max();
        //cout<<"BREAK"<<endl;
        found = true;
        continue;//break;
      }*/
    }
    //todo tohle by melo byt nad tim nahore
    if(eps > 0.5) // ignoring edges, if it gets to this, there was probably a problem with mesh
    {
      eps = 1e-24;
      minEps = numeric_limits<double>::max();
      //
      //       RETURN FALSE Could prevent poly2tri sefault 
      // 
      std::cerr<<"Unable to find a connected edge, mesh might be invalid"<<endl;
      //
      cont = border.back();
      border.pop_back();
      end = border.back();
      border.pop_back();
      if(!polylines[numOfPolylines].empty())
        numOfPolylines++;
      pushToPolylines(polylines[numOfPolylines],end);
      pushToPolylines(polylines[numOfPolylines],cont);
       // if we didnt found it with 0.5 tolerance.. lets just skip it
    }
    for (int i = border.size()-1;i >= 0; i-=2)
    {
      found = false;
      tmp1 = border[i]; 
      tmp2 = border[i-1];

      if(vertexEqual(tmp1,cont) || vertexEqual(tmp2,cont) || vertexEqual(tmp1,end) || vertexEqual(tmp2,end))// tmp1 == cont||tmp2 == cont) // ve found next vertex in polyline
      {
        cout<<"pokracovani nalezeno"<<endl;
        found = true;
          //eps = 1e-24;
        if(vertexEqual(tmp1,cont)) //tmp1 == cont) 
        {            
          pushToPolylines(polylines[numOfPolylines],tmp2);
          cont = tmp2;
          border.erase(border.begin()+(i-1),border.begin()+i+1); // erase doesnt delete the last element IT DOES with set
        }  
        else if(vertexEqual(tmp2,cont))
        {
          pushToPolylines(polylines[numOfPolylines],tmp1);
          cont = tmp1;
          border.erase(border.begin()+(i-1),border.begin()+i+1); // again, different erase with set
        }
        else if(vertexEqual(tmp1,end))
        {
          pushToPolylinesFront(polylines[numOfPolylines],tmp2);
          end=tmp2;
          cout<<border.size()<<" size,i: "<<i<<endl;
          border.erase(border.begin()+(i-1),border.begin()+i+1); 
        }
        else if(vertexEqual(tmp2,end))
        {
          pushToPolylinesFront(polylines[numOfPolylines],tmp1);
          end = tmp1;
          cout<<border.size()<<" size,i: "<<i<<endl;
          border.erase(border.begin()+(i-1),border.begin()+i+1); 
        }
        if(vertexEqual(cont,end)) //(cont == end) )//&& border.size()>1)
        {
          //cout << "END NALEZEN"<<endl;
          if(!polylines[numOfPolylines].empty()) polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
          
          if(border.size()>0)
          {
            numOfPolylines++;
            if(numOfPolylines+5 > polylines.size())
              polylines.resize(numOfPolylines*10);
            end = border.back();
            border.pop_back();
            cont = border.back();
            border.pop_back();
            pushToPolylines(polylines[numOfPolylines],end);
            pushToPolylines(polylines[numOfPolylines],cont);
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
  }
  polylines.resize(numOfPolylines+1);

  cout<<"////////////////////////////"<<endl;
  cout<<"Vypis polylines: "<<"numOfPolylines je "<< numOfPolylines<<endl;

  for (int k = 0; k < polylines.size(); ++k)
  { 
    cout<<"Polyline cislo: "<<k<<" size je "<<polylines[k].size()<<endl;
    for (int i = 0; i < polylines[k].size(); i++)
    {
      cout<<"translate(["<<polylines[k][i]->x<<", "<<polylines[k][i]->y<<",0]) sphere(r = 0.3);"<<endl;//" //"<<polylines[k][i+1]->x<<" "<<polylines[k][i+1]->y<<endl;
    }
  }
  sortPolylines();
  if (polylines.size()==0)
    return false;
  return true;
}

/*
bool Mesh::createBorderPolylines(bool processOnFac)
{
  numOfPolylines = 0;
  if(processOnFac == true) 
    processOnFacets();
  if(border.size() == 0) 
  {
    if(!(bot_facets.size() != 0 && top_facets.size() != 0))//(bot_facets.size()!=0 || top_facets.size()!=0 ))
    {
      cerr<<"Nothing to cut"<<endl;
      return false;
    }
    else
    {
      this->save();
      return false;
    }
  }

  {
    //TODO rozmyslet, jestli to po zmene zachazeni s hranama NA potrebuju
  checkDuplicity(border);
  

  cout<<"Jsem v createborderpolylines"<<endl;
  cout<<"V border je: "<<endl;
  for (int i = 0; i < border.size(); ++i)
  {
    cout<<border[i].x << " "<<border[i].y<<" "<<border[i].z<<endl;
  }
  stl_vertex cont = border.back();
  border.pop_back();
  stl_vertex end = border.back();
  border.pop_back();
  polylines.resize(5);
  pushToPolylines(polylines[numOfPolylines],end);
  pushToPolylines(polylines[numOfPolylines],cont);
  eps = 1e-24;
  minEps = numeric_limits<double>::max(); 
  bool found = true;
  stl_vertex tmp1,tmp2;

  while(border.size()!=0)
  { 


    if(!found) 
    {
      eps = minEps*1.005;
      if( eps <=0.5 && vertexEqual(cont,end))//(cont == end) )//&& border.size()>1)
      {
        cout << "END NALEZEN, eps= "<<eps<<endl;
        cout << "end: "<<end.x <<" "<<end.y   <<" "<<end.z<<endl; 
        cout << "cont: "<<cont.x <<" "<<cont.y<<" "<<cont.z<<endl; 
        if(!polylines[numOfPolylines].empty()) polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
        
        if(border.size()>0)
        {
          numOfPolylines++;
          if(numOfPolylines+5 > polylines.size()) 
            polylines.resize(numOfPolylines*10);
          end = border.back();
          border.pop_back();
          cont = border.back();
          border.pop_back();
          pushToPolylines(polylines[numOfPolylines],end);
          pushToPolylines(polylines[numOfPolylines],cont);
        }
        //nevim jestli jsou ty posledni 4 radky spravne.. byl tu break a nic, ale neungovalo to na duplu 0.1 0.1 1 0.9
        eps = 1e-24;
        minEps = numeric_limits<double>::max();
        //cout<<"BREAK"<<endl;
        found = true;
        continue;//break;
      }
    }
    //todo tohle by melo byt nad tim nahore
    if(eps > 0.5) // ignoring edges, if it gets to this, there was probably a problem with mesh
    {
      eps = 1e-24;
      minEps = numeric_limits<double>::max();
      //
      //       RETURN FALSE Could prevent poly2tri sefault 
      // 
      std::cerr<<"Unable to find a connected edge, mesh might be invalid"<<endl;
      //
      cont = border.back();
      border.pop_back();
      end = border.back();
      border.pop_back();
      if(!polylines[numOfPolylines].empty())
        numOfPolylines++;
      pushToPolylines(polylines[numOfPolylines],end);
      pushToPolylines(polylines[numOfPolylines],cont);
       // if we didnt found it with 0.5 tolerance.. lets just skip it
    }


    for (int i = border.size()-1;i >= 0; i-=2)
    {
     
      found = false;
      tmp1 = border[i]; 
      tmp2 = border[i-1];
      
      cout << "tmp1: "<<tmp1.x <<" "<<tmp1.y<<" "<<tmp1.z<<endl; 
      cout << "tmp2: "<<tmp2.x <<" "<<tmp2.y<<" "<<tmp2.z<<endl; 
      cout << "end: "<<end.x <<" "<<end.y   <<" "<<end.z<<endl; 
      cout << "cont: "<<cont.x <<" "<<cont.y<<" "<<cont.z<<endl; 
      
      if(vertexEqual(tmp1,cont) || vertexEqual(tmp2,cont))// tmp1 == cont||tmp2 == cont) // ve found next vertex in polyline
      {
        cout<<"pokracovani nalezeno"<<endl;
          found = true;
          //eps = 1e-24;
        if(vertexEqual(tmp1,cont)) //tmp1 == cont) 
        {            
          pushToPolylines(polylines[numOfPolylines],tmp2);
          cont = tmp2;
          border.erase(border.begin()+(i-1),border.begin()+i+1); // erase doesnt delete the last element IT DOES with set
        }  
        else
        {
          pushToPolylines(polylines[numOfPolylines],tmp1);
          cont = tmp1;
          border.erase(border.begin()+(i-1),border.begin()+i+1); // again, different erase with set
        }
        if(vertexEqual(cont,end)) //(cont == end) )//&& border.size()>1)
        {
          //cout << "END NALEZEN"<<endl;
          if(!polylines[numOfPolylines].empty()) polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
          
          if(border.size()>0)
          {
            numOfPolylines++;
            if(numOfPolylines+5 > polylines.size())
              polylines.resize(numOfPolylines*10);
            end = border.back();
            border.pop_back();
            cont = border.back();
            border.pop_back();
            pushToPolylines(polylines[numOfPolylines],end);
            pushToPolylines(polylines[numOfPolylines],cont);
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
  }
  polylines.resize(numOfPolylines+1);

  cout<<"////////////////////////////"<<endl;
  cout<<"Vypis polylines: "<<"numOfPolylines je "<< numOfPolylines<<endl;

  for (int k = 0; k < polylines.size(); ++k)
  { 
    cout<<"Polyline cislo: "<<k<<" size je "<<polylines[k].size()<<endl;
    for (int i = 0; i < polylines[k].size(); i++)
    {
      cout<<"translate(["<<polylines[k][i]->x<<", "<<polylines[k][i]->y<<",0]) sphere(r = 0.3);"<<endl;//" //"<<polylines[k][i+1]->x<<" "<<polylines[k][i+1]->y<<endl;
    }
  }
  

  return true;
}
*/
stl_facet Mesh::createFacet(stl_facet facet, int s, int i, stl_vertex intersect)
{
  stl_facet tmp_facet = facet;
  tmp_facet.vertex[0] = facet.vertex[s];
  tmp_facet.vertex[1] = facet.vertex[(s+i)%3];
  tmp_facet.vertex[2] = intersect;
  return tmp_facet;
}

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
void Mesh::setRemovedAxis()
{
  removedAxis = 'z';
  if(ABS(plane.x) >= ABS(plane.y) && ABS(plane.x) >= ABS(plane.z) )
    removedAxis = 'x';
  if(ABS(plane.y) >= ABS(plane.x) && ABS(plane.y) >= ABS(plane.z) )
    removedAxis = 'y';
}

void Mesh::setVertex(stl_vertex& a, stl_vertex& b,const int &s,const stl_facet & facet)
{
  a = facet.vertex[(s+1)%3];
  b = facet.vertex[(s+2)%3];
}
bool Mesh::haveEqualEdges(tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet1, tuple<stl_facet,stl_position,stl_vertex,stl_vertex>& facet2)
{
  /*
  int samePoints=0;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if((get<0>(facet)).vertex[i] == (get<0>(facet)).vertex[j] )//facet.first.vertex[i] == facet.first.vertex[j])
      {
        samePoints++;
      }
    }
  }*/
  //ifs its on, all 3 points were on and we have to try to match all 3 edges of facet1 to 3 edges of facet2 
  if( (get<1>(facet1) == on && get<1>(facet2)!= on))// || (get<1>(facet2) == on && get<1>(facet1)!= on) )//|| get<1>(facet2) == on)
  {
    pair<stl_vertex,stl_vertex> a,b;
    for (int i = 0; i < 3; ++i)
    {
      a = make_pair(get<0>(facet1).vertex[i] ,       get<0>(facet1).vertex[(i+1)%3]); 
      b = make_pair(get<0>(facet1).vertex[(i+1)%3] , get<0>(facet1).vertex[i]); 
     if(a == (make_pair(get<2>(facet2),get<3>(facet2))) || b == (make_pair(get<2>(facet2),get<3>(facet2))))
      {cout<<"joo"<<endl;return true;}
    }
  }
  if( (get<1>(facet2) == on && get<1>(facet1)!= on))// || (get<1>(facet2) == on && get<1>(facet1)!= on) )//|| get<1>(facet2) == on)
  {
    pair<stl_vertex,stl_vertex> a,b;
    for (int i = 0; i < 3; ++i)
    {
      a = make_pair(get<0>(facet2).vertex[i] ,       get<0>(facet2).vertex[(i+1)%3]); 
      b = make_pair(get<0>(facet2).vertex[(i+1)%3] , get<0>(facet2).vertex[i]); 
     if(a == make_pair(get<2>(facet1),get<3>(facet1)) || b == make_pair(get<2>(facet1),get<3>(facet1)))
      {cout<<"joo2"<<endl;return true;}
    }
  }
  // tests if its a same edge
  if( ((get<2>(facet1) == get<2>(facet2)) && (get<3>(facet1) == get<3>(facet2))) || ((get<3>(facet1) == get<2>(facet2) && get<2>(facet1) == get<3>(facet2))))
     {
      cout<<(get<2> (facet1)).x<<" "<<(get<2> (facet1)).y<<" "<<(get<2> (facet1)).z <<" | "<<(get<3> (facet1)).x<<" "<<(get<3> (facet1)).y<<" "<<(get<3> (facet1)).z<<endl;
            cout<<(get<2> (facet2)).x<<" "<<(get<2> (facet2)).y<<" "<<(get<2> (facet2)).z <<" | "<<(get<3> (facet2)).x<<" "<<(get<3> (facet2)).y<<" "<<(get<3> (facet2)).z<<endl;

      cout<<"joo3"<<endl;return true;
    }
  return false;
}
void Mesh::insertTo(stl_vertex x, stl_vertex y, vector<stl_vertex>& a, set<stl_vertex,setVertComp> & b)
{
  a.push_back(x);
  a.push_back(y);
  b.insert(x);
  b.insert(y);
}
void Mesh::processOnFacets()
{
  vector<stl_vertex>botBorder,topBorder;
  for (int i = 0; i < facetsOnPlane.size(); ++i)
  {
    for (int j = i+1; j < facetsOnPlane.size(); ++j)
    {
     if( haveEqualEdges(facetsOnPlane[i],facetsOnPlane[j]) )
     // if(result)
      {
        //auto k=facetsOnPlane[i].second; auto l=facetsOnPlane[j].second;
        auto k = get<1>(facetsOnPlane[i]); auto l = get<1>(facetsOnPlane[j]);
        if( (k == below && l == above) || (k == above && l == below) )//|| (k == on && l != on) || (l == on && k != on) )
        {
          cout<<(get<2> (facetsOnPlane[i])).x<<" "<<(get<2> (facetsOnPlane[i])).y<<" "<<(get<2> (facetsOnPlane[i])).z <<" | "<<(get<3> (facetsOnPlane[i])).x<<" "<<(get<3> (facetsOnPlane[i])).y<<" "<<(get<3> (facetsOnPlane[i])).z<<endl;
          cout<<(get<2> (facetsOnPlane[j])).x<<" "<<(get<2> (facetsOnPlane[j])).y<<" "<<(get<2> (facetsOnPlane[j])).z <<" | "<<(get<3> (facetsOnPlane[j])).x<<" "<<(get<3> (facetsOnPlane[j])).y<<" "<<(get<3> (facetsOnPlane[j])).z<<endl;
          cout <<"above/below do borderu"<<endl;
          insertTo(get<2> (facetsOnPlane[i]), get<3> (facetsOnPlane[i]), border, border2);
        }
        if((k == on && l == below) || (l == on && k == below))
        {
          /*potrbeuju volat samostatne createborderpolyline findholes a triangulate
            to znamena zapamatovat si border
          */
          cout<<"on/below na border"<<endl;
          int pos = (k == on) ? pos=j:pos=i;
          insertTo(get<2>(facetsOnPlane[pos]), get<3>(facetsOnPlane[pos]), botBorder, border2);
          //cout<<"on a neon"<<endl;
        }
         if((k == on && l == above) || (l == on && k == above))
        {
          cout<<"on/above na border"<<endl;
          int pos = (k == on) ? pos=j:pos=i;
          insertTo(get<2>(facetsOnPlane[pos]), get<3>(facetsOnPlane[pos]), topBorder, border2);
          //cout<<"neon a on"<<endl;
        }
        //if(facetsOnPlane[i].second == above && facetsOnPlanep[j] == below || facetsOnPlane[i].second == below && facetsOnPlanep[j] == on ||) //if(facetsOnPlane[i].second == facetsOnPlane[j].second  && facetsOnPlane[i] != on )

      }
    }
  }
  cout<<"Border size: "<<border.size()<<endl;
  vector<stl_vertex> borderBackUp = border;
  if(botBorder.size()>0)
  {
    border = botBorder;
    cout<<"Bot size: "<<border.size()<<endl;
    if(createBorderPolylines(false))
    {
      findHoles();
      triangulateCut(-1);
    }
  }
  if(topBorder.size()>0)
  {
    border = topBorder;
    cout<<"top size: "<<border.size()<<endl;
    if(createBorderPolylines(false))
    {
      findHoles();
      triangulateCut(1);
    }
  }
  border = borderBackUp;
  cout<<"Aorder size: "<<border.size()<<endl;
}

void Mesh::cut(stl_plane plane)
{
  setPlane(plane);
  setRemovedAxis();
  size_t aboves = 0;
  size_t belows = 0;
  size_t ons = 0;
  stl_position pos[3];
  stl_vertex a,b;
  for (size_t k = 0; k < mesh_file.stats.number_of_facets; k++)
  {
    stl_facet facet = mesh_file.facet_start[k];
    aboves = belows = ons = 0;
    for (int i = 0; i < 3; ++i)
    {
      pos[i] = vertex_position(facet.vertex[i]);
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
    { // last 2 vertexes are not important in case on ons == 3
      facetsOnPlane.push_back(make_tuple(facet,pos[0],facet.vertex[0],facet.vertex[1]));
      continue;
    }
    //if(ons == 3)
    //{
      //continue;     
    //}
    if(ons == 1)
    {
      if(belows == 1 && aboves == 1)
      {
        for (int s = 0; s < 3; ++s)
        {
          if(pos[s] == on)
          {
            setVertex(a,b,s,facet);
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
        continue;
      }
    }
    
    if(ons == 2)
    { 
      cout<<"ONS"<<endl;
      for (int n = 0; n < 3; ++n)
      { 
        if(pos[n] == below)  
        {cout<<"below"<<endl;
          facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));
          bot_facets.push_back(facet); break;   
        }
        if(pos[n] == above) 
        {cout<<"above"<<endl;
          facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));      
          top_facets.push_back(facet); break;
        }
        /*
        if(pos[n] == on)
        {
          facetsOnPlane.push_back(make_tuple(facet,pos[n],facet.vertex[(n+1)%3],facet.vertex[(n+2)%3]));
          break;
        }*/
      }
      continue;
    }
    if(aboves == 1 && belows == 2) // last possibility... the plane cuts the triangle and doesnt intersect with any vertex,
    {
      for (int s = 0; s < 3; ++s)
      {
        if(pos[s] == above) 
        {
          setVertex(a,b,s,facet);
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
    }
    else // below == 1, above == 2
    {
      for (int s = 0; s < 3; ++s)
      {
        if(pos[s] == below) 
        {
          /*stl_vertex a,b;
          a = facet.vertex[(s+1)%3];
          b = facet.vertex[(s+2)%3];*/
          setVertex(a,b,s,facet);
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
  }

  for (unsigned int i = 0; i < border.size(); ++i)
  {
    border2.insert(border[i]);
  }
  //cout<<"BORDER SIZE: "<<border.size()<<endl;
}

 // returns the position of the vertex related to the plane
  stl_position Mesh::vertex_position(stl_vertex vertex) 
  {
    double result = plane.x*vertex.x + plane.y*vertex.y + plane.z*vertex.z + plane.d;
    //cout<<"x: "<<vertex.x<<" y: "<<vertex.y<<" z: "<<vertex.z<<" res: "<<result<<endl;
    //cout<<"px: "<<plane.x<<" py: "<<plane.y<<" pz: "<<plane.z<<" pd: "<<plane.d<<endl;
    if (result > 0) return above;
    if (result < 0) return below;
    return on;
  }

  /*
  *Calculates vertex where edge intersects with plane.
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
          {b = a; a = tmp;}
    stl_vector ab; // vector from A to B
    ab.x = b.x-a.x;
    ab.y = b.y-a.y;
    ab.z = b.z-a.z;
    double t = - (a.x*plane.x + a.y*plane.y + a.z*plane.z + plane.d) / (ab.x*plane.x + ab.y*plane.y + ab.z*plane.z);
    stl_vertex result;
    result.x = a.x + ab.x*t;
    result.y = a.y + ab.y*t;
    result.z = a.z + ab.z*t;
    //cout<<result.x<<" "<<result.y<<" "<<result.z<<endl;
    return result;
  }
void Mesh::open( stl_file file)
{
  mesh_file = file;
}

void Mesh::open( char * name)
{
  stl_open(&mesh_file, name);
  stl_exit_on_error(&mesh_file);
}

stl_file* Mesh::export_stl2(deque<stl_facet> facets) 
{
  stl_file* stl_out = new stl_file;
  initializeStl(stl_out,facets.size());
  /*
  stl_out->stats.type = inmemory;
  stl_out->stats.number_of_facets = facets.size();
  stl_out->stats.original_num_facets = (*stl_out).stats.number_of_facets;
  stl_out->v_indices = NULL;
  stl_out->v_shared = NULL;
  stl_out->neighbors_start = NULL;
  stl_clear_error(stl_out);
  stl_allocate(stl_out);
  */
  int first = 1;
  for (deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++) 
  {
    stl_out->facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(stl_out, *facet, first);
    first = 0;
  }
  /*
   stl_out->stats.degenerate_facets = 0;
   stl_out->stats.edges_fixed = 0;
   stl_out->stats.facets_removed = 0;
   stl_out->stats.facets_added = 0;
   stl_out->stats.facets_reversed = 0;
   stl_out->stats.backwards_edges = 0;
   stl_out->stats.normals_fixed = 0;
   */
  stl_repair(stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//reverse normal directions
  return stl_out;
}
void Mesh::save()
{
  string name="";
  if(!acquireSaveName(name))
    name="Cut_Mesh";
  //c_name+"_1.stl"
  export_stl(top_facets,(name+"_1.stl").c_str());//"Cut_Mesh_1.stl");//"pokus1.stl");
  export_stl(bot_facets,(name+"_2.stl").c_str());//"Cut_Mesh_2.stl");//"pokus2.stl");
  cout<<"Files saved to "<<name<<"_1.stl and "<<name<<"_2.stl"<<endl;
}

std::array<stl_file*,2> Mesh::save2()
{
  return{export_stl2(top_facets),export_stl2(bot_facets)};//"pokus1.stl");

}

bool isStringValid(const std::string &str)
{
    return find_if(str.begin(), str.end(), 
        [](char c) { return !(isalnum(c) || (c == ' ') || (c == '_' )); }) == str.end();
}

bool Mesh::acquireSaveName(string& name)
{
  cout<<"Please enter the new name of meshes ( \"_1.stl\" and \"_2.stl\" will be added)"<<endl;
  cout<<"Or pres Enter for default name (Cut_Mesh_1/2.stl)"<<endl;
  while(true)
  {
    getline(cin,name);//cin>>name;
    if(!isStringValid(name))
    {
      cout<<"Plese use only alphanumeric characters, space and underscore"<<endl;
      name="";
      continue;
    }
    if(name != "")
      return true;
    if(name == "")
      return false;
  }
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

void Mesh::export_stl(deque<stl_facet> facets, const char* name) 
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
  
  stl_repair(&stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//fix normal directions
  stl_write_ascii(&stl_out, name, "stlcut");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

//used by ADMeshGUI
std::array<stl_file*,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes)
{
	stl_plane plane = stl_plane(a,b,c,d);
	Mesh mesh;
  mesh.open(*stlMesh);
  mesh.cut(plane); 
  std::array<stl_file*,2> cutMesh;                              
  if(mesh.createBorderPolylines())
    {
      mesh.findHoles();
      mesh.triangulateCut();
      cutMesh = mesh.save2();
      succes = true;
    }
    else succes = false;
  
  return{cutMesh[0],cutMesh[1]};

}

bool Mesh::t_setRemovedAxis()
{ 
  vector<int> fails;
  int num=1;
  bool succes = true;
  cout<<"Testing setRemovedAxis"<<endl;
  stl_plane plane = stl_plane(1 , 0 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
    {succes = false; fails.push_back(num);}
  num++;
  plane = stl_plane(0 , 1 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'y')
    {succes = false; fails.push_back(num);}
  num++;
  plane = stl_plane(0 , 0 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'z')
    {succes = false; fails.push_back(num);}
  num++;
  plane = stl_plane(1 , 0.5 , 0 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
    {succes = false; fails.push_back(num);}
  num++;
  plane = stl_plane(1 , 1 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'y')
    {succes = false; fails.push_back(num);}
  num++;
  plane = stl_plane(1 , 0 , 1 , 0);
  setPlane(plane);
  setRemovedAxis();
  if(removedAxis != 'x')
    {succes = false; fails.push_back(num);}

  if(fails.size()>0)
  {
    cout<<"   Test/s number: ";
    for (int i = 0; i < fails.size(); ++i)
    {
      cout<<fails[i]<<" ";
    }
    cout<<"failed."<<endl;
  }
return succes;
}


bool Mesh::t_minMaxPointsSame()
{
  stl_get_size(&mesh_file);
  double maxX=mesh_file.stats.max.x;
  cout<<maxX;

}


bool Mesh::runUnitTests()
{
  bool succes = true;
  if( !(t_setRemovedAxis()) ) succes = false;
  return succes; 
}
  
