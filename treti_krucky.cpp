#include <iostream>
#include <deque>
#include <vector>
//#include <utility>
#include <set>
#include <math.h>
#include <map>
#include <algorithm>
#include <admesh/stl.h>
 //#include "/home/debian/admesh-0.98.2/src/stl.h"
#include "poly2tri/poly2tri.h"
#include <string>
 //#include "common/shapes.h"

using namespace p2t;
using namespace std;
double EPS=1e-25;
double minEPS=999999999;
char removedAxis='z';

// vertex position related to the plane
enum stl_position { above, on, below };

// makes more sense in this program
typedef stl_vertex stl_vector;

struct stl_plane
{
  float x;
  float y;
  float z;
  float d;
  stl_plane(float x, float y, float z,float d) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = d;
  }
};

// "=="
bool operator==(stl_vertex a, stl_vertex b) {
  float tolerance=EPS;//1e-25;//1e-15;
  double tol1= ABS(a.x-b.x);
  double tol2= ABS(a.y-b.y);
  double tol3= ABS(a.z-b.z);
  if(removedAxis=='x') { tol1=tol3;a.x=a.z; b.x=b.z; }
  if(removedAxis=='y') { tol2=tol3;a.y=a.z; b.y=b.z; }
  double tmp=tol1>tol2?tol1:tol2;
  if(tmp < minEPS) minEPS=tmp;
  return (ABS((float)(a.x-b.x))<tolerance && ABS((float)(a.y-b.y))<tolerance  );//&& ABS((float)(a.z-b.z))<tolerance);


}


void vertToCout(stl_vertex a)
{
  cout<<a.x<<" "<<a.y<<" "<<a.z;
}

/*
*Checks duplicit vertexes in border
*/
void checkDuplicity(vector<stl_vertex> &border)
{
  //xcout<<"CHECK DUPLICITY"<<endl;
  for (int i = 0; i < border.size(); i+=2)
  {
    for (int j = i+2; j < border.size(); j+=2)
    {
      stl_vertex tmp1 = border[j];
      stl_vertex tmp2 = border[j+1];
      stl_vertex tmp3 = border[i];
      stl_vertex tmp4 = border[i+1];
      if ((tmp1.x==tmp3.x && tmp1.y==tmp3.y &&tmp1.z==tmp3.z && tmp2.x==tmp4.x && tmp2.y==tmp4.y &&  tmp2.z==tmp4.z) || ( tmp2.x==tmp3.x && tmp2.y==tmp3.y &&tmp2.z==tmp3.z && tmp1.x==tmp4.x && tmp1.y==tmp4.y && tmp1.z==tmp4.z ) )
      {
        border.erase(border.begin()+j,border.begin()+(j+2));
        break;
      }
    }
  }
}

/*
* Compares vertexes based of removed axis
*/
struct setVertComp 
 {
    bool operator() (const stl_vertex& lhs, const stl_vertex& rhs) const
    {

      if(removedAxis=='z')
      {
        if(lhs.x<rhs.x) return true;
        else 
          {
            if(lhs.x!=rhs.x) return false;
            if(lhs.y<rhs.y) return true;
             else 
             {
                if(lhs.y!=rhs.y) return false;
                if(lhs.z<rhs.z) return true;
             }
          } 
        return false;
      }

      if(removedAxis=='x')
      {
        if(lhs.y<rhs.y) return true;
        else 
          {
            if(lhs.y!=rhs.y) return false;
            if(lhs.z<rhs.z) return true;
             else 
             {
                if(lhs.z!=rhs.z) return false;
                if(lhs.x<rhs.x) return true;
             }
          } 
        return false;
      }

      if(removedAxis=='y')
      {
        if(lhs.x<rhs.x) return true;
        else 
          {
            if(lhs.x!=rhs.x) return false;
            if(lhs.z<rhs.z) return true;
             else 
             {
                if(lhs.z!=rhs.z) return false;
                if(lhs.y<rhs.y) return true;
             }
          } 
        return false;
      }

      //return lhs<rhs;
      
    }
 };




class Mesh
{
public:
   //Mesh();
  ~Mesh();
  void cut(stl_plane plane);
  stl_position vertex_position(stl_vertex vertex);
  stl_vertex intersection(stl_vertex a, stl_vertex b);
  //bool vertexInPolygon(int nvert, vector<Point* >& vertex,  double testx, double testy);
  void open(char* name);
  void export_stl(deque<stl_facet> facets, const char* name);
  void save();
  void close(){stl_close(&mesh_file);}
  void text();
  bool createBorderPolylines();
  void tmpBorder();
  void findHoles();
  void triangulateCut();
  void remaining();
private:
  double calculatePolygonArea(vector<p2t::Point*> polygon);
  bool vertexInPolygon(int nvert, const vector<Point* >& vertex,  const double &testx, const double testy);
  void createFaces(vector<p2t::Triangle*> &triangles);
  void setMissingCoordinate(const p2t::Point* a,stl_vertex & b);
  void pushToPolylines(vector<p2t::Point*> &vec,stl_vertex vert);
  void volumeTest();
  
  //void calculateVertexCoordinate(stl_vertex &b);
  stl_file mesh_file;
  stl_plane plane=stl_plane(0,0,0,0);
  deque<stl_facet>top_facets,bot_facets;
  vector<stl_vertex>border;
  set<stl_vertex,setVertComp> border2;
  deque<stl_vertex>remainingBorder;
  //set< pair<stl_vertex,stl_vertex> > border2;
  //asi misto toho potrebuju hashmapu/set pairu.. a davat tam rovnou dvojice tech bodu
  //vector<p2t::Point*> polylines;
  vector<vector <p2t::Point*> > polylines;
  vector< vector< pair <vector<p2t::Point*>,int> > > polygonsWithHoles; 
  // this vector contains vectors which contains pair<polygon, -1 for polygon and positive number representing in which polygon is this hole>

  //vector< vector< vector<p2t::Point*> > > polygonsWithHoles;
  int numOfPolylines=0;
  float zCoord;
  //char removedAxis='z';
  
};
/*
void Mesh::calculateVertexCoordinate(stl_vertex &b)
{
  //double coord;
  if(removedAxis=='x')
  {
    b.x= (plane.y*b.y+plane.z*b.z+plane.d) / ((-1.0)* plane.x);
  }
  if(removedAxis=='y')
  {
    b.y= (plane.x*b.x+plane.z*b.z+plane.d) / ((-1.0)* plane.y);  
  }
  if(removedAxis=='z')
  {
    b.z= (plane.x*b.x+plane.y*b.y+plane.d) / ((-1.0)* plane.z);
  }
  
  set<stl_vertex>::iterator itlow,itup;
  itlow=border2.lower_bound (b);                
  itup=border2.upper_bound (b);
  itup--;
  if (itlow==itup) // we found it
  {
    b=(*itlow);
  }   
  
}
*/
Mesh::~Mesh()
{
  for (int i = 0; i < polylines.size(); ++i)
  {
    for (int j = 0; j< polylines[i].size(); ++j)
    {
      //delete (p2t::Point*)polylines[i][j];
    }
  }
}
void Mesh::remaining()
{
  cout<<"REMAINING in border"<<endl;
  for (int i = 0; i < remainingBorder.size(); ++i)
  {
    vertToCout(remainingBorder[i]); //xcout<<endl;
  }

}
/*
*Calculates missing coordinate and then tryes to find it in the border to remove possible inaccuracy
*/
void Mesh::setMissingCoordinate(const p2t::Point* a,stl_vertex &b)
{
  /*
  kdyz vynecham X, prohodil jsem ho se Z
  kdyz vynecham Y, prohodil jsem ho se Z
  kdyz vynecham Z.. vynechal jsem Z
  */
  if(removedAxis=='x')
  {
    b.y=a->y;
    b.z=a->x;
    b.x= (plane.y*b.y+plane.z*b.z+plane.d) / ((-1.0)* plane.x);
    //return true;
  }
  if(removedAxis=='y')
  {
    b.x=a->x;
    b.z=a->y;
    b.y= (plane.x*b.x+plane.z*b.z+plane.d) / ((-1.0)* plane.y);
    //return true;
  }
  if(removedAxis=='z')
  {
    b.x=a->x;
    b.y=a->y;
    b.z= (plane.x*b.x+plane.y*b.y+plane.d) / ((-1.0)* plane.z);
    //return true;
  }
  //return false;
  
  set<stl_vertex>::iterator itlow,itup;
  itlow=border2.lower_bound (b);   //ukazuje na nej(je tam) nebo za nej             
  itup=border2.upper_bound (b);   //vzdy za nej
  itup--;
  if (itlow==itup||itlow==border2.begin()) // we found it
  {
    b=(*itlow);
  }
  else
  {
    itup++;
    //xcout<<"Hledany bod ";vertToCout(b);//xcout<<endl;
    //xcout<<"Nalezeny bod1: ";vertToCout(*itlow);//xcout<<endl;
    //xcout<<"Nalezeny bod2: ";vertToCout(*itup);//xcout<<endl;
    itlow--;
    itup++; 
    //xcout<<"PRED     bod1: ";vertToCout(*itlow);//xcout<<endl;
    //xcout<<"PO     bod2: ";vertToCout(*itup);//xcout<<endl;
    ////xcout<<"POUZIVAM VYPOCITANOU HODNOTU"<<endl;
    itup--;
    double dif1= ABS((*itlow).x-b.x + (*itlow).y-b.y + (*itlow).z-b.z);
    double dif2= ABS((*itup).x -b.x + (*itup).y -b.y + (*itup).z -b.z);
    if(dif1<dif2 && itlow!=border2.end() && dif1<0.1) // find the closest
    {
      b=(*itlow); 
      //xcout<<"Pouzivam    ";vertToCout(*itlow);//xcout<<endl;
      //if(itlow==border2.end())//xcout<<"Protoze to ukazovalo na end"<<endl;
    }
    else if(dif2 <0.1) 
    {
      b=(*itup);
      //xcout<<"Pouzivam    ";vertToCout(*itup);//xcout<<endl;
    }

  }
}


void Mesh::volumeTest()
{
    
    stl_calculate_volume(&mesh_file);
    double org_volume=mesh_file.stats.volume;
    double volume=0.0;
    stl_file cut_mesh;
    stl_open(&cut_mesh, (char*)"Cut_Mesh_1.stl");
    stl_calculate_volume(&cut_mesh);
    //stl_stats_out(&stl_in, stdout, input_file);
    //xcerr<<"OBJEM JE: "<<mesh_file.stats.volume<<endl;
    volume+=mesh_file.stats.volume;
    stl_close(&cut_mesh);
    stl_open(&cut_mesh, (char*)"Cut_Mesh_2.stl");
    stl_calculate_volume(&cut_mesh);
    //stl_stats_out(&stl_in, stdout, input_file);
    //xcerr<<"OBJEM JE: "<<mesh_file.stats.volume<<endl;
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

  if(removedAxis=='x')
  {
    vert.x=vert.z;
  }
  if(removedAxis=='y')
  {
    vert.y=vert.z;
  }

vec.push_back(new p2t::Point(vert.x,vert.y));
}


/*
*Self-explenatory
*/
double Mesh::calculatePolygonArea(vector<p2t::Point*> polygon)
{
    int n = polygon.size();
    int j=0;
    double area = 0.0;
    //shoelace algorithm to calculate polygon area
    for (int i = 0; i < n; ++i)
    {
      j=(i+1)%n;
      area += polygon[i]->x * polygon[j]->y;
      area -= polygon[j]->x * polygon[i]->y;
      //xcout<<"tmpArea: "<<area<<endl;
      //xcout<<"polyhon[i]->x: "<<polygon[i]->x  << "polyhon[i]->y: "<<polygon[i]->y<<endl; 
      //xcout<<"polyhon[j]->x: "<<polygon[j]->x  << "polyhon[j]->y: "<<polygon[j]->y<<endl; 

    }
    area = abs(area) / 2.0;
    return area;
}
/*
*Tests if given vertex is inside given polygon 
* used to find holes
*/
bool Mesh::vertexInPolygon(int nvert, const vector<Point* >& vertex,  const double &testx, const double testy)
{
  int i, j;
  bool in = false;
  //xcout<<"Point k otestovani: "<<testx<<" "<<testy<<endl;
  //xcout<<"Proti pointu      : "<<vertex.at(0)->x<<" "<<vertex.at(0)->y<<endl;
  //vector from point to the right(to infinity) vs edges
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((vertex.at(i)->y > testy) != (vertex.at(j)->y > testy)) &&
     (testx < (vertex.at(j)->x - vertex.at(i)->x) * (testy-vertex.at(i)->y) / (vertex.at(j)->y - vertex.at(i)->y) + vertex.at(i)->x) )
    {
       in = !in;
       //xcout<<"Otacim IN"<<endl;
    }
  }
  return in;
}


void Mesh::findHoles()
{
  //xcout<<"FIND HOLES"<<endl;
  vector<double> polygonArea;
  polygonArea.resize(polylines.size());
  for (int i = 0; i < polylines.size(); ++i)
  {
    polygonArea[i]=calculatePolygonArea(polylines[i]);    
  }

  vector<int> areaOrder;
  for (int i = 0; i < polylines.size(); ++i)
  {
    areaOrder.push_back(i);
  }
  sort(areaOrder.begin(),areaOrder.end(),[&polygonArea](const int &a, const int &b)->bool{return polygonArea[a]<polygonArea[b];});
  //AreaOrdernow contains sorted indexes of polygons based on their area
  
  vector<vector <p2t::Point*> > tmpPolylines;
  tmpPolylines.resize(polylines.size());
  for (int i = 0; i < polylines.size(); ++i)
  {
    //if(polygonArea[areaOrder[i]]!=0)
    {
      tmpPolylines[i]=polylines[areaOrder[i]];
    }  
  }
  vector<int> indexes;
  for (int i = 0; i < tmpPolylines.size(); ++i)
  {
    
    if(polygonArea[areaOrder[i]]==0)
    {
      //xcout<<"PRED VYMAZEM SIZE: "<<tmpPolylines.size();
      //vector<stl_vertex> brokenPolygon = (tmpPolylines[i]);
      //remainingBorder.insert(remainingBorder.end(),tmpPolylines[i].begin(),tmpPolylines[i].end());   
      //tmpPolylines.erase(tmpPolylines.begin()+1); 
      indexes.push_back(i);
      //xcout<<"PO VYMAZU SIZE: "<<tmpPolylines.size();
    }
  }
  sort(indexes.begin(),indexes.end());
  for (int i = indexes.size()-1; i >= 0; --i)
  {
     for (int n = 0; n < tmpPolylines[indexes[i]].size(); ++n)
     {
       delete tmpPolylines[indexes[i]][n];
     }
     //erasing wrong polygons with area=0
     tmpPolylines.erase(tmpPolylines.begin()+indexes[i]);
  }
  //tmpPolylines contains polylines/polygons sorted based of their area.
  polylines=tmpPolylines; 
  /*
  for (int i = 0; i < areaOrder.size(); ++i)
  {
    //xcout<<"area "<<i<<": " <<polygonArea[areaOrder[i]]<<endl;
  }
  */

  vector<p2t::Point*> tmpPolygon=polylines.back();
  polylines.pop_back();
  pair<vector<p2t::Point*>,int >tmpPair=make_pair(tmpPolygon,-1);
  vector<pair<vector<p2t::Point*> , int> > tmpVecPair;
  tmpVecPair.push_back(tmpPair);
  polygonsWithHoles.push_back(tmpVecPair);//tmplist);
  int in;
  int out;
  //bool hole=false;
  int holeIn=-1; //-1 for polygon
  bool placeFound=false;
  int pos=0;
  //
  for(;polylines.size()>0;)
  {
    tmpPolygon = polylines.back(); 
    polylines.pop_back();
    pos=0;
    placeFound=false;
    holeIn=-1; 
    for(int k=0;k<polygonsWithHoles.size();k++)
    {
      in=0;
      out=0;
       
      //test 3 points and based on majority(just in case) decice if this polygons is inside biggerPolygon
      for(int j=0;j<3;j++)
      {
        vector<p2t::Point*>  biggerPolygon=polygonsWithHoles[k][pos].first;//.front();
        double a=tmpPolygon[j]->x; //(*tmpPolygon).at(j)->x;
        double b=tmpPolygon[j]->y;//tmpPolygon->at(j)->y;
        if(vertexInPolygon((int)biggerPolygon.size(),biggerPolygon,a,b))//vertexInPolygon((int)tmpPolygon[j].size(),tmpPolygon,a,b))
          in++;
        else
          out++;
      }
      if(in>out || placeFound)
      {
         ////xcout<<"POINT BYL IN"<<endl;
         ////xcout<<"in: "<<in<<" out: "<<out<<endl;
         if(in > out) 
         {
            if(holeIn!=-1)
              holeIn=-1;  //swaping between hole and not a hole
            else
              holeIn=pos;//pos; // asi tam patri 0, jeden pripad to opravilo, druhy ale zacyklilo
          }
         
         if(!placeFound)  
            {
              //polygonsWithHoles[k].push_back(tmpPolygon);
              //xcout<<"POINT BYL IN"<<endl;
            }
            
         placeFound=true;
         //as long as we didnt get to the end of the vector, continue - basicaly
         if(pos < polygonsWithHoles[k].size()-1) // we dont want to test it vs itself
         {
            pos++;
            k--; // we need to stay at the same index, dont want to use goto
         }
         else //we get to the end, make a pair and push it to the end
         {
            tmpPair=make_pair(tmpPolygon,holeIn);
            polygonsWithHoles[k].push_back(tmpPair);
            break;
         }

      }
    }
    if(!placeFound)//out>in) // if we get here and point wasnt it, condition should be always true
    {                        // and we found a new polygon
        //xcout<<"POINT BYL out"<<endl;
        //xcout<<"in: "<<in<<" out: "<<out<<endl;
        //tmpPair=make_pair(tmpPolygon,false);
        //vector<vector<p2t::Point*>  > list2;
        //list2.push_back(tmpPolygon);
        tmpPair=make_pair(tmpPolygon,-1);
        tmpVecPair.erase(tmpVecPair.begin(),tmpVecPair.end());
        tmpVecPair.push_back(tmpPair);
        polygonsWithHoles.push_back(tmpVecPair);//list2);
    }
  }
  
  
  }


void Mesh::triangulateCut()
{
  //cout<<"SNAZIM SE O TRIANGULACI"<<endl;
  //int numOfBasePolygons=0;
  /*projizdej polygonswithholes po jednotlivych listech*/
  //vector<p2t::CDT> polygons;
  map<int, p2t::CDT*> polygons;//vector<p2t::CDT>> polygons;

  for (int i = 0; i < polygonsWithHoles.size(); ++i)
  {
      for (int j = 0; j < polygonsWithHoles[i].size(); ++j)
      {
          if(polygonsWithHoles[i][j].second==-1) // we found polygon, not a hole
          {   
              p2t::CDT* tmp=new  p2t::CDT(polygonsWithHoles[i][j].first);
              polygons.insert ( pair<int,p2t::CDT*>(j,tmp) );
          }
          else // we found a hole
          {
              int holeIn=polygonsWithHoles[i][j].second;
              map<int, p2t::CDT*>::iterator it;
              it = polygons.find(holeIn);
              it->second->AddHole(polygonsWithHoles[i][j].first);
          }
      }

      if(polygons.size()>0)
        for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
        {
            k->second->Triangulate();
            vector<p2t::Triangle*> triangles = k->second->GetTriangles();
            createFaces(triangles);   
        }

      for (map<int, p2t::CDT*>::iterator k = polygons.begin(); k != polygons.end(); ++k)
      {
        delete (*k).second;
      }

      polygons.erase(polygons.begin(),polygons.end());       
  }
}

void Mesh::createFaces(vector<p2t::Triangle*> &triangles)
{
  
    // for each triangle, create facet
    for (vector<p2t::Triangle*>::iterator i = triangles.begin(); i != triangles.end(); i++) {
    stl_vertex vertex;
    stl_facet facet;
    for (size_t j = 0; j < 3; j++) 
    {
      p2t::Point* p = (*i)->GetPoint(j);
      setMissingCoordinate(p,vertex);
      facet.vertex[j]=vertex;
      //facet.vertex[j].x = vertex.x;
      //facet.vertex[j].y = vertex.y;
      //facet.vertex[j].z = vertex.z;
    }

    // normal goes out of the object, for lower part, it is identical to plane normal
    float test_norm[3];
    test_norm[0] = plane.x;
    test_norm[1] = plane.y;
    test_norm[2] = plane.z;
    stl_normalize_vector(test_norm);
    facet.normal.x = test_norm[0];
    facet.normal.y = test_norm[1];
    facet.normal.z = test_norm[2];
    bot_facets.push_back(facet);
    /*
    facet.normal.x = plane.x;
    facet.normal.y = plane.y;
    facet.normal.z = plane.z;
    bot_facets.push_back(facet);
    */
    
    // for the upper part, we need to invert the normal...
    /*facet.normal.x = -plane.x;
    facet.normal.y = -plane.y;
    facet.normal.z = -plane.z;
    */
    facet.normal.x *= -1.0;
    facet.normal.y *= -1.0;
    facet.normal.z *= -1.0;
    // ...and reverse the order of the vertices
    // TODO check if the order of vertices from poly2tri in fact depends on orientation of the first used edge
    // .. and the order might be reversed anyway
    vertex = facet.vertex[1];
    facet.vertex[1] = facet.vertex[2];
    facet.vertex[2] = vertex;
    top_facets.push_back(facet);
    }
   
}

void Mesh::tmpBorder()
{
  for (int i = 0; i <= 0;i++)//numOfPolylines; ++i) vratit, testoval jsem add hole
  {
    
  
    p2t::CDT cdt(polylines[0]);// = p2t::CDT(polygon);  /*cdt polylines i, addhole smazat, testuju
    //cdt.AddHole(polylines[0]);
    cdt.Triangulate();
    vector<p2t::Triangle*> triangles = cdt.GetTriangles();
  
    // for each triangle, create facet
    for (vector<p2t::Triangle*>::iterator i = triangles.begin(); i != triangles.end(); i++) {
    stl_vertex vertex;
    stl_facet facet;
    for (size_t j = 0; j < 3; j++) 
    {
      p2t::Point* p = (*i)->GetPoint(j);
      vertex.x = p->x;
      vertex.y = p->y;
      vertex.z =zCoord;
      facet.vertex[j].x = vertex.x;
      facet.vertex[j].y = vertex.y;
      facet.vertex[j].z = vertex.z;
    }

    // normal goes out of the object, for lower part, it is identical to plane normal
    facet.normal.x = plane.x;
    facet.normal.y = plane.y;
    facet.normal.z = plane.z;
    bot_facets.push_back(facet);
    
    // for the upper part, we need to invert the normal...
    facet.normal.x = -plane.x;
    facet.normal.y = -plane.y;
    facet.normal.z = -plane.z;
    // ...and reverse the order of the vertices
    // TODO check if the order of vertices from poly2tri in fact depends on orientation of the first used edge
    // .. and the order might be reversed anyway
    vertex = facet.vertex[1];
    facet.vertex[1] = facet.vertex[2];
    facet.vertex[2] = vertex;
    top_facets.push_back(facet);
    }
  
  }
}






bool Mesh::createBorderPolylines()
{
  /*
  if (border.size()%2 !=0)    
  {
    //xcout<<"Mesh had incorrect triangle/s with all vertices in one line on the cutting border";
    return false;
  }
  */
  
  if(border.size()==0 /*|| top_facets.size()==0 || bot_facets.size()==0*/) 
  {
    cerr<<"Nothing to cut"<<endl;
    return false;
  }
  /*
  else
  */
  {
    //xcout<<"Mame tady "<<border.size()<<" bodu ve fronte"<<endl;




    int x;
   
    int tmp=0;
    //float eps=10e-13;
    //int numOfPolylines=0;
    //xcout<<"BEFORE REMOVE DUPLICITY: "<<border.size()<<endl;
      for (int m = 0; m < border.size(); m+=2)
    {
      //xcout<<"x:"<<border[m].x<< " "<<border[m].y<< " "<<border[m].z<< " // "<<border[m+1].x<<" "<<border[m+1].y<< " "<<border[m+1].z<<endl;
    }
    checkDuplicity(border);
    //xcout<<"AFTER REMOVE DUPLICITY: "<<border.size()<<endl;
      for (int m = 0; m < border.size(); m+=2)
    {
      //xcout<<"x:"<<border[m].x<< " "<<border[m].y<< " "<<border[m].z<< " // "<<border[m+1].x<<" "<<border[m+1].y<< " "<<border[m+1].z<<endl;
    }


    stl_vertex cont=border.back();
    border.pop_back();
    stl_vertex end=border.back();
    border.pop_back();
    polylines.resize(5);
    //p2t::Point* xx=new p2t::Point(end.x,end.y);
    //vector<p2t::Point*> neco;
    //neco.push_back(xx);
    //cin>>x;
    pushToPolylines(polylines[numOfPolylines],end);
    pushToPolylines(polylines[numOfPolylines],cont);
    //polylines[numOfPolylines].push_back(new p2t::Point(end.x,end.y));
    //polylines[numOfPolylines].push_back(new p2t::Point(cont.x,cont.y));
    /*
    vezmi dva body vraz je do polylinu
    zapamatuj si prvni jako konec
    zapamatuj si druhy jak opokracovni
    hledej dal dokud nenajdes pokracovani
    kdyz ho najdes, vraz tam pokracovani a k nemu nalezejici bod
    ten druhy z nich oznac za pokracovani
    opakuj dokud nenajdes dvojici, ktera ma jako jeden z bodu konec
    */
    
    //xcout<<"Mame tady "<<border.size()<<" bodu ve fronte"<<endl;
    EPS=1e-24;
    minEPS=999999999; //dont want to add limits just for this
    bool found=true;
    stl_vertex tmp1,tmp2;
    while(border.size()!=0)
    {
      
      //xcout<<border.size()<<endl;
      if(!found) EPS=minEPS*1.005;
      if(EPS > 0.5) // ignoring edges, if it gets to this, there was probably a problem with mesh
      {
        //xcout<<"IGNORUJU EDGE, numOfPolylines: "<<numOfPolylines+1<<endl;
        //tmp1=border[border.size()-1];
        //tmp2=border[border.size()-2];
        //border.erase(border.end()-1,border.end()+1);
        EPS=1e-24;
        minEPS=999999999;
        //remainingBorder.insert(remainingBorder.end(),polylines[numOfPolylines].begin(),polylines[numOfPolylines].end());
        //remainingBorder.push_back(tmp2);
        numOfPolylines++; // if we didnt found it with 0.5 tolerance.. lets just skip it
        cont=border.back();
        border.pop_back();
        end=border.back();
        border.pop_back();
        pushToPolylines(polylines[numOfPolylines],end);
        pushToPolylines(polylines[numOfPolylines],cont);
        //polylines[numOfPolylines].push_back(new p2t::Point(end.x,end.y));
        //polylines[numOfPolylines].push_back(new p2t::Point(cont.x,cont.y));


      }

      //cin>>x;
      for (int i = border.size()-1;i >= 0; i-=2)
      {
       
        found=false;
        tmp1=border[i]; 
        tmp2=border[i-1];

        /*
        // tohle predstavuje jeden par ktery pak porovnavam s cont/end
        misto toho musim vzit prez iterator pair, a do tmp1 a tmp2 z nej dostat vertex.. a pak pokracovat stejne?
        */
       //xcout << "tmp1: "<<tmp1.x <<" "<<tmp1.y<<" "<<tmp1.z<<endl; 
       //xcout << "tmp2: "<<tmp2.x <<" "<<tmp2.y<<" "<<tmp2.z<<endl; 
       //xcout << "end: "<<end.x <<" "<<end.y   <<" "<<end.z<<endl; 
       //xcout << "cont: "<<cont.x <<" "<<cont.y<<" "<<cont.z<<endl; 
        if(tmp1==cont||tmp2==cont) // ve found next vertex in polyline
        {

           /*ASI POTREBUJU, prvky co jsou ON tam budou dvakrat, ale kdzy jsem to tam dal, smazalo to neco v obyc box.stl*/
          /*
          if(alreadyIn(polylines,numOfPolylines,tmp1,tmp2,cont,end))
            {
              //xcout<<"mazu duplicitu"<<endl;
              border.erase(border.begin()+(i-1),border.begin()+i+1);
              continue;
            }
          */
            found=true;
            //EPS=1e-24;

          if(tmp1==cont) 
          {            
            pushToPolylines(polylines[numOfPolylines],tmp2);
            //polylines[numOfPolylines].push_back(new p2t::Point(tmp2.x,tmp2.y));
            //xcout<<"Cont se meni z:";
            vertToCout(cont);//xcout<<" na ";vertToCout(tmp2);//xcout<<endl;
            cont=tmp2;
            border.erase(border.begin()+(i-1),border.begin()+i+1); // erase doesnt delete the last element IT DOES with set
            //xcout<<"zapsal jsem1 ";
           
            //xcout<<"size je "<<polylines[numOfPolylines].size()<<endl;
          }  
          else
          {
            pushToPolylines(polylines[numOfPolylines],tmp1);
            //polylines[numOfPolylines].push_back(new p2t::Point(tmp1.x,tmp1.y));
             //xcout<<"Cont se meni z:";vertToCout(cont);//xcout<<" na ";vertToCout(tmp1);//xcout<<endl;
            cont=tmp1;
            border.erase(border.begin()+(i-1),border.begin()+i+1); // again, different erase with set
            
            //xcout<<"zapsal jsem2 ";
           
            //xcout<<"size je "<<polylines[numOfPolylines].size()<<endl;
          }


          if((cont==end) )//&& border.size()>1)
          {
            //xcout<<"nalezen cont==end"<<endl;
            //xcout<<"cont= "<<cont.x<<" "<<cont.y<<" "<<cont.z<<endl;
            //xcout<<"end= "<<end.x<<" "<<end.y<<" "<<end.z<<endl;
            polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
            
            if(border.size()>0)
            {
            numOfPolylines++;
            if(numOfPolylines+5 > polylines.size()) polylines.resize(numOfPolylines*10);
            end=border.back();
            border.pop_back();
            cont=border.back();
            border.pop_back();
            pushToPolylines(polylines[numOfPolylines],end);
            pushToPolylines(polylines[numOfPolylines],cont);
            //polylines[numOfPolylines].push_back(new p2t::Point(end.x,end.y));
            //polylines[numOfPolylines].push_back(new p2t::Point(cont.x,cont.y));
            //xcout<<"Nastavil jsem nove end a cont na: "<<endl;
            //xcout<<"cont= "<<cont.x<<" "<<cont.y<<" "<<cont.z<<endl;
            //xcout<<"end= "<<end.x<<" "<<end.y<<" "<<end.z<<endl;
            //xcout<<"/////////////////////"<<endl;
            }

            //polylines[numOfPolylines].push_back(new p2t::Point(end.x,end.y));
            //polylines[numOfPolylines].push_back(new p2t::Point(cont.x,cont.y));
            minEPS=999999999;
            break;
          }

          else 
          {
            minEPS=999999999.0;
            EPS=1e-24;
            break;
          }
          
        }
        

      }

    }


  }
  //xcout<<"////////////////////////////"<<endl;
  //xcout<<"Vypis polylines: "<<"numOfPolylines je "<< numOfPolylines<<endl;
  for (int k = 0; k < polylines.size(); ++k)
  { 
  //xcout<<"Polyline cislo: "<<k<<" size je "<<polylines[k].size()<<endl;
  for (int i = 0; i < polylines[k].size(); i++)
  {
    //xcout<<polylines[k][i]->x<<" "<<polylines[k][i]->y<<endl;//" //"<<polylines[k][i+1]->x<<" "<<polylines[k][i+1]->y<<endl;
  }
  }
  polylines.resize(numOfPolylines+1);

  return true;
}



void Mesh::text()
{
  //xcout<<"Top Facets:"<<top_facets.size()<<endl;
  for (int i = 0; i < top_facets.size(); ++i)
  {
    vertToCout(top_facets[i].vertex[0]);
    //xcout<<" / ";
    vertToCout(top_facets[i].vertex[1]);
    //xcout<<" /";
    vertToCout(top_facets[i].vertex[2]);
    //xcout<<endl;
  }
  //xcout<<"------------"<<endl;
   //xcout<<"Bot Facets:"<<bot_facets.size()<<endl;
    for (int i = 0; i < bot_facets.size(); ++i)
  {
    vertToCout(bot_facets[i].vertex[0]);
    //xcout<<" / ";
    vertToCout(bot_facets[i].vertex[1]);
    //xcout<<" /";
    vertToCout(bot_facets[i].vertex[2]);
    //xcout<<endl;
  }
  /*
  for (int i = 0; i < bot_facets.size(); ++i)
  {
    //xcout<<vertToCout(bot_facets[i].vertex[0])<<" "<<vertToCout(bot_facets[i].vertex[1])<<" "<<vertToCout(bot_facets[i].vertex[2])<<endl;
  }*/


}
void Mesh::cut(stl_plane plane)
{
  this->plane=plane;
  size_t aboves = 0;
  size_t belows = 0;
  size_t ons = 0;
  if(ABS(plane.x)>=ABS(plane.y) && ABS(plane.x)>=ABS(plane.z) ) removedAxis='x';
  if(ABS(plane.y)>=ABS(plane.x) && ABS(plane.y)>=ABS(plane.z) ) removedAxis='y';
  //xcout<<"REMOVED AXIS"<<removedAxis<<endl;

  stl_position pos[3];

  //xcout<<"Number of Facets: "<<mesh_file.stats.number_of_facets<<endl;
  for (size_t k = 0; k < mesh_file.stats.number_of_facets; k++)
  {
    stl_facet facet=mesh_file.facet_start[k];
    aboves=belows=ons=0;
    for (int i = 0; i < 3; ++i)
    {
      pos[i]=vertex_position(facet.vertex[i]);
      if(pos[i]==above) aboves++;
      if(pos[i]==on)    ons++;
      if(pos[i]==below) belows++;

      /*porovnej bod s rovinnou
      zjisti jak se to sakra pocita
      a pak to vrat do que nad/pod/na
      */
    }
    //xcout<<"A: "<<aboves<<" B:"<<belows<<" O:"<<ons<<endl;
    if(aboves==3)
      {
        top_facets.push_back(facet);
        continue;
      }
    if(belows==3)
      {
        bot_facets.push_back(facet);
        continue;
      }
      //this should never happen
    if(ons==3)
      {
       for (int j = 0; j < 3; ++j)
       {
          //xcout<<"Warning - Some of models triangles are exactly on the cutting plane ,all of those triangles will be deleted"<<endl;//added to the \"bottom\" model"<<endl;
          //bot_facets.push_back(facet);
          //border.push_back(facet.vertex[j]);
       }
       continue;     
      }
      //1/1/1, jednoduchy rozdeleni trojuhelniku na dva
      /*zjistim intersekci vezmu bod na, bod nad a interesekcni a dam to do nad
      pak bod na pod a intersekcni a dam to do pod fronty
      */
      if(ons==1)
      {
        //xcout<<"JEDEN BOD NA "<<endl;
        if(belows==1 && aboves==1)
          {
            //xcout<<"A JEDEN POD A NAD "<<endl;
          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==on)
            {
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect=intersection(a,b);//ziskam ten novy bod
              zCoord=intersect.z;

              if(pos[(s+1)%3]==below) 
              {
                stl_facet tmp_facet=facet;
                tmp_facet.vertex[0]=facet.vertex[s];
                tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
                tmp_facet.vertex[2]=intersect;
                bot_facets.push_back(tmp_facet); // new triangle added to the bot

                tmp_facet.vertex[1]=facet.vertex[(s+2)%3];
                top_facets.push_back(tmp_facet);
                //xcout<<"Do border davam: "<<endl;
                //xcout<<facet.vertex[s].x<<" "<<facet.vertex[s].y<<" // "<<intersect.x<<" "<<intersect.y<<endl;
                border.push_back(facet.vertex[s]);
                border.push_back(intersect);
              }
              
              else
              {
                stl_facet tmp_facet=facet;
                tmp_facet.vertex[0]=facet.vertex[s];
                tmp_facet.vertex[1]=facet.vertex[(s+2)%3];
                tmp_facet.vertex[2]=intersect;
                bot_facets.push_back(tmp_facet); // new triangle added to the bot

                tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
                top_facets.push_back(tmp_facet);
                //xcout<<"Do border davam: "<<endl;
                //xcout<<facet.vertex[s].x<<" "<<facet.vertex[s].y<<" // "<<intersect.x<<" "<<intersect.y<<endl;
                border.push_back(facet.vertex[s]);
                border.push_back(intersect);

              }
              
              break;
            }
          }
          
          



        continue;
        }
      }
      if(ons==2)
      { //bool flag=false;
        //if(belows==1) flag=true;
        for (int n = 0; n < 3; ++n)
        { 
          if(pos[n]==below) {bot_facets.push_back(facet);}
          if(pos[n]==above) {top_facets.push_back(facet);}
          if(pos[n]==on)    {border.push_back(facet.vertex[n]);} // when 2 are on, we add it just once.. it should appear twice during bordering process
        }

        continue;
      }
    /*ted se zacne rezat a vytvaret nove trojuhelniky*/
      /*pokud jsou 2 nad a 1 pod nebo naopak
      tak musim vytvorit dva nove body
      a 3 trojuhelniky*/
      if(ons==1)
      {
        if(belows==2)
        {
          bot_facets.push_back(facet);
          continue;
        }
        if(aboves==2)
        {
          top_facets.push_back(facet);
          continue;
        }
      }
      else // onst==0 --> need to make 3 new triangles
      {
            //xcout<<"Delam 3 triangly"<<endl;
       if(aboves==1&&belows==2) // last possibility... the plane cuts the triangle and doesnt intersect with any vertex, this if is not neccesary
       {
          //xcout<<"prosel jsem podminkou"<<endl;
          /*zjistit ktery je ten horni
          a pak to poslat na rozdeleni
          v rozdeleni musim zjistit kde se to protina mezi prvni a druhym a prvnim a tretim bodem
          kdyz ty body zjistim, musim vytvorit nove trojuhelniky
          a potom je dat do spravne deque*/
         // while(pos[tmp]!=on)
          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==above) // tohle je moemntalne spatne... funguje to jen na 1 nad 2 pod, ne naopak//snad uz dobre
            {
              //xcout<<"nasel jsem above vertex"<<endl;
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect1=intersection(facet.vertex[s],a);//ziskam ten novy bod
              stl_vertex intersect2=intersection(facet.vertex[s],b);//ziskam ten novy bod
              zCoord=intersect1.z;
               //hotfix to failed floating point operations
              if(intersect1.x==intersect2.x && intersect2.x==facet.vertex[s].x && intersect1.y==intersect2.y && intersect2.y==facet.vertex[s].y && intersect1.z==intersect2.z && intersect2.z==facet.vertex[s].z)
              {
                //xcout<<"HOTFIX IN ACTION"<<endl;
                bot_facets.push_back(facet);
                break;
              }
              /*if(hotfixIntersectionEqualsVertex(a,b,pos[s],bot_facets) )
              {
                break;
              }
              */
              stl_facet tmp_facet=facet;
              tmp_facet.vertex[0]=facet.vertex[s];
              tmp_facet.vertex[1]=intersect1;//facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=intersect2;
              top_facets.push_back(tmp_facet); // facet with above triangle and intersecting vertexes added
              tmp_facet.vertex[0]=intersect1;
              tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=intersect2;             
              bot_facets.push_back(tmp_facet); // facet with intersecting vertexes and below vertex 1
              tmp_facet.vertex[0]=intersect2;
              tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=facet.vertex[(s+2)%3];
              bot_facets.push_back(tmp_facet); // facet with intersecting vertexes and below vertex 2
              //xcout<<"Do border davam: "<<endl;
              //xcout<<intersect1.x<<" "<<intersect1.y <<" "<<intersect1.z<<" // "<<intersect2.x<<" "<<intersect2.y<<" "<<intersect2.z<<endl;
              border.push_back(intersect1);
              border.push_back(intersect2);
              break;
            }
           
          }
           
       }
       else // below==1, above == 2
       {
          //xcout<<"prosel jsem podminkou"<<endl;
          /*zjistit ktery je ten horni
          a pak to poslat na rozdeleni
          v rozdeleni musim zjistit kde se to protina mezi prvni a druhym a prvnim a tretim bodem
          kdyz ty body zjistim, musim vytvorit nove trojuhelniky
          a potom je dat do spravne deque*/
         // while(pos[tmp]!=on)
          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==below) // tohle je moemntalne spatne... funguje to jen na 1 nad 2 pod, ne naopak
            {
              //xcout<<"nasel jsem below vertex"<<endl;
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect1=intersection(facet.vertex[s],a);//ziskam ten novy bod
              stl_vertex intersect2=intersection(facet.vertex[s],b);//ziskam ten novy bod
              zCoord=intersect1.z;
              //hotfix to failed floating point operations
             if(intersect1.x==intersect2.x && intersect2.x==facet.vertex[s].x && intersect1.y==intersect2.y && intersect2.y==facet.vertex[s].y && intersect1.z==intersect2.z && intersect2.z==facet.vertex[s].z)
              {
                //xcout<<"HOTFIX IN ACTION"<<endl;
                top_facets.push_back(facet);
                break;
              }
              /*
              if(hotfixIntersectionEqualsVertex(a,b,pos[s],top_facets))
              {
                  break;
              }
              */
              stl_facet tmp_facet=facet;
              tmp_facet.vertex[0]=facet.vertex[s];
              tmp_facet.vertex[1]=intersect1;//facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=intersect2;
              bot_facets.push_back(tmp_facet); // facet with above triangle and intersecting vertexes added
              tmp_facet.vertex[0]=intersect1;
              tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=intersect2;
              top_facets.push_back(tmp_facet); // facet with intersecting vertexes and above vertex 1
              tmp_facet.vertex[0]=intersect2;
              tmp_facet.vertex[1]=facet.vertex[(s+1)%3];
              tmp_facet.vertex[2]=facet.vertex[(s+2)%3];
              top_facets.push_back(tmp_facet); // facet with intersecting vertexes and above vertex 2
              //xcout<<"Do border davam: "<<endl;
              //xcout<<intersect1.x<<" "<<intersect1.y<<" // "<<intersect2.x<<" "<<intersect2.y<<endl;
              border.push_back(intersect1);
              border.push_back(intersect2);
              break;
            }
           
          }
           
       }
      }


  }

  for (int i = 0; i < border.size(); ++i)
  {
    border2.insert(border[i]);
  }
  //xcout<<"KONEC MESH::CUT"<<endl;
}

 // returns the position of the vertex related to the plane
  stl_position Mesh::vertex_position(stl_vertex vertex) 
  {
    //double result = (double)plane.x*vertex.x + (double)plane.y*vertex.y + (double)plane.z*vertex.z + plane.d;
    double result = plane.x*vertex.x + plane.y*vertex.y + plane.z*vertex.z + plane.d;
    //xcout<<"Vert: ";vertToCout(vertex);//xcout<<" result: "<<result<<endl;
    if (result > 0) return above;
    if (result < 0) return below;
    return on;
  }

  stl_vertex Mesh::intersection(stl_vertex a, stl_vertex b) 
  {
    //if(a.x>b.x && a.y>b.y && a.z>b.z){stl_vertex tmp=b; b=a; a=tmp;} tried to fix non-connected edges
    stl_vector ab; // vector from A to B
    ab.x = b.x-a.x;
    ab.y = b.y-a.y;
    ab.z = b.z-a.z;
    double t = - (a.x*plane.x + a.y*plane.y + a.z*plane.z + plane.d) / (ab.x*plane.x + ab.y*plane.y + ab.z*plane.z);
    
    stl_vertex result;
    result.x = a.x + ab.x*t;
    result.y = a.y + ab.y*t;
    result.z = a.z + ab.z*t;

    return result;
  }

void Mesh::open( char * name)
{
  //mesh_file = file;
  stl_open(&mesh_file, name);
  //xcout<<"Error int: "<<stl_get_error(&mesh_file)<<endl;
  stl_exit_on_error(&mesh_file);
}

void Mesh::save()
{
  export_stl(top_facets,"Cut_Mesh_1.stl");//"pokus1.stl");
  export_stl(bot_facets,"Cut_Mesh_2.stl");//"pokus2.stl");
  cout<<"Files saved to Cut_Mesh_1 and Cut_Mesh_2"<<endl;
}

void Mesh::export_stl(deque<stl_facet> facets, const char* name) 
{
  stl_file stl_out;
  stl_out.stats.type = inmemory;
  stl_out.stats.number_of_facets = facets.size();
  stl_out.stats.original_num_facets = stl_out.stats.number_of_facets;
  stl_out.v_indices = NULL;
  stl_out.v_shared = NULL;
  stl_out.neighbors_start = NULL;
  stl_clear_error(&stl_out);
  stl_allocate(&stl_out);
  
  int first = 1;
  for (deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++) {
    stl_out.facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(&stl_out, *facet, first);
    first = 0;
  }
  
  // check nearby in 2 iterations
  // remove unconnected facets
  // fill holes
  //stl_repair(&stl_out, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 0);

  //stl_check_facets_exact(&stl_out);
  stl_write_ascii(&stl_out, name, "stlcut");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

bool tests()
{
  string names[]={"sphere.stl","trubka2.stl","trubka3.stl","trubka4.stl","trubka5.stl","sphere2.stl","hrana2.stl","boxsphere.stl"};
  //char* name=(char*)names[0].c_str();
  //otevrit prvni mesh, testnout volume
  //a pak testnout volume po rezu
  float org_volume;
  ////xcerr<<"TESTUJU"<<endl;
  stl_file mesh_file;
  int test_ok=0;
  int num_of_tests=8;
 
  for (int i = 0; i < num_of_tests;i++)// num_of_tests; ++i)
  {
     org_volume=0.0;
     stl_open(&mesh_file, (char*)names[i].c_str());
     stl_calculate_volume(&mesh_file);
     //stl_stats_out(&stl_in, stdout, input_file);
     //xcerr<<"OBJEM JE: "<<mesh_file.stats.volume<<endl;
     org_volume=mesh_file.stats.volume;
     stl_close(&mesh_file);
  

    Mesh mesh;
    mesh.open((char*)names[i].c_str());
    mesh.cut(stl_plane(1,1,0,0));
    if(mesh.createBorderPolylines())
      {
        ////xcout<<"CreateBorderPolylines proslo"<<endl;
        mesh.findHoles();
        mesh.triangulateCut();
        //mesh.tmpBorder();
      }
    //mesh.findHoles();
    //mesh.text();
    mesh.save();
    mesh.close();

    float volume=0.0;
    stl_open(&mesh_file, (char*)"Cut_Mesh_1.stl");
    stl_calculate_volume(&mesh_file);
    //stl_stats_out(&stl_in, stdout, input_file);
    //xcerr<<"OBJEM JE: "<<mesh_file.stats.volume<<endl;
    volume+=mesh_file.stats.volume;
    stl_close(&mesh_file);
    stl_open(&mesh_file, (char*)"Cut_Mesh_2.stl");
    stl_calculate_volume(&mesh_file);
    //stl_stats_out(&stl_in, stdout, input_file);
    //xcerr<<"OBJEM JE: "<<mesh_file.stats.volume<<endl;
    volume+=mesh_file.stats.volume;
    stl_close(&mesh_file);
    //xcerr<<"Puvodni OBJEM: "<<org_volume<<" SJEDNOCENY OBJEM: "<<volume<<endl;
    if(abs(org_volume-volume) < abs(org_volume/1000.0)) {cerr<<"Test OBJEMU uspesny"<<endl;test_ok++;};
  }
  cerr<<"Testu probehlo: "<<num_of_tests<<" z toho uspesne: "<<test_ok<<endl;
 

  
  
}


int main(int argc, char **argv) {
  if ((argc != 2 && argc!=6 && argc!=1) ||argv[1]==string("help")) {
    cerr << "Usage: " << argv[0] << " file.stl a b c d (plane-optional or 0 0 1 -1 will be used) " << endl<<"or"<<endl<< argv[0]<<" tests"<<endl;
    return 1;
  }

  //xcerr<<"ODKOMENTUJ REPAIR AZ TO BUDES DODELAVAT"<<endl;
  ////xcerr<<argv[1]<<endl;
  if(argv[1] == string("tests"))
  {
    tests();
    return 0;
  }
   stl_plane plane=stl_plane(0,0,1,-1);
  if(argc==6)
  {
    plane=stl_plane(atof( argv[2]),atof(argv[3]),atof(argv[4]),atof(argv[5]));
  }
 
  Mesh mesh;
  mesh.open(argv[1]);
  mesh.cut(plane); //na tomhle zlobila mrizka(1,-6,0.00,2) //(1,1,0,0)) na tomhle boxsphere
                                //0 0 1 -5 problem pro bo
  if(mesh.createBorderPolylines())
    {
      //xcout<<"CreateBorderPolylines proslo"<<endl;
      mesh.findHoles();
      mesh.triangulateCut();
      mesh.text();
      mesh.save();
    }
  //mesh.findHoles();
  mesh.close();
  mesh.remaining();

  

return 0;
}