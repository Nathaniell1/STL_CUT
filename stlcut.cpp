#include "stlcut.h"
using namespace p2t;
using namespace std;

double EPS=1e-25;
double minEPS=999999999;
char removedAxis='z';


  stl_plane::stl_plane(float x, float y, float z,float d) 
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = -1*d;
  }


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
      if ((tmp1.x==tmp3.x && tmp1.y==tmp3.y &&tmp1.z==tmp3.z && tmp2.x==tmp4.x && tmp2.y==tmp4.y &&  tmp2.z==tmp4.z) || ( tmp2.x==tmp3.x && tmp2.y==tmp3.y &&tmp2.z==tmp3.z && tmp1.x==tmp4.x && tmp1.y==tmp4.y && tmp1.z==tmp4.z ) )
      {
        border.erase(border.begin()+j,border.begin()+(j+2));
        j-=2;
       // break;
      }
    }
  }
}

/*
* Compares vertexes based of removed axis
*/

    bool setVertComp::operator() (const stl_vertex& lhs, const stl_vertex& rhs) const
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

    }



Mesh::~Mesh()
{
  for (unsigned  int i = 0; i < polylines.size(); ++i)
  {
    for (unsigned  int j = 0; j< polylines[i].size(); ++j)
    {
      //delete (p2t::Point*)polylines[i][j];
    }
  }
}

/*
*Calculates missing coordinate and then tryes to find it in the border to remove possible inaccuracy
*/
void Mesh::setMissingCoordinate(const p2t::Point* a,stl_vertex &b)
{
  if(removedAxis=='x')
  {
    b.y=a->y;
    b.z=a->x;
    b.x= (plane.y*b.y+plane.z*b.z+plane.d) / ((-1.0)* plane.x);
  }
  if(removedAxis=='y')
  {
    b.x=a->x;
    b.z=a->y;
    b.y= (plane.x*b.x+plane.z*b.z+plane.d) / ((-1.0)* plane.y);
  }
  if(removedAxis=='z')
  {
    b.x=a->x;
    b.y=a->y;
    b.z= (plane.x*b.x+plane.y*b.y+plane.d) / ((-1.0)* plane.z);
  }


  set<stl_vertex>::iterator itlow,itup;
  itlow=border2.lower_bound (b);              
  itup=border2.upper_bound (b);   
  itup--;
  if (itlow==itup||itlow==border2.begin()) // we found it
  {
    b=(*itlow);
  }
  else
  {
    itup++;
    itlow--;
    double dif1=ABS(max(  (max(ABS((*itlow).x-b.x),ABS((*itlow).y-b.y))) , ABS(((*itlow).z-b.z) )));
    double dif2=ABS(max(  (max(ABS((*itup).x-b.x) ,ABS((*itup).y-b.y ))) , ABS(((*itup).z-b.z )) ));
    if(dif1<dif2 && itlow!=border2.end() && dif1<0.00001) // find the closest
    {
      b=(*itlow); 
     
    }
    else 
    {
      if(dif2 <0.00001)     
        {
          b=(*itup);

        }
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

  if(removedAxis=='x')
  {
    vert.x=vert.z;
  }
  if(removedAxis=='y')
  {
    vert.y=vert.z;
  }
  if(vec.size()==0 ||  !(vec.back()->x==vert.x && vec.back()->y==vert.y))
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
    }
  }
  return in;
}


void Mesh::findHoles()
{
  //xcout<<"FIND HOLES"<<endl;
  vector<double> polygonArea;
  polygonArea.resize(polylines.size());
  for (unsigned  int i = 0; i < polylines.size(); ++i)
  {
    polygonArea[i]=calculatePolygonArea(polylines[i]);    
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
      tmpPolylines[i]=polylines[areaOrder[i]];
    }  
  }
  vector<int> indexes;
  for (unsigned int i = 0; i < tmpPolylines.size(); ++i)
  {
    
    if(polygonArea[areaOrder[i]]==0)
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
     //erasing wrong polygons with area=0
     tmpPolylines.erase(tmpPolylines.begin()+indexes[i]);
  }
  polylines=tmpPolylines; 


  vector<p2t::Point*> tmpPolygon=polylines.back();
  polylines.pop_back();
  pair<vector<p2t::Point*>,int >tmpPair=make_pair(tmpPolygon,-1);
  vector<pair<vector<p2t::Point*> , int> > tmpVecPair;
  tmpVecPair.push_back(tmpPair);
  polygonsWithHoles.push_back(tmpVecPair);//tmplist);
  int in;
  int out;
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
    for(unsigned  int k=0;k<polygonsWithHoles.size();k++)
    {
      in=0;
      out=0;
       
      //test 3 points and based on majority(just in case) decide if this polygons is inside biggerPolygon
      for(int j=0;j<3;j++)
      {
        vector<p2t::Point*>  biggerPolygon=polygonsWithHoles[k][pos].first;//.front();
        double a=tmpPolygon[j]->x; 
        double b=tmpPolygon[j]->y;
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
              holeIn=-1;  //swaping between hole and not a hole
            else
              holeIn=pos;//pos; 
          }
         
            
         placeFound=true;
         //as long as we didnt get to the end of the vector, continue - basicaly
         if(pos < polygonsWithHoles[k].size()-1) // we dont want to test it vs itself
         {
            pos++;
            k--; // we need to stay at the same index, 
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

        tmpPair=make_pair(tmpPolygon,-1);
        tmpVecPair.erase(tmpVecPair.begin(),tmpVecPair.end());
        tmpVecPair.push_back(tmpPair);
        polygonsWithHoles.push_back(tmpVecPair);//list2);
    }
  }
  
  
  }


void Mesh::triangulateCut()
{
  map<int, p2t::CDT*> polygons;//vector<p2t::CDT>> polygons;

  for (unsigned int i = 0; i < polygonsWithHoles.size(); ++i)
  {
      for (unsigned  int j = 0; j < polygonsWithHoles[i].size(); ++j)
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
    }
    float test_norm[3];
    test_norm[0] = plane.x;
    test_norm[1] = plane.y;
    test_norm[2] = plane.z;
    stl_normalize_vector(test_norm);
    facet.normal.x = test_norm[0];
    facet.normal.y = test_norm[1];
    facet.normal.z = test_norm[2];
    bot_facets.push_back(facet);
    //reverse normals
    facet.normal.x *= -1.0;
    facet.normal.y *= -1.0;
    facet.normal.z *= -1.0;
    vertex = facet.vertex[1];
    facet.vertex[1] = facet.vertex[2];
    facet.vertex[2] = vertex;
    top_facets.push_back(facet);
    }
   
}

bool Mesh::createBorderPolylines()
{

  if(border.size()==0 /*|| top_facets.size()==0 || bot_facets.size()==0*/) 
  {
    cerr<<"Nothing to cut"<<endl;
    return false;
  }

  {

    checkDuplicity(border);

    stl_vertex cont=border.back();
    border.pop_back();
    stl_vertex end=border.back();
    border.pop_back();
    polylines.resize(5);
    pushToPolylines(polylines[numOfPolylines],end);
    pushToPolylines(polylines[numOfPolylines],cont);
    EPS=1e-24;
    minEPS=999999999; //dont want to add limits just for this
    bool found=true;
    stl_vertex tmp1,tmp2;
    while(border.size()!=0)
    {
      
      if(!found) EPS=minEPS*1.005;
      if(EPS > 0.5) // ignoring edges, if it gets to this, there was probably a problem with mesh
      {
        EPS=1e-24;
        minEPS=999999999;
        //
        //       RETURN FALSE Could prevent poly2tri sefault 
        // 

        std::cerr<<"Unable to find a connected edge, mesh might be invalid"<<endl;
        //
        cont=border.back();
        border.pop_back();
        end=border.back();
        border.pop_back();
        if(!polylines[numOfPolylines].empty())
          numOfPolylines++;
        pushToPolylines(polylines[numOfPolylines],end);
        pushToPolylines(polylines[numOfPolylines],cont);
         // if we didnt found it with 0.5 tolerance.. lets just skip it
      }

 
      for (int i = border.size()-1;i >= 0; i-=2)
      {
       
        found=false;
        tmp1=border[i]; 
        tmp2=border[i-1];

        if(tmp1==cont||tmp2==cont) // ve found next vertex in polyline
        {
            found=true;
            //EPS=1e-24;

          if(tmp1==cont) 
          {            
            pushToPolylines(polylines[numOfPolylines],tmp2);
            cont=tmp2;
            border.erase(border.begin()+(i-1),border.begin()+i+1); // erase doesnt delete the last element IT DOES with set

          }  
          else
          {
            pushToPolylines(polylines[numOfPolylines],tmp1);
            cont=tmp1;
            border.erase(border.begin()+(i-1),border.begin()+i+1); // again, different erase with set

          }


          if((cont==end) )//&& border.size()>1)
          {

            if(!polylines[numOfPolylines].empty()) polylines[numOfPolylines].pop_back(); // delete last one, we dont want it twice
            
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
            }

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
  polylines.resize(numOfPolylines+1);

  return true;
}


void Mesh::cut(stl_plane plane)
{
  this->plane=plane;
  float norm[3];
  norm[0]=plane.x;
  norm[1]=plane.y;
  norm[2]=plane.z;
  stl_normalize_vector(norm);
  this->plane.x=norm[0];
  this->plane.y=norm[1];
  this->plane.z=norm[2];
  size_t aboves = 0;
  size_t belows = 0;
  size_t ons = 0;
  if(ABS(plane.x)>=ABS(plane.y) && ABS(plane.x)>=ABS(plane.z) ) removedAxis='x';
  if(ABS(plane.y)>=ABS(plane.x) && ABS(plane.y)>=ABS(plane.z) ) removedAxis='y';
  stl_position pos[3];


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

    }

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
    if(ons==3)
      {
       continue;     
      }
      if(ons==1)
      {
        if(belows==1 && aboves==1)
          {

          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==on)
            {
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect=intersection(a,b);
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
      { 
        for (int n = 0; n < 3; ++n)
        { 
          if(pos[n]==below) {bot_facets.push_back(facet);}
          if(pos[n]==above) {top_facets.push_back(facet);}
          if(pos[n]==on)    {border.push_back(facet.vertex[n]);} 
        }

        continue;
      }
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
      else 
      {
       if(aboves==1&&belows==2) // last possibility... the plane cuts the triangle and doesnt intersect with any vertex, this if is not neccesary
       {

          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==above) 
            {
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect1=intersection(facet.vertex[s],a);
              stl_vertex intersect2=intersection(facet.vertex[s],b);
              zCoord=intersect1.z;
               //hotfix to failed floating point operations
              if(intersect1.x==intersect2.x && intersect2.x==facet.vertex[s].x && intersect1.y==intersect2.y && intersect2.y==facet.vertex[s].y && intersect1.z==intersect2.z && intersect2.z==facet.vertex[s].z)
              {
                bot_facets.push_back(facet);
                break;
              }
              stl_facet tmp_facet=facet;
              tmp_facet.vertex[0]=facet.vertex[s];
              tmp_facet.vertex[1]=intersect1;
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
              border.push_back(intersect1);
              border.push_back(intersect2);
              break;
            }
           
          }
           
       }
       else // below==1, above == 2
       {
          for (int s = 0; s < 3; ++s)
          {
            if(pos[s]==below) 
            {
              stl_vertex a,b;
              a=facet.vertex[(s+1)%3];
              b=facet.vertex[(s+2)%3];
              stl_vertex intersect1=intersection(facet.vertex[s],a);
              stl_vertex intersect2=intersection(facet.vertex[s],b);
              zCoord=intersect1.z;
             if(intersect1.x==intersect2.x && intersect2.x==facet.vertex[s].x && intersect1.y==intersect2.y && intersect2.y==facet.vertex[s].y && intersect1.z==intersect2.z && intersect2.z==facet.vertex[s].z)
              {
                top_facets.push_back(facet);
                break;
              }

              stl_facet tmp_facet=facet;
              tmp_facet.vertex[0]=facet.vertex[s];
              tmp_facet.vertex[1]=intersect1;
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
              border.push_back(intersect1);
              border.push_back(intersect2);
              break;
            }
           
          }
           
       }
      }


  }

  for (unsigned int i = 0; i < border.size(); ++i)
  {
    border2.insert(border[i]);
  }
}

 // returns the position of the vertex related to the plane
  stl_position Mesh::vertex_position(stl_vertex vertex) 
  {
    
    double result = plane.x*vertex.x + plane.y*vertex.y + plane.z*vertex.z + plane.d;
    if (result > 0) return above;
    if (result < 0) return below;
    return on;
  }

  stl_vertex Mesh::intersection(stl_vertex a, stl_vertex b) 
  {
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
void Mesh::open( stl_file file)
{
  mesh_file = file;
}

void Mesh::open( char * name)
{
  stl_open(&mesh_file, name);
  stl_exit_on_error(&mesh_file);
}

stl_file* Mesh::export_stl2(deque<stl_facet> facets, const char* name) 
{
  stl_file* stl_out=new stl_file;
  stl_out->stats.type = inmemory;
  stl_out->stats.number_of_facets = facets.size();
  stl_out->stats.original_num_facets = (*stl_out).stats.number_of_facets;
  stl_out->v_indices = NULL;
  stl_out->v_shared = NULL;
  stl_out->neighbors_start = NULL;
  stl_clear_error(stl_out);
  stl_allocate(stl_out);
  
  int first = 1;
  for (deque<stl_facet>::const_iterator facet = facets.begin(); facet != facets.end(); facet++) {
    stl_out->facet_start[facet - facets.begin()] = *facet;
    stl_facet_stats(stl_out, *facet, first);
    first = 0;
  }
   stl_out->stats.degenerate_facets=0;
   stl_out->stats.edges_fixed=0;
   stl_out->stats.facets_removed=0;
   stl_out->stats.facets_added=0;
   stl_out->stats.facets_reversed=0;
   stl_out->stats.backwards_edges=0;
   stl_out->stats.normals_fixed=0;
  stl_repair(stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//reverse normal directions
  return stl_out;
}
void Mesh::save()
{
  export_stl(top_facets,"Cut_Mesh_1.stl");//"pokus1.stl");
  export_stl(bot_facets,"Cut_Mesh_2.stl");//"pokus2.stl");
  cout<<"Files saved to Cut_Mesh_1 and Cut_Mesh_2"<<endl;
}

std::array<stl_file*,2> Mesh::save2()
{
  return{export_stl2(top_facets,"Cut_Mesh_1.stl"),export_stl2(bot_facets,"Cut_Mesh_2.stl")};//"pokus1.stl");

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
  
  stl_repair(&stl_out,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,0,0);//reverse normal directions
  stl_write_ascii(&stl_out, name, "stlcut");
  stl_clear_error(&stl_out);
  stl_close(&stl_out);
}

bool tests()
{
  string names[]={"sphere.stl","trubka2.stl","trubka3.stl","trubka4.stl","trubka5.stl","sphere2.stl","hrana2.stl","boxsphere.stl"};
  float org_volume;

  stl_file mesh_file;
  int test_ok=0;
  int num_of_tests=8;
 
  for (int i = 0; i < num_of_tests;i++)
  {
     org_volume=0.0;
     stl_open(&mesh_file, (char*)names[i].c_str());
     stl_calculate_volume(&mesh_file);
     org_volume=mesh_file.stats.volume;
     stl_close(&mesh_file);
  

    Mesh mesh;
    mesh.open((char*)names[i].c_str());
    mesh.cut(stl_plane(1,1,0,0));
    if(mesh.createBorderPolylines())
      {
        mesh.findHoles();
        mesh.triangulateCut();
      }
    mesh.save();
    mesh.close();

    float volume=0.0;
    stl_open(&mesh_file, (char*)"Cut_Mesh_1.stl");
    stl_calculate_volume(&mesh_file);
    volume+=mesh_file.stats.volume;
    stl_close(&mesh_file);
    stl_open(&mesh_file, (char*)"Cut_Mesh_2.stl");
    stl_calculate_volume(&mesh_file);
    volume+=mesh_file.stats.volume;
    stl_close(&mesh_file);
    if(abs(org_volume-volume) < abs(org_volume/1000.0)) {cerr<<"Test OBJEMU uspesny"<<endl;test_ok++;};
  }
  cerr<<"Testu probehlo: "<<num_of_tests<<" z toho uspesne: "<<test_ok<<endl;
 

  
  
}

std::array<stl_file*,2> stlCut(stl_file* stlMesh,double a, double b, double c, double d,bool & succes)
{

	
	stl_plane plane=stl_plane(a,b,c,d);

	Mesh mesh;
  mesh.open(*stlMesh);
  mesh.cut(plane); 
  std::array<stl_file*,2> cutMesh;
                                
  if(mesh.createBorderPolylines())
    {
      mesh.findHoles();
      mesh.triangulateCut();
      cutMesh=mesh.save2();
      succes=true;
    }
    else succes=false;
  
  return{cutMesh[0],cutMesh[1]};

}
  
