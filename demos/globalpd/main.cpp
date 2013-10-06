#include <tinyxml.h>
#include "CKL.h"
#include "MatVec.h"
#include <vector>
#include <string>

using namespace CKL;


static void loadSceneFile(const std::string& filename,
                          std::vector<std::vector<CKL_REAL> >& points_array,
                          std::vector<std::vector<int> >& triangles_array,
                          std::vector<std::pair<Transform, Transform> >& motions)
{
  TiXmlDocument doc(filename.c_str());
  if(doc.LoadFile())
  {
    TiXmlHandle hdoc(&doc);
    TiXmlElement* model, *motion;
    model = doc.FirstChildElement("MODEL");
    if(model)
    {
      TiXmlElement* object = model->FirstChildElement("OBJECT");
      int i = 0;
      while(object)
      {
        int object_id = -1;
        object->Attribute("id", &object_id);
        // std::cout << object_id << std::endl;

        std::vector<CKL_REAL> points;
        std::vector<int> triangles;

        TiXmlElement* grid = object->FirstChildElement("GRID");
        int n_vertices = 0;
        while(grid)
        {
          int grid_id;
          double grid_x, grid_y, grid_z;

          grid->Attribute("id", &grid_id);
          grid->Attribute("x", &grid_x);
          grid->Attribute("y", &grid_y);
          grid->Attribute("z", &grid_z);

          if(3 * (grid_id - 1) == (int)points.size())
          {
            points.push_back(grid_x);
            points.push_back(grid_y);
            points.push_back(grid_z);
          }
          else if(3 * (grid_id - 1) < (int)points.size())
          {
            points[3 * (grid_id - 1)] = grid_x;
            points[3 * (grid_id - 1) + 1] = grid_y;
            points[3 * (grid_id - 1) + 2] = grid_z;
          }
          else // if(3 * (grid_id - 1) > points.size())
          {
            points.resize(3 * grid_id);
            points[3 * (grid_id - 1)] = grid_x;
            points[3 * (grid_id - 1) + 1] = grid_y;
            points[3 * (grid_id - 1) + 2] = grid_z;
          }
                          
          n_vertices++;
          grid = grid->NextSiblingElement("GRID");
        }

        // std::cout << "#vertices " << n_vertices << std::endl;

        TiXmlElement* tri = object->FirstChildElement("TRIA");
        int n_tris = 0;
        while(tri)
        {
          int tri_id;
          int v1, v2, v3;

          tri->Attribute("id", &tri_id);
          tri->Attribute("g1", &v1);
          tri->Attribute("g2", &v2);
          tri->Attribute("g3", &v3);

          if(3 * (tri_id - 1) == (int)triangles.size())
          {
            triangles.push_back(v1-1);
            triangles.push_back(v2-1);
            triangles.push_back(v3-1);
          }
          else if(3 * (tri_id - 1) < (int)triangles.size())
          {
            triangles[3* (tri_id - 1)] = v1-1;
            triangles[3* (tri_id - 1) + 1] = v2-1;
            triangles[3* (tri_id - 1) + 2] = v3-1;            
          }
          else
          {
            triangles.resize(3 * tri_id);
            triangles[3* (tri_id - 1)] = v1-1;
            triangles[3* (tri_id - 1) + 1] = v2-1;
            triangles[3* (tri_id - 1) + 2] = v3-1;
          }
          
          n_tris++;
          tri = tri->NextSiblingElement("TRIA");
        }

        // std::cout << "#triangles " << n_tris << std::endl;

        if(object_id - 1 == (int)points_array.size())
        {
          points_array.push_back(points);
          triangles_array.push_back(triangles);
        }
        else if(object_id - 1 < (int)points_array.size())
        {
          points_array[object_id - 1] = points;
          triangles_array[object_id - 1] = triangles;
        }
        else
        {
          points_array.resize(object_id);
          triangles_array.resize(object_id);
          points_array.back() = points;
          triangles_array.back() = triangles;
        }

        object = object->NextSiblingElement("OBJECT");
        i++;
      }

      // std::cout << "#objects " << i << std::endl;
    }

    motion = doc.FirstChildElement("MOTION");
    if(motion)
    {
      TiXmlElement* frame = motion->FirstChildElement("FRAME");
      int n_frame = 0;
      while(frame)
      {
        int frame_id;
        double frame_time;
        frame->Attribute("id", &frame_id);
        frame->Attribute("time", &frame_time);


        CKL_REAL T1[3], T2[3];
        CKL_REAL R1[3][3], R2[3][3];
        const char* obj1_pos_string = frame->Attribute("obj1_pos");
        const char* obj2_pos_string = frame->Attribute("obj2_pos");
        const char* obj1_dc_string = frame->Attribute("obj1_dc");
        const char* obj2_dc_string = frame->Attribute("obj2_dc");

        std::stringstream s1_pos(obj1_pos_string);
        s1_pos >> T1[0] >> T1[1] >> T1[2];
        std::stringstream s2_pos(obj2_pos_string);
        s2_pos >> T2[0] >> T2[1] >> T2[2];
        std::stringstream s1_mat(obj1_dc_string);
        for(int j = 0; j < 3; ++j)
          for(int k = 0; k < 3; ++k)
            s1_mat >> R1[j][k];
        std::stringstream s2_mat(obj2_dc_string);
        for(int j = 0; j < 3; ++j)
          for(int k = 0; k < 3; ++k)
            s2_mat >> R2[j][k];


        std::pair<Transform, Transform> motion = std::make_pair(Transform(R1, T1), Transform(R2, T2));
        if(frame_id - 1 == (int)motions.size())
          motions.push_back(motion);
        else if(frame_id - 1 < (int)motions.size())
          motions[frame_id - 1] = motion;
        else
        {
          motions.resize(frame_id);
          motions.back() = motion;
        }
        
        
        frame = frame->NextSiblingElement("FRAME");
        n_frame++;
      }

      // std::cout << "#frames " << n_frame << std::endl;
    }
  }
  else
    std::cerr << "Failed to load file " << filename << std::endl;    
}

void setModel(CKL_Model* o,
              const std::vector<CKL_REAL>& points,
              const std::vector<int>& triangles)
{
  o->BeginModel();
  for(std::size_t i = 0; i < triangles.size() / 3; ++i)
  {
    CKL_REAL p1[3], p2[3], p3[3];
    int tri_id;

    tri_id = triangles[3 * i];
    p1[0] = points[3 * tri_id + 0];
    p1[1] = points[3 * tri_id + 1];
    p1[2] = points[3 * tri_id + 2];

    tri_id = triangles[3 * i + 1];
    p2[0] = points[3 * tri_id + 0];
    p2[1] = points[3 * tri_id + 1];
    p2[2] = points[3 * tri_id + 2];

    tri_id = triangles[3 * i + 2];
    p3[0] = points[3 * tri_id + 0];
    p3[1] = points[3 * tri_id + 1];
    p3[2] = points[3 * tri_id + 2];

    o->AddTri(p1, p2, p3, i);
  }

  o->EndModel();
}

void scenePenetrationTest(const std::string& filename)
{
  std::vector<std::vector<CKL_REAL> > points_array;
  std::vector<std::vector<int> > triangles_array;
  std::vector<std::pair<Transform, Transform> > motions;

  loadSceneFile(filename, points_array, triangles_array, motions);

  CKL_Model* m1 = new CKL_Model();
  setModel(m1, points_array[0], triangles_array[0]);
  CKL_Model* m2 = new CKL_Model();
  setModel(m2, points_array[1], triangles_array[1]);
  
  std::size_t KNN_K = 10;

  std::vector<Transform> contact_vectors = CKL_PenetrationDepthModelLearning(m1, m2, 100000, KNN_K);

  for(std::size_t frame_id = 0; frame_id < motions.size(); ++frame_id)
  {
    CKL_PenetrationDepthResult result;

    CKL_PenetrationDepth(&result,
                         motions[frame_id].first.R, motions[frame_id].first.T, m1,
                         motions[frame_id].second.R, motions[frame_id].second.T, m2,
                         contact_vectors);
    
    std::cout << "pd value " << frame_id << " " << result.pd_value << std::endl;
    std::cout << "resolved tf translation " << result.T << std::endl;
    std::cout << "resolved tf rotation " << result.R << std::endl;

    CKL_CollideResult cresult;
    CKL_Collide(&cresult,
                motions[frame_id].first.R, motions[frame_id].first.T, m1,
                result.R, result.T, m2,
                CKL_FIRST_CONTACT);

     if(cresult.NumPairs() > 0)
    {
      std::cout << "contact pos " << cresult.pairs[0].point << std::endl;
      std::cout << "contact normal " << cresult.pairs[0].normal << std::endl;
    }
  }
}



int main()
{

  std::string filename1("Model_1.xml");
  scenePenetrationTest(filename1);

  std::string filename2("Model_4.xml");
  scenePenetrationTest(filename2);

  std::string filename3("Model_5.xml");
  scenePenetrationTest(filename3);

  
}
