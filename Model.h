//
//  Model.h
//  mcaq
//
//  Created by 赵子扬 on 14-7-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#ifndef __mcaq__Model__
#define __mcaq__Model__

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "Structures.h"

#if defined(GLUI_FREEGLUT)
#include <GL/freeglut.h>

#elif defined(GLUI_OPENGLUT)
#include <GL/openglut.h>

#else

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <GLUI/glui.h>// Header File For The GLUI Library
#include <OpenGL/gl.h> // Header File For The OpenGL Library
#include <OpenGL/glu.h> // Header File For The GLu Library
#else
#include <GL/glut.h>
#include <GL/glui.h>// Header File For The GLUI Library
#include <GL/gl.h> // Header File For The OpenGL Library
#include <GL/glu.h> // Header File For The GLu Library
#endif

#endif

using namespace std;


class Model
{
public:
    Model();//build model
    Model(int vertexNum, int faceNum);//build model
    ~Model();
    
    Vertex* AddVertex(int, GLfloat, GLfloat, GLfloat);//add a vertex to the model
    void AddFace(int, int, int);//add a face to the model
    
    void SaveModel();//save the model to a smf file
    
    void DrawModel(bool);//draw the model
    
    void DecimateModel(int);//decimate the model
    
private:
    list<Vertex*> vertexList;//list of vertexs in the model
    list<Face*> faceList;//list of edges in the model
    int vertexNum;//number of vertex
    int faceNum;//number of face
    
    
    GLfloat xMiddle, yMiddle, zMiddle;//middle point's coordinate
    GLfloat range;//model maximum range
    
    void ComputeVertexNormalAndModelPosition();
    void FindEdgeToDecimate();
    
    void DeleteFace(Face * face);//delete a face
    void DeleteVertex(Vertex * vertex);//delete a vertex
    
    int GetEdgeNum();
    /*User the max/min x/y/z coordinate to compute the model's middle point and range*/
    void ComputePosition(GLfloat, GLfloat, GLfloat, GLfloat, GLfloat, GLfloat);
};

#endif /* defined(__mcaq__Model__) */
