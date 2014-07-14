//
//  Structures.h
//  mcaq
//
//  Created by 赵子扬 on 14-7-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#ifndef mcaq_Structures_h
#define mcaq_Structures_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <list>

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

struct Vertex;
struct Edge;
struct Face;

struct Vertex
{
    GLfloat coordinate[4];
    GLfloat normal[3];
    GLfloat quadric[10];
    int order;
    list<Face *> faceList;
};

struct Face
{
    GLfloat panel[4];
    Vertex *vertex[3];
};
#endif
