//
//  Model.cpp
//  mcaq
//
//  Created by 赵子扬 on 14-7-1.
//  Copyright (c) 2014年 ziyang zhao. All rights reserved.
//

#include "Model.h"
/**************float round(float, float)**********************/
float round(float orig, float bit)
{
    for (int i = 0; i < bit; i++) {
        orig *= 10;
    }
    
    int temp1 = orig;
    float temp2 = temp1;
    for (int i = 0; i < bit; i++) {
        temp2 /= 10;
    }
    return temp2;
}

/**************ComputeCoordinate(float*, float *)**************************/
/***Compute the optimal Coordinate of the new vertex**************************/
bool ComputeCoordinate(float* Q, float* coordinate)
{
     float d, x, y, z;
     d = Q[0]*Q[4]*Q[7] + Q[1]*Q[5]*Q[2] + Q[2] * Q[1] * Q[5] - Q[2]*Q[4]*Q[2] - Q[0]*Q[5]*Q[5] - Q[1]*Q[1]*Q[7];
     if (d == 0) {
         return false;
     }
    x = 0-Q[3] * Q[4] * Q[7] - Q[1]*Q[5]*Q[8] - Q[2] * Q[6] * Q[5] + Q[2]*Q[4]*Q[8] + Q[3]*Q[5]*Q[5] + Q[1]*Q[6]*Q[7];
    y = 0- Q[0]*Q[6]*Q[7] - Q[3]*Q[5]*Q[2] - Q[2]*Q[1]*Q[8] + Q[2]*Q[6]*Q[2] + Q[3]*Q[1]*Q[7] + Q[0]*Q[5]*Q[8];
    z = 0-Q[0]*Q[4]*Q[8] - Q[1]*Q[6]*Q[2] -Q[3]*Q[5]*Q[1] +Q[3]*Q[4]*Q[2] + Q[1]*Q[1]*Q[8] + Q[0]*Q[6]*Q[5];
    
    coordinate[0] = x/d;
    coordinate[1] = y/d;
    coordinate[2] = z/d;
    coordinate[3] = 1;
    return true;
}

/***************GLfloat determinant(GLfloat[][], int)***********************************************/
GLfloat determinant(GLfloat b[][4],int m)
{
    int i,j;
    GLfloat sum = 0, c[4][4];
    if(m==2)
    {
        sum = b[0][0]*b[1][1] - b[0][1]*b[1][0];
        return sum;
    }
    for(int p=0;p<m;p++)
    {
        int h = 0,k = 0;
        for(i=1;i<m;i++)
        {
            for( j=0;j<m;j++)
            {
                if(j==p)
                    continue;
                c[h][k] = b[i][j];
                k++;
                if(k == m-1)
                {
                    h++;
                    k = 0;
                }
            }
        }
        sum = sum + b[0][p]*pow(-1,p)*determinant(c,m-1);
    }
    return sum;
}

/********************ComputeError()***********************/
float ComputeError(float coordinate[4], float Q[][4])
{
    float error = 0, errorTemp[4] = {0};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            errorTemp[i] += coordinate[j] * Q[j][i];
        }
        error += errorTemp[i] * coordinate[i];
    }
    return error;
}



/*********GLfloat ComputeVertexAndError(Vertex*, Vertex*, GLfloat*)********/
GLfloat ComputeVertexAndError(Vertex *v1, Vertex *v2, GLfloat * newVertex)
{
    GLfloat error = 0;
    GLfloat Q[10];
    GLfloat QTemp[4][4], matric[4][4];
    
    //compute the new quadric
    for (int i = 0; i < 10; i++) {
            Q[i] = v1->quadric[i] + v2->quadric[i];
    }
    QTemp[0][0] = Q[0];
    QTemp[0][1] = Q[1];
    QTemp[0][2] = Q[2];
    QTemp[0][3] = Q[3];
    QTemp[1][0] = Q[1];
    QTemp[1][1] = Q[4];
    QTemp[1][2] = Q[5];
    QTemp[1][3] = Q[6];
    QTemp[2][0] = Q[2];
    QTemp[2][1] = Q[5];
    QTemp[2][2] = Q[7];
    QTemp[2][3] = Q[8];
    QTemp[3][0] = Q[3];
    QTemp[3][1] = Q[6];
    QTemp[3][2] = Q[8];
    QTemp[3][3] = Q[9];
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            matric[i][j] = QTemp[i][j];
        }
    }
    matric[3][0] = matric[3][1] = matric[3][2] = 0;
    matric[3][3] = 1;
    
    GLfloat det = round(determinant(matric, 4),8);
    if ((0 != det) && ComputeCoordinate(Q, newVertex)) {
        error = ComputeError(newVertex, QTemp);
        
    }else{//if the inverse of the matrix does not exist
        GLfloat errors[3]= {0};
        GLfloat **coordinates = new GLfloat *[3];
        int mark = 0;
        
        //get the coordinate of v1, v2 and the midpoint of v1and v2
        coordinates[0] = v1->coordinate;
        coordinates[1] = v2->coordinate;
        coordinates[2] = new GLfloat[4];
        for (int i = 0; i < 4; i++) {
            coordinates[2][i] = (coordinates[0][i] + coordinates[1][i]) / 2;
        }
        
        //compute the error of this three new vertex
        for (int i = 0; i < 3; i++) {
            errors[i] = ComputeError(coordinates[i], QTemp);
        }
        
        //choose the one with minimum error
        if (errors[1] < errors[0]) {
            mark = 1;
            if (errors[2] < errors[1]) {
                mark = 2;
            }
        }else{
            if (errors[2] < errors[0]) {
                mark = 2;
            }
        }
        for (int i = 0; i < 4; i++) {
            newVertex[i] = coordinates[mark][i];
        }
        error = errors[mark];
        delete [] coordinates[2];
        delete [] coordinates;
    }
    return error;
}

/********************ComputeNormal()***********************/
/*Compute the normal vecter and panel function of a triangle face*/
void ComputeFacePanel(Face *face)
{
/*    GLfloat temp1[3];
    GLfloat temp2[3];
    GLfloat *coordinate;
    GLfloat *panel = face->panel;
    
    
    for (int i = 0; i < 3; i++) {
        temp1[i] = face->vertex[0]->coordinate[i] - face->vertex[1]->coordinate[i];
        temp2[i] = face->vertex[1]->coordinate[i] - face->vertex[2]->coordinate[i];
    }
    
    panel[0] = temp1[1] * temp2 [2] - temp1[2] * temp2[1];
    panel[1] = temp1[0] * temp2 [2] - temp1[2] * temp2[0];
    panel[2] = temp1[0] * temp2 [1] - temp1[1] * temp2[0];
    
    /*convert this normal vector to a unit vector*/
/*    GLfloat length = sqrt(panel[0] * panel[0] + panel[1] * panel[1] + panel[2] * panel[2]);
    if (length == 0) {
        length = 1.0f;
    }
    for (int i = 0; i < 3; i++) {
        panel[i] = panel[i] / length;
    }
    
    coordinate = face->vertex[0]->coordinate;
    panel[3] = 0-(panel[0] * coordinate[0] +
                              panel[1] * coordinate[1] +
                              panel[2] * coordinate[2]) ;*/
    GLfloat *plane = face->panel;
    GLfloat *A = face->vertex[0]->coordinate;
    GLfloat *B = face->vertex[1]->coordinate;
    GLfloat *C = face->vertex[2]->coordinate;
    plane[0] = (B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]);
    plane[1] = (B[2] - A[2]) * (C[0] - A[0]) - (C[2] - A[2]) * (B[0] - A[0]);
    plane[2] = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    // normalize it with a^2 + b^2 + c^2 = 1
    double normalizer = sqrt( plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2] );
    plane[0] /= normalizer;
    plane[1] /= normalizer;
    plane[2] /= normalizer;
    plane[3] = - (plane[0] * A[0] + plane[1] * A[1] + plane[2] * A[2]);

    return;
}

/********************ComputeVertexNormal()***********************/
/*Compute the normal vecter of all the vertexes*/
void Model::ComputeVertexNormalAndModelPosition()
 {
     GLfloat length;
     GLfloat normalTemp = 0,
                    *normalTempFace,
                    *normalTempVertex;
     
     GLfloat xMin = 10000, xMax = -10000, yMin = 10000, yMax = -10000, zMin = 10000, zMax = -10000, x, y, z;
     
     for (list<Vertex*>::iterator v_iter = vertexList.begin(); v_iter != vertexList.end(); v_iter++)
     {
         length = 0;
         normalTempVertex = (*v_iter)->normal;
         for (int i = 0; i < 3; i++) {
             normalTempVertex[i] = 0;
         }
         for (list<Face *>::iterator f_iter = (*v_iter)->faceList.begin(); f_iter != (*v_iter)->faceList.end(); f_iter++) {
             normalTempFace = (*f_iter)->panel;
             for (int i = 0; i < 3; i++)
             {
                 normalTemp = normalTempFace[i];
                 normalTempVertex[i] += normalTemp;
                 length += normalTemp * normalTemp;
             }
         }
         
         length == 0?length = 1.0f : length = sqrt( length);
         for (int j = 0;  j < 3;  j++)
         {
             normalTempVertex[j]  /= length;
         }
         
         x = (*v_iter)->coordinate[0];
         y = (*v_iter)->coordinate[1];
         z = (*v_iter)->coordinate[2];
         
         if (x >xMax) {
             xMax = x;
         }else if (x < xMin) {
             xMin = x;
         }
         if (y >yMax) {
             yMax = y;
         }else if (y < yMin) {
             yMin = y;
         }
         if (z >zMax) {
             zMax = z;
         }else if (z < zMin) {
             zMin = z;
         }
     }
     
     ComputePosition(xMin, xMax, yMin, yMax, zMin, zMax);
 }



/****************Model()************************/
Model::Model()
{
    vertexNum = 0;
    faceNum = 0;
}

/****************Model()************************/
/***use the number of vertex and face to build model************/
Model::Model(int vn, int fn)
{
    vertexNum = 0;
    faceNum = 0;
}

/****************~Model()************************/
Model::~Model()
{
    vertexList.clear();
    faceList.clear();
}

/****************void AddVertex(int, GLfloat, GLfloat, GLfloat)************************/
Vertex* Model::AddVertex(int order, GLfloat x, GLfloat y, GLfloat z)
{
    Vertex *vertex = new Vertex;//create a new vertex
    GLfloat *coordinate = vertex->coordinate;
    
    vertex->order = order;//set the order
    //set the coordinate to x, y, z, 1
    coordinate[0] = x;
    coordinate[1] = y;
    coordinate[2] = z;
    coordinate[3] = 1;
    
    //initial the quadric
    for (int i = 0; i < 10; i++) {
        vertex->quadric[i] = 0;
    }
    
    //add this vertex to the list of vertex
    vertexList.push_back(vertex);
    
    //number of vertex increased
    vertexNum++;
    return vertex;
}

/****************void DeleteVertex(Face *face)*********************************/
void Model::DeleteVertex(Vertex *vertex)
{
    //clear the list in the vertex
    vertex->faceList.clear();
    //remove this vertex from the vertex list
    vertexList.remove(vertex);
    //decrease the number of vertex
    vertexNum--;
    delete vertex;//delete the vertex
}

/****************void AddFace(int, int, int)************************/
void Model::AddFace(int v1, int v2, int v3)
{
    Face * face = new Face;//create a new face
    
    Vertex*vertexTemp = NULL;
    
    int vertexOrder = 0;
    int found = 0;
    
    GLfloat quadricTemp[10], *panelTemp;
    
    //find vertex associated to this face
    for(list<Vertex *>::iterator list_iter = vertexList.begin(); list_iter != vertexList.end(); list_iter++)
    {
        vertexOrder = (*list_iter)->order;
        if (vertexOrder == v1) {
           face->vertex[0] = *list_iter;
            found ++;
        }else if ( vertexOrder == v2){
           face->vertex[1] = *list_iter;
            found ++;
        }else if(vertexOrder  == v3){
            face->vertex[2] = *list_iter;
            found ++;
        }
        if (found == 3) {
            break;
        }
    }
    
    //compute the panel and normal of the face
    ComputeFacePanel(face);
    panelTemp = face->panel;
    
    quadricTemp[0] = panelTemp[0] * panelTemp[0];
    quadricTemp[1] = panelTemp[0] * panelTemp[1];
    quadricTemp[2] = panelTemp[0] * panelTemp[2];
    quadricTemp[3] = panelTemp[0] * panelTemp[3];
    quadricTemp[4] = panelTemp[1] * panelTemp[1];
    quadricTemp[5] = panelTemp[1] * panelTemp[2];
    quadricTemp[6] = panelTemp[1] * panelTemp[3];
    quadricTemp[7] = panelTemp[2] * panelTemp[2];
    quadricTemp[8] = panelTemp[2] * panelTemp[3];
    quadricTemp[9] = panelTemp[3] * panelTemp[3];
    //compute the quadric of this face
    for (int i = 0; i< 3; i++)
    {
        vertexTemp = face->vertex[i];
        vertexTemp->faceList.push_back(face);
        for (int i = 0; i < 10; i++) {
            vertexTemp->quadric[i] += quadricTemp[i];
        }
    }
    
    //put the face into the list of face
    faceList.push_back(face);
    
    //increase the number of face
    faceNum++;
}

/**************void DeleteFace(Face*)**********************/
void Model::DeleteFace(Face *face)
{
    //remove this face from all the vertex's list
    for (int j = 0; j < 3; j++) {
        face->vertex[j]->faceList.remove(face);
    }
    
    //remove this face from the list of face
    faceList.remove(face);
    
    delete face;//delete the face
    
    //decrease the number of face
    faceNum --;
}

/****************void AddPosition()************************/
void Model::ComputePosition(GLfloat xMin, GLfloat xMax, GLfloat yMin, GLfloat yMax, GLfloat zMin, GLfloat zMax)
{
    /*compute the middle point of the model*/
    xMiddle = (xMin + xMax) / 2;
    yMiddle = (yMin + yMax) / 2;
    zMiddle = (zMin + zMax) / 2;
    
    /*compute the range of the model: find the maximum among  x-range, y-range and z-range*/
    range = xMax - xMin;
    if ((xMax - xMin) < (yMax - yMin)){
        range = yMax - yMin;
    }
    if (range < (zMax - zMin)) {
        range = zMax - zMin;
    }
}

/****************void DrawModel()************************/
void Model::DrawModel(bool isFlat)
{
    //compute the normal of vertex, and compute the position and scale of this model
    ComputeVertexNormalAndModelPosition();
    
    GLfloat scale = 1 / range;
    Vertex *v1, *v2, *v3;
    
    /*put the middle point of the model at the middle of the screen*/
    glTranslatef(-xMiddle, -yMiddle, -zMiddle);
    
    /*adjust the model to the best scale*/
    glScalef( scale, scale, scale );
    glEnable(GL_NORMALIZE);
    
    glBegin(GL_TRIANGLES);//draw triangles
    
    /*for each face, get its vertexs' index, use the indexes to get the vertex's coordinate, connect the vertex*/
    for (list<Face *>::iterator f_iter = faceList.begin(); f_iter != faceList.end(); f_iter ++)
    {
        //get the three vertex of a face
        v1 = (*f_iter)->vertex[0];
        v2 = (*f_iter)->vertex[1];
        v3 = (*f_iter)->vertex[2];
        
        if (isFlat) {
            glNormal3fv((*f_iter)->panel);
            glVertex3fv(v1->coordinate);
            glVertex3fv(v2->coordinate);
            glVertex3fv(v3->coordinate);
        }else{
            glNormal3fv(v1->normal);//get the normal of each vertex
            glVertex3fv(v1->coordinate);
            glNormal3fv(v2->normal);
            glVertex3fv(v2->coordinate);
            glNormal3fv(v3->normal);
            glVertex3fv(v3->coordinate);
        }
    }
    
    glScalef( range, range, range );
    glEnd();
}

/****************void SaveModel()************************/
void Model::SaveModel()
{
    string save_path;//the path to save the file
    string format;//format of the file
    ofstream file;//the file to save
    
    /*ask to user to input a usable path*/
    while (true)
    {
        cout<<"Please enter the path to save:"<<endl;
        getline(cin, save_path);
        
        /*if the path does not ends with '.smf', ask the user to input again*/
        format = save_path.substr(save_path.length()-4, 4);
        if (0 != strcmp(format.c_str(), ".smf"))
        {
            cout<<"This path does not end with '.smf', please enter again"<<endl;
            continue;
        }//end if
        
        file.open(save_path.c_str());
        if (file.is_open())
        {
            break;
        }//end outer else
        else//if the file does not exist, ask the user to input another path
        {
            cout<<"Failed to open the path!"<<endl;
            continue;
        }//end outer else
    }//end while
    

    /*write the first line to the file*/
    file<<"# "<<vertexNum<<" "<<faceNum<<"\n";
    
    int vOrder= 1;
    GLfloat *vertexCoor = NULL;
    Vertex ** vertexTemp = NULL;
    /*write every vertex to the file*/
    for(list<Vertex *>::iterator list_iter = vertexList.begin(); list_iter != vertexList.end(); list_iter++)
    {
        (*list_iter)->order = vOrder;
        vertexCoor = (*list_iter)->coordinate;
        file<<"v "<<vertexCoor[0]<<" "<<vertexCoor[1]<<" "<<vertexCoor[2]<<"\n";
        vOrder++;
    }
    
    /*write every face to the file*/
    for(list<Face *>::iterator list_iter = faceList.begin(); list_iter != faceList.end(); list_iter++)
    {
        vertexTemp = (*list_iter)->vertex;
        file<<"f "<<vertexTemp[0]->order
        <<" "<<vertexTemp[1]->order
        <<" "<<vertexTemp[2]->order<<"\n";
    }
    
    cout<<"successfully saved!"<<endl;
    file.close();
} 

/****************int GetEdgeNum()************************/
int Model::GetEdgeNum()
{
    //get the number of edge
    int edgeNum = faceNum * 3/ 2;
    return edgeNum;
}

/****************void DecimateModel()************************/
void Model::DecimateModel(int number)
{
    //whether  the number is usable or not
    if (number < 0 || number > GetEdgeNum()) {
        cout<<"Wrong number! Too big or too small, please try another one!"<<endl;
        return;
    }
    /* initialize random seed: */
    srand (time(0));

    //we need to decimate "number" of edges
    for (int i = 0; i < number; i++) {
        FindEdgeToDecimate();
    }
    
    cout<<"Decimate finished."<<endl;
    return;
}



/****************void DecimateOneEdge()************************/
void Model::FindEdgeToDecimate()
{
    int number = 8;
    Vertex* vArray[number][2],//store 8 vertexs
                 * vertex1, *vertex2, *newVertex;//the two vertexs that has been choose to decimate
    
    int vOrder[number][2];//the order of the face, and the order of vertexs in the face
    
    GLfloat newVertexs[number][4];//store the coordinate of the 8 new vertexs
    
    int count = 0, found = 0, mark = 0;
    GLfloat minError = 10000, error = 0;
    
    //get the order of the vertexs by random
    for (int i = 0; i < number; i++) {
        vOrder[i][0] = rand()%faceNum;
        vOrder[i][1] = rand()%3;
    }
    
    //get the 8*2 vertexs
    for (list<Face *>::iterator iter_face = faceList.begin(); iter_face != faceList.end(); iter_face++, count++) {
        for (int i = 0; i < number; i++) {
            //if this is the vertex we are looking for
            if (count == vOrder[i][0]) {
                //get the two vertex(the edge)
                vArray[i][0] = (*iter_face)->vertex[vOrder[i][1]];
                vArray[i][1] = (*iter_face)->vertex[(vOrder[i][1] + 1)%3];
                //compute the coordinate of the new vertex and of the error
                error = ComputeVertexAndError(vArray[i][0], vArray[i][1], newVertexs[i]);
                //find the smallest error
                if (error < minError) {
                    minError = error;
                    mark = i;
                }
                found++;
            }
        }
        if (found == number) {
            break;
        }
    }
    
    //we find the two vertex(the edge to decimate) by the minimum error
    vertex1 = vArray[mark][0];
    vertex2 = vArray[mark][1];
    
    //create the new vertex
    newVertex = AddVertex(0, newVertexs[mark][0], newVertexs[mark][1], newVertexs[mark][2]);
    for (int i = 0; i < 10; i++) {
        newVertex->quadric[i] = vertex1->quadric[i] + vertex2->quadric[i];
    }
    
    //for all the face associated with vertex1
    for (list<Face *>::iterator iter_face = vertex1->faceList.begin(); iter_face !=  vertex1->faceList.end();  iter_face++) {
        for (int i = 0; i < 3; i++) {
            //if this vertex is vertex1, update this face
            if ((*iter_face)->vertex[i] == vertex1) {
                (*iter_face)->vertex[i] = newVertex;
                newVertex->faceList.push_back(*iter_face);
                ComputeFacePanel(*iter_face);
                break;
            }
        }
    }
    
    Face *faces[2];
    int ffound = 0;
    bool isNew;
    //for all the face associated with vertex2
    for (list<Face *>::iterator iter_face = vertex2->faceList.begin(); iter_face !=  vertex2->faceList.end();  iter_face++) {
        isNew = false;
        for (int i = 0; i < 3; i++) {
            //if this vertex is newVertex, mark this face
            if ((*iter_face)->vertex[i] == newVertex) {
                faces[ffound] = *iter_face;
                ffound++;
                isNew = true;
                //DeleteFace(*iter_face);
                break;
            }
        }
        if (isNew == true) {
            continue;
        }
        for (int i = 0; i < 3; i++) {
            //if this vertex is vertex2, update thsi face
            if ((*iter_face)->vertex[i] == vertex2) {
                (*iter_face)->vertex[i] = newVertex;
                newVertex->faceList.push_back(*iter_face);
                ComputeFacePanel(*iter_face);
                break;
            }
        }
    }
    
    //delete the faces and vertexs
    DeleteFace(faces[0]);
    DeleteFace(faces[1]);
    DeleteVertex(vertex1);
    DeleteVertex(vertex2);
}
