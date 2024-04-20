// Main File for OpenGl objectes for project AANG

// Steven Priddy & Colin Craighead


#include "CSCIx229.h"
#include <stdio.h>
#include <string.h>
// Mode Globals
int mainMode = 0; // Modes that the project can be in: 0 map open, 1 looking at a scene, 2 in Appa cut scene
#define sceneNumTotal 3 // total number of scenes
int sceneNum = 0; // Scene index or num that we are in
// View Globals
int axes=0;       //  Display axes
int th=0;         //  Azimuth of view angle
int ph=30;         //  Elevation of view angle
double asp=1;     //  Aspect ratio
double dim[sceneNumTotal] = {80,110,300};   //  Size of Scenes - chinCity = 55 - 
int globalWidth = 600; // Global vars that track the size of the window
int globalHeight = 600;
// Light Globals
int show_ball = 0; // toggles weather the light source position is shown
int emission  =   0;  // Emission intensity (%)
int ambient[sceneNumTotal]   =  {40,25,26};  // Ambient intensity (%)
int diffuse[sceneNumTotal]   = {170,68,80};  // Diffuse intensity (%)
int specular  =   0;  // Specular intensity (%)
int shininess =   0;  // Shininess (power of two)
float shiny   =   1;    // Shininess (value)
int zh        =  90;  // Light azimuth
float light_x[sceneNumTotal] = {102,-48,279};  // x of the secens light
float light_y[sceneNumTotal] = {39,79,180};  // y of the secens light
float light_z[sceneNumTotal] = {-62,54,62};  // z of the secens light
// Texture Globals
unsigned int texture[40]; // Texture names
int sky[6];   //  Sky textures
// Map Globals
char* cityNames[sceneNumTotal] = {"* Chin City", "* Fire Nation Castle", "* Eastern Air Temple"}; // Names for the Map
int cityNamesLen[sceneNumTotal] = {92,155,170}; // Width of the name
int cityHovered[sceneNumTotal] = {0,0,0}; // 0 if city is not being hovered, 1 if it is
double cityPercentX[sceneNumTotal] = {.61,.31,.77}; // location % of the width
double cityPercentY[sceneNumTotal] = {.35,.43,.39}; // location % of the height
// Viewing - First person
int view_mode = 0; // overhead perspective and first-person views
// First person eyes - origin
double ex_1[sceneNumTotal] = {3,0,92};
double ey_1[sceneNumTotal] = {3,-15,-32};
double ez_1[sceneNumTotal] = {-3,150,-19};
// First person eyes - will be edited when the player is moving
double ex_1_0[sceneNumTotal] = {3,0,92};
double ey_1_0[sceneNumTotal] = {3,2,-32};
double ez_1_0[sceneNumTotal] = {-3,150,-19};
// First-person location
double firstLoc_x[sceneNumTotal] = {0.0,0.0,0.0};
double firstLoc_z[sceneNumTotal] = {0.0,0,0.0};
int th_1[sceneNumTotal] = {180,270,0};     //  first person x angle
double walk[sceneNumTotal] = {0.5,0.75,1}; // speed of walking
double walkingRadius[sceneNumTotal] = {18,5,20}; // radius of area the player can walk with respect to the original position of the player
double distance[sceneNumTotal] = {0.0,0.0,0.0}; // distance from the scene's origin to the player
double walkStep_x = 0.0; // x value of the player moving forward/backward
double walkStep_z = 0.0; // z value of the player moving forward/backward
// Viewing - Appa
int th_ap = 0;
int ph_ap = 0;

typedef struct {float x,y,z;} vtx;
// Terrain
// Array needed for an 11x11 vertex - 10 by 10 grid 
vtx TerrVerts10[11][11];
vtx TerrNorms10[10][20]; // 10x10 squares but 2 triangles per square
vtx TerrVertNorms10[11][11];
// Array needed for an 21x21 vertex - 20 by 20 grid 
vtx TerrVerts20[21][21];
vtx TerrNorms20[20][40]; // 10x10 squares but 2 triangles per square
vtx TerrVertNorms20[21][21];
// Array needed for an 11x11 vertex - 10 by 10 grid -- CC for clamped on all edges
vtx TerrVerts10CC[11][11];
vtx TerrNorms10CC[10][20]; 
vtx TerrVertNorms10CC[11][11];
// Array needed for an 11x11 vertex - 10 by 10 grid -- CC for clamped on all edges -- Up becuse all heights are rand >= 0
vtx TerrVerts10CCUp[11][11];
vtx TerrNorms10CCUp[10][20]; 
vtx TerrVertNorms10CCUp[11][11];
// Array needed for an 51x51 vertex - 50 by 50 triangle terrain
vtx TerrVerts50[51][51];
vtx TerrNorms50[50][100];
vtx TerrVertNorms50[51][51];
// Arrays needed for a 11x11 vertex - will be filled with dome heights
vtx TerrVertsDome[11][11];
vtx TerrNormsDome[10][20]; 
vtx TerrVertNormsDome[11][11];
// Arrays needed for a 11x11 vertex - will be filled with peak heights
vtx TerrVertsPeak[11][11];
vtx TerrNormsPeak[10][20]; 
vtx TerrVertNormsPeak[11][11];
// appa call list int
int appa_list_num=0;
// appas position when flying around the scene
int appa_scene_x[sceneNumTotal] = {0,0,16};
int appa_scene_y[sceneNumTotal] = {0,32,0};
int appa_scene_z[sceneNumTotal] = {0,0,-233};
double appa_scene_s[sceneNumTotal] = {.14,.28,.80}; // his scale
int appa_scene_h[sceneNumTotal] = {25,38,250}; // his fly height
int appa_scene_r[sceneNumTotal] = {75,100,340}; // his fly radi








// Uses the call list made from the obj file appa
void appa(double x, double y, double z, double scale, double rx, double ry, double rz){
      glPushMatrix();
      glTranslated(x,y,z);
      glRotated(ry,0,1,0);
      glRotated(rx,1,0,0);
      glRotated(rz,0,0,1);
      glScaled(scale,scale,scale);
      glCallList(appa_list_num);
      glPopMatrix();
}

/* 
 *  Draw sky box
 */
static void Sky(double D,unsigned int sidetex,unsigned int topbottex)
{
   glColor3f(1,1,1);
   glEnable(GL_TEXTURE_2D);

   //  Sides
   glBindTexture(GL_TEXTURE_2D,sidetex);
   glBegin(GL_QUADS);
   glTexCoord2f(0.00,0); glVertex3f(-D,-D,-D);
   glTexCoord2f(0.25,0); glVertex3f(+D,-D,-D);
   glTexCoord2f(0.25,1); glVertex3f(+D,+D,-D);
   glTexCoord2f(0.00,1); glVertex3f(-D,+D,-D);

   glTexCoord2f(0.25,0); glVertex3f(+D,-D,-D);
   glTexCoord2f(0.50,0); glVertex3f(+D,-D,+D);
   glTexCoord2f(0.50,1); glVertex3f(+D,+D,+D);
   glTexCoord2f(0.25,1); glVertex3f(+D,+D,-D);

   glTexCoord2f(0.50,0); glVertex3f(+D,-D,+D);
   glTexCoord2f(0.75,0); glVertex3f(-D,-D,+D);
   glTexCoord2f(0.75,1); glVertex3f(-D,+D,+D);
   glTexCoord2f(0.50,1); glVertex3f(+D,+D,+D);

   glTexCoord2f(0.75,0); glVertex3f(-D,-D,+D);
   glTexCoord2f(1.00,0); glVertex3f(-D,-D,-D);
   glTexCoord2f(1.00,1); glVertex3f(-D,+D,-D);
   glTexCoord2f(0.75,1); glVertex3f(-D,+D,+D);
   glEnd();

   //  Top and bottom
   glBindTexture(GL_TEXTURE_2D,topbottex);
   glBegin(GL_QUADS);
   glTexCoord2f(0.0,0); glVertex3f(+D,+D,-D);
   glTexCoord2f(0.5,0); glVertex3f(+D,+D,+D);
   glTexCoord2f(0.5,1); glVertex3f(-D,+D,+D);
   glTexCoord2f(0.0,1); glVertex3f(-D,+D,-D);

   glTexCoord2f(1.0,1); glVertex3f(-D,-D,+D);
   glTexCoord2f(0.5,1); glVertex3f(+D,-D,+D);
   glTexCoord2f(0.5,0); glVertex3f(+D,-D,-D);
   glTexCoord2f(1.0,0); glVertex3f(-D,-D,-D);
   glEnd();

   glDisable(GL_TEXTURE_2D);
}

/*
 *  Draw map
 *  The map is a 2d image that opens over the whole scene
 *  it has the city names printed to the window location
 */
void map()
{
   //  Save transform attributes (Matrix Mode and Enabled Modes)
   glPushAttrib(GL_TRANSFORM_BIT|GL_ENABLE_BIT);
   // setup map texture
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texture[4]);
   //  Save projection matrix and set unit transform
   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   glOrtho(-asp,+asp,-1,1,-1,1);
   //  Save model view matrix and set to indentity
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();
   //  Draw instrument panel with texture
   glColor3f(1,1,1);
   glBegin(GL_QUADS);
   glTexCoord2d(0,0);glVertex2f(-1*asp,-1);
   glTexCoord2d(1,0);glVertex2f(+1*asp,-1);
   glTexCoord2d(1,1);glVertex2f(+1*asp, 1);
   glTexCoord2d(0,1);glVertex2f(-1*asp, 1);
   glEnd();
   glColor3f(0,0,0);
   for(int i=0;i<sceneNumTotal;i++){
      if(cityHovered[i] == 1){
         glColor3f(.78,.74,.38);
      }
      glWindowPos2d(globalWidth*cityPercentX[i],globalHeight*cityPercentY[i]);
      Print(cityNames[i]);
      if(cityHovered[i] == 1){
         glColor3f(0,0,0);
      }
   }
   glDisable(GL_TEXTURE_2D);
   
   //  Reset model view matrix
   glPopMatrix();
   //  Reset projection matrix
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();
   //  Pop transform attributes (Matrix Mode and Enabled Modes)
   glPopAttrib();
}


/*
 *  Draw vertex in polar coordinates - from ex18 
 */
static void Vertex(int th,int ph)
{
   double x = Cos(th)*Cos(ph);
   double y = Sin(th)*Cos(ph);
   double z =         Sin(ph);
   glNormal3d(x,y,z);
   glTexCoord2d(th/360.0,ph/180.0+0.5);
   glVertex3d(x,y,z);
}

/*
 *  sphere function from ex18, DrawPlanet
 */
void sphere(double x, double y, double z, double dx, double dy, double dz, unsigned int texnum)
{
   int th,ph;

   /*
    *  Draw surface of the planet
    */
   //  Set texture
   glEnable(GL_TEXTURE_2D);
   glBindTexture(GL_TEXTURE_2D,texnum);
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glScaled(dx,dy,dz);
   //  Latitude bands
   glColor3f(1,1,1);
   for (ph=-90;ph<90;ph+=45)
   {
      glBegin(GL_QUAD_STRIP);
      for (th=0;th<=360;th+=45)
      {
         Vertex(th,ph);
         Vertex(th,ph+45);
      }
      glEnd();
   }
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}



/*
 * Curved Bridge
 * Bridge where the curve is based off an ellipse shape
 * a sets the upper curve of the bridge
 * b sets the lower curve of the bridge
 * th is rotation around y
 */
void bridge(double x, double y, double z, double dx, double dy, double dz, double th, double a, double b, unsigned int texnum){
   //  Set specular color to white
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   // Texture setup
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texnum);
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glRotated(th,0,1,0);
   glScaled(dx,dy,dz);
   // base color
   glColor3f(1,1,1);
   glBegin(GL_QUADS);
   // far left cap
   glNormal3f(-1,0,0);
   glTexCoord2f(0,1); glVertex3f(-1.2,.7,-.1);
   glTexCoord2f(1,1); glVertex3f(-1.2,.7,.1);
   glTexCoord2f(1,0); glVertex3f(-1.2,0,.1);
   glTexCoord2f(0,0); glVertex3f(-1.2,0,-.1);
   // far right cap
   glNormal3f(1,0,0);
   glTexCoord2f(1,1); glVertex3f(1.2,.7,-.1);
   glTexCoord2f(0,1); glVertex3f(1.2,.7,.1);
   glTexCoord2f(0,0); glVertex3f(1.2,0,.1);
   glTexCoord2f(1,0); glVertex3f(1.2,0,-.1);
   // Curved part of the Bridge
   int numSteps = 24; // Number of points along the ellipse
   double step = 2.4/numSteps;
   // Loop from left to right side
   for(int i = 0; i<numSteps; i++){
      double x1 = -1.2+(i*step); // x at the current step
      double x2 = -1.2+((i+1)*step); // one step ahead
      double a2 = a*a; // a^2 for the upper ellipse equation
      double b2 = b*b; // b^2 for the lower ellipse equation
      // Y values for the lower curve
      double Ly1;
      double Ly2;
      // Noramls x and y for the lower curve
      double NLx1;
      double NLy1;
      double NLx2;
      double NLy2;
      // Starts and ends flat but matches the ellipse inbetween
      if(i < 2 || i > 21){
         Ly1 = 0;
         Ly2 = 0;
         NLx1 = 0;
         NLy1 = 1; // Points up b/c its flipped later
         NLx2 = 0;
         NLy2 = 1; // Points up b/c its flipped later
      }else{
         Ly1 = sqrt(b2 - b2*x1*x1); // lower ellipse equation 
         Ly2 = sqrt(b2 - b2*x2*x2);
         // Normal to ellipse equation for lower ellipse
         NLx1 = b2*x1;
         NLy1 = b == 0? 1 : Ly1; // Turnaray is for if the ellipse is flat this would be 0 but needs to be 1 -- same for other y normals
         NLx2 = b2*x2;
         NLy2 = b == 0? 1 : Ly2;
      }
      // Upper curve is an ellipse curve the whole time
      // Y values for the upper curve
      double Uy1 = sqrt(a2 - a2*x1*x1/1.44) + .7;
      double Uy2 = sqrt(a2 - a2*x2*x2/1.44) + .7;
      // Normals to upper ellipse
      double NUx1 = a2*x1;
      double NUy1 = a == 0? 1 : 1.44*(Uy1-.7); 
      double NUx2 = a2*x2;
      double NUy2 = a == 0? 1 : 1.44*(Uy2-.7); 
      // Back side of the curve
      glNormal3f(0,0,-1);
      glTexCoord2f(0,1); glVertex3f(x1,Uy1,-.1);
      glTexCoord2f(1,1); glVertex3f(x2,Uy2,-.1);
      glTexCoord2f(1,0); glVertex3f(x2,Ly2,-.1);
      glTexCoord2f(0,0); glVertex3f(x1,Ly1,-.1);
      // Front side of the curve
      glNormal3f(0,0,1);
      glTexCoord2f(0,1); glVertex3f(x1,Uy1,.1);
      glTexCoord2f(1,1); glVertex3f(x2,Uy2,.1);
      glTexCoord2f(1,0); glVertex3f(x2,Ly2,.1);
      glTexCoord2f(0,0); glVertex3f(x1,Ly1,.1);
      // Top of the bridge
      glNormal3f(NUx2,NUy2,0); glTexCoord2f(1,0); glVertex3f(x2,Uy2,.1);
      glNormal3f(NUx1,NUy1,0); glTexCoord2f(0,0); glVertex3f(x1,Uy1,.1);
      glNormal3f(NUx1,NUy1,0); glTexCoord2f(0,1); glVertex3f(x1,Uy1,-.1);
      glNormal3f(NUx2,NUy2,0); glTexCoord2f(1,1); glVertex3f(x2,Uy2,-.1);
      // Bottem of the bridge
      // Normals are flipped to the inside of the ellipse
      glNormal3f(-NLx2,-NLy2,0); glTexCoord2f(1,0); glVertex3f(x2,Ly2,.1);
      glNormal3f(-NLx1,-NLy1,0); glTexCoord2f(0,0); glVertex3f(x1,Ly1,.1);
      glNormal3f(-NLx1,-NLy1,0); glTexCoord2f(0,1); glVertex3f(x1,Ly1,-.1);
      glNormal3f(-NLx2,-NLy2,0); glTexCoord2f(1,1); glVertex3f(x2,Ly2,-.1);
   }
   glEnd();
   //  Undo transformations
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

/*
 * CALCULATE NORMAL VECTOR FOR ONE SIDE OF THE FRUSTUM
 */
static vtx frustumNormals(vtx A,vtx B,vtx C)
{
   //  Planar vector 0
   float dx0 = A.x-B.x;
   float dy0 = A.y-B.y;
   float dz0 = A.z-B.z;
   //  Planar vector 1
   float dx1 = (C.x-A.x);
   float dy1 = (C.y-A.y);
   float dz1 = (C.z-A.z);
   //  Normal
   float Nx = dy0*dz1 - dy1*dz0;
   float Ny = dz0*dx1 - dz1*dx0;
   float Nz = dx0*dy1 - dx1*dy0;

   vtx normVec = { Nx, Ny, Nz};
   //  DRAW STAR SIDES
   //glNormal3f(Nx,Ny,Nz);
   return normVec; 
}


/*
 *  Draw a frustum with a base center at (x,y,z) with a height of 'h'. The base will have a side lengths
    of x_b and z_b and the top will have a side lengths of x_t and z_t
 */
static void frustum(double x,double y,double z, double h, double x_b, double z_b, double x_t, double z_t, double th, const char *textureStyle, unsigned int texNum1, unsigned int texNum2, unsigned int texNum3, unsigned int texNum4,double r, double g, double b)
{ 

   //  frustum coordinates
   vtx xyz_frustum[] =
   {
      { x_b/2, 0, z_b/2} , {-x_b/2, 0, z_b/2} , { x_t/2, h, z_t/2} , //front vertices
      {-x_b/2, 0,-z_b/2} , { x_b/2, 0,-z_b/2} , {-x_t/2, h,-z_t/2} , //back vertices  
      { x_b/2, 0,-z_b/2} , { x_b/2, 0, z_b/2} , { x_t/2, h,-z_t/2} , //right vertices    
      {-x_b/2, 0, z_b/2} , {-x_b/2, 0,-z_b/2},  {-x_t/2, h, z_t/2}   //left vertices
   };

   //  Set specular color to white
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glRotated(th,0,1,0);
   //glScaled(0.5,0.5,0.5);

   //  Enable textures
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   // glColor3f(150.0/255,64/255.0,0); for fire building
   glColor3f(r,g,b);
   //glColor3f(108.0/255,101.0/255,9.0/255); for chin city roofs
   glBindTexture(GL_TEXTURE_2D,texNum1);

   //  Front
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   // Calculating the frustum side's normal vector
   vtx nVec1 = frustumNormals(xyz_frustum[0],xyz_frustum[1],xyz_frustum[2]);
   glNormal3f( nVec1.x, nVec1.y, nVec1.z);
   glTexCoord2f(1,0); glVertex3f(-x_b/2,0, z_b/2);
   glTexCoord2f(0,0); glVertex3f( x_b/2,0, z_b/2);
   // Specify wrapped or cut off texture
   if(strcmp(textureStyle,"wrapped")==0)
   {
      glTexCoord2f(0,1); glVertex3f( x_t/2,h, z_t/2);
      glTexCoord2f(1,1); glVertex3f(-x_t/2,h, z_t/2);
   }
   else if (strcmp(textureStyle,"chopped")==0)
   {
      float a = (x_b/2.0-x_t/2.0)/x_b;
      float b = (x_b/2.0+x_t/2.0)/x_b;
      glTexCoord2f(a,1); glVertex3f( x_t/2,h, z_t/2);
      glTexCoord2f(b,1); glVertex3f(-x_t/2,h, z_t/2);
   }
   glEnd();

   //  Back
   glBindTexture(GL_TEXTURE_2D,texNum2);
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   // Calculating the frustum side's normal vector
   vtx nVec2 = frustumNormals(xyz_frustum[3],xyz_frustum[4],xyz_frustum[5]);
   glNormal3f( nVec2.x, nVec2.y, nVec2.z);
   glTexCoord2f(1,0); glVertex3f( x_b/2,0,-z_b/2);
   glTexCoord2f(0,0); glVertex3f(-x_b/2,0,-z_b/2);
   if(strcmp(textureStyle,"wrapped")==0)
   {
      glTexCoord2f(0,1); glVertex3f(-x_t/2,h,-z_t/2);
      glTexCoord2f(1,1); glVertex3f( x_t/2,h,-z_t/2);
   }
   else if (strcmp(textureStyle,"chopped")==0)
   {
      float a = (x_b/2.0-x_t/2.0)/x_b;
      float b = (x_b/2.0+x_t/2.0)/x_b;
      glTexCoord2f(a,1); glVertex3f(-x_t/2,h,-z_t/2);
      glTexCoord2f(b,1); glVertex3f( x_t/2,h,-z_t/2);
   }
   glEnd();

   //  Right
   glBindTexture(GL_TEXTURE_2D,texNum3);
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   // Calculating the frustum side's normal vector
   vtx nVec3 = frustumNormals(xyz_frustum[6],xyz_frustum[7],xyz_frustum[8]);
   glNormal3f( nVec3.x, nVec3.y, nVec3.z);
   glTexCoord2f(1,0); glVertex3f(x_b/2,0, z_b/2);
   glTexCoord2f(0,0); glVertex3f(x_b/2,0,-z_b/2);
   if(strcmp(textureStyle,"wrapped")==0)
   {
      glTexCoord2f(0,1); glVertex3f(x_t/2,h,-z_t/2);
      glTexCoord2f(1,1); glVertex3f(x_t/2,h, z_t/2);
   }
   else if (strcmp(textureStyle,"chopped")==0)
   {
      float a = (z_b/2.0-z_t/2.0)/z_b;
      float b = (z_b/2.0+z_t/2.0)/z_b;
      glTexCoord2f(a,1); glVertex3f(x_t/2,h,-z_t/2);
      glTexCoord2f(b,1); glVertex3f(x_t/2,h, z_t/2);
   }
   glEnd();

   //  Left
   glBindTexture(GL_TEXTURE_2D,texNum4);
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   // Calculating the frustum side's normal vector
   vtx nVec4 = frustumNormals(xyz_frustum[9],xyz_frustum[10],xyz_frustum[11]);
   glNormal3f( nVec4.x, nVec4.y, nVec4.z);
   glTexCoord2f(1,0); glVertex3f(-x_b/2,0,-z_b/2);
   glTexCoord2f(0,0); glVertex3f(-x_b/2,0, z_b/2);
   if(strcmp(textureStyle,"wrapped")==0)
   {
      glTexCoord2f(0,1); glVertex3f(-x_t/2,h, z_t/2);
      glTexCoord2f(1,1); glVertex3f(-x_t/2,h,-z_t/2);
   }
   else if(strcmp(textureStyle,"chopped")==0)
   {
      float a = (z_b/2.0-z_t/2.0)/z_b;
      float b = (z_b/2.0+z_t/2.0)/z_b;
      glTexCoord2f(a,1); glVertex3f(-x_t/2,h, z_t/2);
      glTexCoord2f(b,1); glVertex3f(-x_t/2,h,-z_t/2);
   }

   glEnd();
   //  Top
   glBindTexture(GL_TEXTURE_2D,texNum1);
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   glNormal3f( 0,+1, 0);
   glTexCoord2f(0,0); glVertex3f(-x_t/2,h, z_t/2);
   glTexCoord2f(1,0); glVertex3f( x_t/2,h, z_t/2);
   glTexCoord2f(1,1); glVertex3f( x_t/2,h,-z_t/2);
   glTexCoord2f(0,1); glVertex3f(-x_t/2,h,-z_t/2);
   glEnd();
   //  Bottom
   glBindTexture(GL_TEXTURE_2D,texNum1);
   //glColor3f(1,1,1);
   glBegin(GL_QUADS);
   glNormal3f( 0,-1, 0);
   glTexCoord2f(0,0); glVertex3f(-x_b/2,0,-z_b/2);
   glTexCoord2f(1,0); glVertex3f( x_b/2,0,-z_b/2);
   glTexCoord2f(1,1); glVertex3f( x_b/2,0, z_b/2);
   glTexCoord2f(0,1); glVertex3f(-x_b/2,0, z_b/2);
   glEnd();
   //  Undo transformations and textures
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}


/*
 * FROM CSCI 4229 EX CODE
 *  Draw a cube
 *     at (x,y,z)
 *     dimensions (dx,dy,dz)
 *     rotated th about the y axis
 */
static void cube(double x,double y,double z,
                 double dx,double dy,double dz,
                 unsigned int texNum1, double r, 
                 double g, double b, int rep)
{
   //  Set specular color to white
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(dx,dy,dz);
   //  Enable textures
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glColor3f(r,g,b);
   glBindTexture(GL_TEXTURE_2D,texNum1);
   //  Front
   glBegin(GL_QUADS);
   glNormal3f( 0, 0, 1);
   glTexCoord2f(0,0); glVertex3f(-1,-1, 1);
   glTexCoord2f(rep,0); glVertex3f(+1,-1, 1);
   glTexCoord2f(rep,rep); glVertex3f(+1,+1, 1);
   glTexCoord2f(0,rep); glVertex3f(-1,+1, 1);
   glEnd();
   //  Back
   glBegin(GL_QUADS);
   glNormal3f( 0, 0,-1);
   glTexCoord2f(0,0); glVertex3f(+1,-1,-1);
   glTexCoord2f(rep,0); glVertex3f(-1,-1,-1);
   glTexCoord2f(rep,rep); glVertex3f(-1,+1,-1);
   glTexCoord2f(0,rep); glVertex3f(+1,+1,-1);
   glEnd();
   //  Right
   glBegin(GL_QUADS);
   glNormal3f(+1, 0, 0);
   glTexCoord2f(0,0); glVertex3f(+1,-1,+1);
   glTexCoord2f(rep,0); glVertex3f(+1,-1,-1);
   glTexCoord2f(rep,rep); glVertex3f(+1,+1,-1);
   glTexCoord2f(0,rep); glVertex3f(+1,+1,+1);
   glEnd();
   //  Left
   glBegin(GL_QUADS);
   glNormal3f(-1, 0, 0);
   glTexCoord2f(0,0); glVertex3f(-1,-1,-1);
   glTexCoord2f(rep,0); glVertex3f(-1,-1,+1);
   glTexCoord2f(rep,rep); glVertex3f(-1,+1,+1);
   glTexCoord2f(0,rep); glVertex3f(-1,+1,-1);
   glEnd();
   //  Top
   glBegin(GL_QUADS);
   glNormal3f( 0,+1, 0);
   glTexCoord2f(0,0); glVertex3f(-1,+1,+1);
   glTexCoord2f(rep,0); glVertex3f(+1,+1,+1);
   glTexCoord2f(rep,rep); glVertex3f(+1,+1,-1);
   glTexCoord2f(0,rep); glVertex3f(-1,+1,-1);
   glEnd();
 
   glBegin(GL_QUADS);
   glNormal3f( 0,-1, 0);
   glTexCoord2f(0,0); glVertex3f(-1,-1,-1);
   glTexCoord2f(rep,0); glVertex3f(+1,-1,-1);
   glTexCoord2f(rep,rep); glVertex3f(+1,-1,+1);
   glTexCoord2f(0,rep); glVertex3f(-1,-1,+1);
   glEnd();
   //  Undo transformations and textures
   glPopMatrix();
   glDisable(GL_TEXTURE_2D);
}

/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius r
 */
static void ball(double x,double y,double z,double r)
{
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball
   glColor3f(1,1,1);
   glutSolidSphere(1.0,16,16);
   //  Undo transofrmations
   glPopMatrix();
}

/*
   Draws the cylander cap for the cylanders in the easter air temple
 */
static void cylinderCap(double x,double y,double z,double dx,double dy,double dz,double ry,unsigned int texnum_roof,unsigned int texnum_gold)
{
	
   //  Save transformation
   glPushMatrix();
   //  Offset and scale
   glTranslated(x,y,z);
   glRotatef(ry,0,1,0);
   glScaled(dx,dy,dz);

   // gold bottem frusts
   frustum(0,0,0,.4,2,2,2,2,0,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // bottom frust
   frustum(0,.4,0,.6,2.4,2.4,2.4,2.4,0,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // box above bottem
   // rotate frusts upside down to be the wings out
   glPushMatrix();
   glRotatef(180,1,0,0);
   frustum(-1,-1,1,1,2.1,.1,.1,.1,45,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // front left corrner
   frustum(-1,-1,-1,1,2.1,.1,.1,.1,-45,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // back right corrner
   frustum(1,-1,-1,1,2.1,.1,.1,.1,45,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // back right corrner
   frustum(1,-1,1,1,2.1,.1,.1,.1,-45,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // front right corrner
   glPopMatrix();
   // Roof frusts
   frustum(0,1,0,1.2,2.6,2.6,.8,.8,0,"wrapped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // roof trapaziod
   frustum(0,2.2,0,1.2,.8,.8,.8,.8,0,"wrapped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // box collumn from traps
   frustum(0,3.4,0,.5,1,1,.5,.5,0,"wrapped",texnum_gold,texnum_gold,texnum_gold,texnum_gold,1.0,1.0,1.0); // gold cap
   
   glPopMatrix();
}

/*
   DRAWS A CYLINDER AT (x,y,z) OF RADIUS r AND HEIGHT h
 */
static void cylinder(double x,double y,double z,double r,double h,double th,unsigned int texnum)
{
   float d=10.0;
   float theta = 0.0;

   //  Save transformation
   glPushMatrix();
   //  Offset and scale
   glTranslated(x,y,z);
   glRotatef(th,0,0,1);
   glScaled(r,h,r);

   // TEXTURES 
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texnum);

   glBegin(GL_QUAD_STRIP);
   //  Latitude bands
   glColor3f(1, 1,1); 
   for (theta=0;theta<=360;theta+=d)
   {
      glNormal3d(Sin(theta),0,Cos(theta));
      glTexCoord2f( theta/360.0, 1); glVertex3f(Sin(theta), 1, Cos(theta));
      glTexCoord2f( theta/360.0, 0);glVertex3f(Sin(theta), 0, Cos(theta));
   
   }
   glEnd();

   // BOTTOM OF CYLINDER
   glBegin(GL_POLYGON);
   for (theta=0;theta<=360;theta+=d)
   {
      glNormal3d(0,-1,0);
      glTexCoord2f(2/2*Cos(theta)+0.5,2/2*Sin(theta)+0.5);
      glVertex3f(Sin(theta), 0, Cos(theta));
   }
   glEnd();

   // TOP OF CYLINDER
   glBegin(GL_POLYGON);
   for (theta=0;theta<=360;theta+=d)
   {
      glNormal3d(0,1,0);
      glTexCoord2f(3/2*Cos(theta)+0.5,3/2*Sin(theta)+0.5);
      glVertex3f(Sin(theta), 1, Cos(theta));
   }
   glEnd();
   glDisable(GL_TEXTURE_2D);
   //  Undo transformations
   glPopMatrix();
}

// Draws a tree made of a cylander trunk and sphere leaves
void tree_1(double x, double y, double z, double dx, double dy, double dz, double alpha, unsigned int textNum_bark,unsigned int texNum_leave)
{
	//  Save transformation
   	glPushMatrix();
   	//  Offset and scale
   	glTranslated(x,y,z);
   	glRotatef(alpha,0,1,0);
   	glScaled(dx,dy,dz);

	cylinder(0,0,0,0.5,5,0,textNum_bark);
	sphere(0, 6, 0, 3, 3, 3, texNum_leave);
	sphere(-2, 5, 1, 2, 2, 2, texNum_leave);
	sphere(-2, 5, 1, 2, 2, 2, texNum_leave);
	sphere(3, 6, 0, 2, 2, 2, texNum_leave);
	sphere(1, 6, 3, 2, 2, 2, texNum_leave);
	sphere(2, 6, 2, 2.5, 2.5, 2.5, texNum_leave);
	sphere(-2, 6, -2, 2.5, 2.5, 2.5, texNum_leave);
	sphere(0, 6, -3, 2, 2, 2, texNum_leave);
	sphere(0, 4, 2, 2, 1.5, 1.5, texNum_leave);
	sphere(-1, 5, 2, 2, 1.5, 1.5, texNum_leave);

	glPopMatrix();
}


/* 
 * Fills in a square Heights array with semi random but Peak shaped floats
*/
void PeakTerrHeights(int sqr_size, float heights[sqr_size][sqr_size], unsigned seed){
   srand(seed);
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         heights[i][j] = (rand() % 1000)/500.0 - 1; // semi random floats -1 < x < 1
      }
   }
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         int x = i-sqr_size/2;
         int y = j-sqr_size/2;
         heights[i][j] += sqrt((-2*abs(x)+(sqr_size-1))*(-2*abs(y)+(sqr_size-1))); // absloute value functions fliped over make a peak
      }
   }
   // clamp left right col then clamp top and bot row
   for(int i=0;i<sqr_size;i++){
      heights[i][0] = 0;
      heights[i][sqr_size-1] = 0;
   }
   for(int j=0;j<sqr_size;j++){
      heights[0][j] = 0;
      heights[sqr_size-1][j] = 0;
   }
}


/* 
 * Fills in a square Heights array with semi random but dome shaped floats
*/
void DomeTerrHeights(int sqr_size, float heights[sqr_size][sqr_size], unsigned seed){
   srand(seed);
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         heights[i][j] = (rand() % 1000)/500.0 - 1; // semi random floats -1 < x < 1
      }
   }
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         int x = i-sqr_size/2;
         int y = j-sqr_size/2;
         heights[i][j] += -1*(.5/(sqr_size-1))*(x*x + y*y) + sqr_size/2.0;
      }
   }
   // clamp left right col then clamp top and bot row
   for(int i=0;i<sqr_size;i++){
      heights[i][0] = 0;
      heights[i][sqr_size-1] = 0;
   }
   for(int j=0;j<sqr_size;j++){
      heights[0][j] = 0;
      heights[sqr_size-1][j] = 0;
   }
}


/* 
 * Fills in a square Heights array with random floats from -1 to 1
 * clamp mode: 0 no clamp, 1 left right clamped, 2 all sides clamped - clamped means that the said edge is all at 0
*/
void RandTerrHeights(int clampMode, int sqr_size, float heights[sqr_size][sqr_size], unsigned seed){
   srand(seed);
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         heights[i][j] = (rand() % 1000)/500.0 - 1; // semi random floats -1 < x < 1
      }
   }
   // clap edges based on mode
   // clamp left and right col
   if(clampMode == 1){
      for(int i=0;i<sqr_size;i++){
         heights[i][0] = 0;
         heights[i][sqr_size-1] = 0;
      }
   }
   // clamp left right col then clamp top and bot row
   if(clampMode == 2){
      for(int i=0;i<sqr_size;i++){
         heights[i][0] = 0;
         heights[i][sqr_size-1] = 0;
      }
      for(int j=0;j<sqr_size;j++){
         heights[0][j] = 0;
         heights[sqr_size-1][j] = 0;
      }
   }
}

/* 
 * Fills in a square Heights array with random floats from 0 to 2
*/
void RandTerrHeightsUp(int clampMode, int sqr_size, float heights[sqr_size][sqr_size], unsigned seed){
   srand(seed);
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         heights[i][j] = (rand() % 1000)/500.0; // semi random floats 0 < x < 2
      }
   }
   // clap edges based on mode
   // clamp left and right col
   if(clampMode == 1){
      for(int i=0;i<sqr_size;i++){
         heights[i][0] = 0;
         heights[i][sqr_size-1] = 0;
      }
   }
   // clamp left right col then clamp top and bot row
   if(clampMode == 2){
      for(int i=0;i<sqr_size;i++){
         heights[i][0] = 0;
         heights[i][sqr_size-1] = 0;
      }
      for(int j=0;j<sqr_size;j++){
         heights[0][j] = 0;
         heights[sqr_size-1][j] = 0;
      }
   }
}


/* 
 * Creates 3 2d Vertex arrays all needed for the terrain
 * sqr_size is the number of vertices in the terrain - !!! square size needs to be odd
 * TerrVerts holds the x,y,z values of the vertices
 * TerrNorms holds the face normals to every triangle
 * TerrVertNorms holds the avg face norms for every vertex
 * seed realates to the radomization of the heights
*/
void createTerrMats(int heightsMode, int sqr_size, vtx TerrVerts[sqr_size][sqr_size], vtx TerrNorms[sqr_size-1][2*(sqr_size-1)], vtx TerrVertNorms[sqr_size][sqr_size], unsigned seed){
   // Set the vertex grid up based on random heights
   float heights[sqr_size][sqr_size];
   if(heightsMode == 0){ // No clamped edges -- clamped means the edges heights are held at 0
      RandTerrHeights(0,sqr_size,heights,seed);
   }else if(heightsMode == 1){ // LR edges clamped
      RandTerrHeights(1,sqr_size,heights,seed);
   }else if(heightsMode == 2){ // all four edges clamped
      RandTerrHeights(2,sqr_size,heights,seed);
   }else if(heightsMode == 3){ // Dome with all four edges clamped
      DomeTerrHeights(sqr_size,heights,seed);
   }else if(heightsMode == 4){ // Peak with all four edges clamped
      PeakTerrHeights(sqr_size,heights,seed);
   }else if(heightsMode == 5){ // all 4 sides clamped and terain never dips below 0
      RandTerrHeightsUp(2,sqr_size,heights,seed);
   }

   int topLeftCoord =  -1*(sqr_size/2); // top left corner x,z value 
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         TerrVerts[i][j].x = topLeftCoord + j;
         TerrVerts[i][j].y = heights[i][j];
         TerrVerts[i][j].z = topLeftCoord + i;
      }
   }
   // Create the FACE Normals
   for(int i=0;i<sqr_size-1;i++){
      for(int j=0;j<sqr_size-1;j++){
         // upper triangle normal
         vtx V; // vector going left to right along the top of the triangle
         V.x = TerrVerts[i][j+1].x - TerrVerts[i][j].x; 
         V.y = TerrVerts[i][j+1].y - TerrVerts[i][j].y;
         V.z = TerrVerts[i][j+1].z - TerrVerts[i][j].z;
         vtx U; // vector going down along the left of the triangle
         U.x = TerrVerts[i+1][j].x - TerrVerts[i][j].x; 
         U.y = TerrVerts[i+1][j].y - TerrVerts[i][j].y;
         U.z = TerrVerts[i+1][j].z - TerrVerts[i][j].z;
         vtx N; // N = UxV
         N.x = (U.y*V.z) - (U.z*V.y);
         N.y = (U.z*V.x) - (U.x*V.z);
         N.z = (U.x*V.y) - (U.y*V.x);
         // Make the Vector of length one
         float mag = sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
         N.x = N.x/mag; N.y = N.y/mag; N.z = N.z/mag;
         TerrNorms[i][2*j] = N;

         // Lower triangle normal
         // new V
         V.x = TerrVerts[i+1][j].x - TerrVerts[i+1][j+1].x; 
         V.y = TerrVerts[i+1][j].y - TerrVerts[i+1][j+1].y;
         V.z = TerrVerts[i+1][j].z - TerrVerts[i+1][j+1].z;
         // new U
         U.x = TerrVerts[i][j+1].x - TerrVerts[i+1][j+1].x; 
         U.y = TerrVerts[i][j+1].y - TerrVerts[i+1][j+1].y;
         U.z = TerrVerts[i][j+1].z - TerrVerts[i+1][j+1].z;
         // N = UxV
         N.x = (U.y*V.z) - (U.z*V.y);
         N.y = (U.z*V.x) - (U.x*V.z);
         N.z = (U.x*V.y) - (U.y*V.x);
         // Make the Vector of length one
         mag = sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
         N.x = N.x/mag; N.y = N.y/mag; N.z = N.z/mag;
         TerrNorms[i][2*j + 1] = N;
      }
   }
   // Create Vert Norms by Avg Face Norms
   int faceNormArrRows = sqr_size-1;
   int faceNormArrCols = (sqr_size-1)*2;
   for(int i=0;i<sqr_size;i++){
      for(int j=0;j<sqr_size;j++){
         if(i == 0 && j == 0){ // Top left Vertex
            TerrVertNorms[i][j] = TerrNorms[0][0];
         }else if(i == 0 && j == sqr_size-1){ // Top Right Vertex
            TerrVertNorms[i][j].x = (TerrNorms[0][faceNormArrCols-2].x + TerrNorms[0][faceNormArrCols-1].x)/2;
            TerrVertNorms[i][j].y = (TerrNorms[0][faceNormArrCols-2].y + TerrNorms[0][faceNormArrCols-1].y)/2;
            TerrVertNorms[i][j].z = (TerrNorms[0][faceNormArrCols-2].z + TerrNorms[0][faceNormArrCols-1].z)/2;
         }else if(i == sqr_size-1 && j == 0){ // Bottem left
            TerrVertNorms[i][j].x = (TerrNorms[faceNormArrRows-1][0].x + TerrNorms[faceNormArrRows-1][1].x)/2;
            TerrVertNorms[i][j].y = (TerrNorms[faceNormArrRows-1][0].y + TerrNorms[faceNormArrRows-1][1].y)/2;
            TerrVertNorms[i][j].z = (TerrNorms[faceNormArrRows-1][0].z + TerrNorms[faceNormArrRows-1][1].z)/2;
         }else if(i == sqr_size-1 && j == sqr_size-1){ // Bottem Right Vertex
            TerrVertNorms[i][j] = TerrNorms[faceNormArrRows-1][faceNormArrCols-1];
         }else if(i == 0){ // Top Edge Vertices
            TerrVertNorms[i][j].x = (TerrNorms[0][2*j-2].x + TerrNorms[0][2*j-1].x + TerrNorms[0][2*j].x)/3;
            TerrVertNorms[i][j].y = (TerrNorms[0][2*j-2].y + TerrNorms[0][2*j-1].y + TerrNorms[0][2*j].y)/3;
            TerrVertNorms[i][j].z = (TerrNorms[0][2*j-2].z + TerrNorms[0][2*j-1].z + TerrNorms[0][2*j].z)/3;
         }else if(i == sqr_size-1){ // Bottem Edge Vertices
            TerrVertNorms[i][j].x = (TerrNorms[faceNormArrRows-1][2*j-1].x + TerrNorms[faceNormArrRows-1][2*j].x + TerrNorms[faceNormArrRows-1][2*j+1].x)/3;
            TerrVertNorms[i][j].y = (TerrNorms[faceNormArrRows-1][2*j-1].y + TerrNorms[faceNormArrRows-1][2*j].y + TerrNorms[faceNormArrRows-1][2*j+1].y)/3;
            TerrVertNorms[i][j].z = (TerrNorms[faceNormArrRows-1][2*j-1].z + TerrNorms[faceNormArrRows-1][2*j].z + TerrNorms[faceNormArrRows-1][2*j+1].z)/3;
         }else if(j == 0){ // Left Edge Vertices
            TerrVertNorms[i][j].x = (TerrNorms[i-1][0].x + TerrNorms[i-1][0].x + TerrNorms[i][0].x)/3;
            TerrVertNorms[i][j].y = (TerrNorms[i-1][0].y + TerrNorms[i-1][0].y + TerrNorms[i][0].y)/3;
            TerrVertNorms[i][j].z = (TerrNorms[i-1][0].z + TerrNorms[i-1][0].z + TerrNorms[i][0].z)/3;
         }else if(j == sqr_size-1){ // Right Edge Vertices
            TerrVertNorms[i][j].x = (TerrNorms[i-1][faceNormArrCols-1].x + TerrNorms[i][faceNormArrCols-2].x + TerrNorms[i][faceNormArrCols-1].x)/3;
            TerrVertNorms[i][j].y = (TerrNorms[i-1][faceNormArrCols-1].y + TerrNorms[i][faceNormArrCols-2].y + TerrNorms[i][faceNormArrCols-1].y)/3;
            TerrVertNorms[i][j].z = (TerrNorms[i-1][faceNormArrCols-1].z + TerrNorms[i][faceNormArrCols-2].z + TerrNorms[i][faceNormArrCols-1].z)/3;
         }else{ // All other Vertices are not on a edge or corrner
            TerrVertNorms[i][j].x = (TerrNorms[i-1][2*j-1].x + TerrNorms[i-1][2*j].x + TerrNorms[i-1][2*j+1].x + TerrNorms[i][2*j-2].x + TerrNorms[i][2*j-1].x + TerrNorms[i][2*j].x)/6;
            TerrVertNorms[i][j].y = (TerrNorms[i-1][2*j-1].y + TerrNorms[i-1][2*j].y + TerrNorms[i-1][2*j+1].y + TerrNorms[i][2*j-2].y + TerrNorms[i][2*j-1].y + TerrNorms[i][2*j].y)/6;
            TerrVertNorms[i][j].z = (TerrNorms[i-1][2*j-1].z + TerrNorms[i-1][2*j].z + TerrNorms[i-1][2*j+1].z + TerrNorms[i][2*j-2].z + TerrNorms[i][2*j-1].z + TerrNorms[i][2*j].z)/6;
         }
      }
   }
}


/* 
 * Uses a vertex and vertex normal mat created by createTerrMats in main to draw a terrain
 * sqr_size needs to match the arrays from createTerrMats
 * TerrVerts holds the x,y,z values of the vertices
 * TerrVertNorms holds the avg face norms for every vertex
 * rest params are same for any object
*/
void terrain(int sqr_size, vtx TerrVerts[sqr_size][sqr_size], vtx TerrVertNorms[sqr_size][sqr_size], double x, double y, double z, double dx, double dy, double dz, double rx, double ry, double rz, unsigned int texnum){
   //  Set specular color to white
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glRotated(rx,1,0,0);
   glRotated(ry,0,1,0);
   glRotated(rz,0,0,1);
   glScaled(dx,dy,dz);
   // Texture setup
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texnum);
   // Draw trianlges of the terrain
   glBegin(GL_TRIANGLES);
   glColor3f(.5,.5,.5);
   for(int i=0;i<sqr_size-1;i++){
      for(int j=0;j<sqr_size-1;j++){
         // square's upper trianlge
         glNormal3f(TerrVertNorms[i][j].x,TerrVertNorms[i][j].y,TerrVertNorms[i][j].z); glTexCoord2f(0,1); glVertex3f(TerrVerts[i][j].x,TerrVerts[i][j].y,TerrVerts[i][j].z);
         glNormal3f(TerrVertNorms[i][j+1].x,TerrVertNorms[i][j+1].y,TerrVertNorms[i][j+1].z); glTexCoord2f(1,1); glVertex3f(TerrVerts[i][j+1].x,TerrVerts[i][j+1].y,TerrVerts[i][j+1].z);
         glNormal3f(TerrVertNorms[i+1][j].x,TerrVertNorms[i+1][j].y,TerrVertNorms[i+1][j].z); glTexCoord2f(0,0); glVertex3f(TerrVerts[i+1][j].x,TerrVerts[i+1][j].y,TerrVerts[i+1][j].z);
         // square's lower trianlge
         glNormal3f(TerrVertNorms[i][j+1].x,TerrVertNorms[i][j+1].y,TerrVertNorms[i][j+1].z); glTexCoord2f(1,1); glVertex3f(TerrVerts[i][j+1].x,TerrVerts[i][j+1].y,TerrVerts[i][j+1].z);
         glNormal3f(TerrVertNorms[i+1][j].x,TerrVertNorms[i+1][j].y,TerrVertNorms[i+1][j].z); glTexCoord2f(0,0); glVertex3f(TerrVerts[i+1][j].x,TerrVerts[i+1][j].y,TerrVerts[i+1][j].z);
         glNormal3f(TerrVertNorms[i+1][j+1].x,TerrVertNorms[i+1][j+1].y,TerrVertNorms[i+1][j+1].z); glTexCoord2f(1,0); glVertex3f(TerrVerts[i+1][j+1].x,TerrVerts[i+1][j+1].y,TerrVerts[i+1][j+1].z);
      }
   }
   glEnd();
   //  Undo transformations
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();
}

/* 
 * This is a helper function for cliff - it creates a trianglar silce for the top of the cliff
 * Will create a triangluar slice that fits into the top of the passed Terrain
 * sqr_size needs to match the arrays from createTerrMats
 * TerrVerts holds the x,y,z values of the vertices
 * rest params are same for any object
*/
void terrainSlice(int side_num, int sqr_size, vtx TerrVerts[sqr_size][sqr_size], double x, double y, double z, double dx, double dy, double dz, double rx, double ry, double rz, unsigned int texnum){
   // double r - this is the radius of the circle instribe inside the regular poly. this means r is the distance from the center to any side at a perp
   double r = .5*(sqr_size-1)*(cos(M_PI/side_num)/sin(M_PI/side_num));
   //  Set specular color to white
   float white[] = {1,1,1,1};
   float black[] = {0,0,0,1};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glRotated(rx,1,0,0);
   glRotated(ry,0,1,0);
   glRotated(rz,0,0,1);
   glScaled(dx,dy,dz);
   // Texture setup
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV , GL_TEXTURE_ENV_MODE , GL_MODULATE);
   glBindTexture(GL_TEXTURE_2D,texnum);
   // Draw trianlges of the terrain
   glBegin(GL_TRIANGLES);
   glColor3f(.75,.75,.75);
   glNormal3f(0,0,-1);
   for(int i=0;i<sqr_size-1;i++){
      glTexCoord2f(0,0); glVertex3f(TerrVerts[0][i].x,TerrVerts[0][i].y,TerrVerts[0][i].z);
      glTexCoord2f(0,1); glVertex3f(TerrVerts[0][i+1].x,TerrVerts[0][i+1].y,TerrVerts[0][i+1].z);
      glTexCoord2f(r,0); glVertex3f(0,r*-1,(sqr_size/2)*-1);
   }
   glEnd();
   //  Undo transformations
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();
}

/* 
 * This function will create a cliff from terrains fliped vertically
 * side_num is how many faces the cliff will have 3 or greater
 * sqr_size needs to match the arrays from createTerrMats
 * TerrVerts holds the x,y,z values of the vertices
 * rest params are same for any object
*/
void cliff(int side_num, int sqr_size, vtx TerrVerts[sqr_size][sqr_size], vtx TerrVertNorms[sqr_size][sqr_size], double x, double y, double z, double dx, double dy, double dz, double rx, double ry, double rz, unsigned int texnum_sides, unsigned int texnum_top){
   // double r - this is the radius of the circle instribe inside the regular poly. this means r is the distance from the center to any side at a perp
   int side_width = sqr_size - 1;
   double r = .5*side_width*(cos(M_PI/side_num)/sin(M_PI/side_num));
   double interior_angle = 360/side_num;
   //  Save transformation
   glPushMatrix();
   //  Offset
   glTranslated(x,y,z);
   glRotated(rx,1,0,0);
   glRotated(ry,0,1,0);
   glRotated(rz,0,0,1);
   glScaled(dx,dy,dz);
   // rotate and translate each terrain to be a side of a cliff
   for(int i=0;i<side_num;i++){
      glPushMatrix();
      glTranslated(r*Sin(interior_angle*i),0,r*Cos(interior_angle*i));
      glRotated(interior_angle*i,0,1,0);
      glRotated(90,1,0,0);
      glScaled(1,1,1);
      terrain(sqr_size,TerrVerts,TerrVertNorms,0,0,0,1,1,1,0,0,0,texnum_sides);
      terrainSlice(side_num,sqr_size,TerrVerts,0,0,0,1,1,1,0,0,0,texnum_top);
      glPopMatrix();
   }
   //  Undo transformations
   glPopMatrix();
}

// Draws the eastern air temple scene
void easternAirTemple(double x, double y, double z){
      // center lower cliff with temple
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y-195,z-160,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y+-375,z-160,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      // rock and gran frustums on this cliff
      frustum(x+15,y+-128,z+-151,40,135,135,135,135,20,"wrapped",texture[38],texture[38],texture[38],texture[38],1.0,1.0,1.0);
      frustum(x+15,y+-88,z+-152,10,125,125,110,110,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      frustum(x+15,y+-78,z+-152,10,100,100,85,85,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      // buildings on these^ frustums
      frustum(x+13,y+-68,z+-169,15,30,50,30,50,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // temple block
      frustum(x+41,y+-68,z+-162,15,30,40,30,40,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // cylander block
      frustum(x+18,y+-68,z+-143,15,5,20,5,2,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // temple stairs
      frustum(x+11,y+-53,z+-172,20,25,45,25,45,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // main part of temple
      frustum(x+11,y+-33,z+-172,20,25,45,0,45,20,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // temple roof
      frustum(x+13,y+-53,z+-167,15,20,40,20,40,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // inner temple face
      frustum(x+13,y+-38,z+-167,15,20,40,0,40,20,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+15,y+-53,z+-164,10,16,36,16,36,20,"wrapped",texture[31],texture[2],texture[2],texture[2],1.0,1.0,1.0); // inner inner temple face
      frustum(x+15,y+-43,z+-164,10,16,36,0,36,20,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+15,y+-53,z+-172,42,8,8,8,8,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // chimney ^
      frustum(x+15,y+-11,z+-172,6,8,8,8,8,20,"wrapped",texture[28],texture[28],texture[28],texture[28],1.0,1.0,1.0); // chimney topping
      cylinder(x+47,y+-53,z+-155,10,35,0,texture[2]); // cylander next to temple
      cylinderCap(x+47,y+-18,z+-155,8,12,8,15,texture[28],texture[29]); // cap
      // Bridge conecting small forward cliff to to center cliff
      bridge(x+46,y+-95,z+-68,40,40,60,110,0,.5,texture[2]);
      // small center forward cliff
      cliff(5,11,TerrVerts10,TerrVertNorms10,x+102,y+-162,z+-6,8,20,8,0,0,0,texture[38],texture[2]);
      cliff(5,11,TerrVerts10,TerrVertNorms10,x+102,y+-340,z+-6,8,20,8,180,180,0,texture[38],texture[2]); // extension of cliff
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+93,y+-82,z+-18,2,8,2,0,0,0,texture[38],texture[2]); // small cliff in center
      cliff(3,11,TerrVerts10,TerrVertNorms10,x+89,y+-88,z+3,2,8,2,0,0,0,texture[38],texture[2]); // small path off small cliff
      // large center cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-75,z+-363,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-270,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-452,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+39,y+11,z+-365,19,18,17,0,0,0,texture[38]);
      // low left cliff from/in center cliff
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+-34,y+-218,z+-267,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+-34,y+-398,z+-267,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      // mid right cliff from/in center cliff
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+113,y+-142,z+-292,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+113,y+-322,z+-292,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      // stairs on center cliffs
      frustum(x+92,y+-96,z+-203,55,55,8,10,8,100,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // lowwer
      frustum(x+100,y+-44,z+-254,70,60,8,10,8,0,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // middle
      frustum(x+46,y+24,z+-297,75,65,8,10,8,35,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // upper
      frustum(x+65,y+82,z+-342,85,65,8,10,8,90,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // upper to peak
      frustum(x+-40,y+-149,z+-225,55,55,8,10,8,110,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // backside lower left
      // Cylander buildings on central cliff
      cylinder(x+59,y+25,z+-274,10,35,0,texture[2]); // middle
      cylinderCap(x+59,y+60,z+-274,9,11,9,60,texture[28],texture[29]); // cap
      cylinder(x+35,y+158,z+-360,20,82,0,texture[2]); // upper center
      cylinderCap(x+35,y+240,z+-360,18,22,18,30,texture[28],texture[29]); // cap
      cylinder(x+10,y+126,z+-386,12,82,0,texture[2]); // upper left
      cylinderCap(x+10,y+208,z+-386,10,16,10,25,texture[28],texture[29]); // cap
      cylinder(x+78,y+95,z+-401,12,82,0,texture[2]); // upper right
      cylinderCap(x+78,y+177,z+-401,10,16,10,30,texture[28],texture[29]); // cap
      // houses on center cliff
      frustum(x+141,y+24,z+-369,20,25,45,25,45,-38,"wrapped",texture[34],texture[34],texture[32],texture[2],1.0,1.0,1.0); // middle right
      frustum(x+141,y+44,z+-369,20,25,45,0,45,-38,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+-20,y+24,z+-284,20,25,45,25,45,50,"wrapped",texture[31],texture[2],texture[2],texture[32],1.0,1.0,1.0); // middle left
      frustum(x+-20,y+44,z+-284,20,25,45,0,45,50,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+170,y+-42,z+-279,24,35,45,35,45,-45,"wrapped",texture[34],texture[34],texture[33],texture[2],1.0,1.0,1.0); // lower right
      frustum(x+170,y+-18,z+-279,20,35,45,0,45,-45,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+-76,y+-118,z+-273,24,35,45,35,45,15,"wrapped",texture[34],texture[2],texture[2],texture[33],1.0,1.0,1.0); // very low back left
      frustum(x+-76,y+-94,z+-273,24,35,45,0,45,15,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof

      // bridge from center to right cliff
      bridge(x+230,y+-67,z+-330, 42, 30, 50, 20, .3, .9, texture[2]);

      // Right cliff
      glPushMatrix();
      glTranslated(x+240,y+-62,z+-81);
      glRotated(-20,0,1,0);
      glScaled(.8,.8,.8);
      // center lower cliff with temple
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y-195,z-160,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y+-375,z-160,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      // rock and gran frustums on this cliff
      frustum(x+15,y+-128,z+-151,40,135,135,135,135,20,"wrapped",texture[38],texture[38],texture[38],texture[38],1.0,1.0,1.0);
      frustum(x+15,y+-88,z+-152,10,125,125,110,110,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      frustum(x+15,y+-78,z+-152,10,100,100,85,85,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      // temple on these frustums
      frustum(x+11,y+-68,z+-172,20,25,45,25,45,-70,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // main part of temple
      frustum(x+11,y+-48,z+-172,20,25,45,0,45,-70,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // temple roof
      frustum(x+6,y+-68,z+-170,15,20,40,20,40,-70,"wrapped",texture[31],texture[2],texture[2],texture[2],1.0,1.0,1.0); // inner temple face
      frustum(x+6,y+-53,z+-170,15,20,40,0,40,-70,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+15,y+-68,z+-164,10,16,36,16,36,-70,"wrapped",texture[2],texture[2],texture[32],texture[2],1.0,1.0,1.0); // inner inner temple face
      frustum(x+15,y+-58,z+-164,10,16,36,0,36,-70,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+15,y+-69,z+-172,42,8,8,8,8,-70,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // chimney ^
      frustum(x+15,y+-28,z+-172,6,8,8,8,8,-70,"wrapped",texture[28],texture[28],texture[28],texture[28],1.0,1.0,1.0); // chimney topping
      cylinder(x+47,y+-68,z+-155,10,35,0,texture[2]); // cylander next to temple
      cylinderCap(x+47,y+-33,z+-155,8,12,8,15,texture[28],texture[29]); // cap
      // large center cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-75,z+-363,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-270,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-452,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+39,y+11,z+-365,19,18,17,0,0,0,texture[38]); // center peak
      // Cylander on top of peak
      cylinder(x+35,y+158,z+-360,20,82,0,texture[2]); // upper center
      cylinderCap(x+35,y+240,z+-360,18,22,18,30,texture[28],texture[29]); // cap
      // Cylander in the middle, right of peak
      cylinder(x+135,y+25,z+-331,15,36,0,texture[2]);
      cylinderCap(x+135,y+61,z+-331,13,12,13,30,texture[28],texture[29]); // cap
      // Cylander in the middle, left of peak
      cylinder(x+-55,y+25,z+-382,12,62,0,texture[2]);
      cylinderCap(x+-55,y+87,z+-382,10,12,10,30,texture[28],texture[29]); // 
      // houses
      frustum(x+-15,y+25,z+-281,40,25,30,25,30,-35,"wrapped",texture[34],texture[2],texture[32],texture[32],1.0,1.0,1.0); // middle tall just left of stairs
      frustum(x+-15,y+65,z+-281,28,25,30,0,30,-35,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      // stairs
      frustum(x+46,y+24,z+-297,75,65,8,10,8,35,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // upper
      frustum(x+65,y+82,z+-342,85,65,8,10,8,90,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // upper to peak

      glPopMatrix();

      // bridge from lower center to left cliff
      bridge(x+-171,y+-141,z+-218, 60, 30, 50, 30, .3, .9, texture[2]);
      //  bridge from upper center to left cliff
      bridge(x+-124,y+-59,z+-392, 48, 30, 50, 20, .3, .9, texture[2]);

      // left cliff
      glPushMatrix();
      glTranslated(x+-252,y+-59,z+-48);
      glRotated(10,0,1,0);
      glScaled(.85,.85,.85);
      // center lower cliff with temple
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y-195,z-160,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+20,y+-375,z-160,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      // rock and gran frustums on this cliff
      frustum(x+15,y+-128,z+-151,40,135,135,135,135,20,"wrapped",texture[38],texture[38],texture[38],texture[38],1.0,1.0,1.0);
      frustum(x+15,y+-88,z+-152,10,125,125,110,110,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      frustum(x+15,y+-78,z+-152,10,100,100,85,85,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
      frustum(x+-32,y+-88,z+-154,20,50,50,50,50,20,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // left brick
      // temples on frustums
      cylinder(x+-34,y+-68,z+-156,14,48,0,texture[2]);// cylander
      cylinderCap(x+-34,y+-20,z+-156,12,12,12,20,texture[28],texture[29]); // cap
      frustum(x+14,y+-68,z+-165,40,50,35,50,35,20,"wrapped",texture[34],texture[2],texture[32],texture[2],1.0,1.0,1.0); 
      frustum(x+14,y+-28,z+-165,38,50,35,0,35,20,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // roof
      // large center cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-75,z+-363,10,20,10,0,0,0,texture[38],texture[2]);
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-270,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      cliff(8,11,TerrVerts10,TerrVertNorms10,x+41,y+-452,z+-363,10,20,10,180,0,0,texture[38],texture[2]); // extension of cliff
      cliff(6,11,TerrVerts10,TerrVertNorms10,x+8,y+94,z+-458,6,20,6,0,45,0,texture[38],texture[2]); // top cliff behind peak
      terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+39,y+11,z+-365,19,18,17,0,0,0,texture[38]);
      // cylander one peak
      cylinder(x+8,y+193,z+-449,30,70,0,texture[2]);// one on peak
      cylinderCap(x+8,y+263,z+-449,27,25,27,30,texture[28],texture[29]); // cap
      // houses
      frustum(x+49,y+25,z+-266,30,35,45,35,45,90,"wrapped",texture[2],texture[31],texture[2],texture[32],1.0,1.0,1.0); // middle left
      frustum(x+49,y+55,z+-266,28,35,45,0,45,90,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      frustum(x+109,y+25,z+-297,30,35,45,35,45,115,"wrapped",texture[31],texture[2],texture[2],texture[33],1.0,1.0,1.0); // middle right
      frustum(x+109,y+55,z+-297,25,35,45,0,45,115,"wrapped",texture[2],texture[2],texture[28],texture[28],1.0,1.0,1.0); // ^ roof
      // stairs
      frustum(x+20,y+106,z+-354,90,95,10,12,10,90,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // to peak
      frustum(x+-1,y+17,z+-310,90,125,10,12,10,0,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // upper
      frustum(x+16,y+77,z+-320,30,30,30,30,30,0,"wrapped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // connect these two stairs with this box
      glPopMatrix();

      // trees
      // center cliffs
      tree_1(x+69,y+-42,z+-221,1.6,2,1.6,0,texture[7],texture[6]);
      tree_1(x+120,y+-42,z+-221,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+140,y+-42,z+-235,1.6,1.7,1.6,0,texture[7],texture[6]);
      tree_1(x+153,y+25,z+-333,1.6,1.6,1.6,0,texture[7],texture[6]);
      tree_1(x+77,y+25,z+-263,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+45,y+25,z+-263,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-50,y+25,z+-297,1.6,1.9,1.6,0,texture[7],texture[6]);
      tree_1(x+-67,y+25,z+-325,1.6,2,1.6,0,texture[7],texture[6]);
      tree_1(x+-60,y+25,z+-403,1.6,1.7,1.6,0,texture[7],texture[6]);
      tree_1(x+71,y+109,z+-313,1.6,1.6,1.6,0,texture[7],texture[6]);
      tree_1(x+83,y+127,z+-353,1.6,2,1.6,0,texture[7],texture[6]);
      tree_1(x+-16,y+127,z+-353,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+1,y+144,z+-353,1.6,1.7,1.6,0,texture[7],texture[6]);
      tree_1(x+-13,y+96,z+-334,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-28,y+79,z+-334,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-109,y+-118,z+-297,1.6,1.9,1.6,0,texture[7],texture[6]);
      tree_1(x+-81,y+-118,z+-220,1.6,1.7,1.6,0,texture[7],texture[6]);
      tree_1(x+26,y+120,z+-314,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+107,y+88,z+-337,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+115,y+72,z+-362,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-29,y+96,z+-362,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+18,y+134,z+-332,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+33,y+106,z+-307,1.6,1.8,1.6,0,texture[7],texture[6]);
      // right cliff
      tree_1(x+350,y+-42,z+-263,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+279,y+-42,z+-323,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+308,y+-42,z+-266,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+418,y+30,z+-339,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+411,y+47,z+-339,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+338,y+47,z+-339,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+325,y+28,z+-339,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+319,y+12,z+-339,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+343,y+21,z+-312,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+382,y+18,z+-293,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+396,y+-15,z+-283,1.6,1.8,1.6,0,texture[7],texture[6]);
      // left cliff
      tree_1(x+-293,y+-38,z+-276,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-351,y+-38,z+-324,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-375,y+-38,z+-354,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-252,y+25,z+-312,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-256,y+69,z+-378,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-226,y+53,z+-371,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-274,y+106,z+-440,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-279,y+106,z+-412,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-356,y+77,z+-412,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-357,y+37,z+-412,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-363,y+1,z+-440,1.6,1.8,1.6,0,texture[7],texture[6]);
      tree_1(x+-315,y+58,z+-363,1.6,1.8,1.6,0,texture[7],texture[6]);
}

// draws the fire castle scene
void fireNationCastle(double x, double y, double z)
{
	//Base
    frustum(x-2,y,z-75,10,109,50,109,50,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-2,y,z-7,10,109,109,109,109,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);

    //Front Border
    frustum(x-32,y+10,z+45,10,49,5,48,3,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-32,y+20,z+45,1,48,3,48,3,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+28,y+10,z+45,10,49,5,48,3,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+28,y+20,z+45,1,48,3,48,3,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    // front gate
    frustum(x-7.5,y+10,z+46,10,3,7,3,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+3.5,y+10,z+46,10,3,7,3,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-2,y+20,z+46,3,8,7,8,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-7.5,y+20,z+46,2,5,9,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+3.5,y+20,z+46,2,5,9,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-2,y+23,z+46,2,10,9,10,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-2,y+10,z+44,8,8,1,8,1,0,"chopped",texture[5],texture[5],texture[5],texture[5],1.0,1.0,1.0);


    //Right Border
    frustum(x+50,y+10,z,10,5,90,3,85,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+50,y+20,z,1,3,85,3,85,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+50,y+10,z-55,10,5,90,3,85,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+50,y+20,z-55,1,3,85,3,85,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Left Border
    frustum(x-55,y+10,z,10,5,90,3,85,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-55,y+20,z,1,3,85,3,85,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-55,y+10,z-55,10,5,90,3,85,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-55,y+20,z-55,1,3,85,3,85,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Entrance Buildings
    frustum(x,y+10,z+17,12,20,15,20,15,0,"chopped",texture[25],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x,y+22,z+17,.1,20,15,20,15,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // cover top tex
    frustum(x,y+22,z+17,4,18,13,18,13,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x,y+26,z+17,3,25,17,18,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x,y+25,z+17,7,18,10,18,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-3,y+32,z+17,4,16,14,14,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+5,y+32,z+17,4,8,14,0,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Top Right
    frustum(x+23,y+10,z+17,10,13,19,13,19,0,"chopped",texture[25],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+23,y+20,z+17,7,16,22,0,16,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Middle right
    frustum(x+47.5,y+10,z-4.5,18,10,23,10,23,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+47.5,y+28,z-4.5,5,11,23,0,23,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Top Left
    frustum(x-33,y+10,z-3,12,15,38,13,30,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-33,y+10,z-10,12,15,38,13,30,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-33,y+22,z+1,5,11,20,11,20,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-33,y+22,z-15,5,11,20,11,20,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-33,y+27,z-5,2,13,35,5,27,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-33,y+27,z-7,5,13,32,0,30,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    // Middle Buildings
    frustum(x-5,y+10,z-30,20,35,35,35,35,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-23,y+10,z-40,20,35,35,35,35,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+15,y+10,z-38,20,20,10,15,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x,y+30,z-27,7,21,25,21,25,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // a
    frustum(x-1,y+37,z-27,4,26,29,8,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0); // a

    frustum(x-17,y+10,z-10,35,30,27,26,23,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-17,y+45,z-10,10,23,20,23,20,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-17,y+50,z-10,5,30,27,17,15,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-17,y+55,z-10,8,30,27,17,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    frustum(x-5,y+30,z-35,50,20,30,14,24,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); //tower
    frustum(x-5,y+80,z-35,3,14,24,14,24,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-5,y+80,z-35,6,13,23,13,23,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-5,y+86,z-35,3,15,25,10,15,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-5,y+89,z-35,2,10,15,10,15,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-5,y+91,z-35,3,12,17,8,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-5,y+94,z-35,1.5,8,10,8,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-5,y+95.5,z-35,7,12,17,0,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);


    // Right Tower Wing
    frustum(x+12,y+30,z-38,13,21,7,21,7,0,"chopped",texture[2],texture[2],texture[23],texture[2],1.0,1.0,1.0);
    frustum(x+7,y+43,z-38,8,15,7,15,7,0,"chopped",texture[2],texture[2],texture[23],texture[2],1.0,1.0,1.0); 
    frustum(x+7,y+51,z-38,8,7,7,7,7,0,"chopped",texture[2],texture[2],texture[23],texture[2],1.0,1.0,1.0); frustum(x+6,y+51,z-38,8,7,7,7,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);  
    frustum(x+3,y+59,z-38,8,7,7,7,7,0,"chopped",texture[2],texture[2],texture[23],texture[2],1.0,1.0,1.0); 
    // RTW Roofs
    frustum(x+19,y+43,z-38,3,9,9,9,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+13,y+51,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+13,y+51,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+9,y+59,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+5,y+67,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);


    // Left Tower Wing
    frustum(x-23,y+30,z-38,13,21,7,21,7,0,"chopped",texture[2],texture[2],texture[2],texture[23],1.0,1.0,1.0);
    frustum(x-18,y+43,z-38,8,15,7,15,7,0,"chopped",texture[2],texture[2],texture[2],texture[23],1.0,1.0,1.0); 
    frustum(x-17.5,y+51,z-38,8,8,7,8,7,0,"chopped",texture[2],texture[2],texture[2],texture[23],1.0,1.0,1.0); 
    frustum(x-14,y+59,z-38,8,7,7,7,7,0,"chopped",texture[2],texture[2],texture[2],texture[23],1.0,1.0,1.0);
    // LTW Roofs
    frustum(x-30,y+43,z-38,3,9,9,9,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-24,y+51,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-20,y+59,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-15,y+67,z-38,3,5,8,5,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0); 


    // Front Walkway
    frustum(x-13,y+10,z-3,10,40,20,40,20,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); //walkway
    frustum(x+4.5,y+20,z-16,10,5,25,5,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); // stairs
    frustum(x+15,y+10,z-38,20,5,25,5,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);

    // Back Buildings
    frustum(x-5.5,y+20,z-52,10,15,11,0,11,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); //stairs - upper
    frustum(x-7,y+10,z-59,10,15,11,0,11,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); //stairs - lower
    frustum(x+8,y+10,z-55,10,30,30,30,30,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+15,y+20,z-45,10,15,7,15,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+15,y+30,z-45,5,20,7,13,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+20,y+20,z-50,10,7,15,7,15,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+20,y+30,z-50,5,10,15,0,15,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    frustum(x-23.5,y+30,z-45,10,20,7,20,7,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-24,y+40,z-45,5,21,9,20,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);

    frustum(x-30,y+12,z-80,4,18,13,18,13,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-30,y+16,z-80,3,25,17,18,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-30,y+15,z-80,7,18,10,18,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x-33,y+22,z-80,4,16,14,14,0,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x-25,y+22,z-80,4,8,14,0,10,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);


    // Details
    cube(x+39,y+10.001,z+16.2, 10,0.001,8.75, texture[37], 1,1,1, 1);
    cube(x-1,y+9,z+185, 53,.1,38, texture[37], 1,1,1, 1); // grass at the start of the path
    cube(x+33,y+10.001,z-70, 15,0.001,30, texture[37], 1,1,1, 1);
    frustum(x+39,y+10,z+25,6,20,3,20,2,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); frustum(x+39,y+16,z+27,1,5,2,5,2,0,"chopped",texture[3],texture[3],texture[3],texture[3],1.0,1.0,1.0);
    frustum(x+41,y+10,z+27,6,1,1,1,1,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); frustum(x+37,y+10,z+27,6,1,1,1,1,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0); 
    frustum(x+39,y+10,z+9,6,20,3,20,2,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+34,y+10,z-40,6,30,3,30,2,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    tree_1(x+33,y+10,z-70,1,1,1,0,texture[7],texture[6]);
    tree_1(x+35,y+10,z-75,1,1,1,45,texture[7],texture[6]);
    tree_1(x+30,y+10,z-60,1,1,1,60,texture[7],texture[6]);
    tree_1(x+42,y+10,z-53,1,1,1,90,texture[7],texture[6]);
    tree_1(x+27,y+10,z-45,1,1,1,205,texture[7],texture[6]);
    tree_1(x+30,y+10,z-50,1,1,1,70,texture[7],texture[6]);
    tree_1(x+40,y+10,z-60,1,1,1,10,texture[7],texture[6]);
    tree_1(x+25,y+10,z-85,1,1,1,60,texture[7],texture[6]);
    tree_1(x+43,y+10,z-90,1,1,1,300,texture[7],texture[6]);
    tree_1(x+33,y+10,z-97,1,1,1,250,texture[7],texture[6]);
    tree_1(x+34,y+10,z-87,1,1,1,150,texture[7],texture[6]);
    tree_1(x+42,y+10,z-70,1,1,1,0,texture[7],texture[6]);
    tree_1(x+42,y+10,z-65,1,1,1,90,texture[7],texture[6]);
    tree_1(x+42,y+10,z-75,1,1,1,330,texture[7],texture[6]);



    // Terrain
    // Main square cliff the castle is on
    cliff(4,21,TerrVerts20,TerrVertNorms20,x+-2,y+-37,z+-10,5.8,4,9,0,0,0,texture[38],texture[38]);
    cliff(4,21,TerrVerts20,TerrVertNorms20,x+-1,y+-31,z+129,5.8,4,9,0,0,0,texture[38],texture[38]);
    // terrain by front gate
    terrain(11,TerrVerts10CCUp,TerrVertNorms10CCUp,x+-34,y+9.1,z+97,10,4,5,0,90,0,texture[38]);
    terrain(11,TerrVerts10CCUp,TerrVertNorms10CCUp,x+30,y+9.1,z+97,10,4,5,0,90,0,texture[38]);
    // path from gate
    frustum(x+-2,y+9,z+97,1,14,100,14,100,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    // back side cliff
    cliff(3,11,TerrVerts10,TerrVertNorms10,x+-2,y+-17,z+-123,11,4,10,0,0,0,texture[38],texture[2]);
    // patch of rocks on back cliff
    // left side
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+35,y+3,z+-104,1,1,1,0,0,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+32,y+3,z+-115,1,1,1,0,0,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+23,y+3,z+-126,1,1,1,0,45,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+14,y+3,z+-140,1,1,1,0,5,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+10,y+2,z+-150,.8,.8,.8,0,45,0,texture[2]);
    // right side
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+-35,y+3,z+-104,1,1,1,0,0,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+-32,y+3,z+-115,1,1,1,0,0,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+-23,y+3,z+-126,1,1,1,0,45,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+-14,y+3,z+-140,1,1,1,0,5,0,texture[2]);
    terrain(11,TerrVertsPeak,TerrVertNormsPeak,x+-8,y+2,z+-154,.8,.8,.8,0,45,0,texture[2]);
    // Ramp to back side cliff
    frustum(x+-2,y,z+-104,10,30,10,11,10,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
    frustum(x+-2,y,z+-108,15,39.5,2,11,2,0,"chopped",texture[2],texture[2],texture[2],texture[2],1.0,1.0,1.0);
}

// Filler House - type 1
void ccHouse_1(double x, double y, double z, double xi, double zi, double h, double theta)
{
   frustum(x+xi,y,z+zi,h,6.81,3.46,6.81,3.46,theta,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+xi,y+h,z+zi,3,7,3.7,7,0,theta,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
}

// Filler House - type 2
void ccHouse_2(double x, double y, double z, double xi, double zi, double h, double theta)
{
   frustum(x+xi,y,z+zi,h,6.81,3.46,6.81,3.46,theta,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+xi,y+h,z+zi,3,7,3.7,7,0,theta,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
}

// Filler House - type 3
void ccHouse_3(double x, double y, double z, double xi, double zi, double h, double theta)
{
   frustum(x+xi,y,z+zi,h,6.81,3.46,6.81,3.46,theta,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x+xi,y+h,z+zi,3,7,3.7,7,0,theta,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
}

// draws trees inbetween an inner and outer circle -- to surround chin city with trees
void treeFinder(double x, double y, double z, double xi, double yi, double zi, double r_i, double r_o, int step)
{
   // Increment for spacing the grid
   double incrementGrid = 2.0*r_o/(step-1.0);

   double outerGrid_x[step];
   double outerGrid_z[step];

   // x and z will be the current element in question
   double x_c = 0;
   double z_c = 0;
   double distance = 0;

   // number of correct matches
   int n = 0;

   // creating grid elements
   for(int i=0;i<step;i++){
      outerGrid_x[i] = -r_o+incrementGrid*i+xi;
      outerGrid_z[i] = -r_o+incrementGrid*i+zi;
   }

   // Based on each z value, is the corresponding x value out of the inner circle and in the outer circle
   for(int j=0;j<step;j++)
   {
      z_c = outerGrid_z[j];
      for(int k=0;k<step;k++)
      {
         x_c = outerGrid_x[k];
         // finding the distance from (x,z) to center of circle
         distance = sqrt((x_c-xi)*(x_c-xi)+(z_c-zi)*(z_c-zi));
         if(distance>r_i && distance<r_o){
            n++;
            tree_1(x+x_c,y,z+z_c,0.75,1,0.75,n%360,texture[7],texture[6]);
         }
      }
   }


}

/*
   Draws the temple in chin city using frustums
 */
static void chinCityTemple(double x,double y,double z,double dx,double dy,double dz,double ry,unsigned int texnum_base,unsigned int texnum_wall,unsigned int texnum_roof,unsigned int texnum_door)
{
	
   //  Save transformation
   glPushMatrix();
   //  Offset and scale
   glTranslated(x,y,z);
   glRotatef(ry,0,1,0);
   glScaled(dx,dy,dz);

   // Base blocks
   frustum(-30,0,0,20,40,100,40,100,0,"chopped",texnum_base,texnum_base,texnum_base,texnum_base,1.0,1.0,1.0); // left
   frustum(+30,0,0,20,40,100,40,100,0,"chopped",texnum_base,texnum_base,texnum_base,texnum_base,1.0,1.0,1.0); // right
   frustum(0,0,-5,20,20,90,20,90,0,"chopped",texnum_base,texnum_base,texnum_base,texnum_base,1.0,1.0,1.0); // center behind stairs
   frustum(0,0,40,20,20,20,20,0,0,"chopped",texnum_base,texnum_base,texnum_base,texnum_base,1.0,1.0,1.0); // stairs
   //first floor
   frustum(0,20,0,25,60,60,60,60,0,"chopped",texnum_door,texnum_wall,texnum_wall,texnum_wall,1.0,1.0,1.0); // base box
   frustum(0,50,0,15,73,73,54,54,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // right side up roof
   // second floor
   frustum(0,65,0,12,54,54,54,54,0,"chopped",texnum_wall,texnum_wall,texnum_wall,texnum_wall,1.0,1.0,1.0); // base box
   frustum(0,82,0,10,66,66,40,40,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // right side up roof
   // third floor
   frustum(0,92,0,6,40,40,40,40,0,"chopped",texnum_wall,texnum_wall,texnum_wall,texnum_wall,1.0,1.0,1.0); // base box
   frustum(0,110,0,15,58,58,30,10,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // right side up roof

   // top roof cap
   frustum(0,125,0,10,30,10,15,0,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0);


   // Upside down frustums
   glPushMatrix();
   glRotatef(180,1,0,0);
   frustum(0,-50,0,5,73,73,60,60,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // first floor
   frustum(0,-82,0,5,66,66,54,54,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // second floor
   frustum(0,-110,0,12,58,58,40,40,0,"chopped",texnum_roof,texnum_roof,texnum_roof,texnum_roof,1.0,1.0,1.0); // third floor
   glPopMatrix();
   
   glPopMatrix();
}

// draws the cin chity scene
void chinCity(double x, double y, double z)
{
   treeFinder(x,y,z,2.66,0.0,-6.62, 70, 110, 40);

	tree_1(x-0.6,y,z+38.55,0.75,0.75,0.75,0,texture[7],texture[6]);
	tree_1(x-7.72,y,z+45.7,0.5,0.5,0.5,45,texture[7],texture[6]);
	tree_1(x-11.30,y,z+52.43,0.6,0.6,0.6,-45,texture[7],texture[6]);
	tree_1(x-0.16,y,z+47.59,0.5,0.5,0.5,0,texture[7],texture[6]);
	tree_1(x+15,y,z+48.69,0.6,0.75,0.3,90,texture[7],texture[6]);
	tree_1(x+18,y,z+54.27,0.75,0.9,0.75,0,texture[7],texture[6]);
	tree_1(x-28.21,y,z+34.18,0.5,0.7,0.5,30,texture[7],texture[6]);

	cylinder(x,y+0.001,z,25,0.1,0,texture[36]);
   cube(0,-0.001,0, 115,0.001,115, texture[2], 240.0/255.0,172.0/255.0,121.0/255.0, 5);

	// First radii of houses
	frustum(x-22.88,y,z-17.52,5,6.81,3.46,6.81,3.46,60,"chopped",texture[8],texture[9],texture[13],texture[12],1.0,1.0,1.0); frustum(x-22.88,y+5,z-17.52,3,7,3.7,7,0,60,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-16.82,y,z-23.37,5,6.81,3.46,6.81,3.46,35,"chopped",texture[10],texture[11],texture[12],texture[13],1.0,1.0,1.0); frustum(x-16.82,y+5,z-23.37,3,7,3.7,7,0,35,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-8.97,y,z-26.81,5,6.81,3.46,6.81,3.46,15,"chopped",texture[16],texture[16],texture[12],texture[12],1.0,1.0,1.0); frustum(x-8.97,y+5,z-26.81,3,7,3.7,7,0,15,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x-1.99,y,z-29.38,5,6.81,3.46,6.81,3.46,90,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x-1.99,y+5,z-29.38,3,7,3.7,7,0,90,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
	frustum(x+4.5,y,z-27.85,5,6.81,3.46,6.81,3.46,-5,"chopped",texture[10],texture[10],texture[13],texture[12],1.0,1.0,1.0); frustum(x+4.5,y+5,z-27.85,3,7,3.7,7,0,-5,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+12.47,y,z-25.95,5,6.81,3.46,6.81,3.46,-30,"chopped",texture[17],texture[19],texture[18],texture[19],1.0,1.0,1.0); frustum(x+12.47,y+5,z-25.95,3,7,3.7,7,0,-30,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
	frustum(x+18.92,y,z-22.83,5,5.86,2.98,5.86,2.98,-45,"chopped",texture[16],texture[16],texture[12],texture[12],1.0,1.0,1.0); frustum(x+18.92,y+5,z-22.83,3,7,3.7,7,0,-45,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+23.31,y,z-17.96,5,5.99,3.05,5.99,3.05,-45,"chopped",texture[11],texture[13],texture[13],texture[12],1.0,1.0,1.0); frustum(x+23.31,y+5,z-17.96,3,7,3.7,7,0,-45,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	//
	frustum(x+29.39,y,z-7.01,5,6.81,3.46,6.81,3.46,0,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+29.39,y+5,z-7.01,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	//
	frustum(x+29.39,y,z+7.01,5,6.81,3.46,6.81,3.46,0,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x+29.39,y+5,z+7.01,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+27.64,y,z+13.37,5,6.81,3.46,6.81,3.46,-20,"chopped",texture[8],texture[8],texture[12],texture[12],1.0,1.0,1.0); frustum(x+27.64,y+5,z+13.37,3,7,3.7,7,0,-20,"chopped",texture[5],texture[5],texture[5],texture[5],1,1,1);
	frustum(x+22.43,y,z+18.63,5,6.81,3.46,6.81,3.46,50,"chopped",texture[19],texture[17],texture[18],texture[19],1.0,1.0,1.0); frustum(x+22.43,y+5,z+18.63,3,7,3.7,7,0,50,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
	frustum(x+18.06,y,z+23.56,5,5.86,2.98,5.86,2.98,-50,"chopped",texture[17],texture[17],texture[19],texture[18],1.0,1.0,1.0); frustum(x+18.06,y+5,z+23.56,3,7,3.7,7,0,-50,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
	frustum(x+13.05,y,z+26.94,5,6.81,3.46,6.81,3.46,-70,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+13.05,y+5,z+26.94,3,7,3.7,7,0,-70,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
	frustum(x+5.44,y,z+28.07,5,6.81,3.46,6.81,3.46,10,"chopped",texture[11],texture[11],texture[12],texture[13],1.0,1.0,1.0); frustum(x+5.44,y+5,z+28.07,3,7,3.7,7,0,10,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x-3.74,y,z+28.74,5,6.81,3.46,6.81,3.46,0,"chopped",texture[10],texture[11],texture[13],texture[12],1.0,1.0,1.0); frustum(x-3.74,y+5,z+28.74,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-11.36,y,z+28.07,5,6.81,3.46,6.81,3.46,0,"chopped",texture[8],texture[8],texture[12],texture[12],1.0,1.0,1.0); frustum(x-11.36,y+5,z+28.07,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-22.99,y,z+20.69,5,6.81,3.46,6.81,3.46,-45,"chopped",texture[9],texture[8],texture[12],texture[13],1.0,1.0,1.0); frustum(x-22.99,y+5,z+20.69,3,7,3.7,7,0,-45,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	// temple
   chinCityTemple(x-34.22,y,z,.15,.15,.15,90,texture[1],texture[21],texture[5],texture[22]);

	// Rest of houses
	frustum(x-28.33,y,z-26.81,5,6.81,3.46,6.81,3.46,0,"chopped",texture[8],texture[9],texture[12],texture[12],1.0,1.0,1.0); frustum(x-28.33,y+5,z-26.81,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-19.73,y,z-27.85,5,6.81,3.46,6.81,3.46,15,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x-19.73,y+5,z-27.85,3,7,3.7,7,0,15,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-19.73,y,z-33.25,5,6.81,3.46,6.81,3.46,0,"chopped",texture[10],texture[10],texture[12],texture[10],1.0,1.0,1.0); frustum(x-19.73,y+5,z-33.25,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-12.64,y,z-33.25,5,6.81,3.46,6.81,3.46,90,"chopped",texture[9],texture[9],texture[15],texture[12],1.0,1.0,1.0); frustum(x-12.64,y+5,z-33.25,3,7,3.7,7,0,90,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-5.97,y,z-36.08,5,6.81,3.46,6.81,3.46,-45,"chopped",texture[8],texture[8],texture[12],texture[15],1.0,1.0,1.0); frustum(x-5.97,y+5,z-36.08,3,7,3.7,7,0,-45,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+3.78,y,z-34.29,5,6.81,3.46,6.81,3.46,90,"chopped",texture[9],texture[9],texture[12],texture[12],1.0,1.0,1.0); frustum(x+3.78,y+5,z-34.29,3,7,3.7,7,0,90,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x+10.27,y,z-31.92,5,6.81,3.46,6.81,3.46,0,"chopped",texture[16],texture[11],texture[12],texture[15],1.0,1.0,1.0); frustum(x+10.27,y+5,z-31.92,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+13.05,y,z-37.11,5,6.81,3.46,6.81,3.46,-25,"chopped",texture[17],texture[19],texture[18],texture[19],1.0,1.0,1.0); frustum(x+13.05,y+5,z-37.11,3,7,3.7,7,0,-25,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	frustum(x+18.92,y,z-29.07,5,6.81,3.46,6.81,3.46,-40,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+18.92,y+5,z-29.07,3,7,3.7,7,0,-40,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	// Bottom
	frustum(x-34.9,y,z+22.75,5,6.81,3.46,6.81,3.46,50,"chopped",texture[17],texture[17],texture[19],texture[18],1.0,1.0,1.0); frustum(x-34.9,y+5,z+22.75,3,7,3.7,7,0,50,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	frustum(x-32.83,y,z+24.53,5,4.25,2.02,4.25,2.02,50,"chopped",texture[8],texture[9],texture[12],texture[12],1.0,1.0,1.0); frustum(x-32.83,y+5,z+24.53,3,4.25,2.02,4.25,0,50,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-31.9,y,z+32.55,5,5.06,2.69,5.06,2.69,50,"chopped",texture[10],texture[10],texture[12],texture[10],1.0,1.0,1.0); frustum(x-31.9,y+5,z+32.55,3,7,3.7,7,0,50,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x-27.49,y,z+38.89,5,6.81,3.46,6.81,3.46,-35,"chopped",texture[8],texture[8],texture[12],texture[15],1.0,1.0,1.0); frustum(x-27.49,y+5,z+38.89,3,7,3.7,7,0,-35,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-22.13,y,z+37.91,5,3.69,3.96,3.69,3.96,-35,"chopped",texture[13],texture[12],texture[9],texture[13],1.0,1.0,1.0);
	frustum(x-22.28,y,z+29.80,5,6.81,3.46,6.81,3.46,-20,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x-22.28,y+5,z+29.80,3,7,3.7,7,0,-20,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-16.05,y,z+33.54,5,6.81,3.46,6.81,3.46,-20,"chopped",texture[9],texture[9],texture[15],texture[12],1.0,1.0,1.0); frustum(x-16.05,y+5,z+33.54,3,7,3.7,7,0,-20,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-12.13,y,z+40.19,5,6.81,3.46,6.81,3.46,-65,"chopped",texture[9],texture[9],texture[12],texture[12],1.0,1.0,1.0); frustum(x-12.13,y+5,z+40.19,3,7,3.7,7,0,-65,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-5.34,y,z+51.10,5,6.81,3.46,6.81,3.46,40,"chopped",texture[17],texture[17],texture[19],texture[18],1.0,1.0,1.0); frustum(x-5.34,y+5,z+51.10,3,7,3.7,7,0,40,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	frustum(x-4.37,y,z+43.15,5,6.81,3.46,6.81,3.46,-25,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x-4.37,y+5,z+43.15,3,7,3.7,7,0,-25,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	frustum(x-0.51,y,z+32.58,5,6.17,3.14,6.17,3.14,0,"chopped",texture[17],texture[19],texture[18],texture[19],1.0,1.0,1.0); frustum(x-0.51,y+5,z+32.58,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);//a
	frustum(x+3.83,y,z+47.06,5,6.81,3.46,6.81,3.46,-95,"chopped",texture[8],texture[9],texture[12],texture[12],1.0,1.0,1.0); frustum(x+3.83,y+5,z+47.06,3,7,3.7,7,0,-95,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+5.51,y,z+36.20,5,6.81,3.46,6.81,3.46,15,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x+5.51,y+5,z+36.20,3,7,3.7,7,0,15,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x+12.05,y,z+46.06,5,6.81,3.46,6.81,3.46,-75,"chopped",texture[9],texture[9],texture[15],texture[12],1.0,1.0,1.0); frustum(x+12.05,y+5,z+46.06,3,7,3.7,7,0,-75,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x+19.16,y,z+42.00,5,6.81,3.46,6.81,3.46,-75,"chopped",texture[8],texture[8],texture[12],texture[15],1.0,1.0,1.0); frustum(x+19.16,y+5,z+42.00,3,7,3.7,7,0,-75,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x+9.24,y,z+31.51,3,3.96,2.02,3.96,2.02,-15,"chopped",texture[10],texture[10],texture[12],texture[10],1.0,1.0,1.0); frustum(x+9.24,y+3,z+31.51,2,3.96,2.02,3.96,0,-15,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);


	// More houses near look-out point
	frustum(x-41.27,y,z+18,5,6.81,3.46,6.81,3.46,30,"chopped",texture[9],texture[9],texture[15],texture[12],1.0,1.0,1.0); frustum(x-41.27,y+5,z+18,3,7,3.7,7,0,30,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-45.62,y,z+26.48,5,6.81,3.46,6.81,3.46,50,"chopped",texture[10],texture[10],texture[12],texture[10],1.0,1.0,1.0); frustum(x-45.62,y+5,z+26.48,3,7,3.7,7,0,50,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x-40.49,y,z+32.55,5,6.81,3.46,6.81,3.46,30,"chopped",texture[8],texture[8],texture[12],texture[15],1.0,1.0,1.0); frustum(x-40.49,y+5,z+32.55,3,7,3.7,7,0,30,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-29.38,y,z+45.89,5,6.81,3.46,6.81,3.46,-15,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x-29.38,y+5,z+45.89,3,7,3.7,7,0,-15,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
	frustum(x-21.71,y,z+49.83,5,6.81,3.46,6.81,3.46,60,"chopped",texture[11],texture[11],texture[12],texture[13],1.0,1.0,1.0); frustum(x-21.71,y+5,z+49.83,3,7,3.7,7,0,60,"chopped",texture[5],texture[5],texture[13],texture[13],1,1,1);
	frustum(x-12.13,y,z+48.72,3,3.96,2.02,3.96,2.02,-10,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x-12.13,y+3,z+48.72,2,3.96,2.02,3.96,0,-10,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);

   frustum(x+22.94,y,z+32.84,3,3.96,2.02,3.96,2.02,-10,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+22.94,y+3,z+32.84,2,3.96,2.02,3.96,0,-10,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
   frustum(x+28.26,y,z+36.20,5,3.96,2.02,3.96,2.02,170,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+28.26,y+5,z+36.20,2,3.96,2.02,3.96,0,170,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);

   frustum(x+27.64,y,z+24.92,5,6.81,3.46,6.81,3.46,-20,"chopped",texture[10],texture[11],texture[12],texture[12],1.0,1.0,1.0); frustum(x+27.64,y+5,z+24.92,3,7,3.7,7,0,-20,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
   frustum(x+32.05,y,z+20.69,5,6.81,3.46,6.81,3.46,-35,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+32.05,y+5,z+20.69,3,7,3.7,7,0,-35,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
   frustum(x+36.27,y,z+15.53,5,6.81,3.46,6.81,3.46,-15,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+36.27,y+5,z+15.53,3,7,3.7,7,0,-15,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
   frustum(x+41.83,y,z+22.00,5,6.81,3.46,6.81,3.46,0,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+41.83,y+5,z+22.00,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);
   frustum(x+45.21,y,z+16.32,5,6.81,3.46,6.81,3.46,0,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+45.21,y+5,z+16.32,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);

   frustum(x+37.61,y,z+7.01,5,6.81,3.46,6.81,3.46,0,"chopped",texture[17],texture[19],texture[19],texture[18],1.0,1.0,1.0); frustum(x+37.61,y+5,z+7.01,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[20],texture[20],1,1,1);
   frustum(x+45.83,y,z+7.01,5,6.81,3.46,6.81,3.46,0,"chopped",texture[9],texture[13],texture[13],texture[13],1.0,1.0,1.0); frustum(x+45.83,y+5,z+7.01,3,7,3.7,7,0,0,"chopped",texture[5],texture[5],texture[14],texture[14],1,1,1);


   // Filler Houses - Inputted from Left to Right
   ccHouse_1(x,y,z,-44.99,-30.81, 5, 90);
   ccHouse_3(x,y,z,-44.63,-38.98, 5, 90);
   ccHouse_2(x,y,z,-39.15,-25.78, 6, 15);
   ccHouse_1(x,y,z,-36.77,-31.20, 5, 15);
   ccHouse_2(x,y,z,-37.68,-37.91, 6, 40);
   ccHouse_1(x,y,z,-37.90,-43.94, 4, 40);
   ccHouse_1(x,y,z,-30.53,-48.95, 5, 40);
   ccHouse_3(x,y,z,-29.84,-41.79, 5, 65);
   ccHouse_1(x,y,z,-28.33,-33.12, 5, -30);
   ccHouse_1(x,y,z,-23.40,-41.79, 6, 40);
   ccHouse_2(x,y,z,-16.33,-40.23, 5, 40);
   ccHouse_1(x,y,z,-23.65,-55.28, 5, 0);
   ccHouse_3(x,y,z,-22.11,-49.99, 5, 0);
   ccHouse_2(x,y,z,-16.33,-47.38, 4, 40);
   ccHouse_1(x,y,z,-13.57,-53.91, 5, 0);
   ccHouse_1(x,y,z,-8.97,-48.17, 5, -70);
   ccHouse_3(x,y,z,-7.13,-41.79, 5, 0);
   ccHouse_2(x,y,z,-2.56,-48.20, 6, 90);
   ccHouse_1(x,y,z,-2.35,-55.78, 5, 0);
   ccHouse_3(x,y,z,1.80,-41.79, 6, 0);
   ccHouse_1(x,y,z,4.19,-51.60, 5, 0);
   ccHouse_2(x,y,z,7.90,-45.43, 4, -45);
   ccHouse_2(x,y,z,16.28,-41.56, 4, -30);
   ccHouse_1(x,y,z,22.63,-33.62, 5, -30);
   ccHouse_3(x,y,z,14.96,-47.68, 5, 0);
   ccHouse_1(x,y,z,14.30,-55.76, 5, 0);
   ccHouse_3(x,y,z,23.31,-43.52, 6, -35);
   ccHouse_2(x,y,z,27.69,-37.22, 5, -35);
   ccHouse_1(x,y,z,25.22,-23.37, 5, -40);
   ccHouse_1(x,y,z,28.93,-27.85, 4, -65);
   ccHouse_2(x,y,z,34.79,-30.88, 6, 30);
   ccHouse_3(x,y,z,41.80,-34.67, 5, 15);

   // Back Row
   ccHouse_1(x,y,z,22.86,-54.66, 5, 0);
   ccHouse_2(x,y,z,30.86,-56.32, 5, 0);
   ccHouse_2(x,y,z,30.86,-54.57, 4, 0);
   ccHouse_1(x,y,z,23.23,-49.10, 5, 0);
   ccHouse_3(x,y,z,31.64,-50.37, 5, 0);
   ccHouse_1(x,y,z,41.71,-50.36, 6, 0);
   ccHouse_3(x,y,z,31.96,-43.84, 5, 0);
   ccHouse_2(x,y,z,36.04,-38.83, 5, 0);
   ccHouse_3(x,y,z,40.31,-45.35, 4, 0);
   ccHouse_1(x,y,z,45.32,-40.56, 5, 0);
   ccHouse_3(x,y,z,45.32,-40.56, 4, 0);
   ccHouse_2(x,y,z,49.12,-45.77, 5, 0);
   ccHouse_2(x,y,z,50.18,-36.08, 6, 180);

   //Very Back Row
   ccHouse_1(x,y,z,-17.22,-61.6, 5, 0);
   ccHouse_1(x,y,z,-12.45,-66.68, 5, 180);
   ccHouse_2(x,y,z,-6.26,-62.07, 5, 0);
   ccHouse_3(x,y,z,-0.99,-67.16, 6, 0);
   ccHouse_1(x,y,z,4.5,-61.98, 5, 0);
   ccHouse_2(x,y,z,10.15,-67.08, 5, 0);
   ccHouse_2(x,y,z,15.93,-62.07, 4, 180);
   ccHouse_1(x,y,z,21.45,-67.16, 5, 0);
   ccHouse_3(x,y,z,26.97,-61.87, 4, 0);
   ccHouse_1(x,y,z,58.03,-46.34, 6, 0);

   // Right upper-middle houses
   ccHouse_1(x,y,z,37.62,-7.01, 5, 0);
   ccHouse_3(x,y,z,34.22,-12.52, 5, 30);
   ccHouse_1(x,y,z,42.43,-15.32, 5, 180);
   ccHouse_3(x,y,z,45.83,-7.01, 5, 0);
   ccHouse_2(x,y,z,48.61,-15.36, 5, 90);
   ccHouse_2(x,y,z,55.73,-7.54, 5, 0);
   ccHouse_1(x,y,z,56.52,-13.46, 5, 0);
   ccHouse_3(x,y,z,55.73,-19.11, 5, 30);
   ccHouse_1(x,y,z,63.02,-19.75, 5, 90);
   ccHouse_1(x,y,z,63.80,-10.79, 5, 270);
   ccHouse_2(x,y,z,68.60,-26.64, 5, 90);
   ccHouse_3(x,y,z,68.99,-17.42, 5, 90);
   ccHouse_1(x,y,z,68.99,-8.75, 5, 90);

   // Bottom-outer-right ring
   ccHouse_1(x,y,z,3.82,54.10, 5, 0);
   ccHouse_1(x,y,z,21.92,48.23, 5, 180);
   ccHouse_2(x,y,z,24.45,53.08, 5, 0);
   ccHouse_3(x,y,z,31.22,48.43, 5, 0);
   ccHouse_2(x,y,z,33.55,52.78, 5, 0);
   ccHouse_1(x,y,z,28.26,41.24, 5, -20);
   ccHouse_2(x,y,z,32.58,30.41, 5, 20);
   ccHouse_3(x,y,z,35.42,36.01, 5, 80);
   ccHouse_1(x,y,z,40.49,41.79, 5, 80);
   ccHouse_2(x,y,z,45.83,41.79, 5, 90);
   ccHouse_1(x,y,z,43.81,34.54, 5, 20);
   ccHouse_3(x,y,z,41.48,27.97, 5, 80);
   ccHouse_1(x,y,z,50.16,36.09, 5, 90);
   ccHouse_2(x,y,z,49.79,25.37, 5, 90);
   ccHouse_2(x,y,z,55.79,25.02, 5, 80);
   ccHouse_3(x,y,z,56.61,34.86, 5, 80);
   ccHouse_1(x,y,z,62.39,26.54, 5, 80);
   ccHouse_2(x,y,z,41.59,11.38, 5, 0);
   ccHouse_3(x,y,z,52.49,11.02, 5, 0);
   ccHouse_1(x,y,z,62.18,10.46, 5, 0);
   ccHouse_2(x,y,z,65.86,14.63, 5, 0);
   ccHouse_1(x,y,z,55.11,16.26, 5, -20);

   ccHouse_1(x,y,z,-43.26,39.99, 5, 90);
   ccHouse_2(x,y,z,-36.15,38.47, 5, 0);
   ccHouse_3(x,y,z,-36.99,43.49, 5, 0);



}

/*
 *  Call back function for glutDisplayFunc
 */
void display()
{
   //  Length of axes
   const double len=25.0;
   //  Over head Eye position
   double Ex = -2*dim[sceneNum]*Sin(th)*Cos(ph);
   double Ey = +2*dim[sceneNum]        *Sin(ph);
   double Ez = +2*dim[sceneNum]*Cos(th)*Cos(ph);
   //  Erase the window and the depth buffer
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   //  Enable Z-buffering in OpenGL
   glEnable(GL_DEPTH_TEST);
   //  Set perspective
   glLoadIdentity();
   if(view_mode==0){
      gluLookAt(Ex,Ey,Ez , 0,0,0 , 0,Cos(ph),0);
   } else if(view_mode==1){
      firstLoc_x[sceneNum] = Cos(th_1[sceneNum])+ex_1[sceneNum];
      firstLoc_z[sceneNum] = Sin(th_1[sceneNum])+ez_1[sceneNum];
      gluLookAt(ex_1[sceneNum],ey_1[sceneNum],ez_1[sceneNum], firstLoc_x[sceneNum],ey_1[sceneNum], firstLoc_z[sceneNum] , 0,1,0); // 2nd 3 -> where youre looking
   } else if(view_mode==2){
      gluLookAt(appa_scene_x[sceneNum]+appa_scene_r[sceneNum]*Cos(zh), appa_scene_y[sceneNum]+appa_scene_h[sceneNum], appa_scene_z[sceneNum]+appa_scene_r[sceneNum]*Sin(zh),appa_scene_x[sceneNum],appa_scene_y[sceneNum],appa_scene_z[sceneNum],0,1,0);
   }
   //  smooth shading
   glShadeModel(GL_SMOOTH);
   //  Translate intensity to color vectors
   float Ambient[]   = {0.01*ambient[sceneNum],0.01*ambient[sceneNum],0.01*ambient[sceneNum],1.0};
   float Diffuse[]   = {0.01*diffuse[sceneNum],0.01*diffuse[sceneNum],0.01*diffuse[sceneNum],1.0};
   float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
   //  Light postion
   float Position[]  = {light_x[sceneNum],light_y[sceneNum],light_z[sceneNum],1};
   //  Draw light position as ball
   glColor3f(1,1,1);
   if(show_ball == 1){
      ball(Position[0],Position[1],Position[2] , 1);
   }
   // normalize all normal vectors
   glEnable(GL_NORMALIZE);
   //  Enable lighting
   glEnable(GL_LIGHTING);
   //  glColor sets ambient and diffuse color materials
   glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);
   //  Enable light 0
   glEnable(GL_LIGHT0);      //  Set ambient, diffuse, specular components and position of light 0
   glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
   glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
   glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
   glLightfv(GL_LIGHT0,GL_POSITION,Position);

   // appa in every scene flying around
   appa(appa_scene_x[sceneNum]+appa_scene_r[sceneNum]*Cos(zh),appa_scene_y[sceneNum]+appa_scene_h[sceneNum],appa_scene_z[sceneNum]+appa_scene_r[sceneNum]*Sin(zh),appa_scene_s[sceneNum],90,-1*zh-90,-90);
   
   //  Draw scenes
   if(sceneNum == 0){
   	  chinCity(0,0,0);
        glDisable(GL_LIGHTING);
        Sky(512,sky[2],sky[3]);
        glEnable(GL_LIGHTING);
   }else if(sceneNum == 1){
      fireNationCastle(0,-30,0);
      glDisable(GL_LIGHTING);
      Sky(512,sky[4],sky[5]);
      glEnable(GL_LIGHTING);
   }else if(sceneNum == 2){
      easternAirTemple(0,0,0);
      glDisable(GL_LIGHTING);
      Sky(750,sky[0],sky[1]);
      glEnable(GL_LIGHTING);
   }
   
   // no lighting on axis
   glDisable(GL_LIGHTING);
   // axes at origin
   glColor3f(1,1,1);
   if (axes)
   {
      glBegin(GL_LINES);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(len,0.0,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,len,0.0);
      glVertex3d(0.0,0.0,0.0);
      glVertex3d(0.0,0.0,len);
      glEnd();
      //  Label axes
      glRasterPos3d(len,0.0,0.0);
      Print("X");
      glRasterPos3d(0.0,len,0.0);
      Print("Y");
      glRasterPos3d(0.0,0.0,len);
      Print("Z");
   }
   //  Disable Z-buffering to print to sceen
   glDisable(GL_DEPTH_TEST);
   //  Display scene name
   glWindowPos2i(5,5);
   Print("%s", cityNames[sceneNum]+2);
   
   //  Draw Map
   if (mainMode == 0)
      map();
   //  Render the scene and make it visible
   ErrCheck("display");
   glFlush();
   glutSwapBuffers();
}

/*
 *  callback function when no events are occuring
 */
void idle()
{
   //  Elapsed time in seconds
   double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
   zh = fmod(15*t,360.0);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  callback function when an arrow key is pressed
 */
void special(int key,int x,int y)
{
   if(mainMode == 1) // Key Binds When viewing a scene
   {
      if(view_mode==0){
         //  Right arrow key - increase angle by 5 degrees
         if (key == GLUT_KEY_RIGHT )
            th += 5;
         //  Left arrow key - decrease angle by 5 degrees
         else if (key == GLUT_KEY_LEFT)
            th -= 5;
         //  Up arrow key - increase elevation by 5 degrees
         else if (key == GLUT_KEY_UP)
            ph += 5;
         //  Down arrow key - decrease elevation by 5 degrees
         else if (key == GLUT_KEY_DOWN)
            ph -= 5;
         //  Keep angles to +/-360 degrees
         th %= 360;
         ph %= 360;
      } else if (view_mode==1) {
         //  Right arrow key - increase angle by 5 degrees
         if (key == GLUT_KEY_RIGHT)
         {
            th_1[sceneNum] += 5;
         }
         //  Left arrow key - decrease angle by 5 degrees
         else if (key == GLUT_KEY_LEFT) 
         {
            th_1[sceneNum] -= 5;
         }
         //  Up arrow key - increase elevation by 5 degrees
         else if (key == GLUT_KEY_UP)
         { 
            walkStep_x = walk[sceneNum]*Cos(th_1[sceneNum]);
            walkStep_z = walk[sceneNum]*Sin(th_1[sceneNum]);

            // Calculate distance between player and scene origin
            distance[sceneNum] = sqrt((ex_1_0[sceneNum]-(ex_1[sceneNum]+walkStep_x))*(ex_1_0[sceneNum]-(ex_1[sceneNum]+walkStep_x))+(ez_1_0[sceneNum]-(ez_1[sceneNum]+walkStep_z))*(ez_1_0[sceneNum]-(ez_1[sceneNum]+walkStep_z)));

            // Check to see if player is out of walking range
            if(distance[sceneNum]<walkingRadius[sceneNum]){
               ex_1[sceneNum] += walkStep_x;
               ez_1[sceneNum] += walkStep_z;
            }
         }
         //  Down arrow key - decrease elevation by 5 degrees
         else if (key == GLUT_KEY_DOWN)
         {
            walkStep_x = walk[sceneNum]*Cos(th_1[sceneNum]);
            walkStep_z = walk[sceneNum]*Sin(th_1[sceneNum]);

            // Calculate distance between player and scene origin
            distance[sceneNum] = sqrt((ex_1_0[sceneNum]-(ex_1[sceneNum]-walkStep_x))*(ex_1_0[sceneNum]-(ex_1[sceneNum]-walkStep_x))+(ez_1_0[sceneNum]-(ez_1[sceneNum]-walkStep_z))*(ez_1_0[sceneNum]-(ez_1[sceneNum]-walkStep_z)));
            // Check to see if player is out of walking range
            if(distance[sceneNum]<walkingRadius[sceneNum]){
               ex_1[sceneNum] -= walkStep_x;
               ez_1[sceneNum] -= walkStep_z;
            }
         }
      } else if(view_mode==2){
         //  Right arrow key - increase angle by 5 degrees
         if (key == GLUT_KEY_RIGHT)
            th_ap += 5;
         //  Left arrow key - decrease angle by 5 degrees
         else if (key == GLUT_KEY_LEFT)
            th_ap -= 5;
         //  Up arrow key - increase elevation by 5 degrees
         else if (key == GLUT_KEY_UP)
            ph_ap += 5;
         //  Down arrow key - decrease elevation by 5 degrees
         else if (key == GLUT_KEY_DOWN)
            ph_ap -= 5;
         //  Keep angles to +/-360 degrees
         th_ap %= 360;
         ph_ap %= 360;
      }
   }
   else if(mainMode == 0) // Key Binds when viewing the map
   {
   }
   //  Update projection
   Project(45,asp,dim[sceneNum]);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  callback function when a normal keyboard key is pressed
 */
void key(unsigned char ch,int x,int y)
{
   if(mainMode == 1) // Key Binds When viewing a scene
   {
      if (ch == 27)
         exit(0);
      //  Reset view angle
      else if (ch == '0')
         th = ph = 0;
      // open and close map
      else if (ch == 'm'){
         mainMode = 0;
      }
      //  Toggle axes
      else if (ch == 'x' || ch == 'X'){
         axes = 1-axes;
      }
      if(ch=='1'){ // overhead perspective view 
         view_mode = 0;
      } else if(ch=='2'){ // first person view 
         view_mode = 1;
      } else if (ch=='3'){ // appa view
         view_mode = 2;
      }
      else if(ch == 'b'){ // shows and hides the light source ball
         show_ball = 1-show_ball;
      }
   }
   else if(mainMode == 0) // Key Binds when viewing the map
   {
      //  Exit on ESC
      if (ch == 27)
         exit(0);
      else if (ch == 'm'){
         mainMode = 1;
      }
   }
   
   //  Reproject
   Project(45,asp,dim[sceneNum]);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  callback function when a mouse button is clicked
 */
void mouseClick(int button, int state,int x, int y){
   if(mainMode == 1) // Key Binds When viewing a scene
   {
      
   }
   else if(mainMode == 0) // Key Binds when viewing the map
   {
      if(button == GLUT_LEFT_BUTTON && state == GLUT_UP){
         for(int i =0;i<sceneNumTotal;i++){
            if(cityHovered[i] == 1){
               sceneNum = i;
               mainMode = 1;
            }
         }
      }
   }
   //  Reproject
   Project(45,asp,dim[sceneNum]);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}
/*
 *  callback function when the mouse moves with no buttons held down
 */
void mouseMove(int x, int y){
   y = globalHeight - y;
   if(mainMode == 1) // Key Binds When viewing a scene
   {
      
   }
   else if(mainMode == 0) // Key Binds when viewing the map
   {
      for(int i =0;i<sceneNumTotal;i++){
         if((globalWidth*cityPercentX[i] <= x && globalWidth*cityPercentX[i]+cityNamesLen[i] >= x) && (globalHeight*cityPercentY[i] <= y && globalHeight*cityPercentY[i]+14 >= y)){
            cityHovered[i] = 1;
         }else{
            cityHovered[i] = 0;
         }
      }
   }
   //  Reproject
   Project(45,asp,dim[sceneNum]);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  callback function when the window is resized
 */
void reshape(int width,int height)
{
   globalHeight = height;
   globalWidth = width;
   //  Ratio of the width to the height of the window
   asp = (height>0) ? (double)width/height : 1;
   //  Set the viewport to the entire window
   glViewport(0,0, width,height);
   //  Set projection
   Project(45,asp,dim[sceneNum]);
}

/*
 *  Start up GLUT and loads in images and .objs and loads up terrain arrays
 */
int main(int argc,char* argv[])
{
   //  Initialize GLUT
   glutInit(&argc,argv);

   //  Request double buffered, true color window with Z buffering at 600x600
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(600,600);
   glutCreateWindow("Final Project: Colin Craighead & Steven Priddy");
   //  Set callbacks
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   glutPassiveMotionFunc(mouseMove);
   glutMouseFunc(mouseClick);
   glutIdleFunc(idle);
   //  Load textures
   sky[0] = LoadTexBMP("eatskybox1.bmp");
   sky[1] = LoadTexBMP("eatskybox2.bmp");
   sky[2] = LoadTexBMP("chincityskybox1.bmp");
   sky[3] = LoadTexBMP("chincityskybox2.bmp");
   sky[4] = LoadTexBMP("firecastleskybox1.bmp");
   sky[5] = LoadTexBMP("firecastleskybox2.bmp");
   texture[0] = LoadTexBMP("roof.bmp");
   texture[1] = LoadTexBMP("test2.bmp"); // stone wall
   texture[2] = LoadTexBMP("limestone2.bmp"); // white stone
   texture[3] = LoadTexBMP("redMetal4.bmp");
   texture[4] = LoadTexBMP("map.bmp"); // map image
   texture[5] = LoadTexBMP("wood2.bmp"); // roof wood
   texture[6] = LoadTexBMP("leaves.bmp"); // tree
   texture[7] = LoadTexBMP("bark.bmp"); // tree
   texture[8] = LoadTexBMP("threewindowedwall.bmp"); //v chin city house texs v
   texture[9] = LoadTexBMP("onewindowedwall.bmp");
   texture[10] = LoadTexBMP("stripedwall.bmp");
   texture[11] = LoadTexBMP("thickstripedwall.bmp");
   texture[12] = LoadTexBMP("doorwall.bmp");
   texture[13] = LoadTexBMP("plainwall.bmp");
   texture[14] = LoadTexBMP("roofsidewall.bmp");
   texture[15] = LoadTexBMP("doorwallR.bmp");
   texture[16] = LoadTexBMP("onewindowedwallR.bmp");
   texture[17] = LoadTexBMP("altonewindowedwall.bmp");
   texture[18] = LoadTexBMP("altdoorwall.bmp");
   texture[19] = LoadTexBMP("altplainwall.bmp");
   texture[20] = LoadTexBMP("altroofsidewall.bmp");
   texture[21] = LoadTexBMP("chincitytemplewall.bmp");
   texture[22] = LoadTexBMP("chincitytempledoor.bmp");
   texture[23] = LoadTexBMP("singlewindowhigh.bmp"); //v fire castle building texs v
   texture[24] = LoadTexBMP("singlewindowhighright.bmp");
   texture[25] = LoadTexBMP("mediumdoor.bmp");
   texture[26] = LoadTexBMP("doublewindowdoor.bmp");
   texture[27] = LoadTexBMP("threewindowedwallfc.bmp");
   texture[28] = LoadTexBMP("eatroof.bmp"); //v Eastern air temple texs v
   texture[29] = LoadTexBMP("gold.bmp");
   texture[31] = LoadTexBMP("eatdoor.bmp");
   texture[32] = LoadTexBMP("eatthreewindowedwall.bmp");
   texture[33] = LoadTexBMP("eatonewindowedwall.bmp");
   texture[34] = LoadTexBMP("eatdoorwindows.bmp");
   texture[36] = LoadTexBMP("dirt.bmp");
   texture[37] = LoadTexBMP("grass.bmp");
   texture[38] = LoadTexBMP("smoothstone.bmp");

   // Build Terr Arrays
   // Params: heightsMode, square size, 3 global arrays pre made,,, seed for random
   createTerrMats(1,11,TerrVerts10,TerrNorms10,TerrVertNorms10,100);
   createTerrMats(1,21,TerrVerts20,TerrNorms20,TerrVertNorms20,100);
   createTerrMats(2,11,TerrVerts10CC,TerrNorms10CC,TerrVertNorms10CC,100);
   createTerrMats(1,51,TerrVerts50,TerrNorms50,TerrVertNorms50,100);
   createTerrMats(3,11,TerrVertsDome,TerrNormsDome,TerrVertNormsDome,100);
   createTerrMats(4,11,TerrVertsPeak,TerrNormsPeak,TerrVertNormsPeak,100);
   createTerrMats(5,11,TerrVerts10CCUp,TerrNorms10CCUp,TerrVertNorms10CCUp,100);
   // Read a obj file with appa
   appa_list_num = LoadOBJ("appa_color.obj");
   //  Pass control to GLUT so it can interact with the user
   ErrCheck("init");
   glutMainLoop();
   return 0;
}
