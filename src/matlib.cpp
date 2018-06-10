/*

Dummy routines for matrix transformations.

These are for you to write!

*/


#include <stdio.h>
#include <math.h>
#include <fstream>
#include <Windows.h>
#include <ostream>
#include <sstream>
#include <string>

#include "osuGraphics.h"
#include "lines.h"

//-------------------------------------------------
float xnorm;
float ynorm;
float znorm;

class RGB { 
public: 
	float r_diffuse;
	float g_diffuse;
	float b_diffuse;
	float r_specular;
	float g_specular;
	float b_specular;
	double phong_exp;
};

class LightSrc {
public:
	int size;
	float x;
	float y;
	float z;
	LightSrc() {
		x = 0;
		y = 0;
		z = 0;
	}
};



class Illumination {
public: 
	LightSrc pointsources[4];
	LightSrc dirsources[4];
	LightSrc ambientsource;
	int ps_size;
	int ds_size;
	double ambient;
	RGB color;
	Illumination() {
		ps_size = 0;
		ds_size = 0;
	}
};


class osuVertex {
	double x;
	double y;
	double z;
	double xnormal;
	double ynormal;
	double znormal;
	osuVertex(double a, double b, double c) {
		x = a;
		y = b;
		z = c;
		xnormal = xnorm;
		ynormal = ynorm;
		znormal = znorm;
	}
};

class Matr4X4 {
public:
	float x1y1;
	float x1y2;
	float x1y3;
	float x1y4;
	float x2y1;
	float x2y2;
	float x2y3;
	float x2y4;
	float x3y1;
	float x3y2;
	float x3y3;
	float x3y4;
	float x4y1;
	float x4y2;
	float x4y3;
	float x4y4;
	int size;
	Matr4X4() {
		x1y1 = 1;
		x1y2 = 0;
		x1y3 = 0;
		x1y4 = 0;
		x2y1 = 0;
		x2y2 = 1;
		x2y3 = 0;
		x2y4 = 0;
		x3y1 = 0;
		x3y2 = 0;
		x3y3 = 1;
		x3y4 = 0;
		x4y1 = 0;
		x4y2 = 0;
		x4y3 = 0;
		x4y4 = 1;
		size = 0;
	};
};

class Matr4X1 {
public:
	double x;
	double y;
	double z;
	double w;
	double canonz;
	int size;
	Matr4X1() {
		x = 0;
		y = 0;
		z = 0;
		w = 1;
		canonz = 0;
	};
	Matr4X1(float a, float b, float c, float d);
};

Matr4X1::Matr4X1(float a, float b, float c, float d) {
	x = a;
	y = b;
	z = c;
	w = d;
}

class Pixel {
public: 	
	double z;
	double r;
	double g;
	double b;
};

int drawmode;
int xy = 0;
float nearclip;
float farclip;
Matr4X1 endpoints[4];

Matr4X4 *matrix = new Matr4X4[32];

Pixel zbuffer[360][360];

float r;
float g;
float b;

LightSrc source = LightSrc();

Matr4X4 orthmatrix;

Matr4X4 perspective;

Matr4X4 camera = Matr4X4();

Matr4X4 windowing_transform;

Illumination room;


Matr4X4 matrix_Multiply4x4(Matr4X4 left, Matr4X4 right) {
	Matr4X4 matrix1;
	matrix1.x1y1 = (left.x1y1 * right.x1y1) + (left.x2y1 * right.x1y2) + (left.x3y1 * right.x1y3) + (left.x4y1 * right.x1y4);
	matrix1.x2y1 = (left.x1y1 * right.x2y1) + (left.x2y1 * right.x2y2) + (left.x3y1 * right.x2y3) + (left.x4y1 * right.x2y4);
	matrix1.x3y1 = (left.x1y1 * right.x3y1) + (left.x2y1 * right.x3y2) + (left.x3y1 * right.x3y3) + (left.x4y1 * right.x3y4);
	matrix1.x4y1 = (left.x1y1 * right.x4y1) + (left.x2y1 * right.x4y2) + (left.x3y1 * right.x4y3) + (left.x4y1 * right.x4y4);

	matrix1.x1y2 = (left.x1y2 * right.x1y1) + (left.x2y2 * right.x1y2) + (left.x3y2 * right.x1y3) + (left.x4y2 * right.x1y4);
	matrix1.x2y2 = (left.x1y2 * right.x2y1) + (left.x2y2 * right.x2y2) + (left.x3y2 * right.x2y3) + (left.x4y2 * right.x2y4);
	matrix1.x3y2 = (left.x1y2 * right.x3y1) + (left.x2y2 * right.x3y2) + (left.x3y2 * right.x3y3) + (left.x4y2 * right.x3y4);
	matrix1.x4y2 = (left.x1y2 * right.x4y1) + (left.x2y2 * right.x4y2) + (left.x3y2 * right.x4y3) + (left.x4y2 * right.x4y4);


	matrix1.x1y3 = (left.x1y3 * right.x1y1) + (left.x2y3 * right.x1y2) + (left.x3y3 * right.x1y3) + (left.x4y3 * right.x1y4);
	matrix1.x2y3 = (left.x1y3 * right.x2y1) + (left.x2y3 * right.x2y2) + (left.x3y3 * right.x2y3) + (left.x4y3 * right.x2y4);
	matrix1.x3y3 = (left.x1y3 * right.x3y1) + (left.x2y3 * right.x3y2) + (left.x3y3 * right.x3y3) + (left.x4y3 * right.x3y4);
	matrix1.x4y3 = (left.x1y3 * right.x4y1) + (left.x2y3 * right.x4y2) + (left.x3y3 * right.x4y3) + (left.x4y3 * right.x4y4);

	matrix1.x1y4 = (left.x1y4 * right.x1y1) + (left.x2y4 * right.x1y2) + (left.x3y4 * right.x1y3) + (left.x4y4 * right.x1y4);
	matrix1.x2y4 = (left.x1y4 * right.x2y1) + (left.x2y4 * right.x2y2) + (left.x3y4 * right.x2y3) + (left.x4y4 * right.x2y4);
	matrix1.x3y4 = (left.x1y4 * right.x3y1) + (left.x2y4 * right.x3y2) + (left.x3y4 * right.x3y3) + (left.x4y4 * right.x3y4);
	matrix1.x4y4 = (left.x1y4 * right.x4y1) + (left.x2y4 * right.x4y2) + (left.x3y4 * right.x4y3) + (left.x4y4 * right.x4y4);

	return matrix1;
}

Matr4X1 matrix_Multiply4x1(Matr4X4 left, Matr4X1 right) {
	Matr4X1 matrix1;
	matrix1.x = (left.x1y1 * right.x) + (left.x2y1 * right.y) + (left.x3y1 * right.z) + (left.x4y1 * right.w);
	matrix1.y = (left.x1y2 * right.x) + (left.x2y2 * right.y) + (left.x3y2 * right.z) + (left.x4y2 * right.w);
	matrix1.z = (left.x1y3 * right.x) + (left.x2y3 * right.y) + (left.x3y3 * right.z) + (left.x4y3 * right.w);
	matrix1.w = (left.x1y4 * right.x) + (left.x2y4 * right.y) + (left.x3y4 * right.z) + (left.x4y4 * right.w);
	matrix1.canonz = right.canonz;
	return matrix1;
}

void osuWindowTransform() {
	int windowWidth, windowHeight;
	osuGetFramebufferSize(&windowWidth, &windowHeight);

	windowing_transform.x1y1 = windowWidth / 2.0f;
	windowing_transform.x1y2 = 0;
	windowing_transform.x1y3 = 0;
	windowing_transform.x1y4 = 0;
	windowing_transform.x2y1 = 0;
	windowing_transform.x2y2 = windowHeight / 2.0f;
	windowing_transform.x2y3 = 0;
	windowing_transform.x2y4 = 0;
	windowing_transform.x3y1 = 0;
	windowing_transform.x3y2 = 0;
	windowing_transform.x3y3 = 1;
	windowing_transform.x3y4 = 0;
	windowing_transform.x4y1 = (windowWidth / 2.0f) - 0.5f;
	windowing_transform.x4y2 = (windowHeight / 2.0f) - 0.5f;
	windowing_transform.x4y3 = 0;
	windowing_transform.x4y4 = 1;
}

void osuOrtho(double left, double right, double bottom, double top, double nearp, double farp)
{ 
		Matr4X4 matrix2 = Matr4X4();
		//nearp = -abs(nearp);
		matrix2.x1y1 = (2.0f / (right - left));
		matrix2.x1y2 = 0;
		matrix2.x1y3 = 0;
		matrix2.x1y4 = 0;
		matrix2.x2y1 = 0;
		matrix2.x2y2 = (2.0f / (top - bottom));
		matrix2.x2y3 = 0;
		matrix2.x2y4 = 0;
		matrix2.x3y1 = 0;
		matrix2.x3y2 = 0;
		matrix2.x3y3 = -(2.0f / (farp - nearp));
		matrix2.x4y1 = -(left + right) / (right - left);
		matrix2.x4y2 = -(top + bottom) / (top - bottom);
		matrix2.x4y3 = -(nearp + farp) / (farp - nearp);
		matrix2.x4y4 = 1;

		orthmatrix = matrix2;
		nearclip = nearp;
		farclip = farp;

		//Z-buffer init
		for (int i = 0; i < 360; i++) {
			for (int j = 0; j < 360; j++) {
				zbuffer[i][j].r = 0;
				zbuffer[i][j].g = 0;
				zbuffer[i][j].b = 0;
				zbuffer[i][j].z = -farp;
			}
		}

}

void osuPerspective(double fovy, double nearp, double farp) 
{  	
	double h = tanf(fovy * 3.14159 / 360.0f);
	double top = h * nearp;
	double bottom = -top;
	double right = top;
	double left = -top;
	osuOrtho(left, right, bottom, top, nearp, farp);

	float temp, temp2, temp3, temp4;
	temp = 2.0 * nearp;
	temp2 = right - left;
	temp3 = top - bottom;
	temp4 = farp - nearp;
	
	perspective.x1y1 = temp /temp2;
	perspective.x1y2 = 0;
	perspective.x1y3 = 0;
	perspective.x1y4 = 0;
	perspective.x2y1 = 0;
	perspective.x2y2 = temp/temp3;
	perspective.x2y3 = 0;
	perspective.x2y4 = 0;
	perspective.x3y1 = (right + left) / temp2;
	perspective.x3y2 = (top + bottom) /temp3;
	perspective.x3y3 = (-farp - nearp) / temp4;
	perspective.x3y4 = -1.0f;
	perspective.x4y1 = 0;
	perspective.x4y2 = 0;
	perspective.x4y3 = (-temp * farp) /temp4;
	perspective.x4y4 = 0;
	orthmatrix = perspective;
	nearclip = nearp;
	farclip = farp;
}




void osuBegin(int mode)
{
	osuWindowTransform();

	
	switch (mode) {
		case OSU_TRIANGLE:
			drawmode = 1;
			break;

		case OSU_POLYGON:
			drawmode = 2;
			break;
	}



}
void osuEnd()
{

}

double* cross(double ax, double ay, double az, double bx, double by, double bz) {
	double ret[3];
	ret[0] = (ay * bz) - (az * by);
	ret[1] = (az * bx) - (ax * bz);
	ret[2] = (ax * by) - (ay * bx);
	return ret;
}

double dot(double ax, double ay, double az, double bx, double by, double bz) {
	return ((ax*bx) + (ay*by) + (az*bz));
}

void osuColor3f(double red, double green, double blue)
{
	r = red;
	g = green;
	b = blue;
}

void osuVertex2f(double x, double y)
{
}
Matr4X1 homogenize(Matr4X1 m) {
	m.x = m.x / m.w;
	m.y = m.y / m.w;
	m.z = m.z / m.w;
	m.w = m.w / m.w;
	return m;
}


void triangularize(Matr4X1 vert1, Matr4X1 vert2, Matr4X1 vert3)
{
	float xmaxbound = max(vert1.x, max(vert2.x, vert3.x));
	float ymaxbound = max(vert1.y, max(vert2.y, vert3.y));
	float xminbound = min(vert1.x, min(vert2.x, vert3.x));
	float yminbound = min(vert1.y, min(vert2.y, vert3.y));
	if (xmaxbound > 360) { xmaxbound = 360; }
	if (ymaxbound > 360) { ymaxbound = 360; }
	if (xminbound < 0) { xminbound = 0; }
	if (yminbound < 0) { yminbound = 0; }

	//i=x, j=y
	int f;
	int g;
	for (double i = xminbound; i < xmaxbound; i++) {
		for (double j = yminbound; j < ymaxbound; j++) {
			f = i;
			g = j;
			double det = ((vert2.y - vert3.y) * (vert1.x - vert3.x)) + ((vert3.x - vert2.x) * (vert1.y - vert3.y));
			double bary1 = (((vert2.y - vert3.y) * (i - vert3.x)) + ((vert3.x - vert2.x) * (j - vert3.y))) / det;
			double bary2 = (((vert3.y - vert1.y) * (i - vert3.x)) + ((vert1.x - vert3.x) * (j - vert3.y))) / det;
			double bary3 = 1 - bary1 - bary2;
			if (bary1 >= 0 && bary2 >= 0 && bary3 >= 0) {
				double zzz = (vert1.canonz * bary1) + (vert2.canonz * bary2) + (vert3.canonz * bary3);
				if (zzz > zbuffer[f][g].z) {
					double norm[3];
					double angle;
					double angle2;

					double* temp;

					double edge1[3];

					double edge2[3];
					double eye[3];
					double reflect[3];


					edge1[0] = vert1.x - vert2.x;
					edge1[1] = vert1.y - vert2.y;
					edge1[2] = vert1.z - vert2.z;
					edge2[0] = vert2.x - vert3.x;
					edge2[1] = vert2.y - vert3.y;
					edge2[2] = vert2.z - vert3.z;
					eye[0] = -(vert1.x + vert2.x + vert3.x) / 3.0f;
					eye[1] = -(vert1.y + vert2.y + vert3.y) / 3.0f;
					eye[2] = -(vert1.z + vert2.z + vert3.z) / 3.0f;



					temp = cross(edge1[0], edge1[1], edge1[2], edge2[0], edge2[1], edge2[2]);
					norm[0] = *temp;
					temp++;
					norm[1] = *temp;
					temp++;
					norm[2] = *temp;

					double mag = sqrt((norm[0] * norm[0]) + (norm[1] * norm[1]) + (norm[2] * norm[2]));
					norm[0] = (norm[0] / mag);
					norm[1] = (norm[1] / mag);
					norm[2] = (norm[2] / mag);


					mag = sqrt((room.pointsources[0].x * room.pointsources[0].x) + (room.pointsources[0].y * room.pointsources[0].y) + (room.pointsources[0].z * room.pointsources[0].z));
					room.pointsources[0].x = (room.pointsources[0].x / mag);
					room.pointsources[0].y = (room.pointsources[0].y / mag);
					room.pointsources[0].z = (room.pointsources[0].z / mag);

					angle = dot(norm[0], norm[1], norm[2], room.pointsources[0].x, room.pointsources[0].y, room.pointsources[0].z);

					mag = sqrt((eye[0] * eye[0]) + (eye[1] * eye[1]) + (eye[2] * eye[2]));
					eye[0] = (eye[0] / mag);
					eye[1] = (eye[1] / mag);
					eye[2] = (eye[2] / mag);
					
					reflect[0] = (2.0f * norm[0] * angle) - room.pointsources[0].x;
					reflect[1] = (2.0f * norm[1] * angle) - room.pointsources[0].y;
					reflect[2] = (2.0f * norm[2] * angle) - room.pointsources[0].z;

					angle2 = dot(eye[0], eye[1], eye[2], reflect[0], reflect[1], reflect[2]);
					double val = pow(max(0, angle2), room.color.phong_exp);
					//angle = acos(angle);
					//norm holds the normal vector
					if (angle > 0) {
						r = (room.ambient * 256.0f) + (room.color.r_diffuse * 256.0f * angle) + (256.0f  * room.color.r_specular * val);
						g = (room.ambient * 256.0f) + (room.color.g_diffuse * 256.0f * angle) + (256.0f * room.color.g_specular * val);
						b = (room.ambient * 256.0f) + (room.color.b_diffuse * 256.0f * angle) + (256.0f * room.color.b_specular * val);
					}

					osuWritePixel(i, j, r, g, b);
					zbuffer[f][g].z = zzz;
				}
			}
		}
	}
}

void osuNormal(double x, double y, double z) {
	xnorm = x;
	ynorm = y;
	znorm = z;
}

void osuVertex3f(double x, double y, double z)
{
	Matr4X1 temp(x, y, z, 1);
	Matr4X1 temp2;
	Matr4X4 transform = matrix[matrix[0].size];


	temp = matrix_Multiply4x1(transform, temp);
	temp = homogenize(temp);


	temp = matrix_Multiply4x1(camera, temp);
	temp = homogenize(temp);
	temp.canonz = temp.z;  


	temp = matrix_Multiply4x1(orthmatrix, temp);
	temp = homogenize(temp);

	temp = matrix_Multiply4x1(windowing_transform, temp);
	temp = homogenize(temp);
	endpoints[xy] = temp;
	xy++;
	
	if (xy == 3 && drawmode == 1) {
		triangularize(endpoints[0], endpoints[1], endpoints[2]);
		//draw_line((endpoints[0].x), (endpoints[0].y), (endpoints[1].x), (endpoints[1].y));
		//draw_line((endpoints[1].x), (endpoints[1].y), (endpoints[2].x), (endpoints[2].y));
		//draw_line((endpoints[2].x), (endpoints[2].y), (endpoints[0].x), (endpoints[0].y));
		xy = 0;
	}
	else if (xy == 4 && drawmode == 2) {
		triangularize(endpoints[0], endpoints[1], endpoints[2]);
		triangularize(endpoints[0], endpoints[2], endpoints[3]);
		draw_line((endpoints[0].x), (endpoints[0].y), (endpoints[1].x), (endpoints[1].y));
		draw_line((endpoints[1].x), (endpoints[1].y), (endpoints[2].x), (endpoints[2].y));
		draw_line((endpoints[2].x), (endpoints[2].y), (endpoints[3].x), (endpoints[3].y));
		draw_line((endpoints[3].x), (endpoints[3].y), (endpoints[0].x), (endpoints[0].y));
		xy = 0;
	}

		
	}
	





void osuInitialize() 
{ 
}


void osuPushMatrix()
{
	matrix[0].size = matrix[0].size + 1;
	int size = matrix[0].size;
	matrix[size] = matrix[size - 1];
}

void osuPopMatrix()
{
	int size = matrix[0].size;
	matrix[size].x1y1 = 1;
	matrix[size].x1y2 = 0;
	matrix[size].x1y3 = 0;
	matrix[size].x1y4 = 0;
	matrix[size].x2y1 = 0;
	matrix[size].x2y2 = 1;
	matrix[size].x2y3 = 0;
	matrix[size].x2y4 = 0;
	matrix[size].x3y1 = 0;
	matrix[size].x3y2 = 0;
	matrix[size].x3y3 = 1;
	matrix[size].x3y4 = 0;
	matrix[size].x4y1 = 0;
	matrix[size].x4y2 = 0;
	matrix[size].x4y3 = 0;
	matrix[size].x4y4 = 1;
	matrix[0].size = matrix[0].size - 1;
}

void osuLoadIdentityMatrix()
{
}
void osuTranslate(double tx, double ty, double tz)
{
	int size = matrix[0].size;

	Matr4X4 matrix1;


	matrix1.x1y1 = 1;
	matrix1.x2y1 = 0;
	matrix1.x3y1 = 0;
	matrix1.x4y1 = tx;

	matrix1.x1y2 = 0;
	matrix1.x2y2 = 1;
	matrix1.x3y2 = 0;
	matrix1.x4y2 = ty;

	matrix1.x1y3 = 0;
	matrix1.x2y3 = 0;
	matrix1.x3y3 = 1;
	matrix1.x4y3 = tz;

	matrix1.x1y4 = 0;
	matrix1.x2y4 = 0;
	matrix1.x3y4 = 0;
	matrix1.x4y4 = 1;
	matrix[size] = matrix_Multiply4x4(matrix[size], matrix1);
}

void osuScale(double sx, double sy, double sz)
{
	int size = matrix[0].size;
	Matr4X4 matrix1;

	matrix1.x1y1 = sx;
	matrix1.x1y2 = 0;
	matrix1.x1y3 = 0;
	matrix1.x1y4 = 0;

	matrix1.x2y1 = 0;
	matrix1.x2y2 = sy;
	matrix1.x2y3 = 0;
	matrix1.x2y4 = 0;

	matrix1.x3y1 = 0;
	matrix1.x3y2 = 0;
	matrix1.x3y3 = sz;
	matrix1.x3y4 = 0;

	matrix1.x4y1 = 0;
	matrix1.x4y2 = 0;
	matrix1.x4y3 = 0;
	matrix1.x4y4 = 1;
	matrix[size] = matrix_Multiply4x4(matrix[size], matrix1);
}



void osuRotate(double angle, double ax, double ay, double az) {

	int size = matrix[0].size;
	Matr4X4 matrix1;
	Matr4X4 matrix2;
	Matr4X4 matrix3;

	if (ax != 0) {
		matrix1.x1y1 = 1;
		matrix1.x1y2 = 0;
		matrix1.x1y3 = 0;
		matrix1.x1y4 = 0;

		matrix1.x2y1 = 0;
		matrix1.x2y2 = cos(angle * 3.14159 / 180.0f);
		matrix1.x2y3 = sin(angle * 3.14159 / 180.0f) * ax;
		matrix1.x2y4 = 0;

		matrix1.x3y1 = 0;
		matrix1.x3y2 = -sin(angle * 3.14159 / 180.0f) * ax;
		matrix1.x3y3 = cos(angle * 3.14159 / 180.0f);
		matrix1.x3y4 = 0;

		matrix1.x4y1 = 0;
		matrix1.x4y2 = 0;
		matrix1.x4y3 = 0;
		matrix1.x4y4 = 1;
	}
	if (ay != 0) {
		matrix2.x1y1 = cos(angle * 3.14159 / 180.0f);
		matrix2.x1y2 = 0;
		matrix2.x1y3 = -sin(angle * 3.14159 / 180.0f) * ay;
		matrix2.x1y4 = 0;

		matrix2.x2y1 = 0;
		matrix2.x2y2 = 1;
		matrix2.x2y3 = 0;
		matrix2.x2y4 = 0;

		matrix2.x3y1 = sin(angle * 3.14159 / 180.0f) * ay;
		matrix2.x3y2 = 0;
		matrix2.x3y3 = cos(angle * 3.14159 / 180.0f);
		matrix2.x3y4 = 0;

		matrix2.x4y1 = 0;
		matrix2.x4y2 = 0;
		matrix2.x4y3 = 0;
		matrix2.x4y4 = 1;
	}
	matrix2 = matrix_Multiply4x4(matrix2, matrix1);
	if (az != 0) {
		matrix3.x1y1 = cos(angle * 3.14159 / 180.0f);
		matrix3.x1y2 = -sin(angle * 3.14159 / 180.0f) * az;
		matrix3.x1y3 = 0;
		matrix3.x1y4 = 0;

		matrix3.x2y1 = sin(angle * 3.14159 / 180.0f) * az;
		matrix3.x2y2 = cos(angle * 3.14159 / 180.0f);
		matrix3.x2y3 = 0;
		matrix3.x2y4 = 0;

		matrix3.x3y1 = 0;
		matrix3.x3y2 = 0;
		matrix3.x3y3 = 1;
		matrix3.x3y4 = 0;

		matrix3.x4y1 = 0;
		matrix3.x4y2 = 0;
		matrix3.x4y3 = 0;
		matrix3.x4y4 = 1;
	}
	matrix3 = matrix_Multiply4x4(matrix3, matrix2);


	matrix[size] = matrix_Multiply4x4(matrix[size], matrix3);

}


void osuLookat(double from[3], double at[3], double up[3])
{
	double forward[3]; // new zaxis
	double side[3]; //new xaxis
	double temp[3];
	double* gaze;

	forward[0] = from[0] - at[0];
	forward[1] = from[1] - at[1];
	forward[2] = from[2] - at[2];

	// w = -g/||g||

	double forwardmag = sqrt((forward[0] * forward[0]) + (forward[1] * forward[1]) + (forward[2] * forward[2]));
	forward[0] = (forward[0] / forwardmag);
	forward[1] = (forward[1] / forwardmag);
	forward[2] = (forward[2] / forwardmag);

	gaze = cross(forward[0], forward[1], forward[2], up[0], up[1], up[2]);
	//gaze = cross(up[0], up[1], up[2], forward[0], forward[1], forward[2]);
	side[0] = *gaze;
	gaze++;
	side[1] = *gaze;
	gaze++;
	side[2] = *gaze;

	//t x w

	double umag = sqrt((side[0] * side[0]) + (side[1] * side[1]) + (side[2] * side[2]));
	if (umag == 0) {
		umag = 1;
	}
	side[0] = (side[0] / umag);
	side[1] = (side[1] / umag);
	side[2] = (side[2] / umag);

	gaze = cross(side[0], side[1], side[2], forward[0], forward[1], forward[2]);
	up[0] = *gaze;
	gaze++;
	up[1] = *gaze;
	gaze++;
	up[2] = *gaze;

	Matr4X4 matrix2 = Matr4X4();
	Matr4X4 matrix3 = Matr4X4();

	matrix3.x1y1 = side[0];
	matrix3.x1y2 = side[1];
	matrix3.x1y3 = side[2];
	matrix3.x1y4 = 0;
	matrix3.x2y1 = up[0];
	matrix3.x2y2 = up[1];
	matrix3.x2y3 = up[2];
	matrix3.x2y4 = 0;
	matrix3.x3y1 = -forward[0];
	matrix3.x3y2 = -forward[1];
	matrix3.x3y3 = -forward[2];
	matrix3.x3y4 = 0;
	matrix3.x4y1 = from[0];
	matrix3.x4y2 = from[1];
	matrix3.x4y3 = from[2];
	matrix3.x4y4 = 1;

	camera = matrix3;
}




void osuDiffuse(double r, double g, double b) {

	room.color.r_diffuse = r;
	room.color.g_diffuse = g;
	room.color.b_diffuse = b;


}

void osuSpecular(double r, double g, double b, double s) {
	room.color.r_specular = r;
	room.color.g_specular = g;
	room.color.b_specular = b;
	room.color.phong_exp = s;
}

void osuPointLight(double* lpos, double i) {
	room.pointsources[0].x = *lpos;
	lpos++;
	room.pointsources[0].y = *lpos;
	lpos++;
	room.pointsources[0].z = *lpos;
	room.ps_size++;
}


void osuAmbientLight(double i) {
	room.ambient = i;
}

void osuDirectionalLight(double* pos, double i) {
	room.dirsources[0].x = *pos;
	pos++;
	room.dirsources[0].y = *pos;
	pos++;
	room.dirsources[0].z = *pos;
	room.ds_size++;
}

void osuClearZ() {}


/*
void osuVertex3f(double x, double y, double z)
{
	Matr4X1 temp(x, y, z, 1);
	Matr4X1 temp2;
	Matr4X4 transform = matrix[matrix[0].size];


	temp = matrix_Multiply4x1(transform, temp);
	temp = homogenize(temp);

	temp = matrix_Multiply4x1(perspective, temp);
	temp = homogenize(temp);

	temp = matrix_Multiply4x1(orthmatrix, temp);
	temp = homogenize(temp);

	temp = matrix_Multiply4x1(windowing_transform, temp);
	temp = homogenize(temp);
	if (xy == 0) {

		endpoints[0] = temp;
		xy++;
	}
	else if (xy == 1) {

		endpoints[1] = temp;
		if (1 == near_far_clip(-1.0f, 1.0f, &temp.x, &temp.y, &temp.z, &temp2.x, &temp2.y, &temp2.z)) {
			draw_line((endpoints[0].x), (endpoints[0].y), (endpoints[1].x), (endpoints[1].y));
		}
		else {
			//draw_line((endpoints[0].x), (endpoints[0].y), (endpoints[1].x), (endpoints[1].y));
		}
		xy = 0;
	}

}*/