#include <stdio.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"
#include <Windows.h>
#include "ObjLoader.h"
#include "GL/glut.h"

#define PI_10 3.1415926535897
#define INF 1000000
#define DEPTH 3

class Vect {
public:
	void Normalize();
	double x, y, z;
	Vect(); //constructor

	Vect(double x,double y, double z); //constructor
};

Vect::Vect() { x = 0;	y = 0;	z = 0; }

class Color {
public:
	double r, g, b;
	Color(); //constructor
	Color(double r, double g, double b); //constructor
};
Color::Color() { r = 0;	g = 0;	b = 0; }


Vect::Vect(double xcoord, double ycoord, double zcoord) {
	x = xcoord;
	y = ycoord;
	z = zcoord;
}

Color::Color(double rcolor, double gcolor, double bcolor) {
	r = rcolor;
	g = gcolor;
	b = bcolor;
}

void Vect::Normalize() {
	double mag = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	if (mag != 0) {
		x = x / mag;
		y = y / mag;
		z = z / mag;
	}
}

//Pair of sphere-ray intersections of the form ((ray_orig + ray_dir * i1), (ray_orig + ray_dir * i1))

struct Ipair {
	float i1, i2;
	bool valid;
	Ipair();

	Ipair(float v1, float v2, bool valid);
};
Ipair::Ipair() { i1 = 0; i2 = 0; valid = false; }

Ipair::Ipair(float v1, float v2, bool val) {
	i1 = v1;
	i2 = v2;
	valid = val;
}

class Sphere {
public:
	Vect center;
	Color surface;  //color of the surface
	Color emission; //color of ...?
	double radius;	//sphere radius
	double transparency;  
	double reflection;
	Sphere(Vect c, Color s, Color e, double rad, double trans, double ref);
	Ipair Intersect(Vect o, Vect d);
};

Sphere::Sphere(Vect c, Color s, Color e, double rad, double trans, double ref) {
	center = c;
	surface = s;
	emission = e;
	radius = rad;
	transparency = trans;
	reflection = ref;
}

Vect Vect_Subtract(Vect v1, Vect v2) {
	Vect new_Vect((v1.x - v2.x), (v1.y - v2.y), (v1.z - v2.z));
	return new_Vect;
}

double dot(Vect a, Vect b) {
	return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

Ipair Sphere::Intersect(Vect ray_orig, Vect ray_dir) {
	Vect d = Vect_Subtract(center, ray_orig);
	float g = dot(d, ray_dir);
	if (g < 0) {
		return Ipair(0, 0, false);
	}
	float h = dot(d, d) - (g * g);
	if (h > (radius * radius)) {
		return Ipair(0, 0, false);
	}
	float i = sqrt((radius * radius) - h);
	h = g - i;
	i = g + i;
	return Ipair(h, i, true);
}

double rand01()
{
	return ((double)rand() / RAND_MAX);
}

bool trace(std::vector<Sphere> &s, Vect &ray_orig, Vect &ray_dir, Color& color) {
	float tnear = INF;
	Sphere* sphere = NULL;
	Ipair hitpoints;
	int i = 0;
	while (i < s.size()) {
		sphere = &s[i];
		hitpoints = sphere->Intersect(ray_orig, ray_dir);
		if (hitpoints.valid == true) {
			color = Color(0, 0, 255);
			return true;
		}
		i++;
	}
	return false;
}


int main(int argc , char **argv)
{
    int width = 360;
    int height = 360;
	int OSU_FLAT = 1;

	Vect *rays = new Vect[(width * height)];
  
    int num = atoi(argv[1]);


    /* create a framebuffer */
    osuBeginGraphics(width, height);
    osuFlush();
    osuInitialize();
 
	float invWidth = 1 / float(width), invHeight = 1 / float(height);
	float fov = 30, aspectratio = width / float(height);
	float angle = tan(PI_10 * 0.5 * fov / 180.);

	std::vector<Sphere> spheres;
	spheres.push_back(Sphere(Vect(0, 0, -5), Color(0, 0, 255), Color(), 1, .5, .25));
  
	Color pixel = Color();
	Vect raydir;
	Vect rayorig;

	for (unsigned y = 0; y < height; ++y) {
		for (unsigned x = 0; x < width; ++x) {
			float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle;
			float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
			Vect raydir = Vect(xx, yy, -1);
			Vect rayorig = Vect();
			raydir.Normalize();
			if (trace(spheres, rayorig, raydir, pixel)) {
				osuWritePixel(x, y, pixel.r, pixel.g, pixel.b);
			}
		}
	}

  switch (num) {
    case 1:
		break;
	default: 
		("No such test case\n");
  }

  osuFlush();
  printf ("Press 'escape' to exit.\n");
  osuWaitOnEscape();
  osuEndGraphics();

}