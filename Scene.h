#ifndef _SCENE_H_
#define _SCENE_H_
#include "Vec3.h"
#include "Vec4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"
#include "Matrix4.h"
#include <map>
#include <vector>

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color>> image;
	std::vector<std::vector<double>> depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera *camera);
	void compute_model_transformation_for_meshes(std::vector<std::map<int, Vec4>> &meshes_transformed_vertices);
	void compute_camera_transformation_for_meshes(std::vector<std::map<int, Vec4>> &meshes_transformed_vertices, const Camera *camera);
	void compute_projection_transformation_for_meshes(const Camera *camera, std::vector<std::map<int, Vec4>> &meshes_transformed_vertices);
};

struct Line
{
	double x0, y0, z0, x1, y1, z1;
	int c0, c1;
	Line(Vec3 v0, Vec3 v1, int c0, int c1)
	{
		x0 = v0.x;
		y0 = v0.y;
		z0 = v0.z;
		x1 = v1.x;
		y1 = v1.y;
		z1 = v1.z;
		this->c0 = c0;
		this->c1 = c1;
	};
};
Matrix4 calculate_viewport_trans_matrix(Camera *camera);

bool is_backfaced(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
void print_matrix4(Matrix4 matrix);

#endif
