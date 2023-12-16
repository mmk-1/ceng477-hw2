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
	void compute_model_transformation(std::vector<std::map<int, Vec3>> &meshes_transformed_vertices, const Matrix4 &matrix_camera, const Matrix4 &matrix_projection);
};

struct Line
{
	double x0, y0, x1, y1, z0, z1;
	Color c0, c1;
	Line(double x0, double y0, double x1, double y1, double z0, double z1, Color c0, Color c1)
	{
		this->x0 = x0;
		this->y0 = y0;
		this->x1 = x1;
		this->y1 = y1;
		this->z0 = z0;
		this->z1 = z1;
		this->c0 = c0;
		this->c1 = c1;
	}
};

Matrix4 calculate_rotation_transformation(const Rotation *rotation);
Matrix4 calculate_model_trans_matrix(const Mesh *mesh, const Scene *scene);
Matrix4 calculate_camera_trans_matrix(const Camera *camera);
Matrix4 calculate_projection_trans_matrix(const Camera *camera, bool type);
Matrix4 calculate_viewport_trans_matrix(Camera *camera);

bool is_backfaced(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2);
void print_matrix4(Matrix4 matrix);
bool is_visible(double den, double num, double &te, double &tl);
Line clip_line(Line &line, bool &visible, double x_min, double x_max, double y_min, double y_max);

#endif
