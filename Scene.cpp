#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <map>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"
#include "utils.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
 *********************Our Implementation starts here*********************************
 */

/*
 * Find rotation matrix based on orthonormal basis method.
 */
Matrix4 calculate_rotation_transformation(const Rotation *rotation)
{
	// calculate unit vector u
	// find orthoronotmal basis
	// - find min component of u, set it to 0
	// - set v to be perpendicular to u (revers the other 2 while negating one of them)
	// M matrix put the orthonamal basis as rows
	// calculate rotation matrix
	// M-1 matrix put the orthonamal basis as columns
	// multiply M-1 * R * M

	// Rotation *rotation = this->rotations[mesh->transformationIds[i] - 1];
	Vec3 u = Vec3(rotation->ux, rotation->uy, rotation->uz);
	Vec3 v;
	double min_component = std::min(fabs(rotation->ux), fabs(rotation->uy));
	min_component = std::min(min_component, fabs(rotation->uz));

	// Set V vector
	if (min_component == fabs(rotation->ux))
		v = Vec3(0, -1 * rotation->uz, rotation->uy);
	else if (min_component == fabs(rotation->uy))
		v = Vec3(-1 * rotation->uz, 0, rotation->ux);
	else if (min_component == fabs(rotation->uz))
		v = Vec3(-1 * rotation->uy, rotation->ux, 0);

	v = normalizeVec3(v);
	Vec3 w = crossProductVec3(u, v);
	w = normalizeVec3(w);

	// M
	double m[4][4] = {{u.x, u.y, u.z, 0}, {v.x, v.y, v.z, 0}, {w.x, w.y, w.z, 0}, {0, 0, 0, 1}};
	Matrix4 m_matrix(m);
	// M^-1
	double m_inverse[4][4] = {{u.x, v.x, w.x, 0}, {u.y, v.y, w.y, 0}, {u.z, v.z, w.z, 0}, {0, 0, 0, 1}};
	Matrix4 m_matrix_inverse(m_inverse);
	// R
	double cosine_res = cos(rotation->angle * M_PI / 180);
	double sine_res = sin(rotation->angle * M_PI / 180);
	double r[4][4] = {{1, 0, 0, 0}, {0, cosine_res, (-1) * sine_res, 0}, {0, sine_res, cosine_res, 0}, {0, 0, 0, 1}};
	Matrix4 r_matrix(r);
	// M^-1 * R * M
	return multiplyMatrixWithMatrix(m_matrix_inverse, multiplyMatrixWithMatrix(r_matrix, m_matrix));
};

/*
This function calculates the model transformation matrix for a given mesh.
It takes the scene and the mesh as parameters.
It returns the model transformation matrix.
*/
Matrix4 calculate_model_trans_matrix(const Scene *scene, const Mesh *mesh)
{
	// Initalize result matrix to identity matrix
	// multiply instead of assignment to avoid copying
	// TESTED

	Matrix4 result = getIdentityMatrix();

	for (int i = 0; i < mesh->numberOfTransformations; i++)
	{
		// cout << "tranformation type: " << mesh->transformationTypes[i] << endl;
		// cout << "transformatin id: " << mesh->transformationIds[i] << endl;
		// cout << "result matrix before transformation" << endl;
		// print_matrix4(result);
		// Be careful that transormationIds are 1-indexed
		switch (mesh->transformationTypes[i])
		{
		case 'r':
		{
			Rotation *rotation = scene->rotations[mesh->transformationIds[i] - 1];
			Matrix4 rotation_matrix = calculate_rotation_transformation(rotation);
			// cout << "rotation matrix" << endl;
			// print_matrix4(rotation_matrix);
			result = multiplyMatrixWithMatrix(rotation_matrix, result);
		}
		break;
		case 't':
		{
			Translation *translation = scene->translations[mesh->transformationIds[i] - 1];
			double translation_matrix[4][4] = {{1, 0, 0, translation->tx}, {0, 1, 0, translation->ty}, {0, 0, 1, translation->tz}, {0, 0, 0, 1}};
			// cout << "translation matrix" << endl;
			// print_matrix4(translation_matrix);
			result = multiplyMatrixWithMatrix(translation_matrix, result);
		}
		break;
		case 's':
		{
			Scaling *scale = scene->scalings[mesh->transformationIds[i] - 1];
			double scale_matrix[4][4] = {{scale->sx, 0, 0, 0}, {0, scale->sy, 0, 0}, {0, 0, scale->sz, 0}, {0, 0, 0, 1}};
			// cout << "scale matrix" << endl;
			// print_matrix4(scale_matrix);
			result = multiplyMatrixWithMatrix(scale_matrix, result);
		}
		break;
		default:
		{
			cout << "Invalid transformation type" << endl;
			exit(1);
		}
		break;
		}
		// cout << "result matrix after transformation" << endl;
		// print_matrix4(result);
	}
	return result;
}

Matrix4 calculate_camera_trans_matrix(const Camera *camera)
{
	// Translation matrix
	const Vec3 &position = camera->position;
	double t[4][4] = {{1, 0, 0, -(position.x)}, {0, 1, 0, -(position.y)}, {0, 0, 1, -(position.z)}, {0, 0, 0, 1}};

	// Rotation matrix
	const Vec3 u = camera->u;
	const Vec3 v = camera->v;
	const Vec3 w = camera->w;
	double r[4][4] = {{u.x, u.y, u.z, 0}, {v.x, v.y, v.z, 0}, {w.x, w.y, w.z, 0}, {0, 0, 0, 1}};

	// Cast to Matrix4
	Matrix4 r_matrix(r);
	Matrix4 t_matrix(t);

	// cout << "Rotation Matrix" << endl;
	// print_matrix4(r_matrix);
	// cout << "Translation Matrix" << endl;
	// print_matrix4(t_matrix);

	return multiplyMatrixWithMatrix(r_matrix, t_matrix);
}

/*
This function calculates the projection transformation matrix for a given camera.
It takes the camera as parameter.
It returns the projection transformation matrix.
*/
Matrix4 calculate_projection_trans_matrix(const Camera *camera)
{
	Matrix4 result;
	double left = camera->left;
	double right = camera->right;
	double bottom = camera->bottom;
	double top = camera->top;
	double near = camera->near;
	double far = camera->far;
	int projectionType = camera->projectionType;

	if (projectionType == ORTOGRAPHIC_PROJECTION)
	{
		double matrix_orth[4][4] = {
			{2 / (right - left), 0, 0, -((right + left) / (right - left))},
			{0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
			{0, 0, -(2 / (far - near)), -((far + near) / (far - near))},
			{0, 0, 0, 1}};

		result = Matrix4(matrix_orth);
	}
	else if (projectionType == PERSPECTIVE_PROJECTION)
	{
		double matrix_perspective[4][4] = {
			{(2 * near) / (right - left), 0, (right + left) / (right - left), 0},
			{0, (2 * near) / (top - bottom), (top + bottom) / (top - bottom), 0},
			{0, 0, -((far + near) / (far - near)), -((2 * far * near) / (far - near))},
			{0, 0, -1, 0}};

		result = Matrix4(matrix_perspective);
	}
	else
	{
		cout << "Invalid projection type" << endl;
		exit(1);
	}
	return result;
}

Matrix4 calculate_orthographic_projection_matrix(const Camera *camera)
{
	double left = camera->left;
	double right = camera->right;
	double bottom = camera->bottom;
	double top = camera->top;
	double near = camera->near;
	double far = camera->far;
	double matrix_orth[4][4] = {
		{2 / (right - left), 0, 0, -((right + left) / (right - left))},
		{0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
		{0, 0, -(2 / (far - near)), -((far + near) / (far - near))},
		{0, 0, 0, 1}};

	return Matrix4(matrix_orth);
}

Matrix4 calculate_perspective_projection_matrix(const Camera *camera)
{
	double left = camera->left;
	double right = camera->right;
	double bottom = camera->bottom;
	double top = camera->top;
	double near = camera->near;
	double far = camera->far;
	double matrix_perspective[4][4] = {
		{(2 * near) / (right - left), 0, (right + left) / (right - left), 0},
		{0, (2 * near) / (top - bottom), (top + bottom) / (top - bottom), 0},
		{0, 0, -((far + near) / (far - near)), -((2 * far * near) / (far - near))},
		{0, 0, -1, 0}};

	return Matrix4(matrix_perspective);
}

Matrix4 calculate_viewport_trans_matrix(Camera *camera)
{
	double nx_div_2 = camera->horRes / 2.0;
	double ny_div_2 = camera->verRes / 2.0;
	double nx_div_2_minus_1 = (camera->horRes - 1) / 2.0;
	double ny_div_2_minus_1 = (camera->verRes - 1) / 2.0;
	double M_viewport[4][4] = {{nx_div_2, 0, 0, nx_div_2_minus_1},
							   {0, ny_div_2, 0, ny_div_2_minus_1},
							   {0, 0, 0.5, 0.5},
							   {0, 0, 0, 1}};
	return Matrix4(M_viewport);
}

/* Backface culling */
bool is_backfaced(const Vec3 &v0, const Vec3 &v1, const Vec3 &v2)
{
	Vec3 edge0 = subtractVec3(v1, v0);
	Vec3 edge1 = subtractVec3(v2, v0);
	Vec3 normalVector = normalizeVec3(crossProductVec3(edge0, edge1));
	return dotProductVec3(normalVector, v0) < 0;
}

bool is_visible(double den, double num, double &te, double &tl)
{
	double t = num / den;
	if (den > 0)
	{
		if (t > tl)
			return false;

		if (t > te)
			te = t;
	}
	else if (den < 0)
	{
		if (t < te)
			return false;

		if (t < tl)
			tl = t;
	}
	else if (num > 0)
	{
		return false;
	}
	return true;
}
/*



*/
std::vector<Vec4> clip_line(Scene *scene, std::vector<Vec4> &line, std::vector<Color> &line_colors, bool &visible)
{
	// Lian-Barsky Algorithm
	double te = 0;
	double tl = 1;
	double x_min = -1;
	double x_max = 1;
	double y_min = -1;
	double y_max = 1;
	double z_min = -1;
	double z_max = 1;
	double dx = line[1].x - line[0].x;
	double dy = line[1].y - line[0].y;
	double dz = line[1].z - line[0].z;
	visible = false;

	if (is_visible(dx, x_min - line[0].x, te, tl) &&
		is_visible(-dx, line[0].x - x_max, te, tl) &&
		is_visible(dy, y_min - line[0].y, te, tl) &&
		is_visible(-dy, line[0].y - y_max, te, tl) &&
		is_visible(dz, z_min - line[0].z, te, tl) &&
		is_visible(-dz, line[0].z - z_max, te, tl))
	{
		visible = true;
		if (tl < 1)
		{
			line[1].x = line[0].x + tl * dx;
			line[1].y = line[0].y + tl * dy;
			line[1].z = line[0].z + tl * dz;
			line_colors[1].r = line_colors[0].r + tl * (line_colors[1].r - line_colors[0].r);
			line_colors[1].g = line_colors[0].g + tl * (line_colors[1].g - line_colors[0].g);
			line_colors[1].b = line_colors[0].b + tl * (line_colors[1].b - line_colors[0].b);
		}
		if (te > 0)
		{
			line[0].x = line[0].x + te * dx;
			line[0].y = line[0].y + te * dy;
			line[0].z = line[0].z + te * dz;
			line_colors[0].r = line_colors[0].r + te * (line_colors[1].r - line_colors[0].r);
			line_colors[0].g = line_colors[0].g + te * (line_colors[1].g - line_colors[0].g);
			line_colors[0].b = line_colors[0].b + te * (line_colors[1].b - line_colors[0].b);
		}
	}

	return line;
}

/*
This fuction takes Vec4 as input and returns Vec4 after perspective division.
*/
Vec4 perspective_division(Vec4 &v)
{
	v.x = v.x / v.t;
	v.y = v.y / v.t;
	v.z = v.z / v.t;
	v.t = 1;
	return v;
}

void draw(Scene *scene, int x, int y, Color &color, double depth)
{
	if (depth < scene->depth[x][y])
	{
		scene->depth[x][y] = depth;
		scene->image[x][y].r = int(round(color.r));
		scene->image[x][y].g = int(round(color.g));
		scene->image[x][y].b = int(round(color.b));
	}
}

/*
This function takes 2 points representing a line and rasterizes it.
*/
void rasterize_line(Scene *scene, Vec3 v0, Vec3 v1, Color c0, Color c1)
{
	// Midpoint rasterization algorithm
	if (v0.x > v1.x)
	{
		Vec3 temp = v0;
		v0 = v1;
		v1 = temp;
		Color temp_color = c0;
		c0 = c1;
		c1 = temp_color;
	}
	int x0 = int(v0.x);
	int y0 = int(v0.y);
	int z0 = int(v0.z);
	int x1 = int(v1.x);
	int y1 = int(v1.y);
	int z1 = int(v1.z);
	double dx = (x1 - x0);
	double dy = (y1 - y0);
	double slope = dy / dx;

	Color c = c0;

	if (slope >= 0)
	{
		if (slope > 1)
		{
			int x = x0;
			double dz = (z1 - z0) / dy;
			double depth = z0;
			double dr = (c1.r - c0.r) / dy;
			double dg = (c1.g - c0.g) / dy;
			double db = (c1.b - c0.b) / dy;
			double d = (v0.x - v1.x) + .5 * (v1.y - v0.y);
			for (int y = y0; y <= y1; y++)
			{
				draw(scene, x, y, c, depth);
				if (d < 0)
				{
					x = x + 1;
					d += (x0 - x1) + (y1 - y0);
				}
				else
				{
					d += (x0 - x1);
				}
				c.r += dr;
				c.g += dg;
				c.b += db;
				depth += dz;
			}
		}
		else
		{
			int y = y0;
			double dz = (z1 - z0) / dx;
			double depth = z0;
			double dr = (c1.r - c0.r) / dx;
			double dg = (c1.g - c0.g) / dx;
			double db = (c1.b - c0.b) / dx;
			double d = (v0.y - v1.y) + .5 * (v1.x - v0.x);
			for (int x = x0; x <= x1; x++)
			{
				draw(scene, x, y, c, depth);
				if (d < 0)
				{
					y = y + 1;
					d += (y0 - y1) + (x1 - x0);
				}
				else
				{
					d += (y0 - y1);
				}
				c.r += dr;
				c.g += dg;
				c.b += db;
				depth += dz;
			}
		}
	}
	else
	{
		if (slope > -1)
		{
			int y = y0;
			double dz = (z1 - z0) / dx;
			double depth = z0;
			double dr = (c1.r - c0.r) / dx;
			double dg = (c1.g - c0.g) / dx;
			double db = (c1.b - c0.b) / dx;
			double d = (v0.y - v1.y) - .5 * (v1.x - v0.x);
			for (int x = x0; x <= x1; x++)
			{
				draw(scene, x, y, c, depth);
				if (d > 0)
				{
					y = y - 1;
					d += (y0 - y1) - (x1 - x0);
				}
				else
				{
					d += (y0 - y1);
				}
				c.r += dr;
				c.g += dg;
				c.b += db;
				depth += dz;
			}
		}
		else
		{
			int x = x0;
			double dz = (z1 - z0) / dy;
			double depth = z0;
			double dr = (c1.r - c0.r) / dy;
			double dg = (c1.g - c0.g) / dy;
			double db = (c1.b - c0.b) / dy;
			double d = 2 * (v0.x - v1.x) + (v1.y - v0.y);
			for (int y = y0; y >= y1; y--)
			{
				draw(scene, x, y, c, depth);
				if (d < 0)
				{
					x = x + 1;
					d += 2 * ((x0 - x1) + (y0 - y1));
				}
				else
				{
					d += 2 * (x0 - x1);
				}
				c.r -= dr;
				c.g -= dg;
				c.b -= db;
				depth += dz;
			}
		}
	}
}

void clip_and_rasterize_line(Scene *scene, Vec4 v0, Vec4 v1, Matrix4 &viewport_transformation_matrix)
{

	std::vector<Vec4> line = std::vector<Vec4>(2);
	std::vector<Color> line_colors = std::vector<Color>(2);
	bool is_line_visible = false;

	line[0] = v0;
	line[1] = v1;

	line_colors[0] = *(scene->colorsOfVertices[v0.colorId - 1]);
	line_colors[1] = *(scene->colorsOfVertices[v1.colorId - 1]);
	clip_line(scene, line, line_colors, is_line_visible);
	if (is_line_visible)
	{

		// 7. Viewport Transformation for a line
		line[0] = multiplyMatrixWithVec4(viewport_transformation_matrix, line[0]);
		line[1] = multiplyMatrixWithVec4(viewport_transformation_matrix, line[1]);

		Vec3 v0_viewport = Vec3(line[0].x, line[0].y, line[0].z);
		Vec3 v1_viewport = Vec3(line[1].x, line[1].y, line[1].z);

		// 8. Rasterization for a line
		rasterize_line(scene, v0_viewport, v1_viewport, line_colors[0], line_colors[1]);
	}
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// Overall steps for each triangle:
	// 1. Model Transformation [X] Test [X]
	// 2. Camera Transformation [X] Test [X]
	// 3. Projection Transformation [X] Test [X]
	// 4. Perspective Division [X] Test [X]
	// 5. Backface Culling [X] Test [X]
	// 6. Clipping [] Test []
	// 7. Viewport Transformation [X] Test []
	// 8. Rasterization [] Test []
	// 9. Depth Testing [] Test []

	Matrix4 camera_transformation_matrix = calculate_camera_trans_matrix(camera);
	// cout << "Camera Transformation matrix " << endl;
	// print_matrix4(camera_transformation_matrix);
	// cout << "Camera orientation" << endl;
	// print_vec4(Vec4(camera->u.x, camera->u.y, camera->u.z, 1, 0));
	// print_vec4(Vec4(camera->v.x, camera->v.y, camera->v.z, 1, 0));
	// print_vec4(Vec4(camera->w.x, camera->w.y, camera->w.z, 1, 0));
	// print_vec4(Vec4(camera->position.x, camera->position.y, camera->position.z, 1, 0));

	Matrix4 projection_transformation_matrix = calculate_projection_trans_matrix(camera);
	// cout << "Camera stuff" << endl;
	// cout << "Camera left " << camera->left << endl;
	// cout << "Camera right " << camera->right << endl;
	// cout << "Camera bottom " << camera->bottom << endl;
	// cout << "Camera top " << camera->top << endl;
	// cout << "Camera near " << camera->near << endl;
	// cout << "Camera far " << camera->far << endl;
	// cout << "Projection Transformation matrix " << endl;
	// print_matrix4(projection_transformation_matrix);

	Matrix4 viewport_transformation_matrix = calculate_viewport_trans_matrix(camera);
	// cout << "camera horizontal " << camera->horRes << endl;
	// cout << "camera vertical " << camera->verRes << endl;
	// cout << "Viewport Transformation matrix " << endl;
	// print_matrix4(viewport_transformation_matrix);

	for (int m = 0; m < this->meshes.size(); m++)
	{
		Mesh *mesh = this->meshes[m];

		Matrix4 mesh_model_transformation_matrix = calculate_model_trans_matrix(this, mesh);
		// cout << "Model Transformation matrix " << endl;
		// print_matrix4(mesh_model_transformation_matrix);

		// std::cout << "Model Transformation matrix " << std::endl;
		// print_matrix4(mesh_model_transformation_matrix);

		// std::cout << "Camera Transformation matrix " << std::endl;
		// print_matrix4(camera_transformation_matrix);

		// std::cout << "Projection Transformation matrix " << std::endl;
		// print_matrix4(projection_transformation_matrix);

		for (int t = 0; t < mesh->numberOfTriangles; t++)
		{
			Triangle &triangle = mesh->triangles[t];

			Vec3 *vertix0 = this->vertices[triangle.vertexIds[0] - 1];
			Vec3 *vertix1 = this->vertices[triangle.vertexIds[1] - 1];
			Vec3 *vertix2 = this->vertices[triangle.vertexIds[2] - 1];

			Vec4 v0 = Vec4(vertix0->x, vertix0->y, vertix0->z, 1, vertix0->colorId);
			Vec4 v1 = Vec4(vertix1->x, vertix1->y, vertix1->z, 1, vertix1->colorId);
			Vec4 v2 = Vec4(vertix2->x, vertix2->y, vertix2->z, 1, vertix2->colorId);

			// cout << "v0" << endl;
			// print_vec4(v0);
			// cout << "v1" << endl;
			// print_vec4(v1);
			// cout << "v2" << endl;
			// print_vec4(v2);

			// 1. Model Transformation
			v0 = multiplyMatrixWithVec4(mesh_model_transformation_matrix, v0);
			v1 = multiplyMatrixWithVec4(mesh_model_transformation_matrix, v1);
			v2 = multiplyMatrixWithVec4(mesh_model_transformation_matrix, v2);

			// cout << "v0 after multiply with model trans matrix" << endl;
			// print_vec4(v0);
			// cout << "v1 after multiply with model trans matrix" << endl;
			// print_vec4(v1);
			// cout << "v2 after multiply with model trans matrix" << endl;
			// print_vec4(v2);

			// 2. Camera Transformation
			v0 = multiplyMatrixWithVec4(camera_transformation_matrix, v0);
			v1 = multiplyMatrixWithVec4(camera_transformation_matrix, v1);
			v2 = multiplyMatrixWithVec4(camera_transformation_matrix, v2);

			// cout << "v0 after camera trans matrix" << endl;
			// print_vec4(v0);
			// cout << "v1 after camera trans matrix" << endl;
			// print_vec4(v1);
			// cout << "v2 after camera trans matrix" << endl;
			// print_vec4(v2);

			// 3. Projection Transformation
			v0 = multiplyMatrixWithVec4(projection_transformation_matrix, v0);
			v1 = multiplyMatrixWithVec4(projection_transformation_matrix, v1);
			v2 = multiplyMatrixWithVec4(projection_transformation_matrix, v2);

			// cout << "v0 after projection trans matrix" << endl;
			// print_vec4(v0);
			// cout << "v1 after projection trans matrix" << endl;
			// print_vec4(v1);
			// cout << "v2 after projection trans matrix" << endl;
			// print_vec4(v2);

			// 4. Perspective Division
			if (camera->projectionType == PERSPECTIVE_PROJECTION)
			{
				v0 = perspective_division(v0);
				v1 = perspective_division(v1);
				v2 = perspective_division(v2);

				// cout << "v0 after perspective divide" << endl;
				// print_vec4(v0);
				// cout << "v1 after perspective divide" << endl;
				// print_vec4(v1);
				// cout << "v2 after perspective divide" << endl;
				// print_vec4(v2);
			}

			// 5. Backface Culling
			// TODO: Check if backfaced
			if (this->cullingEnabled && is_backfaced(Vec3(v0.x, v0.y, v0.z), Vec3(v1.x, v1.y, v1.z), Vec3(v2.x, v2.y, v2.z)))
				continue;

			// 6. Clipping and rasterization for a line

			if (mesh->type == WIREFRAME_MESH)
			{
				clip_and_rasterize_line(this, v0, v1, viewport_transformation_matrix);
				clip_and_rasterize_line(this, v1, v2, viewport_transformation_matrix);
				clip_and_rasterize_line(this, v2, v0, viewport_transformation_matrix);
			}
			else
			{
				rasterize_triangle(this, camera, viewport_transformation_matrix, v0, v1, v2);
			}
		}
	}
}

/*
 *********************Our Implementation ends here*********************************
 */

double calculate_f(double x, double y, double x0, double y0, double x1, double y1)
{
	return (x * (y0 - y1)) + (y * (x1 - x0)) + (x0 * y1) - (y0 * x1);
}

void rasterize_triangle(Scene *scene, Camera *camera, Matrix4 &viewport_transformation_matrix, Vec4 v0, Vec4 v1, Vec4 v2)
{

	v0 = multiplyMatrixWithVec4(viewport_transformation_matrix, v0);
	v1 = multiplyMatrixWithVec4(viewport_transformation_matrix, v1);
	v2 = multiplyMatrixWithVec4(viewport_transformation_matrix, v2);

	Color *c0 = scene->colorsOfVertices[v0.colorId - 1];
	Color *c1 = scene->colorsOfVertices[v1.colorId - 1];
	Color *c2 = scene->colorsOfVertices[v2.colorId - 1];
	int n_x = camera->horRes;
	int n_y = camera->verRes;

	// x_min and x_max
	int x_min = min(min(v0.x, v1.x), v2.x);
	// Boundary checks
	if (x_min < 0)
	{
		x_min = 0;
	}
	else if (x_min > n_x - 1)
	{
		x_min = n_x - 1;
	}

	int x_max = max(max(v0.x, v1.x), v2.x);
	// Boundary checks
	if (x_max < 0)
	{
		x_max = 0;
	}
	else if (x_max > n_x - 1)
	{
		x_max = n_x - 1;
	}

	// y_min and y_max
	int y_min = min(min(v0.y, v1.y), v2.y);
	// Boundary checks
	if (y_min < 0)
	{
		y_min = 0;
	}
	else if (y_min > n_y - 1)
	{
		y_min = n_y - 1;
	}

	int y_max = max(max(v0.y, v1.y), v2.y);
	// Boundary checks
	if (y_max < 0)
	{
		y_max = 0;
	}
	else if (y_max > n_y - 1)
	{
		y_max = n_y - 1;
	}

	double alpha, beta, gamma;
	Color color;
	for (int y = y_min; y <= y_max; y++)
	{
		for (int x = x_min; x <= x_max; x++)
		{
			// f_01(x, y) / f_01(v2.x, v2.y)
			gamma = calculate_f(x, y, v0.x, v0.y, v1.x, v1.y) / calculate_f(v2.x, v2.y, v0.x, v0.y, v1.x, v1.y);
			// f_12(x, y) / f_12(v0.x, v0.y)
			alpha = calculate_f(x, y, v1.x, v1.y, v2.x, v2.y) / calculate_f(v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
			// f_20(x, y) / f_20(v1.x, v1.y)
			beta = calculate_f(x, y, v2.x, v2.y, v0.x, v0.y) / calculate_f(v1.x, v1.y, v2.x, v2.y, v0.x, v0.y);

			if (gamma >= 0 && alpha >= 0 && beta >= 0)
			{
				color.b = c0->b * alpha + c1->b * beta + c2->b * gamma;
				color.g = c0->g * alpha + c1->g * beta + c2->g * gamma;
				color.r = c0->r * alpha + c1->r * beta + c2->r * gamma;
				double depth = v0.z * alpha + v1.z * beta + v2.z * gamma;
				draw(scene, x, y, color, depth);
			}
		}
	}
}

/*
**********************Test functions starts here*********************************
*/

void print_matrix4(Matrix4 matrix)
{
	for (int i = 0; i < 4; i++)
	{
		cout << "[";
		for (int j = 0; j < 4; j++)
		{
			cout << matrix.values[i][j] << " ";
		}
		cout << "]" << endl;
	}
	std::cout << "---------------------" << std::endl;
}

void print_vec4(Vec4 vec)
{
	cout << "Vec4: " << vec.x << " " << vec.y << " " << vec.z << " " << vec.t << " " << vec.colorId << endl;
}

void print_color(Color color)
{
	cout << "Color: " << color.r << " " << color.g << " " << color.b << endl;
}
