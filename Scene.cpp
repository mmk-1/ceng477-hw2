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

Matrix4 calculate_model_transformation(const Mesh *mesh, const Scene *scene)
{
	// Initalize result matrix to identity matrix
	// multiply instead of assignment

	Matrix4 result = getIdentityMatrix();
	for (int i = 0; i < mesh->numberOfTransformations; i++)
	{
		// Be careful that transormationIds are 1-indexed
		switch (mesh->transformationTypes[i])
		{
		case 'r':
		{
			Rotation *rotation = scene->rotations[mesh->transformationIds[i] - 1];
			result = multiplyMatrixWithMatrix(result, calculate_rotation_transformation(rotation));
		}
		break;
		case 't':
		{
			Translation *translation = scene->translations[mesh->transformationIds[i] - 1];
			double translation_matrix[4][4] = {{1, 0, 0, translation->tx}, {0, 1, 0, translation->ty}, {0, 0, 1, translation->tz}, {0, 0, 0, 1}};
			result = multiplyMatrixWithMatrix(result, translation_matrix);
		}
		break;
		case 's':
		{
			Scaling *scale = scene->scalings[mesh->transformationIds[i] - 1];
			double scale_matrix[4][4] = {{scale->sx, 0, 0, 0}, {0, scale->sy, 0, 0}, {0, 0, scale->sz, 0}, {0, 0, 0, 1}};
			result = multiplyMatrixWithMatrix(result, scale_matrix);
		}
		break;
		default:
		{
			cout << "Invalid transformation type" << endl;
			exit(1);
		}
		break;
		}
	}
	return result;
}

Matrix4 calculate_camera_transformation(const Camera *camera)
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

	return multiplyMatrixWithMatrix(r_matrix, t_matrix);
}

Matrix4 calculate_projection_transformation(const Camera *camera, bool type)
{
	Matrix4 result;
	double left = camera->left;
	double right = camera->right;
	double bottom = camera->bottom;
	double top = camera->top;
	double near = camera->near;
	double far = camera->far;
	switch (type)
	{
	case false:
	{
		// Orthographic
		double matrix_orth[4][4] = {
			{2 / (right - left), 0, 0, -((right + left) / (right - left))},
			{0, 2 / (top - bottom), 0, -((top + bottom) / (top - bottom))},
			{0, 0, -(2 / (far - near)), -((far + near) / (far - near))},
			{0, 0, 0, 1}};

		result = Matrix4(matrix_orth);
	}
	break;
	default:
	{
		// Perspective
		double matrix_perspective[4][4] = {
			{(2 * near) / (right - left), 0, (right + left) / (right - left), 0},
			{0, (2 * near) / (top - bottom), (top + bottom) / (top - bottom), 0},
			{0, 0, -((far + near) / (far - near)), -((2 * far * near) / (far - near))},
			{0, 0, -1, 0}};

		result = Matrix4(matrix_perspective);
	}
	break;
	}
	return result;
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// Overall steps:
	// 1. Model Transformation [X] Test []
	// 2. Camera Transformation [X] Test [X]
	// 3. Projection Transformation [X] Test []
	// 4. Clipping [] Test []
	// 5. Backface Culling [] Test []
	// 6. Viewport Transformation [] Test []
	// 7. Rasterization [] Test []
	// 8. Depth Buffer [] Test []

	Matrix4 matrix_camera = calculate_camera_transformation(camera);
	print_matrix4(matrix_camera);
	Matrix4 matrix_projection = calculate_projection_transformation(camera, camera->projectionType);

	std::vector<std::map<int, Vec3>> meshes_transformed_vertices = std::vector<std::map<int, Vec3>>(meshes.size());

	// Go through all meshes and apply transformations
	for (int m = 0; m < meshes.size(); m++)
	{
		const Mesh *mesh = meshes[m];
		Matrix4 matrix_model = calculate_model_transformation(mesh, this);
		matrix_model = multiplyMatrixWithMatrix(matrix_camera, matrix_model);
		matrix_model = multiplyMatrixWithMatrix(matrix_projection, matrix_model);

		for (int t = 0; t < mesh->triangles.size(); t++)
		{
			const Triangle &triangle = mesh->triangles[t]; // Get triangle
			for (int v = 0; v < 3; v++)
			{

				if (meshes_transformed_vertices[m].find(triangle.vertexIds[v]) != meshes_transformed_vertices[m].end())
					continue;

				// Transform vertex
				Vec3 *vertex = this->vertices[triangle.vertexIds[v] - 1];
				Vec4 vertex_4 = Vec4(vertex->x, vertex->y, vertex->z, 1);
				Vec4 transformed_vertex_4 = multiplyMatrixWithVec4(matrix_model, vertex_4);

				Vec3 transformed_vertex = Vec3(transformed_vertex_4.x, transformed_vertex_4.y, transformed_vertex_4.z);

				// Store transformed vertex
				meshes_transformed_vertices[m][triangle.vertexIds[v]] = transformed_vertex;
			}
		}
	}

	// Clipping (Step 4)
	for (int m = 0; m < meshes.size(); m++)
	{
		const Mesh *mesh = meshes[m];
		if (mesh->type == WIREFRAME_MESH)
			continue;

		// for (int y = 0; y < mesh->triangles.size(); y++)
		// {
		// 	const Triangle &triangle = mesh->triangles[y];
		// }

		// Culling (Step 5)
	}
}

/*
 *********************Our Implementation ends here*********************************
 */

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
}
