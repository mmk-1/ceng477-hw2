#ifndef _UTILS_H_
#define _UTILS_H_
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

Matrix4 calculate_rotation_transformation(const Rotation *rotation);
Matrix4 calculate_model_trans_matrix(const Scene *scene, const Mesh *mesh);
Matrix4 calculate_camera_trans_matrix(const Camera *camera);
Matrix4 calculate_projection_trans_matrix(const Camera *camera);
std::vector<Vec4> clip_line(Scene *scene, std::vector<Vec4> &line, std::vector<Color> &line_colors, bool &visible);
bool is_visible(double den, double num, double &te, double &tl);
void rasterize_line(Scene *scene, Vec3 v0, Vec3 v1, Color c0, Color c1);
void rasterize_triangle(Scene *scene, Camera *camera, Matrix4 &viewport_transformation_matrix, Vec4 v0, Vec4 v1, Vec4 v2);

void print_vec4(Vec4 vec);
void print_matrix4(Matrix4 matrix);
void print_color(Color color);

#endif
