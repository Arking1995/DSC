﻿#include "Camera.h"
#include <glm/gtc/matrix_transform.hpp>
#include "src/glew/include/GL/glew.h"
#include <glm/vec3.hpp>
#include <QVector2D>

#ifndef M_PI
#define M_PI	3.1415926535
#endif

Camera::Camera() {
	xrot = 40.0f;
	yrot = 0.0;
	zrot = 0.0f;
//Liu changed this initialization function on 2019/06/24
	rot = QVector3D(xrot, yrot, zrot);
	//pos = glm::vec3(0, 15, 80);
	pos = QVector3D(0, 15, 80);
	fovy = 45.0f;
}

/**
 * Call this function when the mouse button is pressed.
 */
void Camera::mousePress(int mouse_x, int mouse_y) {
	//Liu changed 'glm::vec2' to 'qvector2d' on 2019/06/25
	mouse_pos = QVector2D(mouse_x, mouse_y);
}

/**
 * Call this function whenever the mouse moves while rotating the model.
 */
void Camera::rotate(int mouse_x, int mouse_y) {
	xrot += mouse_y - mouse_pos.y();
	yrot += mouse_x - mouse_pos.x();
	updateMVPMatrix();
	//Liu changed 'glm::vec2' to 'qvector2d' on 2019/06/25
	mouse_pos = QVector2D(mouse_x, mouse_y);
}

/**
 * Call this function whenever the mouse moves while zooming.
 */
/*void Camera::zoom(int mouse_x, int mouse_y) {
	pos.z() += (mouse_pos.y - mouse_y) * 0.5f;
	updateMVPMatrix();

	mouse_pos = glm::vec2(mouse_x, mouse_y);
}*/

void Camera::zoom(float delta) {
	
	pos.setZ( pos.z() - delta * 0.1f );
	updateMVPMatrix();
}

/**
 * Call this function whenever the mouse moves while moving the model.
 */
void Camera::move(int mouse_x, int mouse_y) {
	pos.setX( pos.x() - (mouse_x - mouse_pos.x()) * 0.1 );
	pos.setY( pos.y() - (mouse_y - mouse_pos.y()) * 0.1 );
	updateMVPMatrix();

	//Liu changed 'glm::vec2' to 'qvector2d' on 2019/06/25
	mouse_pos = QVector2D(mouse_x, mouse_y);
}

/**
 * Update perspective projection matrix, and then, update the model view projection matrix.
 */
void Camera::updatePMatrix(int width,int height) {
	float aspect = (float)width / (float)height;
	float zfar = 3000.0f;
	float znear = 0.1f;
	float f = 1.0f / tan(fovy * M_PI / 360.0f);

	// projection行列
	// ただし、mat4はcolumn majorなので、転置した感じで初期構築する。
	// つまり、下記の一行目は、mat4の一列目に格納されるのだ。
	glm::mat4 pMatrix = glm::mat4(
		 f/aspect,	0,								0,									0,
				0,	f,								0,						 			0,
			    0,	0,		(zfar+znear)/(znear-zfar),		                           -1,
			    0,	0, (2.0f*zfar*znear)/(znear-zfar),									0);

	updateMVPMatrix();
}

/**
 * Update the model view projection matrix
 */
void Camera::updateMVPMatrix() {
	// create model view matrix
	
	//This part was modified by Liu on 2019/06/25	
	mvMatrix = QMatrix4x4();
	mvMatrix.translate(-pos.x(), -pos.y(),-pos.z() );
	mvMatrix.rotate(xrot * (float)M_PI / 180.0f, QVector3D(1, 0, 0));
	mvMatrix.rotate(yrot * (float)M_PI / 180.0f, QVector3D(0, 1, 0));
	mvMatrix.rotate(zrot * (float)M_PI / 180.0f, QVector3D(0, 0, 1));

	// create model view projection matrix
	mvpMatrix = pMatrix * mvMatrix;
}
