#include "RoadVertex.h"

RoadVertex::RoadVertex() {
	this->pt = QVector2D(0.0f, 0.0f);
	this->valid = true;
	//this->onBoundary = false;
	this->eigenType = EIGEN_TYPE_NONE;
	this->new_vertex = false;
}

RoadVertex::RoadVertex(const QVector2D &pt) {
	this->pt.setX(pt.x());
	this->pt.setY(pt.y());
	this->valid = true;
	//this->onBoundary = false;
	this->eigenType = EIGEN_TYPE_NONE;
	this->new_vertex = false;
}

const QVector2D& RoadVertex::getPt() const {
	return pt;
}

