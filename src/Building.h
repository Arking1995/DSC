#pragma once

#include "Polygon3D.h"

class Building {
public:
	Polygon3D buildingFootprint;
	int bldType;
	int numStories;
	int baseLevel;
	QColor color;

public:
	Building() {
		bldType = -1;
		numStories = 0;
	}
};
