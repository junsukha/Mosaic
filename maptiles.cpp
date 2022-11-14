/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */

    //MosaicCanvas* result = new MosaicCanvas(row, column);
    int row = theSource.getRows();
    int column = theSource.getColumns();
    MosaicCanvas* result = new MosaicCanvas(row, column);

    map<Point<3>, int> avgToTilt;
    vector<Point<3>> points;

    for(unsigned int i = 0; i < theTiles.size(); ++i){
      LUVAPixel pix = theTiles[i].getAverageColor();
      avgToTilt[convertToXYZ(pix)] = i;
      points.push_back(convertToXYZ(pix));
    }

    KDTree<3> tree(points);

    LUVAPixel pix;
    Point<3> avgPoint;
    Point<3> nearestPoint;

    for(int r = 0; r < row; r++){
      for(int c = 0; c < column; c++){
        //LUVAPixel pix =
        pix = theSource.getRegionColor(r, c);
        avgPoint = convertToXYZ(pix);
        nearestPoint = tree.findNearestNeighbor(avgPoint);
        int idx = avgToTilt[nearestPoint];
        //result->setTile(r, c, theTiles[idx]);
        result->setTile(r, c, &theTiles[idx]);
      }
    }
    return result;
}
