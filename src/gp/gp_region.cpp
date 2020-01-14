#include "gp_data.h"
#include "gp_region.h"

void findCellRegionDist() {
    for (int i = 0; i < numCells; i++) {
        // for each cell, find nearest region block
        double x = cellX[i];
        double y = cellY[i];
        int region = cellFence[i];
        vector<FRect>::iterator ri = regions[region].rects.begin();
        vector<FRect>::iterator re = regions[region].rects.end();
        double dist = -1;
        double regionx = -1;
        double regiony = -1;
        int regionr = -1;
        for (int r = 0; ri != re; ++ri, r++) {
            double distx = 0;
            double disty = 0;
            double rx = x;
            double ry = y;
            double lx = ri->lx;
            double ly = ri->ly;
            double hx = ri->hx;
            double hy = ri->hy;
            if (x < lx) {
                distx = lx - x;
                rx = lx;
            } else if (x > hx) {
                distx = x - hx;
                rx = hx;
            }
            if (y < ly) {
                disty = ly - y;
                ry = ly;
            } else if (y > hy) {
                disty = y - hy;
                ry = hy;
            }
            if (r == 0 || distx + disty < dist) {
                dist = distx + disty;
                regionx = rx;
                regiony = ry;
                regionr = r;
                if (dist == 0) {
                    break;
                }
            }
        }
        if (regionr < 0) {
            cerr << "no rect defined for region " << region << endl;
            exit(1);
        }
        cellFenceX[i] = regionx;
        cellFenceY[i] = regiony;
        cellFenceRect[i] = regionr;
        cellFenceDist[i] = dist;
    }
}
/*int legal_count=0;
void legalizeRegion(){
    //assume findCellRegionDist() is called already
    double epsilon=circuit.site_w*0.1;
    for(int i=1; i<numCells; i++){
        double x=xCellCoord[i];
        double y=yCellCoord[i];

        x=max(coreLX, min(coreHX, x));
        y=max(coreLY, min(coreHY, y));

        int f=cellFence[i];
        int r=cellFenceRect[i];
        FRegionBlock &block=circuit.regionblocks[f+1][r];
        if(x<block.x){
            xCellCoord[i]=block.x+(x-coreLX)/(block.x-coreLX)*epsilon;
            //it is always x >= coreLX, so when x < block.x, block.x > coreLX , divided by zero avoided
        }else if(x>block.x+block.w){
            xCellCoord[i]=block.x+block.w-epsilon+(x-block.x-block.w)/(coreHX-block.x-block.w)*epsilon;
            //it is always x <= coreHX, so when x > block.x+block.w, block.x+block.w < coreHX , divided by zero avoided
        }else{
            xCellCoord[i]=block.x+epsilon+(x-block.x)/(block.w)*(block.w-epsilon-epsilon);
        }
        if(y<block.y){
            yCellCoord[i]=block.y+(y-coreLY)/(block.y-coreLY)*epsilon;
        }else if(y>block.y+block.h){
            yCellCoord[i]=block.y+block.h-epsilon+(y-block.y-block.h)/(coreHY-block.y-block.h)*epsilon;
        }else{
            yCellCoord[i]=block.y+epsilon+(y-block.y)/(block.h)*(block.h-epsilon-epsilon);
        }
    }
}
*/
