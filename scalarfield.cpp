#include "scalarfield.h"
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

ScalarField::ScalarField(int* _dim, float* _origin, float* _step)
{
    for(int i=0;i<3;i++){
        dim[i] = _dim[i];
        origin[i] = _origin[i];
        step[i] = _step[i];
    }
    data = new float [ dim[0] * dim[1] * dim[2] ];
}
ScalarField::~ScalarField(){
    delete[] data;
}

float ScalarField::getValue(int *indices){
    int x = indices[0];
    int y = indices[1];
    int z = indices[2];
    int index = x*dim[1]*dim[2] + y*dim[2] + z;
    return data[index];
}

float ScalarField::getValue(float *point){
    int indices[3];
    if(point[0]<origin[0]) point[0] = origin[0];
    if(point[1]<origin[1]) point[1] = origin[1];
    if(point[2]<origin[2]) point[2] = origin[2];
    indices[0] = min((int)((point[0]-origin[0])/step[0]), dim[0]-1);
    indices[1] = min((int)((point[1]-origin[1])/step[1]), dim[1]-1);
    indices[2] = min((int)((point[2]-origin[2])/step[2]), dim[2]-1);
    return getValue(indices);
}

float ScalarField::getValueBiLinear(float *point){
    int indices[3];
    if(point[0]<origin[0]) point[0] = origin[0];
    if(point[1]<origin[1]) point[1] = origin[1];
    if(point[2]<origin[2]) point[2] = origin[2];

    float x = (point[0]-origin[0])/step[0];
    int xInt = (int) x;
    float xDiff = x - xInt;

    float y = (point[1]-origin[1])/step[1];
    int yInt = (int) y;
    float yDiff = y - yInt;

    float z = (point[2]-origin[2])/step[2];
    int zInt = (int) z;
    float zDiff = z - zInt;

    indices[0] = min(xInt, dim[0]-1);
    indices[1] = min(yInt, dim[1]-1);
    indices[2] = min(zInt, dim[2]-1);
    float v000 = getValue(indices);

    indices[0] = min(xInt, dim[0]-1);
    indices[1] = min(yInt, dim[1]-1);
    indices[2] = min(zInt + 1, dim[2]-1);
    float v001 = getValue(indices);

    indices[0] = min(xInt, dim[0]-1);
    indices[1] = min(yInt + 1, dim[1]-1);
    indices[2] = min(zInt, dim[2]-1);
    float v010 = getValue(indices);

    indices[0] = min(xInt, dim[0]-1);
    indices[1] = min(yInt + 1, dim[1]-1);
    indices[2] = min(zInt + 1, dim[2]-1);
    float v011 = getValue(indices);

    indices[0] = min(xInt + 1, dim[0]-1);
    indices[1] = min(yInt, dim[1]-1);
    indices[2] = min(zInt, dim[2]-1);
    float v100 = getValue(indices);

    indices[0] = min(xInt + 1, dim[0]-1);
    indices[1] = min(yInt, dim[1]-1);
    indices[2] = min(zInt + 1, dim[2]-1);
    float v101 = getValue(indices);

    indices[0] = min(xInt + 1, dim[0]-1);
    indices[1] = min(yInt + 1, dim[1]-1);
    indices[2] = min(zInt, dim[2]-1);
    float v110 = getValue(indices);

    indices[0] = min(xInt + 1, dim[0]-1);
    indices[1] = min(yInt + 1, dim[1]-1);
    indices[2] = min(zInt + 1, dim[2]-1);
    float v111 = getValue(indices);

    float v00 = v000 * (1-zDiff) + v001 * zDiff;
    float v01 = v010 * (1-zDiff) + v011 * zDiff;
    float v10 = v100 * (1-zDiff) + v101 * zDiff;
    float v11 = v110 * (1-zDiff) + v111 * zDiff;

    float v0 = v00 * (1-yDiff) + v01 * yDiff;
    float v1 = v10 * (1-yDiff) + v11 * yDiff;

    return v0 *(1-xDiff) + v1 * xDiff;
}

Color ScalarField::getColor(Color low, Color high, float val){
    //clamp
    if(val<vmin) val = vmin;
    if(val>vmax) val = vmax;
    float colorVal = (val-vmin)/(vmax-vmin);
    // interpolate
    return Color(low.r + colorVal * (high.r - low.r), low.g + colorVal * (high.g - low.g),
                 low.b + colorVal * (high.b - low.b));
}

Color ScalarField::getColor2(Color low, Color high, float val){
    //clamp
    if(val<vmin) val = vmin;
    if(val>vmax) val = vmax;
    float mid = (vmin + vmax)/2;
    if(val<mid){
        float colorVal = (val-vmin)/(mid-vmin);
        return Color(low.r + colorVal * (1 - low.r), low.g + colorVal * (1 - low.g),
                     low.b + colorVal * (1 - low.b));
    }else{
        float colorVal = (val-mid)/(vmax-mid);
        return Color(1 + colorVal * (high.r - 1), 1 + colorVal * (high.g - 1),
                     1 + colorVal * (high.b - 1));
    }
}

ScalarField* loadField(const char* filename){
    int dim[3];
    float origin[3];
    float step[3];
    FILE* fp = fopen(filename, "r");
    char line[100];
    char temp[40];
    while(true){
        fgets(line, 100, fp);
        char* ptr;
        for(ptr = line; *ptr==' '; ptr++);
        if(*ptr == '#'){
            continue;
        }
        sscanf(line, "%s %s %s %s %s %d %d %d", temp, temp, temp, temp, temp, &dim[0], &dim[1], &dim[2]);
        fgets(line, 100, fp);
        sscanf(line, "%s %f %f %f", temp, &origin[0], &origin[1], &origin[2]);
        fgets(line, 100, fp);
        sscanf(line, "%s %f %s %s", temp, &step[0], temp, temp);
        fgets(line, 100, fp);
        sscanf(line, "%s %s %f %s", temp, temp, &step[1], temp);
        fgets(line, 100, fp);
        sscanf(line, "%s %s %s %f", temp, temp, temp, &step[2]);
        fgets(line, 100, fp);
        fgets(line, 100, fp);
        break;
    }
    ScalarField* field = new ScalarField(dim, origin, step);
    for(int i=0;i<dim[0]*dim[1]*dim[2];i++){
        float val;
        fscanf(fp, "%f", &val);
        if(i==0){
            field->vmin = field->vmax = val;
        }
        field->vmin = min(field->vmin, val);
        field->vmax = max(field->vmax, val);
        field->data[i] = val;
    }
    field->vmin = -2;
    field->vmax = 2;
    fclose(fp);
    return field;
}
