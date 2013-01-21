#ifndef SCALARFIELD_H
#define SCALARFIELD_H

class Color{
public:
    float r, g, b;

    Color(float _r, float _g, float _b):r(_r), g(_g), b(_b){
    }
};

class ScalarField
{
public:

    float* data;
    float vmin, vmax;
    int dim[3];
    float origin[3];
    float step[3];

    float getValue(float* point);
    float getValueBiLinear(float* point);
    float getValue(int* indices);
    Color getColor(Color low, Color high, float val);
    Color getColor2(Color low, Color high, float val);

    ScalarField(int* dim, float* origin, float* step);
    ~ScalarField();
};

ScalarField* loadField(const char* filename);

#endif // SCALARFIELD_H
