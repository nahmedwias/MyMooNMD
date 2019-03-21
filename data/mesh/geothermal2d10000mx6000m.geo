// Geometry (m)

Lx = 10000.; 
Ly = 6000.;
d = 1000.;
r = 5e-4; // well radius 
x_I = Lx/2 - d/2;
x_E = Lx/2 + d/2;
y_W = Ly/2;


h_fine = 5.;
h_max = 200.; // background size


Field[1] = Box;
//Set the option string of the expression-th field. small interior box 1
Field[1].XMin = x_I- h_max/2;
Field[1].YMin = y_W - h_max/2;
Field[1].XMax = x_I + h_max/2;
Field[1].YMax = y_W + h_max/2;
Field[1].ZMin = -1;
Field[1].ZMax = 1;
Field[1].VIn = 0.5 * h_fine;
Field[1].VOut = h_max;

Field[2] = Box;
//Set the option string of the expression-th field.  small interior box 1
Field[2].XMin = x_E-h_max/2;
Field[2].YMin = y_W - h_max/2;
Field[2].XMax = x_E +h_max/2;
Field[2].YMax = y_W + h_max/2;
Field[2].ZMin = -1;
Field[2].ZMax = 1;
Field[2].VIn = 0.5 * h_fine;
Field[2].VOut = h_max;

Field[3] = Box;
//Set the option string of the expression-th field.
Field[3].XMin = x_I - 2*d;
Field[3].YMin = y_W - d;
Field[3].XMax = x_E + 2*d;
Field[3].YMax = y_W + d;
Field[3].ZMin = -1;
Field[3].ZMax = 1;
Field[3].VIn = 4*h_fine;
Field[3].VOut = h_max;


//Field[4] = Box;
//Set the option string of the expression-th field.
//Field[4].XMin = x_I-d/2;
//Field[4].YMin = y_W - d/2;
//Field[4].XMax = x_E +d/2;
//Field[4].YMax = y_W + d/2;
//Field[4].ZMin = -1;
//Field[4].ZMax = 1;
//Field[4].VIn = 12*h_fine;
//Field[4].VOut = h_max;

Field[5] = Min;
//+
//Field[5].FieldsList = {1, 2, 3, 4};
Field[5].FieldsList = {1, 2, 3};

Background Field = 5;


//+
Point(1) = {0, 0, 0, 200.0};
//+
Point(2) = {Lx, 0, 0, 200.0};
//+
Point(3) = {Lx, Ly, 0, 200.0};
//+
Point(4) = {0, Ly, 0, 200.0};

//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};


Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {1};
