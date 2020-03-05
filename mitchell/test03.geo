element_size = 0.20;
Point(1) = {-1.0, -1.0,     0.0, element_size };
Point(2) = {+1.0, -1.0,     0.0, element_size };
Point(3) = {+1.0, -0.00001, 0.0, element_size };
Point(4) = {+0.0,  0.0,     0.0, element_size };
Point(5) = {+1.0, +0.00001, 0.0, element_size };
Point(6) = {+1.0, +1.0,     0.0, element_size };
Point(7) = {-1.0, +1.0,     0.0, element_size };

Line(8) = {1,2};
Line(9) = {2,3};
Line(10) = {3,4};
Line(11) = {4,5};
Line(12) = {5,6};
Line(13) = {6,7};
Line(14) = {7,1};

Line Loop(15) = {8,9,10,11,12,13,14};

Plane Surface(16) = {15};

