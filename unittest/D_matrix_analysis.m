A = [1,2;3;4]
A = [1,2;3,4]
A * diag([2,2])
[1,1,1;2,2,2;3,3,3;4,4,4]
m = [1,1,1;2,2,2;3,3,3;4,4,4];
m
m = [1,1,1;2,2,2;3,3,3;4,4,4]';
m
m = [1,1,1;2,2,2;3,3,3;4,4,4;5,5,5]';
m
diag([0.5,0.5,1/3,1/3,1/3])
D = diag([0.5,0.5,1/3,1/3,1/3]);
m * D
D*m
m*D
V=[1,0;1,0;0,1;0,1;0,1]
m * D
m * D * V
D = sqrt(diag([0.5,0.5,1/3,1/3,1/3]));
D
m * D^2 * V
(m * D^2 * V)*V'
D2 = D;
D(1,2) = D(1,1); D(2,1) = D(1,1);
D
D = sqrt(diag([0.5,0.5,1/3,1/3,1/3]));
D2
D2(1,2) = D2(1,1); D2(2,1) = D2(1,1);
D2
D2(3,4:5) = D2(3,3); D2(4:5,3) = D2(3,3);
D2
D2(4,5) = D2(3,3); D2(5,4) = D2(3,3);
D2
(m * D^2 * V)*V'
(m * D2^2 * V)*V'
D
D2