%% From Class on 30 Jan 2019
clear all
close all
clc

%% Starting matrix
A = [5, 3, -1; -4, 2, 0; -3, 1, -1]

%% Build an augmented matrix
B = [A, eye(3)]

%%
% go through simple Gauss elimination
B(3,:) = B(1,:)*3/5 + B(3,:);
B

B(2,:) = B(1,:)*4/5 + B(2,:);
B

B(3,:) = -2.8/4.4*B(2,:) + B(3,:);
B

B(3,:) = B(3,:)/B(3,3)

B(2,:) = B(2,:) + 0.8*B(3,:);
B

B(1,:) = B(1,:) + B(3,:);
B

B(1,:) = B(1,:) - B(2,:)*3/4.4;
B

B(1,:) = B(1,:)/5;
B(2,:) = B(2,:)/4.4;

B

%%
% We've made the inverse of A in the last 3 colums of B:
C = B(:,4:6)

norm( eye(3) - A*C )

%% Built in function
% The rref function goes through these steps for us
D = [A, eye(3)]
rref(D)


%% Other ways to compute an inverse
inv(A)
A\eye(3)

%% Now with 1 sign change
% the system is singular
A1 = [5, -3, -1; -4, 2, 0; -3, 1, -1]
det(A1)
D = [A1, eye(3)]
rref(D)

%%
% we can see the last row on the left are zeros.  The right side of this is
% the null space of the matrix A1
% A1*v = 0
A1*[1; 2; -1]
