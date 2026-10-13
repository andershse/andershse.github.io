% Big Test: 11 variables, 5x5 system, 3x5 Flat C, 15 Outputs
clear; clc;

% 1. Define exactly 11 source variables
a = sym('a', [1, 5]); % a1...a5
b = sym('b', [1, 5]); % b1...b5
c = sym('c', [1]); % c1
sourceVars = [a, b, c];

% 2. A-matrix (5x5, zero diagonal, sparse cycle)
% Structure: 1->2, 2->3, 3->4, 4->5, 5->1
A = [0, a(1), 0, 0, 0; ...
     0, 0, a(2), 0, 0; ...
     0, 0, 0, a(3), 0; ...
     0, 0, 0, 0, a(4); ...
     a(5), 0, 0, 0, 0];

% 3. B-matrix (5x5, variables mixed with ones and zeros)
% We place variables to ensure every output y_ij is non-trivial
B = [b(1), 1, 0, 0, 0; ...
     0, b(2), 1, 0, 0; ...
     0, 0, b(3), 1, 0; ...
     0, 0, 0, b(4), 1; ...
     b(5), 0, 0, 0, 1];

% 4. C-matrix (3x5, measuring nodes 1, 2, and 4)
C = [c(1), 0, 0, 0, 0; ... % Row 1 picks node 1
     0 1 0 0 0;        % Row 2 picks node 2
     0, 0, 0, 1, 0];   % Row 3 picks node 4
% 5. Compute the expression: y = C * inv(I - A) * B
% 3x5 matrix = 15 rational functions
fprintf('Computing 3x5 rational expression (15 functions)...\n');
expr = C * inv(eye(5) - A) * B;

% 6. Extract Numerators and Denominators
[P_mat, Q_mat] = numden(expr);

% 7. Run the Uddin-based Dimension Analysis
fprintf('Starting image-dimension analysis (15 equations and 11 variables)...\n');
fprintf('This may take a few minutes due to the 5x5 inversion complexity.\n');

tic;
dim = dimMatrixRationalImage(P_mat, Q_mat, sourceVars);
toc;

fprintf('\n------------------------------------------\n');
fprintf('Final Result: The dimension is %d\n', dim);
fprintf('If dim < 11, MATLAB found a hidden algebraic relation in the 5x5 system.\n');


tic
B = [B eye(5)];
C = [C; zeros(5,5)];
D = [zeros(3,10);eye(5) zeros(5,5)];
out = groebner_rank_test_M4(A,B,C,D)
toc

out.case
out.fullRankGeneric
out.Gt
