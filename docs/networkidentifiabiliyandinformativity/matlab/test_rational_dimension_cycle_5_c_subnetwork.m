% Five-node example: identify (G12,Q1,Q2,R1) in the A-subnetwork.
clear; clc;

a = sym('a', [1, 5]);
b = sym('b', [1, 5]);
c = sym('c', [1, 1]);
sourceVars = [a(1), b(1), b(2), c(1)];

A = [0, a(1), 0, 0, 0; ...
     0, 0, a(2), 0, 0; ...
     0, 0, 0, a(3), 0; ...
     0, 0, 0, 0, a(4); ...
     a(5), 0, 0, 0, 0];

B = [b(1), 1, 0, 0, 0; ...
     0, b(2), 1, 0, 0; ...
     0, 0, b(3), 1, 0; ...
     0, 0, 0, b(4), 1; ...
     b(5), 0, 0, 0, 1];

C = [c(1), 0, 0, 0, 0; ...
     0, 1, 0, 0, 0];

% Columns: measured boundary input, r1, r2, r3.
Ffull = C(:,1:2) * inv(eye(2) - A(1:2,1:2)) * ...
        [[0;1], B(1:2,1:3)];
Ffew = Ffull(:,1:3);  % r3 is switched off.

[numFull, denFull] = numden(Ffull);
[numFew, denFew] = numden(Ffew);

tic;
dimFull = dimMatrixRationalImage(numFull, denFull, sourceVars);
fprintf('Image dimension with r1,r2,r3: %d (expected 4)\n', dimFull);
toc;

tic;
dimFew = dimMatrixRationalImage(numFew, denFew, sourceVars);
fprintf('Image dimension with r1,r2: %d (expected 4)\n', dimFew);
toc;

% Informativity uses excitation through the full network.
Bfull = [B, eye(5)];
Cfull = [C; A(2,:); zeros(5,5)];
Dfull = [zeros(3,10); eye(5), zeros(5,5)];

tic;
rankFull = groebner_rank_test_M4(A,Bfull,Cfull,Dfull);
fprintf('Generic full row rank with r1,...,r5: %d (case %d)\n', ...
        rankFull.fullRankGeneric, rankFull.case);
toc;

cols123 = [1:3, 6:10];
rows123 = 1:6;  % Three measured signals and r1,r2,r3.
tic;
rank123 = groebner_rank_test_M4(A,Bfull(:,cols123), ...
          Cfull(rows123,:),Dfull(rows123,cols123));
fprintf('Generic full row rank with r1,r2,r3: %d (case %d)\n', ...
        rank123.fullRankGeneric, rank123.case);
toc;

cols12 = [1:2, 6:10];
rows12 = 1:5;  % Three measured signals and r1,r2.
tic;
rank12 = groebner_rank_test_M4(A,Bfull(:,cols12), ...
         Cfull(rows12,:),Dfull(rows12,cols12));
fprintf('Generic full row rank with r1,r2: %d (case %d)\n', ...
        rank12.fullRankGeneric, rank12.case);
toc;
