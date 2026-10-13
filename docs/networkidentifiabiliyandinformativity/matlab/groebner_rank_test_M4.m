function out = groebner_rank_test_M4(A,B,C,D)
% GROEBNER_RANK_TEST_M4
%
% Implements Theorem 3 of Hansson et al. for k = p.
% Tests generic full row rank of
%
%    Pi = C*(I-A)^(-1)*B + D
%
% INPUTS (symbolic matrices):
%   A : n x n
%   B : n x (mr+me)
%   C : p x n
%   D : p x (mr+me)
%
% OUTPUT (struct):
%   out.G        Gröbner basis
%   out.Gt       Elimination basis (polynomials independent of t)
%   out.case     Which item of Theorem 3 applies
%   out.fullRankGeneric   logical
%
% Requires Symbolic Math Toolbox.

%% sizes
[p,~] = size(C);
n = size(A,1);

%% symbolic variable t
syms t

%% compute Pi
I = eye(n,'sym');
Pi = simplify(C*inv(I-A)*B + D);

%% ---------------------------------------------------------
% Clear denominators row-wise -> Pi_tilde
%% ---------------------------------------------------------

[num,den] = numden(Pi);

Pi_tilde = sym(zeros(size(Pi)));
Drows = sym(ones(p,1));

for i=1:p
    % LCM of denominators in row i
    d = den(i,1);
    for j=2:size(Pi,2)
        d = lcm(d,den(i,j));
    end
    Drows(i) = d;
    Pi_tilde(i,:) = simplify(d*Pi(i,:));
end

% Exclude singular network equations and undefined rational entries.
[detNum,detDen] = numden(det(I-A));
[~,inputDen] = numden([A(:); B(:); C(:); D(:)]);
Dprod = simplify(detNum*detDen*prod(inputDen)*prod(Drows));
if isequal(simplify(Dprod), sym(0))
    error('The network has no admissible point with invertible I-A.');
end

%% ---------------------------------------------------------
% k = p minors
%% ---------------------------------------------------------

k = p;
cols = size(Pi_tilde,2);

if cols < k
    error('Pi has fewer columns than p; cannot test full row rank.');
end

comb = nchoosek(1:cols,k);

minors = sym([]);

for i=1:size(comb,1)
    Msub = Pi_tilde(:,comb(i,:));
    minors(end+1,1) = det(Msub); %#ok<AGROW>
end

%% ---------------------------------------------------------
% Form ideal Ik = < minors , 1 - t*D >
%% ---------------------------------------------------------

polys = [1 - t*Dprod; minors];

% t must precede every model variable for the elimination step.
vars = symvar(polys);
vars = [t, vars(~arrayfun(@(v) isequal(v,t), vars))];

% .' is essential: ' would conjugate the symbolic polynomials.
G = gbasis(polys.',vars,'MonomialOrder','lexicographic');

%% ---------------------------------------------------------
% Analyze Gröbner basis (Theorem 3)
%% ---------------------------------------------------------

out.G = G;
out.Gt = sym([]);

% Case 1: reduced GB is {1}
if numel(G)==1 && isequal(G(1),sym(1))
    out.case = 1;
    out.fullRankGeneric = true;
    return
end

% Remove polynomials containing t
Gt = [];
for i=1:numel(G)
    if ~has(G(i),t)
        Gt = [Gt; G(i)]; %#ok<AGROW>
    end
end
out.Gt = Gt;

% Case 2: basis excluding t is zero ideal
if isempty(Gt)
    out.case = 2;
    out.fullRankGeneric = false;
    return
end

% Case 3: generic rank >= p but not everywhere
out.case = 3;
out.fullRankGeneric = true;

end
