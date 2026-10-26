function compare_interpolation(render_only)
% Reproduce the manuscript's interpolation comparison.
% Run: matlab -batch "compare_interpolation"
% Requires Optimization Toolbox. No third-party packages.
rng(271828,'twister');
out = fullfile(fileparts(mfilename('fullpath')),'comparison_results');
if ~exist(out,'dir'), mkdir(out); end
opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex-highs');
if nargin>0 && render_only
    saved=load(fullfile(out,'results.mat'));
    write_tables(out,saved.comparisons,saved.records,saved.timings);
    plot_results(out,saved.fractions,saved.curve,saved.lower);
    return
end
taus = [.18,.22]; dims = [2,4,8,16]; meshes = [4,16,64];
names = {'linear','one_state','interp4','interp16','interp64','rollout'};
p = 148/9; gamma = .5; tol = 1e-10;
records = {}; timings = {}; construction = {}; checks = struct;
checks.max_lp_objective_difference = 0;
checks.max_interpolation_identity_error = 0;
checks.max_constraint_violation = 0;
checks.max_positive_residual = 0;
checks.max_tail_bound = 0;
% Exact grid-node rollout actions and interval coefficients, shared by blocks.
tables = cell(size(meshes));
pair_rollout([.5;.5]); % Warm up the common evaluator before construction timing.
for k=1:numel(meshes)
    M=meshes(k); nodes=(0:M)/M; build_times=zeros(1,7);
    for rep=1:7
        tic; v=zeros(2,M+1);
        for j=1:M+1, v(:,j)=pair_rollout([nodes(j);1-nodes(j)]); end
        build_times(rep)=toc;
    end
    tables{k}=struct('M',M,'values',v);
    elapsed=median(build_times);
    construction(end+1,:)={M,M+1,numel(v),elapsed}; %#ok<AGROW>
end
for tau=taus
    for n=dims
        [A,B,C,d,K]=data(n,tau);
        nu=d+gamma*B'*(p*ones(n,1)); beta=10*ones(n,1)+gamma*A'*(p*ones(n,1))-p;
        assert(norm((eye(n)-gamma*(A+B*K)')\(10*ones(n,1)+K'*d)-p,inf)<1e-11);
        assert(norm(beta+K'*nu,inf)<1e-11);
        % Three equally sized cohorts; every vector has total inventory one.
        X=rand(n,80); X=X./sum(X,1);
        Y=10.^(-6*rand(n,80)); Y=Y./sum(Y,1);
        Z=rand(n,80);
        for j=1:80, Z(mod(j-1,n)+1,j)=0; end
        Z=Z./sum(Z,1); samples=[X,Y,Z];
        costs=zeros(numel(names),size(samples,2));
        for h=1:numel(names)
            for j=1:size(samples,2)
                x=samples(:,j); discount=1; value=0;
                while discount*p*sum(x)>tol
                    u=policy(x,h,tables);
                    checks.max_constraint_violation=max(checks.max_constraint_violation,max([0;-u;C*u-x]));
                    checks.max_positive_residual=max(checks.max_positive_residual,beta'*x+nu'*u);
                    value=value+discount*(10*sum(x)+d'*u);
                    x=A*x+B*u; assert(min(x)>-1e-12); x=max(x,0);
                    discount=discount*gamma;
                end
                tail=discount*p*sum(x);
                costs(h,j)=value;
                checks.max_tail_bound=max(checks.max_tail_bound,tail);
            end
            savings=100*(1-costs(h,:)/p);
            records(end+1,:)={tau,n,names{h},mean(savings(1:80)),mean(savings(81:160)),mean(savings(161:240)),min(savings),max(savings)}; %#ok<AGROW>
        end
        % Validate exact rollout against independently solved original LPs.
        lp_elapsed=zeros(1,12);
        for j=1:12
            x=samples(:,j); tic;
            [ul,fl,flag]=linprog(nu,C,x,[],[],zeros(2*n,1),[],opts);
            lp_elapsed(j)=toc; assert(flag>0);
            ur=policy(x,6,tables);
            checks.max_lp_objective_difference=max(checks.max_lp_objective_difference,abs(fl-nu'*ur));
            assert(max(C*ul-x)<1e-8);
        end
        for j=1:size(samples,2)
            x=samples(:,j); w=ones(n,1); a=policy(w,6,tables);
            z=x/sum(x); a0=w/sum(w); weight0=min(z./a0);
            weights=z-weight0*a0;
            interpolated=sum(x)*(weight0*a/sum(w)+K*weights);
            checks.max_interpolation_identity_error=max(checks.max_interpolation_identity_error,norm(interpolated-policy(x,2,tables),inf));
        end
        % Median scalar-call latency, including point location, seven repeats.
        for h=1:numel(names)
            policy(samples(:,1),h,tables); times=zeros(1,7);
            for rep=1:7
                tic;
                for j=1:1200, policy(samples(:,mod(j-1,240)+1),h,tables); end
                times(rep)=toc/1200;
            end
            timings(end+1,:)={tau,n,names{h},median(times)*1e6}; %#ok<AGROW>
        end
        timings(end+1,:)={tau,n,'rollout_linprog',median(lp_elapsed)*1e6}; %#ok<AGROW>
    end
end
assert(checks.max_lp_objective_difference<1e-8);
assert(checks.max_interpolation_identity_error<1e-12);
assert(checks.max_constraint_violation<1e-12);
assert(checks.max_positive_residual<1e-11);
% A dense two-dimensional sweep, with finite-horizon LP lower bounds on J*.
fractions=(0:100)/100; curve=cell(1,2); lower=cell(1,2); comparisons={};
checks.max_lp_bound_width=0;
for it=1:2
    tau=taus(it); [A,B,C,d]=data(2,tau);
    J=zeros(6,numel(fractions)); L=zeros(size(fractions));
    for j=1:numel(fractions)
        x0=[fractions(j);1-fractions(j)];
        for h=1:6, J(h,j)=cost(x0,h,tables,A,B,d,p,tol); end
        [L(j),width]=finite_lower(x0,A,B,C,d,24,opts);
        checks.max_lp_bound_width=max(checks.max_lp_bound_width,width);
        assert(all(J(:,j)>=L(j)-1e-8));
    end
    curve{it}=J; lower{it}=L;
    for h=1:6
        comparisons(end+1,:)={tau,names{h},mean(J(h,:)),max(J(h,:)./L-1)*100}; %#ok<AGROW>
    end
    exact_mu=3*p-9/(1-tau/2); exact_roll=3*p-10;
    assert(abs(cost([1;2],2,tables,A,B,d,p,tol)-exact_mu)<tol*2);
    assert(abs(cost([1;2],6,tables,A,B,d,p,tol)-exact_roll)<tol*2);
end
checks.optimal_value_lower_bound_error=p*.4^24;
checks.cost_tolerance=tol;
checks.matlab_version=version;
checks.computer=computer;
checks.seed=271828;
writetable(cell2table(records,'VariableNames',{'tau','n','policy','mean_uniform_saving_pct','mean_log_saving_pct','mean_boundary_saving_pct','min_saving_pct','max_saving_pct'}),fullfile(out,'costs.csv'));
writetable(cell2table(timings,'VariableNames',{'tau','n','policy','median_us'}),fullfile(out,'timing.csv'));
writetable(cell2table(construction,'VariableNames',{'intervals','offline_rollout_evaluations','stored_doubles','construction_seconds'}),fullfile(out,'construction.csv'));
writetable(cell2table(comparisons,'VariableNames',{'tau','policy','mean_grid_cost','max_gap_to_lower_bound_pct'}),fullfile(out,'optimality.csv'));
fid=fopen(fullfile(out,'checks.json'),'w'); fprintf(fid,'%s',jsonencode(checks,PrettyPrint=true)); fclose(fid);
for it=1:2
    writematrix([fractions',curve{it}',lower{it}'],fullfile(out,sprintf('curve_%02d.csv',round(taus(it)*100))));
end
plot_results(out,fractions,curve,lower);
write_tables(out,comparisons,records,timings);
save(fullfile(out,'results.mat'),'records','timings','construction','comparisons','checks','fractions','curve','lower');
disp(checks); disp(cell2table(comparisons));
end

function plot_results(out,fractions,curve,lower)
f=figure('Visible','off','Position',[100,100,1000,400]); taus=[.18,.22];
for it=1:2
    subplot(1,2,it); hold on;
    for h=[2,3,4,6], plot(fractions,100*(curve{it}(h,:)./lower{it}-1),'LineWidth',1.3); end
    xlabel('x_1, with x_1+x_2=1'); ylabel('Gap to optimal-cost lower bound (%)');
    title(sprintf('Nonzero linear feedback, tau = %.2f',taus(it))); grid on;
    legend({'Single correction','Interpolation, 4 intervals','Interpolation, 16 intervals','Exact rollout'},'Location','best');
end
exportgraphics(f,fullfile(out,'cost_comparison.pdf'),'ContentType','vector'); close(f);
end

function write_tables(out,comparisons,records,timings)
fid=fopen(fullfile(out,'grid_rows.tex'),'w');
labels={'Linear','Single correction','Interpolation, $M=4$','Interpolation, $M=16$','Interpolation, $M=64$','Exact rollout'};
for h=1:6
    fprintf(fid,'%s & %.6f & %.4f & %.6f & %.4f \\\\ \n',labels{h},comparisons{h,3},comparisons{h,4},comparisons{h+6,3},comparisons{h+6,4});
end
fprintf(fid,'\\hline\n'); fclose(fid);
fid=fopen(fullfile(out,'scaling_rows.tex'),'w');
for n=[2,4,8,16]
    for h=[2,3,6]
        name={'linear','one_state','interp4','interp16','interp64','rollout'};
        ix=find(cellfun(@(x)x==.22,records(:,1)) & cellfun(@(x)x==n,records(:,2)) & strcmp(records(:,3),name{h}));
        fprintf(fid,'%d & %s & %.3f & %.3f & %.3f \\\\ \n',n,labels{h},records{ix,4},records{ix,5},records{ix,6});
    end
end
fprintf(fid,'\\hline\n'); fclose(fid);
fid=fopen(fullfile(out,'timing_rows.tex'),'w');
for n=[2,4,8,16]
    vals=zeros(1,4); names={'one_state','interp4','rollout','rollout_linprog'};
    for k=1:4
        ix=find(cellfun(@(x)x==.22,timings(:,1)) & cellfun(@(x)x==n,timings(:,2)) & strcmp(timings(:,3),names{k}));
        vals(k)=timings{ix,4};
    end
    fprintf(fid,'%d & %.2f & %.2f & %.2f & %.2f \\\\ \n',n,vals);
end
fprintf(fid,'\\hline\n'); fclose(fid);
end

function [A,B,C,d,K]=data(n,tau)
q=n/2; A=.8*eye(n); C=zeros(n,2*n); D=C; d=zeros(2*n,1); K=zeros(2*n,n);
for b=1:q
    ix=2*b-1:2*b; ju=4*b-3:4*b;
    C(ix,ju)=[eye(2),[1;1],[1;2]];
    D(ix,ju)=[.2*eye(2),[tau;2*tau-.2],[0;0]];
    d(ju)=[24/5;24/5;83/15-(74/3)*tau;28/3];
    K(ju(1:2),ix)=eye(2);
end
B=D-A*C;
assert(all(d>=0)); assert(all(sum(D,1)<=.8*sum(C,1)+1e-12));
end

function v=pair_rollout(x)
if x(2)<=x(1), v=[x(2);0];
elseif x(2)<=2*x(1), v=[2*x(1)-x(2);x(2)-x(1)];
else, v=[0;x(1)]; end
end

function u=policy(x,h,tables)
n=numel(x); u=zeros(2*n,1);
if h==2, theta=min(x); end
for b=1:n/2
    ix=2*b-1:2*b; ju=4*b-3:4*b; z=x(ix);
    if h==1, v=[0;0];
    elseif h==2, v=[theta;0];
    elseif h==6, v=pair_rollout(z);
    elseif sum(z)==0, v=[0;0];
    else
        tab=tables{h-2}; a=tab.M*z(1)/sum(z);
        k=min(floor(a),tab.M-1); lambda=a-k;
        v=sum(z)*((1-lambda)*tab.values(:,k+1)+lambda*tab.values(:,k+2));
    end
    u(ju)=[z-[1,1;1,2]*v;v];
end
end

function value=cost(x,h,tables,A,B,d,p,tol)
discount=1; value=0;
while discount*p*sum(x)>tol
    u=policy(x,h,tables); value=value+discount*(10*sum(x)+d'*u);
    x=max(A*x+B*u,0); discount=discount*.5;
end
end

function [value,width]=finite_lower(x0,A,B,C,d,N,opts)
% States x_0,...,x_N followed by u_0,...,u_(N-1).
n=numel(x0); m=numel(d); nv=n*(N+1)+m*N;
f=zeros(nv,1); E=sparse(n*(N+1),nv); rhs=zeros(n*(N+1),1);
H=sparse(n*N,nv); E(1:n,1:n)=eye(n); rhs(1:n)=x0;
for t=0:N-1
    xi=n*t+(1:n); xp=n*(t+1)+(1:n); ui=n*(N+1)+m*t+(1:m);
    f(xi)=.5^t*10; f(ui)=.5^t*d;
    E(xp,xp)=eye(n); E(xp,xi)=-A; E(xp,ui)=-B;
    H(xi,xi)=-eye(n); H(xi,ui)=C;
end
[~,primal,flag,~,lambda]=linprog(f,H,zeros(n*N,1),E,rhs,zeros(nv,1),[],opts);
assert(flag>0);
% Correct negative dual slacks using known mass bounds for every variable.
% This avoids treating a solver's primal objective as a certified lower bound.
ineq=max(lambda.ineqlin,0); slack=f+H'*ineq+E'*lambda.eqlin;
penalty=0;
for t=0:N
    xi=n*t+(1:n); penalty=penalty+.8^t*max([0;-slack(xi)]);
    if t<N
        ui=n*(N+1)+m*t+(1:m);
        penalty=penalty+.8^t*max([0;-slack(ui)]);
    end
end
value=-lambda.eqlin'*rhs-sum(x0)*penalty-1e-10;
width=primal-value+(148/9)*.4^N*sum(x0);
assert(width>=-1e-8 && width<1e-5);
end
