function assess_polynomial_exactness
% Numerical illustration of the proved polynomial obstruction.
% Run: matlab -batch "assess_polynomial_exactness"
% Requires Optimization Toolbox. Bounds use convex first-order inequalities;
% their floating-point evaluations are not interval-arithmetic certificates.
out = fullfile(fileparts(mfilename('fullpath')), 'polynomial_results');
if ~exist(out, 'dir'), mkdir(out); end
gamma = 0.9; rho = 1; horizon = 18;
kappa = fzero(@(k) rho*k^3*((1-k)^(-3)-gamma^2) ...
    - gamma*(1+gamma), [1e-8, 1-1e-8]);
decay = 1-kappa;
pa = 1+rho*kappa^3; pc = 1+gamma*pa;
Ah = gamma^(1/4)*[0 1; -1 0]; Bh = gamma^(1/4)*[0;1];
S = zeros(2,horizon,horizon+1); T = zeros(2,2,horizon+1);
T(:,:,1) = eye(2);
for t=1:horizon
    S(:,:,t+1) = Ah*S(:,:,t); S(:,t,t+1) = Bh;
    T(:,:,t+1) = Ah*T(:,:,t);
end
options = optimoptions('fminunc','Algorithm','trust-region', ...
    'SpecifyObjectiveGradient',true,'HessianFcn','objective', ...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-14, ...
    'MaxIterations',500,'MaxFunctionEvaluations',2000,'Display','off');
ratios = [0 .05 .1 .2 .3 .5 .75 1 1.5 2 3 5 10 20 Inf];
etas = [0 .2 2]; rows = []; maxGradientCheck = 0; maxHessianCheck = 0;
for eta=etas
    pb = eta/(1-gamma*decay^2);
    for ratio=ratios
        if isinf(ratio), x=[0;1]; else, x=[1;ratio]/sqrt(1+ratio^2); end
        v0 = zeros(horizon,1); xx=x;
        for t=1:horizon
            v0(t)=kappa*xx(1); xx=Ah*xx+Bh*v0(t);
        end
        fun = @(v) objective(v,x,eta,rho,pa,pc,S,T);
        if eta==2 && ratio==1
            trial=v0+.01*sin((1:horizon)');
            [~,gg,hh]=fun(trial); fdg=zeros(horizon,1); fdh=zeros(horizon);
            for j=1:horizon
                ee=zeros(horizon,1); ee(j)=1e-6;
                [fp,gp]=fun(trial+ee); [fm,gm]=fun(trial-ee);
                fdg(j)=(fp-fm)/2e-6; fdh(:,j)=(gp-gm)/2e-6;
            end
            maxGradientCheck=max(abs(fdg-gg));
            maxHessianCheck=max(abs(fdh-hh),[],'all');
        end
        [v,~,exitflag,info] = fminunc(fun,v0,options);
        [f,g,~,last] = fun(v);
        initialUpper=pa*x(1)^4+pb*x(1)^2*x(2)^2+pc*x(2)^4;
        radius=(initialUpper/rho)^(1/4);
        correction=sum(abs(g).*(radius+abs(v)));
        buffer=1e-9;
        lower=f-correction-buffer;
        upper=f+pb*last(1)^2*last(2)^2+buffer;
        baseline=pa*x(1)^4+pc*x(2)^4;
        mixed=x(1)^2*x(2)^2;
        if mixed>0
            effectiveLow=(lower-baseline)/mixed;
            effectiveHigh=(upper-baseline)/mixed;
        else
            effectiveLow=NaN; effectiveHigh=NaN;
        end
        rows(end+1,:)=[eta,ratio,x',lower,upper,upper-lower, ...
            effectiveLow,effectiveHigh,v(1),norm(g,Inf),exitflag,info.iterations]; %#ok<AGROW>
    end
end
assert(maxGradientCheck<1e-7 && maxHessianCheck<1e-6);
names={'eta','ratio','x1','x2','lower','upper','width', ...
    'mixed_coefficient_lower','mixed_coefficient_upper','first_action', ...
    'gradient_norm','exitflag','iterations'};
result=array2table(rows,'VariableNames',names);
writetable(result,fullfile(out,'values.csv'));
baseRows=result.eta==0;
baseExact=pa*result.x1(baseRows).^4+pc*result.x2(baseRows).^4;
assert(all(result.lower(baseRows)<=baseExact+1e-10));
assert(all(result.upper(baseRows)>=baseExact-1e-10));
summary=struct('gamma',gamma,'rho',rho,'horizon',horizon, ...
    'kappa',kappa,'axis_coefficient_x1',pa,'axis_coefficient_x2',pc, ...
    'max_bound_width',max(result.width),'max_gradient_norm',max(result.gradient_norm), ...
    'gradient_check',maxGradientCheck,'hessian_check',maxHessianCheck, ...
    'matlab_version',version,'platform',computer);
fid=fopen(fullfile(out,'checks.json'),'w'); fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true)); fclose(fid);
save(fullfile(out,'results.mat'),'result','summary');
fig=figure('Visible','off','Color','w','Position',[100 100 700 420]);
rr=result.eta==2 & isfinite(result.ratio) & result.ratio>0;
semilogx(result.ratio(rr),.5*(result.mixed_coefficient_lower(rr)+ ...
    result.mixed_coefficient_upper(rr)),'o-','LineWidth',1.3);
xlabel('State ratio x_2 / x_1'); ylabel('Effective mixed coefficient'); grid on;
exportgraphics(fig,fullfile(out,'mixed_coefficient.pdf'),'ContentType','vector'); close(fig);
fid=fopen(fullfile(out,'numerical_rows.tex'),'w');
for ratio=[.1 .3 1 3 10]
    z=result(result.eta==2 & result.ratio==ratio,:);
    fprintf(fid,'%.1f & %.6f & %.6f & %.2e \\\\\n',ratio, ...
        z.mixed_coefficient_lower,z.mixed_coefficient_upper,z.width);
end
fprintf(fid,'\\hline\n'); fclose(fid);
disp(summary); disp(result(result.eta==2,[2 8 9 7 12]));
end

function [f,g,H,last] = objective(v,x,eta,rho,pa,pc,S,T)
N=numel(v); f=rho*sum(v.^4); g=4*rho*v.^3; H=diag(12*rho*v.^2);
for t=1:N+1
    xx=T(:,:,t)*x+S(:,:,t)*v; a=xx(1); b=xx(2);
    if t<=N
        f=f+a^4+b^4+eta*a^2*b^2;
        gx=[4*a^3+2*eta*a*b^2;4*b^3+2*eta*a^2*b];
        hx=[12*a^2+2*eta*b^2,4*eta*a*b;4*eta*a*b,12*b^2+2*eta*a^2];
    else
        f=f+pa*a^4+pc*b^4;
        gx=[4*pa*a^3;4*pc*b^3]; hx=diag([12*pa*a^2,12*pc*b^2]);
    end
    g=g+S(:,:,t)'*gx; H=H+S(:,:,t)'*hx*S(:,:,t);
end
last=xx;
end
