function GreensFunctionOptimizationBulk

%
% code for minimization of p = [\sum_i^N R(x_i;x_i) + sum_{j neq i}^N
% G(x_j;x_i)]/N^2
%
%
close all;

% These two lines activate chunkIE and are only necessary the first time.
%homeDir = getenv('HOME');
%run([homeDir '/chunkie/startup.m'])

cparams = []; cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 0.25; % setting a chunk length helps when the
% frequency is known
pref = [];
pref.k = 16;

cols = get(gca,'colororder');

%o.example = 'disk';
%o.example = 'cassini';
%o.example = 'ellipse';
%o.example = 'barbell';
o.example = 'random';

switch o.example
    case 'barbell'
        verts = chnk.demo.barbell(2.0,2.0,1.0,1.0); % vertices of a barbell
        chnkr = chunkerpoly(verts,cparams,pref);
        filename = 'barbellDomain';

    case 'ellipse'

        %% Ellipse domain
        %a = 1.5;  b=1/a;
        o.a = 1; kap = 0.802; o.b = sqrt(1-kap^2);
        chnkr = chunkerfunc(@(t) ellipse(t,o.a,o.b),cparams,pref);

        filename = 'ellipseDomain';

    case 'disk'

        %% Disk domain
        o.a = 1; o.b = 1;
        chnkr = chunkerfunc(@(t) ellipse(t,o.a,o.b),cparams,pref);
        filename = 'DiskDomain';

    case 'random'

        % % Random domain
        rng(3);
        modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end)));
        ctr = [0.0;0.0];
        chnkr = chunkerfunc(@(t) chnk.curves.bymode(t,modes,ctr));
        filename = 'randomDomain';

        rng('shuffle');
    case 'starfish'
        %% Starfish domain
        narms = 5;
        amp = 0.25;
        chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);
        filename = 'starDomain';

    case 'cassini'
        k = 0.99;
        area = pi;
        Em = ellipticE(k^4);
        b = sqrt(area/Em); a = k*b;
        %a = 1; b = 10;
        chnkr = chunkerfunc(@(t) cassini_oval(t, a, b),cparams,pref);

        hold on;
        plot(chnkr.r(1,:), chnkr.r(2,:));
        hold off;
        axis equal;

        % return
end

export_figs = 0;

N_rep = 5;

o.Sk   = kernel('laplace', 's');
o.Skp  = kernel('laplace', 'sprime');
o.Kgrad = kernel.lap2d('sgrad');

Smat  = chunkermat(chnkr,o.Sk);
Spmat = 1/2*eye(chnkr.npt)+chunkermat(chnkr,o.Skp);

%%% Preliminary functions. These quantities are the same for all solves.

o.xx  = chnkr.r(1,:);   o.yy  = chnkr.r(2,:);
o.nx  = chnkr.n(1,:);   o.ny  = chnkr.n(2,:);

o.vb  = (o.xx.^2+o.yy.^2)/4;    o.vbn = (o.xx.*o.nx+o.yy.*o.ny)/2;

o.area = o.vbn*chnkr.wts(:);
o.v3bn = (o.xx.^3.*o.nx+o.yy.^3.*o.ny)/12;

o.s0 = o.v3bn*chnkr.wts(:)/o.area;
o.s1 = (o.vb.*o.vbn)*chnkr.wts(:)/o.area;

vvec = o.vbn.*(chnkr.wts(:).');
corrmat = ones(chnkr.npt,1)*vvec*Smat;
o.syscorr = Spmat + corrmat/o.area + ones(chnkr.npt,1)*chnkr.wts(:).';
[o.L,o.U,o.P] = lu(o.syscorr);

% Criteria for fmincon

opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',true, ...
    'HessianApproximation','lbfgs', ...'finite-difference'
    'MaxFunctionEvaluations', 5e5, ...
    'MaxIterations', 30, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    'Display','off');
    if (nnz(strcmp(o.example,{'Disk','Ellipse'})>0) )
        optimoptions('fmincon','SpecifyConstraintGradient',true);  % if nonlcon provides grads
        cons = @(x) nonlocDE(x,chnkr,o);
    else
        optimoptions('fmincon','SpecifyConstraintGradient',false);  % if nonlcon provides grads
        cons = @(x) nonloc(x,chnkr,o);
    end

M_range = [1:2:20];
f_min = zeros(1,length(M_range));
x_out = []; y_out = [];

for i = 1:length(M_range)
    o.M = M_range(i);
    Mr = M_range(i);

    out_curr = zeros(N_rep,2*Mr); curr_min = zeros(1,N_rep);

    for j = 1:N_rep
   %parfor j = 1:N_rep

        if (nnz(strcmp(o.example,{'cassini','barbell'})>0))
            x0 = [0.8*(rand(1,Mr)-0.5) 0.2*(rand(1,Mr)-0.5)];
        else
            x0 = [(rand(1,Mr)-0.5) 0.5*(rand(1,Mr)-0.5)];
        end
        %valid = checkGradients(@(x) getpNEW(x,o,chnkr),x0,Display="on")
        out_curr(j,:) = fmincon(@(x) getpNEW(x,o,chnkr),x0,[],[],[],[],[],[],cons,opts);
        curr_min(j) = getpNEW(out_curr(j,:),o,chnkr);

    end
    [~,i_min] = min(curr_min);
    out = out_curr(i_min,:);
    x_out = out(1:Mr); y_out = out(Mr+1:2*Mr);
    f_min(i) = Mr^2 *getpNEW(out,o,chnkr);

    figure(Mr);
    set(gcf,'color','w')
    plot(o.xx,o.yy,'k','linewidth',3)
    hold on;
    plot(x_out,y_out,'.','markersize',32,'markerfacecolor',cols(1,:));
    hold off;
    axis equal;
    xtickformat('%1.1f');
    ytickformat('%1.1f');
    set(gca,'fontsize',16);
    title(['$N = $ ' num2str(Mr,'%g') ' ($p = $ ' num2str(f_min(i),'%5.4f') ')' ],'interpreter','latex','FontSize',24)
    axis off; drawnow; axis tight equal; hold off;
    if (export_figs)
        exportgraphics(gca,[filename 'M=' num2str(o.M) '.png'],'Resolution',300);
    end

    fprintf(1,'N = %g \t p = %5.4f\n', o.M, f_min(i));
end
%
figure('color','w')
plot(M_range, f_min,'-ko')

end


function [out, J_out] =  getpNEW(x_in,o,chnkr)
x = x_in(1:o.M);  y = x_in((o.M+1):2*o.M);

dx = x - x'; dy = y - y';
d2 = dx.^2 + dy.^2;
d2(1:size(d2,1)+1:end) = 1;

d = sqrt( d2 );

c = -1/(2*pi);

gradGx = c*(dx./d2);    gradGy = c*(dy./d2);

out = c*sum(log(d(:)));

targ = [];  targ.r = [x(:)';y(:)'];

rhs_matrix = zeros(size(o.syscorr,1),o.M);

for i = 1:o.M
    src   = {};
    % Opimization points

    src.r = [x(i); y(i)];
    %%%
    t1 = -(src.r(1)^2+src.r(2)^2)/4;
    t2  = (o.vbn*(chnk.lap2d.kern(src,chnkr,'s').*(chnkr.wts(:))));

    rhs_matrix(:,i) = -chnk.lap2d.kern(src,chnkr,'sprime') ...
        - (o.vbn.')/o.area - ones(chnkr.npt,1)*(o.s0+o.s1+t1+t2)/o.area;
end

%sig = o.syscorr\rhs_matrix;
Y = o.L \ (o.P*rhs_matrix);
sig = o.U \ Y;

Jx = zeros(o.M,o.M);
Jy = zeros(o.M,o.M);

for i = 1:o.M
    wbar = sum(sig(:,i).*chnkr.wts(:));
    R = chunkerkerneval(chnkr,o.Sk,sig(:,i),targ) + (x'.^2 + y'.^2)/(4*o.area) + wbar;
    out = out + sum(R);

    % Get Jacobian

    J = chunkerkerneval(chnkr, o.Kgrad, sig(:,i), targ);
    gradU = reshape(J, 2, []);

    %gradU = 2xM;
    Jx(:,i) = gradU(1,:)' + x'/(2*o.area);
    Jy(:,i) = gradU(2,:)' + y'/(2*o.area);
end

out = out/o.M^2;
J_out = 2*[sum(gradGx + Jx') sum(gradGy + Jy')]/o.M^2;

end


function [c,ceq] = nonloc(xin,chnkr,o)
M = length(xin)/2;
x = xin(1:M); y = xin(M+1:end);
pts = [(x(:)).'; (y(:)).'];
in = chunkerinterior(chnkr,pts);

c=[];
if all(in)
    ceq = 0;
else
    ceq = 1e6;
end

end

function [c,ceq,gradc, gradceq] = nonlocDE(xin,chnkr,o)
M = length(xin)/2;
x = xin(1:M); y = xin(M+1:end);
ceq = []; gradceq = [];

c = (x/o.a).^2 + (y/o.b).^2 - 1;
gradc = [2*x/(o.a^2); 2*y/(o.b^2)];

end
