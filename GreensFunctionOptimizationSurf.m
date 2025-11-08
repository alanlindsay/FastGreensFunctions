function GreensFunctionOptimizationSurf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       .   .   .   Green's function optimization 2D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

homeDir = getenv('HOME');
run([homeDir '/chunkie/startup.m'])

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 0.4; % setting a chunk length helps when the
% frequency is known
pref = [];
pref.k = 16;

cols = get(gca,'colororder');

%o.example = 'disk';
%o.example = 'ellipse';
o.example = 'cassini';
%o.example = 'barbell';
o.example = 'random';

switch o.example
    case 'barbell'
        verts = chnk.demo.barbell(2.0,2.0,1.0,1.0); % vertices of a barbell
        chnkr = chunkerpoly(verts,cparams,pref);
        filename = 'barbellDomain';

    case 'ellipse'

        %% Ellipse domain
        a = 2;  b=1/a;
        %a = 1.2; kap = 0.802; b = sqrt(1-kap^2);
        o.shape = @(t) ellipse(t,a,b);

        chnkr = chunkerfunc(o.shape,cparams,pref);

        filename = 'Ellipse_a=2_BDY_';

    case 'disk'

        %% Disk domain
        o.shape = @(t) ellipse(t,1,1);
        chnkr = chunkerfunc(o.shape,cparams,pref);
        filename = 'DiskDomain';

    case 'random'

        % % Random domain
        rng(3);
        modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end)));
        ctr = [0.0;0.0];
        o.shape = @(t) chnk.curves.bymode(t,modes,ctr);
        chnkr = chunkerfunc(o.shape);
        filename = 'randomDomainSurf';

    case 'starfish'
        %% Starfish domain
        narms = 5;
        amp = 0.25;
        chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);
        filename = 'starDomain';

    case 'cassini'
        % Cassini boundary (area scaled to pi)
        k = 0.99; b = sqrt(pi/ellipticE(k^4)); a = k*b;
        o.shape = @(t) cassini_oval(t, a, b);
        chnkr = chunkerfunc(o.shape,cparams,pref);

        filename = 'Cassini';
end

export_figs = 0;

o.Sk   = kernel('laplace', 's');
o.Sk.sing = 'log';
o.Skp  = kernel('laplace', 'sprime');
o.Kgrad = kernel.lap2d('sgrad');

o.Smat  = chunkermat(chnkr,o.Sk);
Spmat = 1/2*eye(chnkr.npt)+chunkermat(chnkr,o.Skp);

%%% Preliminary functions. These quantities are the same for all solves.

o.xx  = chnkr.r(1,:);   o.yy  = chnkr.r(2,:);
o.nx  = chnkr.n(1,:);   o.ny  = chnkr.n(2,:);

%plot(o.xx,o.yy); return;

o.vb  = (o.xx.^2+o.yy.^2)/4;    o.vbn = (o.xx.*o.nx+o.yy.*o.ny)/2;

o.area = o.vbn*chnkr.wts(:);
o.v3bn = (o.xx.^3.*o.nx+o.yy.^3.*o.ny)/12;

%%%
o.s0 = o.v3bn*chnkr.wts(:)/o.area;
o.s1 = (o.vb.*o.vbn)*chnkr.wts(:)/o.area;
%%%
vvec = o.vbn.*(chnkr.wts(:).');
corrmat = ones(chnkr.npt,1)*vvec*o.Smat;
o.syscorr = Spmat + corrmat/o.area + ones(chnkr.npt,1)*chnkr.wts(:).';

opts = optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'SpecifyObjectiveGradient',false, ...
    'SpecifyConstraintGradient',false, ...  %
    'HessianApproximation','lbfgs', ...
    'MaxFunctionEvaluations', 5e5, ...
    'MaxIterations', 20, ...
    'OptimalityTolerance', 1e-5, ...
    'StepTolerance', 1e-5, ...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",...
    'Display','off');

N_rep = 24;
M_range = 1:12;
f_min = zeros(1,length(M_range));

for i = 1:length(M_range)
    o.M = M_range(i);
    Mr = o.M;

    out_curr = zeros(N_rep,o.M); curr_min = zeros(1,N_rep);
    for j = 1:N_rep
    %parfor j = 1:N_rep
        t0= 2*pi*rand(1,Mr);

        out_curr(j,:) = fmincon(@(x) getpNEW(x,o,chnkr),t0,[],[],[],[],zeros(size(t0)),2*pi*ones(size(t0)) ,[],opts);
        curr_min(j) = getpNEW(out_curr(j,:),o,chnkr);
    end

    [~,i_min] = min(curr_min);
    out = out_curr(i_min,:);
    f_min(i) = getpNEW(out,o,chnkr);

    figure(o.M);
    set(gcf,'color','w')
    plot(o.xx,o.yy,'k','linewidth',3)
    hold on;
    r = o.shape(out);
    plot(r(1,:),r(2,:),'.','markersize',64,'markerfacecolor',cols(1,:));
    hold off;
    axis equal;
    title(['$N = $ ' num2str(Mr,'%g') ' ($p = $ ' num2str(f_min(i),'%5.4f') ')' ],'interpreter','latex','FontSize',32)
    xtickformat('%1.1f');
    ytickformat('%1.1f');
    set(gca,'fontsize',16);
    axis off
    drawnow;
    axis tight equal
    if (export_figs)
        exportgraphics(gca,[filename 'M=' num2str(o.M) '.png'],'Resolution',300);
    end
    fprintf(1,'N = %g \t p = %g\n', o.M, f_min(i));

    %pause
end

%figure('color','w')
%plot(M_range,f_min,'-ko')

end


function [out] =  getpNEW(t_in,o,chnkr)

[r,~] = o.shape(t_in);

x = r(1,:);  y = r(2,:);

d2 = (x - x').^2 + (y - y').^2;
d2(1:size(d2,1)+1:end) = 1;

c = -1/pi;

out = c*sum(log(sqrt(d2(:))));

targ = [];  targ.r = [x(:)';y(:)'];

rhs_matrix = zeros(size(o.syscorr,1),o.M);

for i = 1:o.M
    src   = {};
    % Opimization points

    src.r = [x(i); y(i)];
    %%%
    t1 = -(src.r(1)^2+src.r(2)^2)/4;
    t2 = 2*chunkerkerneval(chnkr,o.Sk,o.vbn,src);

    rhs_matrix(:,i) = -2*chnk.lap2d.kern(src,chnkr,'sprime') ...
        - (o.vbn.')/o.area - ones(chnkr.npt,1)*(o.s0+o.s1+t1+t2)/o.area;

end

sig = o.syscorr\rhs_matrix;

for i = 1:o.M
    wbar = sum(sig(:,i).*chnkr.wts(:));
    R = chunkerkerneval(chnkr,o.Sk,sig(:,i),targ) + (x'.^2 + y'.^2)/(4*o.area) + wbar;
    out = out + sum(R);
end

end


