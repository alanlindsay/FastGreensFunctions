function [R,gradR,HessR] = NeumR_2D(x,y,chnkr,opts)

%homeDir = getenv('HOME');
%run([homeDir '/chunkie/startup.m'])

% Regular part of 2D Neumann Green's function for Laplace's equation.
% inputs:
% x - 1 x 2*N_e vector of N_e evaluation points. Stored x = [x_1 x_2]
% y - 1 x 2*N_s vector of N_s source points. Stored y = [y_1 y_2]
% chnkr - a chunkIE object defining the geometry.
% opts.Int - 1 or 0 for interior problem
% opts.Surf - 1 or 0 for surface problem.

% outputs:
% R - a Ne x Ns matrix for the regular part.
% gradR - a 2x1 cell {Rx,Ry} for the partial derivatives of R. Rx and Ry are Ne x Ns
% matricies.
% HessR - a 4x1 cell {Rxx,Rxy,Ryx,Ryy} for the second derivatives matrices
% of R.

% Derivatives may not be accurate for on - surface evaluation (x on ∂Ω ) 

% If code run with no outputs/inputs, a test of each 4 Green's functions is
% produced.

if (nargin == 0)
    close all;
    cparams = [];   cparams.eps = 1.0e-6;
    cparams.nover = 0;  cparams.maxchunklen = 0.5; % setting a chunk length helps when the
    % frequency is known

    pref = [];    pref.k = 16;

    %narms
    % =10;  amp = 0.5;
    %chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);

    % a = 2; b = 1/a;
    % t0 = pi/2;
    % chnkr = chunkerfunc(@(t) ellipse(t,a,b),cparams,pref);
    % [y,yp] = ellipse(t0,a,b); x = y;

    % Plot four examples from Fig 1 of paper

    %Surface sources

    narms = 10;   amp = 0.5; t0 = pi/4; y = starfish(t0,narms,amp);
    start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);

    opts.Surf = 1; opts.Int = 1; NeumR_2D(y,y,chnkr,opts);
    opts.Surf = 1; opts.Int = 0; NeumR_2D(y,y,chnkr,opts);

    % Bulk problems.

    k = 0.99; b = sqrt(pi/ellipticE(k^4)); a = k*b;
    chnkr = chunkerfunc(@(t) cassini_oval(t, a, b),cparams,pref);
    y = [0 0.5]; opts.Surf = 0; opts.Int = 0; NeumR_2D(y,y,chnkr,opts);
    y = [0 0.1]; opts.Surf = 0; opts.Int = 1; NeumR_2D(y,y,chnkr,opts);

    %  a = 1; b = 1/a;
    %  t0 = pi/2;
    %  chnkr = chunkerfunc(@(t) ellipse(t,a,b),cparams,pref);
    %  [y,yp] = ellipse(t0,a,b); x = y;
    %
    % opts.Surf = 1; opts.Int = 1; NeumR_2D(y,y,chnkr,opts);
    % opts.Surf = 1; opts.Int = 0; NeumR_2D(y,y,chnkr,opts);

    return
end

Sk = kernel('laplace', 's');
Sk.sing = 'log';
Skp = kernel('laplace', 'sprime');
Kgrad = kernel.lap2d('sgrad');
Kgrad.sing = 'pv';

% Build the kernel object for Hessian
KH = kernel();                  % empty kernel object
KH.eval   = @hesslap_s_eval;    % point to evaluator
KH.opdims = [4 1];              % 4 components at target, scalar at source
KH.sing   = 'hs';           % off-surface target evals are smooth

xx  = chnkr.r(1,:); yy  = chnkr.r(2,:);
nx  = chnkr.n(1,:); ny  = chnkr.n(2,:);

if (opts.Int == 1)

    Smat  = chunkermat(chnkr,Sk);
    Spmat = 1/2*eye(chnkr.npt)+chunkermat(chnkr,Skp);
    syscorr = Spmat;

    vb  = (xx.^2+yy.^2)/4;  vbn = (xx.*nx+yy.*ny)/2;

    area = vbn*chnkr.wts(:);
    v3bn = (xx.^3.*nx+yy.^3.*ny)/12;

    s0 = v3bn*chnkr.wts(:)/area;
    s1 = (vb.*vbn)*chnkr.wts(:)/area;

    vvec = vbn.*(chnkr.wts(:).');
    corrmat = ones(chnkr.npt,1)*vvec*Smat;
    syscorr = syscorr + corrmat/area + ones(chnkr.npt,1)*chnkr.wts(:).';

else

    %Smat  = chunkermat(chnkr,Sk);
    Spmat = -1/2*eye(chnkr.npt)+chunkermat(chnkr,Skp);
    syscorr = Spmat;

end

%%%%%%%%%%%%%%%

% number of sources (y_in)

Ns = numel(y)/2;  xs = y(1:Ns);  ys = y((Ns+1):2*Ns);

% number of evaluation points (x_in)

Ne = numel(x)/2;  xe = x(1:Ne);  ye = x((Ne+1):2*Ne);

% Solution output;
if (nargout >= 1)
    R = zeros(Ne,Ns);
end

% Derivatives will only accurate for off surface evaluation (x off
% surface).
% Gradient output;
if (nargout >= 2)
    gradR.dx = zeros(Ne,Ns);
    gradR.dy = zeros(Ne,Ns);
end

if (nargout == 3)
    % Hessian output;
    HessR.dxx = zeros(Ne,Ns);
    HessR.dxy = zeros(Ne,Ns);
    HessR.dyx = zeros(Ne,Ns);
    HessR.dyy = zeros(Ne,Ns);
end

rhs_mat = zeros(chnkr.npt,Ns);

for i = 1:Ns
    src.r = [xs(i); ys(i)];

    if (opts.Int==1 && opts.Surf==1)

        t1 = -(src.r(1)^2+src.r(2)^2)/4;
        t2b = 2*chunkerkerneval(chnkr,Sk,vbn,src);

        rhs_mat(:,i) = - (vbn.')/area ...
            - 2*chnk.lap2d.kern(src,chnkr,'sprime') ...
            - ones(chnkr.npt,1)*(s0+s1+t1+t2b)/area;

    elseif (opts.Int==1 && opts.Surf==0)

        t1 = -(src.r(1)^2+src.r(2)^2)/4;       
        t2 = chunkerkerneval(chnkr,Sk,vbn,src);

        rhs_mat(:,i) = - (vbn.')/area ...
            - chnk.lap2d.kern(src,chnkr,'sprime') ...
            - ones(chnkr.npt,1)*(s0+s1+t1+t2)/area;

    elseif (opts.Int==0 && opts.Surf==1)
        rhs_mat(:,i) = -2*chnk.lap2d.kern(src,chnkr,'sprime');

    elseif (opts.Int==0 && opts.Surf==0)
        rhs_mat(:,i) = -chnk.lap2d.kern(src,chnkr,'sprime');
    end

end

sig = syscorr\rhs_mat;

if (nargout == 0)
    % This is just run for the diagnostic plot
    xv = min(xx):0.02:max(xx);
    yv = min(yy):0.02:max(yy);

    [X,Y] = ndgrid(1.5*xv,1.5*yv);
    sz = size(X);   X = X(:);   Y = Y(:);

    targ = [];  targ.r = [X.';Y.'];

    [in] = chunkerinterior(chnkr,targ);
    [fints] = chunkerkerneval(chnkr,Sk,sig,targ);

    if (opts.Surf == 1)
        fints = fints + 2*chnk.lap2d.kern(src,targ,'s');
    else
        fints = fints + chnk.lap2d.kern(src,targ,'s');
    end

    if (opts.Int == 1)
        fints = fints + sum(sig(:,i).*chnkr.wts(:)) + (X.^2+Y.^2)/(4*area);
        fints(~in) = NaN;
    else
        fints(in) = NaN;
    end

    Fmat = reshape(fints,sz);
    Xmat = reshape(X,sz);   Ymat = reshape(Y,sz);

    % if opts.Int == 1
    %     ex = (-1/2 +Xmat.^2 + Ymat.^2 )/(4*pi);
    %     pcolor(Xmat,Ymat,log10(abs(Fmat-ex))); shading interp;
    %     c = colorbar;
    %     c.TickLabelInterpreter = 'latex';
    %     c.Ticks = -13:1:-12; clim([-13 -12]);
    %     c.TickLabels = {'$10^{-13}$','$10^{-12}$'};
    % else
    %     ex = log(sqrt(Xmat.^2 + Ymat.^2))/(2*pi);
    %     pcolor(Xmat,Ymat,log10(abs(Fmat-ex))); shading interp;
    %     c = colorbar;
    %     c.TickLabelInterpreter = 'latex';
    %     c.Ticks = -15:1:-12; clim([-15 -12]);
    %     c.TickLabels = {'$10^{-15}$','$10^{-14}$','$10^{-13}$','$10^{-12}$','$10^{-11}$'};
    % end
    % hold on;

    figure(); hold on;
    pcolor(Xmat,Ymat,log10(abs(Fmat))); shading interp;
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Ticks = -3:1:0; clim([-3 0]);
    c.TickLabels = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'};
    c.FontSize = 24;

    plot(xx,yy,'-k','linewidth',3);

    set(gcf,'color','w');

    plot(y(1),y(2),'r.','markersize',32)
    axis off;   hold off; axis equal;
    drawnow;

    return
end

targ = [];  targ.r = [xe(:)'; ye(:)'];

for i = 1:Ns
    % Evaluate Solution
    if opts.Int == 1

        wbar = sum(sig(:,i).*chnkr.wts(:));
        R(:,i) = chunkerkerneval(chnkr,Sk,sig(:,i),[xe(:)';ye(:)']) + ...
            (xe(:).^2+ye(:).^2)/(4*area) + wbar;

        % Get Gradients
        if (nargout >= 2)
            J = chunkerkerneval(chnkr, Kgrad, sig(:,i), targ);
            gradU = reshape(J, 2, []);
            gradR.dx(:,i) = gradU(1,:)' + xe(:)/(2*area);
            gradR.dy(:,i) = gradU(2,:)' + ye(:)/(2*area);
        end

        if (nargout == 3)
            % Get Hessian
            H = chunkerkerneval(chnkr, KH, sig(:,i), targ);
            HessR.dxx(:,i) = 1/(2*area) + H(1:Ne);
            HessR.dxy(:,i) = H((Ne+1):2*Ne);
            HessR.dyx(:,i) = H((2*Ne+1):3*Ne);
            HessR.dyy(:,i) = 1/(2*area)+ H((3*Ne+1):4*Ne);
        end


    else

        R(:,i) = chunkerkerneval(chnkr,Sk,sig(:,i),targ);

        % Get Jacobian
        if (nargout >= 2)
            J = chunkerkerneval(chnkr, Kgrad, sig(:,i), targ);
            gradU = reshape(J, 2, []);
            gradR.dx(:,i) = gradU(1,:)';
            gradR.dy(:,i) = gradU(2,:)';
        end

        % Get Hessian
        if (nargout == 3)
            H = chunkerkerneval(chnkr, KH, sig(:,i), targ);
            HessR.Rxx(:,i) = H(1:Ne);
            HessR.Rxy(:,i) = H((Ne+1):2*Ne);
            HessR.Ryx(:,i) = H((2*Ne+1):3*Ne);
            HessR.Ryy(:,i) = H((3*Ne+1):4*Ne);
        end

    end

end


end


function out = hesslap_s_eval(s, t)
%HESSLAP_S_EVAL  Hessian of 2D Laplace single-layer kernel wrt target x
%  For sources y=s.r and targets x=t.r:
%     G(x,y) = -(1/(2*pi)) log |x-y|
%     d_i d_j G(x,y) = -(1/(2*pi)) * [ \delta_ij/|r|^2 - 2 r_i r_j / |r|^4 ],
%  where r = x - y.

xs = s.r(1,:);  ys = s.r(2,:);
xt = t.r(1,:);  yt = t.r(2,:);

DX = xt.' - xs;      % nt x ns
DY = yt.' - ys;      % nt x ns
R2 = DX.^2 + DY.^2;  % nt x ns
R4 = R2.^2;

c  = -1/(2*pi);
Hxx = c*( 1./R2 - 2*DX.^2 ./ R4 );
Hyy = c*( 1./R2 - 2*DY.^2 ./ R4 );
Hxy = c*( -2*DX.*DY ./ R4 );   % = Hyx

% Stack rows as [Hxx; Hxy; Hyx; Hyy], blocks of 4x1 per target-source pair
out = [Hxx; Hxy; Hxy; Hyy];    % (4*nt) x ns

end
