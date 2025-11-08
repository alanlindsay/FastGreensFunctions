function ConvergenceAnalysis

%homeDir = getenv('HOME');
%run([homeDir '/chunkie/startup.m'])

close all;

cparams = [];
cparams.eps = 1.0e-8;
%cparams.nover = 0;
cparams.maxchunklen = 0.1; % setting a chunk length helps when the
cparams.ta = 0; cparams.tb = 2*pi;
% frequency is known
pref = [];  pref.k = 16;

export_fig = 0;

%% Ellipse and disk domains. Run for a = 1 and a>1 fo separate results.

o.a = 1.5; o.b = 1/o.a;

% Opimization points
src   = {};
src.r = [o.a/4; o.b/3];

surf_t0 = pi;
surf_scr = ellipse(surf_t0,o.a,o.b);
surf_scr = surf_scr(:)';

for k = 1:4
    pref.k = 2^(k+1);
    j_range = 1:20;

    for i = 1:length(j_range)

        nchk = j_range(i);
        chnkr = chunkerfuncuni(@(t) ellipse(t,o.a,o.b),nchk,cparams,pref);

        opts.Surf = 0; opts.Int = 1; 

        if (o.a>1)

            R_n = NeumR_2D(src.r,src.r,chnkr,opts);
            R_t = getRmEllipse(o,src.r,src.r);
            errR(i) = abs(R_t-R_n)/abs(R_t);

        else

            [R_n,J,H] = NeumR_2D(src.r,src.r,chnkr,opts);
            [R_t,J_t,H_t] = getRxDisk(src.r,src.r);

            J_n = [J.dx J.dy];
            H_n = [H.dxx H.dxy; H.dyx H.dyy];
            errR(i) = abs(R_t-R_n)/abs(R_t);
            errJ(i) = abs(norm(J_t-J_n))/norm(J_t);
            errH(i) = abs(norm(H_t-H_n))/norm(H_t);

        end

        opts.Surf = 1; opts.Int = 1; R_n = NeumR_2D(surf_scr,surf_scr,chnkr,opts);
        opts.Surf = 1; opts.Int = 0; R2_n = NeumR_2D(surf_scr,surf_scr,chnkr,opts);

        if (o.a>1)
            R_t = exactR_Int_surf_ellipse(surf_scr,o.a,o.b);
            R2_t = exactR_Ext_surf_ellipse(surf_scr,o.a,o.b);
            errIntSurfR(i) = abs(R_t-R_n)/abs(R_t);
            errExtSurfR(i) = abs(R2_t-R2_n)/abs(R2_t);
        else
            R_t = 1/(8*pi);
            errSurfR(i) = abs(R_t-R_n)/abs(R_t);
        end
    end

    str = ['$k = $ ' num2str(pref.k)];
    if (o.a == 1)       
        figure(1); hold on; plot(j_range,errR,'-o','linewidth',3,'displayname',str);
        title('Disk: $\mathcal{E}[R^{\textrm{int}}_b]$','interpreter','latex','fontsize',24); hold off;
        figure(2); hold on; plot(j_range,errJ,'-o','linewidth',3,'displayname',str);
        title('Disk: $\mathcal{E}[\nabla R ^{\textrm{int}}_b]$','interpreter','latex','fontsize',24); hold off;
        figure(3); hold on; plot(j_range,errH,'-o','linewidth',3,'displayname',str);
        title('Disk: $\mathcal{E}[\nabla^2 R^{\textrm{int}}_b]$','interpreter','latex','fontsize',24); hold off;
        figure(4); hold on; plot(j_range,errSurfR,'-o','linewidth',3,'displayname',str);
        title('Disk: $\mathcal{E}[R^{\textrm{int}}_s]$','interpreter','latex','fontsize',24); hold off;
    else
        figure(1); hold on; plot(j_range,errR,'-o','linewidth',3,'displayname',str);
        title('Ellipse: $\mathcal{E}[R^{\textrm{int}}_b]$','interpreter','latex','fontsize',24); hold off;
        figure(2); hold on; plot(j_range,errIntSurfR,'-o','linewidth',3,'displayname',str);
        title('Ellipse: $\mathcal{E}[R^{\textrm{int}}_s]$','interpreter','latex','fontsize',24); hold off;
        figure(3); hold on; plot(j_range,errExtSurfR,'-o','linewidth',3,'displayname',str);
        title('Ellipse: $\mathcal{E}[R^{\textrm{ext}}_s]$','interpreter','latex','fontsize',24); hold off;
    end
end

o.a = 1; o.b = 1/o.a;

% Opimization points
src   = {};
src.r = [o.a/4; o.b/3];

surf_t0 = pi;
surf_scr = ellipse(surf_t0,o.a,o.b);
surf_scr = surf_scr(:)';

if (o.a ==1)
    for j = 1:4
        figure(j);  yscale('log');
        set(gcf,'color','w');   set(gca,'fontsize',24);
        xlabel('Number of panels','FontSize',24,'Interpreter','latex')
        ylabel('Relative error','FontSize',24,'Interpreter','latex')
        lgd = legend;   lgd.Interpreter = 'latex';
        lgd.FontSize = 24;
        if (export_fig)
            if j==1, exportgraphics(figure(j),'DiskconvergenceR.png','resolution',300); end;
            if j==2, exportgraphics(figure(j),'DiskconvergenceJ.png','resolution',300); end;
            if j==3, exportgraphics(figure(j),'DiskconvergenceH.png','resolution',300); end;
        end
    end

else
    for j = 1:3
        figure(j)
        yscale('log');
        set(gcf,'color','w');
        set(gca,'fontsize',24);
        xlabel('Number of panels','FontSize',24,'Interpreter','latex')
        ylabel('Relative error','FontSize',24,'Interpreter','latex')
        lgd = legend;   lgd.Interpreter = 'latex';
        lgd.FontSize = 24; ylim([1e-15 1e0])
    end
    if (export_fig)
        exportgraphics(figure(1),['EllipseRIntBulk_a=' num2str(o.a) '.png'],'resolution',300);
        exportgraphics(figure(2),['EllipseRIntSurf_a=' num2str(o.a) '.png'],'resolution',300);
        exportgraphics(figure(3),['EllipseRExtSurf_a=' num2str(o.a) '.png'],'resolution',300);
    end

    %exportgraphics(figure(1),[ 'ellipseconvergence_a=' num2str(o.a)  '.png'],'resolution',300)

end

%R_true
%[R_e,~,~] = getRxDisk(src.r,src.r);

end

function [G,GradG,HessG] = getGxDisk(x,y)

c = 1/(2*pi);

x1 = x(1); x2 = x(2);
y1 = y(1); y2 = y(2);

[R,GradR,HessR] = getRxDisk(x,y);

G = -c*log(norm(x-y)) + R;

Gx = (-x1 + y1)/((x1 - y1)^2 + (x2 - y2)^2);
Gy = (-x2 + y2)/((x1 - y1)^2 + (x2 - y2)^2);

GradG = c*[Gx Gy] + GradR;

Gxx = ((x1 - y1)^2 - (x2 - y2)^2)/((x1 - y1)^2 + (x2 - y2)^2)^2;
Gxy = (2*(x1 - y1)*(x2 - y2))/((x1 - y1)^2 + (x2 - y2)^2)^2;
Gyy = (-(x1 - y1)^2 + (x2 - y2)^2)/((x1 - y1)^2 + (x2 - y2)^2)^2;


HessG = c*[Gxx Gxy; Gxy Gyy] + HessR;

end

function [R,GradR,HessR] = getRxDisk(x,y)

c = -1/(2*pi);

x1 = x(1); x2 = x(2);
y1 = y(1); y2 = y(2);

ax = x1^2 + x2^2; ay = y1^2 + y2^2; xy = 2*(x1*y1+x2*y2);

Z = 1 + ax*ay - xy;

R = 0.5*log(Z) - 0.5*(ax+ay) + 3/4;
R = c*R;

Rx = (ay*x1-y1)/Z -x1;
Ry = (ay*x2-y2)/Z -x2;

GradR = c*[Rx;Ry]';

Rxx = (ay*Z - 2*(ay*x1-y1)^2 )/Z^2 -1;
Rxy = -2*(ay*x1-y1)*(ay*x2-y2)/Z^2;
Ryy = (ay*Z - 2*(ay*x2-y2)^2 )/Z^2 -1;

HessR = c*[Rxx Rxy; Rxy Ryy];

end

function out = exactR_Int_surf_ellipse(y,a,b)

Ny = length(y)/2;
y1 = y(1:Ny); y2 = y(Ny+1:2*Ny);

beta = (a-b)/(a+b);

s2th0 = a^2*y2.^2/(b^2 * y1.^2 + a^2 * y2.^2);
th0 = asin(-sqrt(s2th0));

out = (y1.^2 + y2.^2)/(2*pi*a*b) - 3*(a^2+b^2)/(16*pi*a*b) + ...
    log( b^2 + (a^2-b^2)*s2th0 )/(2*pi);

for n = 1:90
    out = out - (2/pi)*(log(1-beta^(2*n)) + log(abs(1- beta^(2*n-1) *exp(2*1i*th0) ) ));
end

end

function out = exactR_Ext_surf_ellipse(x,a,b)

N = length(x)/2;

aE = sqrt(a^2-b^2);
a0 = x(1:N);
t0 = x(N+1:2*N);
ab = atanh(b/a);

%  x_0=a_e cosh(alpha_0) cos(\theta)
%  y_0=a_e sinh(alpha_0) sin(\theta)

t0 =  asin(x(N+1:2*N)/(aE*sinh(ab)));
a0 = ab;

out = real( log(2*aE) - ab + log( sinh(a0).^2 + sin(t0).^2 ) )/(2*pi);

end
