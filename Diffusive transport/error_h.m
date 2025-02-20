clc 
clear all
set(0,'defaultfigurecolor','w');
button=1; %1for isotropic, 2for anisotropic
count=10; 
dt=3600;
maxtimestep=1000;
sigma=0.5;  %particle irregularity  0.25,0.5,1;
p=4;
if(button==1)
%isotropic
Dxx=1e-6;
Dxy=0;
Dyy=1e-6;
w=2;
else
%anisotropic
Dxx=1e-6;
Dxy=0;
Dyy=1e-7;
w=1.2;
end
x0=30;
y0=30;
c0=1;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;

step = 100;
 t=count*dt*step;
 A1=2*t*Dyy+w^2; 
 A2=2*t*Dxx+w^2;
 A3=t*Dxy;
 C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);
%dx=0.5-----------------------------------------------------------------------------------------------------------------------------------------------
results1=load(['results_dx0.5_var' num2str(sigma) '_h.mat']);
results2=load(['results_dx1.5_dx0.5_var' num2str(sigma) '_GI_q.mat']);
ntotal=10201; %160801,40401,10201,103041
dx=0.5;
hh=[1.5*dx,1.5*dx,2*dx,2.5*dx,3*dx,3.5*dx,4*dx];
if(sigma==1)
   h_GI=1.2384;   %var1
elseif(sigma==0.5)
   h_GI=1.057;  %var0.5
elseif(sigma==0.25)
    h_GI=0.8331; %var0.25
else
    h_GI=1.5*dx;   %var0
end

qoi_s1=results2.s;
qoi_s2=results1.s1;
qoi_s3=results1.s2;
qoi_s4=results1.s3;
qoi_s5=results1.s4;
qoi_s6=results1.s5;
qoi_s7=results1.s6;
qoi_s8=results2.sh;

c=c0*w^2/C4*exp((-(qoi_s1(step,:,2)-x0).^2*A1 -(qoi_s1(step,:,3)-y0).^2*A2  +4.*(qoi_s1(step,:,2)-x0).*(qoi_s1(step,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
e2_s(1)=norm(qoi_s1(step,:,4)-c,2)/ntotal;
e2_s(2)=norm(qoi_s2(step,:,4)-c,2)/ntotal;
e2_s(3)=norm(qoi_s3(step,:,4)-c,2)/ntotal;
e2_s(4)=norm(qoi_s4(step,:,4)-c,2)/ntotal;
e2_s(5)=norm(qoi_s5(step,:,4)-c,2)/ntotal;
e2_s(6)=norm(qoi_s6(step,:,4)-c,2)/ntotal;
e2_s(7)=norm(qoi_s7(step,:,4)-c,2)/ntotal;
e2_s(8)=norm(qoi_s8(step,:,4)-c,2)/ntotal;
p1=semilogy(hh,e2_s(1:7),'--kx','linewidth',1.2,'markersize',7);
hold on
p5=semilogy(h_GI,e2_s(8),'xb','linewidth',1.2,'markersize',7);
hold on
%dx=0.25-------------------------------------------------------------------------------------------------------------------------
results1=load(['results_dx0.25_var' num2str(sigma) '_h.mat']);
results2=load(['results_dx1.5_dx0.25_var' num2str(sigma) '_GI_q.mat']);

ntotal=40401; %160801,103041,40401,10201
dx=0.25;
hh=[1.5*dx,2*dx,3*dx,4*dx,5*dx,6*dx,7*dx];

if(sigma==1)
   h_GI=0.8687;   %var1
elseif(sigma==0.5)
   h_GI=0.7706;  %var0.5
elseif(sigma==0.25)
    h_GI=0.5874; %var0.25
else
    h_GI=1.5*dx;   %var0
end

qoi_s1=results2.s;
qoi_s2=results1.s1;
qoi_s3=results1.s2;
qoi_s4=results1.s3;
qoi_s5=results1.s4;
qoi_s6=results1.s5;
qoi_s7=results1.s6;
qoi_s8=results2.sh;

c=c0*w^2/C4*exp((-(qoi_s1(step,:,2)-x0).^2*A1 -(qoi_s1(step,:,3)-y0).^2*A2  +4.*(qoi_s1(step,:,2)-x0).*(qoi_s1(step,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
e2_s(1)=norm(qoi_s1(step,:,4)-c,2)/ntotal;
e2_s(2)=norm(qoi_s2(step,:,4)-c,2)/ntotal;
e2_s(3)=norm(qoi_s3(step,:,4)-c,2)/ntotal;
e2_s(4)=norm(qoi_s4(step,:,4)-c,2)/ntotal;
e2_s(5)=norm(qoi_s5(step,:,4)-c,2)/ntotal;
e2_s(6)=norm(qoi_s6(step,:,4)-c,2)/ntotal;
e2_s(7)=norm(qoi_s7(step,:,4)-c,2)/ntotal;
e2_s(8)=norm(qoi_s8(step,:,4)-c,2)/ntotal;
p2=semilogy(hh,e2_s(1:7),'--ks','linewidth',1.2,'markersize',7);
hold on
p6=semilogy(h_GI,e2_s(8),'sb','linewidth',1.2,'markersize',7);
hold on
%dx=0.15625---------------------------------------------------------------------------------------------------------------------------
results1=load(['results_dx0.15625_var' num2str(sigma) '_h.mat']);
results2=load(['results_dx1.5_dx0.15625_var' num2str(sigma) '_GI_q.mat']);
ntotal=103041; %160801,40401,10201,2601,103041
dx=0.15625;
hh=[1.5*dx,2*dx,3*dx,4*dx,5*dx,6*dx,7*dx];

if(sigma==1)
   h_GI=0.7026;   %var1
elseif(sigma==0.5)
   h_GI=0.6251;  %var0.5
elseif(sigma==0.25)
    h_GI=0.475; %var0.25
else
    h_GI=1.5*dx;   %var0
end


qoi_s1=results2.s;
qoi_s2=results1.s1;
qoi_s3=results1.s2;
qoi_s4=results1.s3;
qoi_s5=results1.s4;
qoi_s6=results1.s5;
qoi_s7=results1.s6;
qoi_s8=results2.sh;

c=c0*w^2/C4*exp((-(qoi_s1(step,:,2)-x0).^2*A1 -(qoi_s1(step,:,3)-y0).^2*A2  +4.*(qoi_s1(step,:,2)-x0).*(qoi_s1(step,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
e2_s(1)=norm(qoi_s1(step,:,4)-c,2)/ntotal;
e2_s(2)=norm(qoi_s2(step,:,4)-c,2)/ntotal;
e2_s(3)=norm(qoi_s3(step,:,4)-c,2)/ntotal;
e2_s(4)=norm(qoi_s4(step,:,4)-c,2)/ntotal;
e2_s(5)=norm(qoi_s5(step,:,4)-c,2)/ntotal;
e2_s(6)=norm(qoi_s6(step,:,4)-c,2)/ntotal;
e2_s(7)=norm(qoi_s7(step,:,4)-c,2)/ntotal;
e2_s(8)=norm(qoi_s8(step,:,4)-c,2)/ntotal;
p3=semilogy(hh,e2_s(1:7),'--k^','linewidth',1.2,'markersize',5);
hold on
p7=semilogy(h_GI,e2_s(8),'^b','linewidth',1.2,'markersize',5);
hold on
%dx=0.125------------------------------------------------------------------------------------------------------------------------
results1=load(['results_dx0.125_var' num2str(sigma) '_h.mat']);
results2=load(['results_dx1.5_dx0.125_var' num2str(sigma) '_GI_q.mat']);
ntotal=160801; %160801,40401,10201,2601,103041
dx=0.125;
hh=[1.5*dx,2*dx,3*dx,4*dx,5*dx,6*dx,7*dx];

if(sigma==1)
   h_GI=0.6316;   %var1
elseif(sigma==0.5)
   h_GI=0.5593;  %var0.5
elseif(sigma==0.25)
    h_GI=0.4265; %var0.25
else
    h_GI=1.5*dx;   %var0
end

qoi_s1=results2.s;
qoi_s2=results1.s1;
qoi_s3=results1.s2;
qoi_s4=results1.s3;
qoi_s5=results1.s4;
qoi_s6=results1.s5;
qoi_s7=results1.s6;
qoi_s8=results2.sh;

c=c0*w^2/C4*exp((-(qoi_s1(step,:,2)-x0).^2*A1 -(qoi_s1(step,:,3)-y0).^2*A2  +4.*(qoi_s1(step,:,2)-x0).*(qoi_s1(step,:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
e2_s(1)=norm(qoi_s1(step,:,4)-c,2)/ntotal;
e2_s(2)=norm(qoi_s2(step,:,4)-c,2)/ntotal;
e2_s(3)=norm(qoi_s3(step,:,4)-c,2)/ntotal;
e2_s(4)=norm(qoi_s4(step,:,4)-c,2)/ntotal;
e2_s(5)=norm(qoi_s5(step,:,4)-c,2)/ntotal;
e2_s(6)=norm(qoi_s6(step,:,4)-c,2)/ntotal;
e2_s(7)=norm(qoi_s7(step,:,4)-c,2)/ntotal;
e2_s(8)=norm(qoi_s8(step,:,4)-c,2)/ntotal;
p4=semilogy(hh,e2_s(1:7),'--ko','linewidth',1.2,'markersize',5);
hold on
p8=semilogy(h_GI,e2_s(8),'ob','linewidth',1.2,'markersize',5);
hold on
xlim([0,2.5]);
ylim([4e-7,1e-4]);

if(sigma==1)
   lgd1=legend([p1;p2;p3;p4;p5;p6;p7;p8],'N=10201','N=40401','N=103041','N=160801','h_{GI}=1.22','h_{GI}=0.86','h_{GI}=0.69','h_{GI}=0.63');  %var1
elseif(sigma==0.5)
   lgd1=legend([p1;p2;p3;p4;p5;p6;p7;p8],'N=10201','N=40401','N=103041','N=160801','h_{GI}=1.04','h_{GI}=0.76','h_{GI}=0.62','h_{GI}=0.55');  %var0.5
elseif(sigma==0.25)
    lgd1=legend([p1;p2;p3;p4;p5;p6;p7;p8],'N=10201','N=40401','N=103041','N=160801','h_{GI}=0.82','h_{GI}=0.58','h_{GI}=0.47','h_{GI}=0.42');  %var0.25
else
    % lgd1=legend([p1;p2;p3;p4;p5;p6;p7;p8],'N=10201','N=40401','N=103041','N=160801','h_{GI}=0.75','h_{GI}=0.375','h_{GI}=0.234','h_{GI}=0.1875');  %var0.25
end
% set(lgd1,'FontSize',14,'FontName','Times New Roman','orientation','vertical','NumColumns',2,'box','on'); %horizontal

xlabel('Smoothing length h(m)','FontSize',14,'FontName','Times New Roman ');
ylabel('L_2','FontSize',14,'FontName','Times New Roman ');
set(gca, 'FontSize',14,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);
set(gcf,'Units','centimeters','Position',[1 1 14 10]);
print(gcf,'param','-dpng','-r600');