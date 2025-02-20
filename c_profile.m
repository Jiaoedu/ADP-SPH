clc 
clear all
set(0,'defaultfigurecolor','w');
button=2; %1for unity, 2for concentration profile
select=1; %1for isotropic case, 2for anisotropic case
count=10; 
dt=3600;
maxtimestep=1000;
ntotal=160801; %160801
p=6;
if(select==1)
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

if (button==1) %//////////////////////////////////////////////////////////////////////////
  disp('unity')
count_d=1*count;
for i=1:1:maxtimestep/count
    % 
    % s_temp = load(['s' num2str(count*i) 'dx2.0.dat']);
    % s1(i,:,:) = s_temp;
    % s_temp = load(['s' num2str(count*i) 'dx3.0.dat']);
    % s2(i,:,:) = s_temp;
    % s_temp = load(['s' num2str(count*i) 'dx4.0.dat']);
    % s3(i,:,:) = s_temp;
    % s_temp = load(['s' num2str(count*i) 'dx5.0.dat']);
    % s4(i,:,:) = s_temp;
    % s_temp = load(['s' num2str(count*i) 'dx6.0.dat']);
    % s5(i,:,:) = s_temp;
    % s_temp = load(['s' num2str(count*i) 'dx7.0.dat']);
    % s6(i,:,:) = s_temp;
    % 
    % s_temp = load(['s' num2str(count_d*i) 'dx1.5.dat']);
    % s(i,:,:) = s_temp;
    % sh_temp = load(['sh' num2str(count*i) 'dx7.0.dat']);
    % sh(i,:,:) = sh_temp;
end
disp(i)
   % save('results_dx0.125_var1_h.mat','s1','s2','s3','s4','s5','s6')
save('results_dx1.5_dx0.5_var0.25_GI_q.mat','s','sh');
disp('done')

else %(button=2) %////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ntotal=160801;  %160801,103041,40401,10201
disp('plot')
dt=3600;
if(select==1)
  results1=load('results_dx0.125_var1_h.mat');
  results2=load('results_dx1.5_dx0.125_var1_GI_q.mat');
else
  results1=load('results_dx0.125_var1_h_D0.01.mat');
  results2=load('results_dx1.8_dx0.125_var1_GI_q_D0.01.mat');
end

step=100; %timestep
t=count*dt*step;
A1=2*t*Dyy+w^2; 
A2=2*t*Dxx+w^2;
A3=t*Dxy;
C2=Dxx*Dyy-Dxy^2;
C3=Dxx+Dyy;
C4=sqrt(4*t^2*C2+2*t*w^2*C3+w^4);

qoi_s(:,:)=results1.s1(step,:,:);
qoi_f(:,:)=results2.sh(step,:,:);

cmax_num_s=max(qoi_s(:,4)); 
cmax_num_f=max(qoi_f(:,4));
cmin_num_s=min(qoi_s(:,4));
cmin_num_f=min(qoi_f(:,4));
c=c0*w^2/C4*exp((-(qoi_s(:,2)-x0).^2*A1 -(qoi_s(:,3)-y0).^2*A2  +4.*(qoi_s(:,2)-x0).*(qoi_s(:,3)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));
e1_s=norm(qoi_s(:,4)-c(:),1)/ntotal;
e1_f=norm(qoi_f(:,4)-c(:),1)/ntotal; 
e2_s=norm(qoi_s(:,4)-c(:),2)/ntotal;
e2_f=norm(qoi_f(:,4)-c(:),2)/ntotal;
einf_s=norm(qoi_s(:,4)-c(:),inf);
einf_f=norm(qoi_f(:,4)-c(:),inf);
[M,I] = max(c);     
cpeak_num_s=qoi_s(I,4); 
cpeak_num_f=qoi_f(I,4);
cpeak_ana=c0*w^2./C4;
ecpeak_s=cpeak_num_s-M;
ecpeak_f=cpeak_num_f-M;
%--------------------------------------------------------------------------
e=0.02;
yy=30;
k=0;
k2=0;
k3=0;
k4=0;
k5=0;
for n=1:ntotal
  
    if(  abs(results1.s1(step,n,3)-yy) <= e  )  
        k=k+1;
        a(k,1)=results1.s1(step,n,2); %x
        a(k,2)=results1.s1(step,n,3); %y
        a(k,3)=c0*w^2/C4*exp((-( a(k,1)-x0)^2*A1-(a(k,2)-y0)^2*A2+4*( a(k,1)-x0)*( a(k,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%analytical solution
      
        a(k,4)=results1.s1(step,n,4); 
        a(k,5)=results1.s2(step,n,4); 
        a(k,6)=results1.s3(step,n,4);
        a(k,7)=results1.s4(step,n,4); 
        a(k,8)=results1.s5(step,n,4); 
        a(k,9)=results1.s6(step,n,4); 

    end 
    %------------------------------------------------------
    if(  abs(results1.s2(step,n,3)-yy) <= e  )  
         k2=k2+1;
         a2(k2,1)=results2.s(step,n,2); %x
         a2(k2,2)=results2.s(step,n,3); %y
         a2(k2,3)=c0*w^2/C4*exp((-( a2(k2,1)-x0)^2*A1-(a2(k2,2)-y0)^2*A2+4*( a2(k2,1)-x0)*( a2(k2,2)-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%analytical solution
        
         a2(k2,4)=results2.s(step,n,4);     
         a2(k2,5)=results2.sh(step,n,4); 
    end

end
b=sortrows(a,1); 
b2=sortrows(a2,1); 

%analytical solution-------------------------------------
xa=0:0.1:60;
ya=30;
cc=c0*w^2/C4.*exp((-( xa(:)-x0).^2*A1-(ya-y0)^2*A2+4*( xa(:)-x0).*( ya-y0)*A3)/(8*t^2*C2+4*w^2*t*C3+2*w^4));%½âÎö½â
p0=plot(xa,cc,'k-','linewidth',1.5);
hold on
%solution with different h------------------------------------------
p2=plot(b(:,1),b(:,4),'s','color',[50,127,205]/255,'MarkerSize',5,'linewidth',1.2) ; 
hold on
p3=plot(b(:,1),b(:,5),'^','color',[45,141,154]/255,'MarkerSize',3,'linewidth',1.2) ; 
hold on
p4=plot(b(:,1),b(:,6),'o','color',[62,199,157]/255,'MarkerSize',4,'linewidth',1.2) ; 
hold on
p5=plot(b(:,1),b(:,7),'s','color',[249,201,102]/255,'MarkerSize',5,'linewidth',1.2) ;
hold on
p6=plot(b(:,1),b(:,8),'^','color',[247,153,117]/255,'MarkerSize',3,'linewidth',1.2); 
hold on
p7=plot(b(:,1),b(:,9),'o','color',[247,116,146]/255,'MarkerSize',5,'linewidth',1.2);
hold on
%curve of s and sh solution------------------------------------------
p1=plot(b2(:,1),b2(:,4),'x','color',[51,95,114]/255,'MarkerSize',5,'linewidth',1.2) ;
hold on
p8=plot(b2(:,1),b2(:,5), '*b','MarkerSize',5,'linewidth',1.2);
hold on

xlim([15,45]);
ylim([0,0.45]);

ylabel('C/C_0','FontSize',12,'FontName','Times New Roman');
xlabel('x (m)','FontSize',12,'FontName','Times New Roman');
if(select==1)
    lgd1=legend([p0;p1;p2;p3;p4;p5;p6;p7;p8],'Analytical','h = 1.5\Deltax','h = 2\Deltax','h = 3\Deltax','h = 4\Deltax','h = 5\Deltax','h = 6\Deltax','h = 7\Deltax','h_{GI}=5.02\Deltax'); 
else
    lgd1=legend([p0;p1;p2;p3;p4;p5;p6;p7;p8],'Analytical','h = 1.8\Deltax','h = 2\Deltax','h = 3\Deltax','h = 4\Deltax','h = 5\Deltax','h = 6\Deltax','h = 7\Deltax','h_{GI}=5.02\Deltax'); 
end
set(lgd1,'FontSize',16,'FontName','Times New Roman','orientation','horizontal','NumColumns',1,'box','off');
set(gca, 'FontSize',16,'FontName','Times New Roman ','FontWeight','normal','linewidth',1.5);

disp('done')
end  %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////