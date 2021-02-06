function [theta_r,x_f,y_f]=double_slit(xi)

%% Set parameters
mem = 0.9;
Gam = mem*4.21;
Nx = 256; Ny = Nx; 
Lx = 40; Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/8;
plotoption = inf;
p = problem_setup_double_slit(Nx,Ny,Lx,Ly,Gam,dt_desired);

%% Run simulation
p.xi = xi; p.yi = -10; p.ui =0; p.vi = 0.2;
% theta = 0*pi/180;
% p.ui=speed_steady*cos(theta); p.vi = speed_steady*sin(theta);
p.nimpacts = 600;      % Number of impacts
%p.l1 = 14.7/4.75/2; % Half of opening length
p.l1 = 2;
%p.l1 = 1.55;
%p.l1 = 3.1/2;
%p.l2 = 6/4.75;   % Slit "breadth"
p.l2 = 0.6632;
%p.l2 = 6/4.75;
%p.l2 = 1;
%l2 = -inf;
p.l3 = 3.105;   % Center of slits
deep = (p.yy>p.l2|p.yy<0|(abs(p.xx-p.l3)<p.l1)|(abs(p.xx+p.l3)<p.l1));

p.d = p.d0_deep.*deep + p.d0_shallow.*(~deep);
p.a = p.a0_deep.*(p.xx<0|p.yy>0) + p.a0_shallow.*(p.xx>=0 & p.yy<=0);



if p.useGPU == 1
    p.d = gpuArray(p.d);  p.a = gpuArray(p.a);
end

[x_data,y_data,t_data, eta_data]  = trajectory(p,plotoption); 
l = length(x_data);
theta = atan(diff(x_data)./diff(y_data))*180/pi; 
k =length(theta);
theta_r = theta(k);
x_f =x_data(l);
y_f =y_data(l);

 
% save(['single_slit_wavefield_l1_',num2str(p.l1),'_l2_',num2str(p.l2),...
%     '_xin_',num2str(p.xi,'%.4f'),'_droptype_',num2str(p.drop_type),...
%       '_mem_',num2str(mem),'_N_',num2str(p.Nx),'_L_',num2str(p.Lx),'.mat'],...
%         'x_data','y_data','t_data','p','eta_data','theta');
%save(['double_slit_wavefield_l1_',num2str(p.l1),'_l2_',num2str(p.l2),...
 %   '_xin_',num2str(p.xi,'%.4f'),'_droptype_',num2str(p.drop_type),...
 %     '_mem_',num2str(mem),'_N_',num2str(p.Nx),'_L_',num2str(p.Lx),'.mat'],...
  %      'x_data','y_data','t_data','theta');





