clear,clc
% This code is to reproduce the figure 3 of the JFM paper 'energy transfer
% in turbulent channel flows and implications for resolvent modelling'. The
% main part of the code is to calculate the three energy: production,
% dissipation and nonlinear transfer in Fourier space for turbulent channel 
% flow using equation 2.11.
% You may need to upgrade your MATLAB to the lastest version to use
% function 'pagemtimes'
%% main code
% specify the wavenumber you want to consider
kx = [0:2:6,-6:2:-2];
ky = [0:4:12];
kx_array = repelem(kx,length(ky));
ky_array = repmat(ky,1,length(kx));
kx_array = kx_array.';
ky_array = ky_array.';

% load and organise (switch y,z coordinates and scaling) the P4U data 
load('P4U_Retau85.mat');
Ly = Lz; clear Lz;
temp = Ny;
Ny = Nz;
Nz = temp;
dkx = 2;
dky = 4;
temp = y_grid_array;
y_grid_array = z_grid_array;
z_grid_array = temp;            % cheb(80)
U = permute(U, [1 3 2]) .* (Re/Retau);
V = permute(V, [1 3 2]) .* (Re/Retau);
W = permute(W, [1 3 2]) .* (Re/Retau);
temp = V;
V = W;
W = temp; 

% calculate the mean velocity
Uxy = squeeze(mean(mean(U,2),1));
Vxy = squeeze(mean(mean(V,2),1));
Wxy = squeeze(mean(mean(W,2),1));

% calculate the wall-normal derivative of the streamwise mean velocity
[Diff,~] = cheb(Nz-1);
Diff = Diff(2:(Nz-1),2:(Nz-1));
dUdz = Diff * Uxy(2:end-1);

% calculate the fluctuation velocity
Uxy = repmat(Uxy,1,24,24);
U = permute(U,[3 2 1]) - Uxy;
V = permute(V,[3 2 1]);
W = permute(W,[3 2 1]);

% discard the data at the z=1/-1 (because of no-slip boundary condition)
U = U(2:(Nz-1),:,:);
V = V(2:(Nz-1),:,:);
W = W(2:(Nz-1),:,:);

% calculate the three energy
[Pro, Diss, NonT] = get_three_energy(U,V,W,kx_array,ky_array,Nx,Ny,Nz,dkx,dky,Retau,dUdz,Diff);

% calculate the total positive and negative energy
Pro = Pro(1:4,1:4);
Pro(2:4,2:4) = Pro(2:4,2:4) .* 2;
Diss = Diss(1:4,1:4);
Diss(2:4,2:4) = Diss(2:4,2:4) .* 2;
NonT = NonT(1:4,1:4);
NonT(2:4,2:4) = NonT(2:4,2:4) .* 2;
sum_posi = sum(sum(Pro(Pro>0))) + sum(sum(Diss(Diss>0))) + sum(sum(NonT(NonT>0)));
sum_nega = sum(sum(Pro(Pro<0))) + sum(sum(Diss(Diss<0))) + sum(sum(NonT(NonT<0)));

%% plot figures
[x_coor,y_coor] = meshgrid([0:2:6],[0:4:12]);

% for the square edge
Pro_p = Pro;
Pro_p(Pro_p<0) = NaN;
Pro_n = Pro;
Pro_n(Pro_n>0) = NaN;
NonT_p = NonT;
NonT_p(NonT_p<0) = NaN;
NonT_n = NonT;
NonT_n(NonT_n>0) = NaN;

figure
adjust_size = 1500; % adjust this value for the square size (try and error)
for k_row = 1: size(x_coor,2)
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(Pro(:,k_row)),Pro(:,k_row)./(0.5*sum_posi-0.5*sum_nega),'filled','s'), hold on 
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(Pro_p(:,k_row)),[1.0000    0.4118    0.1608],'s')
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(Pro_n(:,k_row)),[0    0.4471    0.7412],'s')
    colormap(redblue)
    caxis([-1 1])
end
cbh = colorbar;
cbh.Ticks = [-1 0 1]; 
cbh.TickLabels = {'$-1$','$0$','$1$'};
cbh.TickLabelInterpreter = 'latex';
xline(1,'Color',[0.6510    0.6510    0.6510])
xline(3,'Color',[0.6510    0.6510    0.6510])
xline(5,'Color',[0.6510    0.6510    0.6510])
yline(2,'Color',[0.6510    0.6510    0.6510])
yline(6,'Color',[0.6510    0.6510    0.6510])
yline(10,'Color',[0.6510    0.6510    0.6510])
axis([-1 7 -2 14])
axis square
set(gca,'xtick',[0 2 4 6])
set(gca,'xticklabel',{'$0$','$2$','$4$','$6$'})
set(gca,'ytick',[0 4 8 12])
set(gca,'yticklabel',{'$0$','$4$','$8$','$12$'})
set(gca,'Fontsize',24)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
xlabel('$k_x$','Interpreter','Latex')
ylabel('$k_y$','Interpreter','Latex')
box on


figure
for k_row = 1: size(x_coor,2)
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(Diss(:,k_row)),Diss(:,k_row)./(0.5*sum_posi-0.5*sum_nega),'filled','s'), hold on 
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(Diss(:,k_row)),[0    0.4471    0.7412],'s')
    colormap(redblue)
    caxis([-1 1])
end
cbh = colorbar;
cbh.Ticks = [-1 0 1]; 
cbh.TickLabels = {'$-1$','$0$','$1$'};
cbh.TickLabelInterpreter = 'latex';
xline(1,'Color',[0.6510    0.6510    0.6510])
xline(3,'Color',[0.6510    0.6510    0.6510])
xline(5,'Color',[0.6510    0.6510    0.6510])
yline(2,'Color',[0.6510    0.6510    0.6510])
yline(6,'Color',[0.6510    0.6510    0.6510])
yline(10,'Color',[0.6510    0.6510    0.6510])
axis([-1 7 -2 14])
axis square
set(gca,'xtick',[0 2 4 6])
set(gca,'xticklabel',{'$0$','$2$','$4$','$6$'})
set(gca,'ytick',[0 4 8 12])
set(gca,'yticklabel',{'$0$','$4$','$8$','$12$'})
set(gca,'Fontsize',24)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
xlabel('$k_x$','Interpreter','Latex')
ylabel('$k_y$','Interpreter','Latex')
box on

figure
for k_row = 1: size(x_coor,2)
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(NonT(:,k_row)),NonT(:,k_row)./(0.5*sum_posi-0.5*sum_nega),'filled','s'), hold on 
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(NonT_p(:,k_row)),[1.0000    0.4118    0.1608],'s')
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(NonT_n(:,k_row)),[0    0.4471    0.7412],'s')
    colormap(redblue)
    caxis([-1 1])
end
cbh = colorbar;
cbh.Ticks = [-1 0 1]; 
cbh.TickLabels = {'$-1$','$0$','$1$'};
cbh.TickLabelInterpreter = 'latex';
xline(1,'Color',[0.6510    0.6510    0.6510])
xline(3,'Color',[0.6510    0.6510    0.6510])
xline(5,'Color',[0.6510    0.6510    0.6510])
yline(2,'Color',[0.6510    0.6510    0.6510])
yline(6,'Color',[0.6510    0.6510    0.6510])
yline(10,'Color',[0.6510    0.6510    0.6510])
axis([-1 7 -2 14])
axis square
set(gca,'xtick',[0 2 4 6])
set(gca,'xticklabel',{'$0$','$2$','$4$','$6$'})
set(gca,'ytick',[0 4 8 12])
set(gca,'yticklabel',{'$0$','$4$','$8$','$12$'})
set(gca,'Fontsize',24)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
xlabel('$k_x$','Interpreter','Latex')
ylabel('$k_y$','Interpreter','Latex')
box on


