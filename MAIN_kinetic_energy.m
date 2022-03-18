clear,clc

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));
ifft_x = @(x) ifft(x, [], 3) * size(x, 3);
ifft_y = @(x) ifft(x, [], 2) * size(x, 2);
ifft_xy = @(x) ifft_x(ifft_y(x));

% load and process data (switch yz coordinates and do scaling)
load('P4U_Retau85.mat');
Ly = Lz; clear Lz;
temp = Ny;
Ny = Nz;
Nz = temp;
temp = y_grid_array;
y_grid_array = z_grid_array;
z_grid_array = temp;            % cheb(80)
U = permute(U, [1 3 2]) .* (Re/Retau);;
V = permute(V, [1 3 2]) .* (Re/Retau);;
W = permute(W, [1 3 2]) .* (Re/Retau);;
temp = V;
V = W;
W = temp; 

% calculate mean velocity
Uxy = squeeze(mean(mean(U,2),1));
Vxy = squeeze(mean(mean(V,2),1));
Wxy = squeeze(mean(mean(W,2),1));

%% plot mean velocity
figure
plot(Uxy,z_grid_array,'LineWidth',2,'Color','k')
xlabel('$U_{mean}$','interpreter','latex')
ylabel('$z$','interpreter','latex')
set(gca,'Fontsize',16)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';

figure
plot(Vxy,z_grid_array,'LineWidth',2,'Color','k')
xlabel('$V_{mean}$','interpreter','latex')
ylabel('$z$','interpreter','latex')
set(gca,'Fontsize',16)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';

figure
plot(Wxy,z_grid_array,'LineWidth',2,'Color','k')
xlabel('$W_{mean}$','interpreter','latex')
ylabel('$z$','interpreter','latex')
set(gca,'Fontsize',16)
ax_obj = gca;
ax_obj.TickLabelInterpreter = 'Latex';
%%

% calculate fluctuation velocity
Uxy = repmat(Uxy,1,24,24);
U = permute(U,[3 2 1]) - Uxy;
V = permute(V,[3 2 1]);
W = permute(W,[3 2 1]);

U_F2 = abs(fft_xy(U)).^2;
V_F2 = abs(fft_xy(V)).^2;
W_F2 = abs(fft_xy(W)).^2;

[~,WEIGHT] = clenCurt(length(z_grid_array)-1);

U_F2 = squeeze(pagemtimes(WEIGHT,U_F2));
V_F2 = squeeze(pagemtimes(WEIGHT,V_F2));
W_F2 = squeeze(pagemtimes(WEIGHT,W_F2));

KE_F2 = U_F2 + V_F2 + W_F2;

KE_F2 = KE_F2(1:4,1:4);
KE_F2(2:4,2:4) = KE_F2(2:4,2:4) .* 2;

%% plot kinetic energy distribution
[x_coor,y_coor] = meshgrid([0:2:6],[0:4:12]);

figure
adjust_size = 4000;
for k_row = 1: size(x_coor,2)
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(KE_F2(:,k_row)),KE_F2(:,k_row)./sum(sum(KE_F2)),'filled','s'), hold on 
    scatter(x_coor(:,k_row),y_coor(:,k_row),adjust_size.*abs(KE_F2(:,k_row)),[0.5020    0.5020    0.5020],'s')
    colormap(pink)
    oldcmap = colormap;
    colormap(flipud(oldcmap));
    caxis([0 1])
end
cbh = colorbar;
cbh.Ticks = [0 0.5 1]; 
cbh.TickLabels = {'$0$','$0.5$','$1$'};
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
