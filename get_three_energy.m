% calculate the three energy: production (Pro), Diss (dissipation), NonT
% (nonlinear transfer) in Fourier space
function [Pro, Diss, NonT] = three_energy(u,v,w,kx_array,ky_array,nx,ny,nz,dkx,dky,Retau,dUdz,Diff)

[kx_m,ky_m] = getkm(kx_array,ky_array,nx,ny,dkx,dky);

fft_x = @(x) fft(x, [], 3) / size(x, 3);
fft_y = @(x) fft(x, [], 2) / size(x, 2);
fft_xy = @(x) fft_x(fft_y(x));
    
uu = u .* u;
vv = v .* v;
ww = w .* w;
uv = u .* v;
uw = u .* w;
vw = v .* w;

u_F  = fft_xy(u);
v_F  = fft_xy(v);
w_F  = fft_xy(w);
uu_F = fft_xy(uu);
vv_F = fft_xy(vv);
ww_F = fft_xy(ww);
uv_F = fft_xy(uv);
uw_F = fft_xy(uw);
vw_F = fft_xy(vw);

[duF_dx,duF_dy,duF_dz]    = get_3d(u_F,Diff,kx_m,ky_m);
[dvF_dx,dvF_dy,dvF_dz]    = get_3d(v_F,Diff,kx_m,ky_m);
[dwF_dx,dwF_dy,dwF_dz]    = get_3d(w_F,Diff,kx_m,ky_m);
[duuF_dx,~,~]             = get_3d(uu_F,Diff,kx_m,ky_m);
[~,dvvF_dy,~]             = get_3d(vv_F,Diff,kx_m,ky_m);
[~,~,dwwF_dz]             = get_3d(ww_F,Diff,kx_m,ky_m);
[duvF_dx,duvF_dy,~]       = get_3d(uv_F,Diff,kx_m,ky_m);
[duwF_dx,~,duwF_dz]       = get_3d(uw_F,Diff,kx_m,ky_m);
[~,dvwF_dy,dvwF_dz]       = get_3d(vw_F,Diff,kx_m,ky_m);

Pro  = -(dUdz.*u_F) .* conj(w_F) - conj(dUdz.*u_F) .* w_F;
Diss = 2.* (-duF_dx.*conj(duF_dx)-duF_dy.*conj(duF_dy)-duF_dz.*conj(duF_dz)-dvF_dx.*conj(dvF_dx)-dvF_dy.*conj(dvF_dy)-dvF_dz.*conj(dvF_dz)-dwF_dx.*conj(dwF_dx)-dwF_dy.*conj(dwF_dy)-dwF_dz.*conj(dwF_dz))./Retau;
NonT = -u_F.*conj(duuF_dx)-u_F.*conj(duvF_dy)-u_F.*conj(duwF_dz)-v_F.*conj(duvF_dx)-v_F.*conj(dvvF_dy)-v_F.*conj(dvwF_dz)-w_F.*conj(duwF_dx)-w_F.*conj(dvwF_dy)-w_F.*conj(dwwF_dz) - conj(u_F).*duuF_dx-conj(u_F).*duvF_dy-conj(u_F).*duwF_dz-conj(v_F).*duvF_dx-conj(v_F).*dvvF_dy-conj(v_F).*dvwF_dz-conj(w_F).*duwF_dx-conj(w_F).*dvwF_dy-conj(w_F).*dwwF_dz;

[~,W] = clenCurt(nz-1);
W = W(2:(nz-1));

Diss = squeeze(pagemtimes(W,Diss));
Pro  = squeeze(pagemtimes(W,Pro));
NonT = squeeze(pagemtimes(W,NonT));

end
