basepath;

%% load model output

load('Ar_diags.mat')
load('PO4mmavg.mat','PO4');
load('PO4premmavg.mat','PO4premmavg');
load('PIC_OU_oxygen_decomp.mat','O2_glob','O2eq_glob');

% O:P Redfield ratio from MOBI 2.0
rOP = -10.6.*16;

%% construct tracers
% all units to mmol m-3
PO4a = 1e-3.*mean(PO4,4);
% preformed PO4 annual
PO4prea = 1e-3.*mean(PO4premmavg,4);
O2 = O2_glob{1};
O2eq = O2eq_glob{1};
Ar_Arsat = 1+mean(TRsatanommm,4);

% negative values for resp
O2_bio = O2 - Ar_Arsat.*O2eq;

% pos values corresponds to preformed P
O2_bio_star = O2_bio - PO4a.*rOP;
% PO4 preformed, converted to O2 units (pos values)
O2prea = -PO4prea.*rOP;


dpre = O2_bio_star - O2prea;

%%
cm = cmocean('curl');
% 30W
iAtl = 92;
crng = [-400:10:400];
clim = [-400 400];

B = PO4a.*rOP;
figure;
subplot(3,1,1);
[C,h] = contourf(y,z,squeeze(1000.*O2_bio(92,:,:))',crng);
h.LineStyle = 'none';
colormap(cm);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('O_2^{bio}');
colorbar;
caxis(clim);


subplot(3,1,2);
[C,h] = contourf(y,z,squeeze(-1000.*B(92,:,:))',crng);
h.LineStyle = 'none';
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('PO_4\timesr_{OP}');
colorbar;
caxis(clim);

subplot(3,1,3);
colormap(cm);
[C,h] = contourf(y,z,squeeze(1000.*O2_bio_star(92,:,:))',crng);
h.LineStyle = 'none';
colormap(cm);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('O2^{bio}^*');
colorbar;
caxis(clim);