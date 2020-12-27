% path to 1atm Ar simlulations
basepath;

%% load model output

load(fullfile(bpath,'Ar_diags.mat'));
load(fullfile(bpath,'Armmavg.mat'));
load(fullfile(bpath,'Areqpremmavg.mat'));
load('PO4mmavg.mat','PO4');
load('PO4premmavg.mat','PO4premmavg');
load('PIC_OU_oxygen_decomp.mat','O2_glob','O2eq_glob');
load('Salt_bc.mat');
load('Theta_bc.mat');
load('PIC_OU_oxygen_decomp.mat','O2_bio_glob','O2_phy_glob','O2dis_glob');

% O:P Redfield ratio from MOBI 2.0
rOP = -10.6.*16;

%% construct tracers
% all units to mmol m-3
PO4a = 1e-3.*mean(PO4,4);
% preformed PO4 annual
PO4prea = 1e-3.*mean(PO4premmavg,4);
O2 = O2_glob{1};
O2eq = O2eq_glob{1};

% Ar saturation using in situ vs preformed equilibrium concentrations
DAr = mean(Armmavg./TReqmm,4);

Areqpre = Areqpremmavg;
% surface sat is treated as passive tracer -- no mixing supesaturation
DArpre = mean(Armmavg./Areqpre,4);



O2_bio_cliff = O2_bio_glob{1};
O2_phy_cliff = O2_phy_glob{1};

%% in situ sat approach
Tm = mean(Tbc,4);
Sm = mean(Sbc,4);

% Oxygen equilibrium at ref 1 atm
O2eqST = mean(o2satfunc(Tbc,Sbc),4);

% Two ways to calculate O2_bio using Ar
% negative values for resp
O2_bio = O2 - DAr.*O2eqST;
%O2_bioST = O2 - DArST.*O2eqST;
%O2_bioST = O2 - DAr.*O2eqST;

O2_biopre = O2 - DArpre.*O2eqST;

% compare difference between O2_bio estimates
del_O2_bio = O2_bio_cliff-O2_bio;
%del_O2_bioST = O2_bio_cliff-O2_bioST;
% pos values corresponds to preformed P
O2_bio_star = O2_bio - PO4a.*rOP;
% PO4 preformed, converted to O2 units (pos values)
O2prea = -PO4prea.*rOP;


dpre = O2_bio_star - O2prea;

%%
cm = cmocean('curl');
% 30W
iAtl = 92;
crng = [-200:10:20];


B = PO4a.*rOP;
figure;
subplot(3,1,1);
[C,h] = contourf(y,z,squeeze(1000.*O2_bio_cliff(92,:,:))',crng);
h.LineStyle = 'none';
ax1 = gca;
ax1.Colormap = cmocean('haline');
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaO_2^{bio} (mmol m^{-3})');
colorbar;
%caxis(clim);

%% This should be the most consistent version w/ Cliff paper
% Includs 
% Units are mmol m-3
% (set to 1 in Cliff paper, but used monthly climatology in Ar sims?)
%
crng = 1000.*[-0.01:0.0001:.01];
subplot(3,1,2);
Del = 1000.*(O2_bio_cliff - O2_bio);
[C,h] = contourf(y,z,squeeze(Del(92,:,:))',crng);
h.LineStyle = 'none';
ax2 = gca;
ax2.Colormap = cm;
caxis([-1 1]);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaO_2^{bio}(cliff) - \DeltaO_2^{bio} (mmol m^{-3})');
colorbar;



%% This version is even further off b/c O2_bio uses the Arpre which doesn't 
% include the mixing-induced supersat
subplot(3,1,3);
Del = 1000.*(O2_bio_cliff - O2_biopre);
[C,h] = contourf(y,z,squeeze(Del(92,:,:))',crng);
h.LineStyle = 'none';
ax3 = gca;
colormap(ax3,cm);
caxis([-3 3]);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaO_2^{bio}(cliff) - \DeltaO_2^{bio}(pref) (mmol m^{-3})');
colorbar;

%% Second figure shows Ar sat from 1atm run
figure;
% Ar % saturation 
subplot(3,1,1);
Del = 100.*(DAr-1);
[C,h] = contourf(y,z,squeeze(Del(92,:,:))',crng);
h.LineStyle = 'none';
ax3 = gca;
colormap(ax3,cm);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaAr');
caxis([-4 4]);
colorbar;

% DAr_pre (%) This version removes the mixing supersat
subplot(3,1,2);
Del = 100.*(DArpre-1);
[C,h] = contourf(y,z,squeeze(Del(92,:,:))',crng);
h.LineStyle = 'none';
ax3 = gca;
colormap(ax3,cm);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaAr_{pref}');
caxis([-4 4]);
colorbar;
% O2_phy vs DAr (%): this is the difference between 'abiotic O2' and Ar
subplot(3,1,3);

Del = 100.*(DAr-1)-100.*(O2_phy_cliff./O2eqST-1);
[C,h] = contourf(y,z,squeeze(Del(92,:,:))',crng);
h.LineStyle = 'none';
ax3 = gca;
colormap(ax3,cm);
set(gca,'YDir','Reverse','TickDir','out','LineWidth',1);
title('\DeltaAr-\DeltaO2_{phy}');
caxis([-0.5 0.5]);
colorbar;
%%
function o2sat=o2satfunc(T,S)

C2K = 273.15;
f1 = log((298.15 - T)./(C2K + T));
f2 = f1.*f1;
f3 = f2.*f1;
f4 = f3.*f1;
f5 = f4.*f1;
o2sat = exp (2.00907 + 3.22014.*f1 + 4.05010*f2 ...
    + 4.94457*f3 - 2.56847e-1*f4 + 3.88767*f5 ...
    + S.*(-6.24523e-3 - 7.37614e-3*f1 - 1.03410e-2*f2 ...
    - 8.17083e-3*f3) - 4.88682e-7.*S.*S);
% Convert from ml/l to mol/m^3
o2sat = o2sat/22391.6*1000.0;
end