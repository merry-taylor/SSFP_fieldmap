%
close all;
clear all;
addpath('lib/phaseUnwrapping');

% % invivo- knee ----------------------------------------------------------------
% 
% pc0=load('field_mapping_2_28_17/meas_MID389_SSFPdjp_TR6_TE3_PC0_FA15_FID3893_Kspace');
% pc90=load('field_mapping_2_28_17/meas_MID391_SSFPdjp_TR6_TE3_PC90_FA15_FID3895_Kspace');
% pc180=load('field_mapping_2_28_17/meas_MID393_SSFPdjp_TR6_TE3_PC180_FA15_FID3897_Kspace');
% pc270=load('field_mapping_2_28_17/meas_MID395_SSFPdjp_TR6_TE3_PC270_FA15_FID3899_Kspace');
% GRE5=load('field_mapping_2_28_17/meas_MID387_gre_TR100_TE5_FID3891_Kspace');
% GRE10=load('field_mapping_2_28_17/meas_MID386_gre_TR100_TE10_FID3890_Kspace');
% 

%phantom- fat/water ---------------------------------------------------------------

pc0=load('field_mapping_2_28_17/meas_MID367_SSFPdjp_TR6_TE3_PC0_FA15_FID3871_Kspace');
pc90=load('field_mapping_2_28_17/meas_MID369_SSFPdjp_TR6_TE3_PC90_FA15_FID3873_Kspace');
pc180=load('field_mapping_2_28_17/meas_MID371_SSFPdjp_TR6_TE3_PC180_FA15_FID3875_Kspace');
pc270=load('field_mapping_2_28_17/meas_MID373_SSFPdjp_TR6_TE3_PC270_FA15_FID3877_Kspace');
GRE5=load('field_mapping_2_28_17/meas_MID365_gre_TR100_TE5_FID3869_Kspace');
GRE10=load('field_mapping_2_28_17/meas_MID364_gre_TR100_TE10_FID3868_Kspace');


pc0c1=double(ifftshift(ifft2(ifftshift(pc0.kSpace(:,:,1))))); 
pc02= pc0c1; %(imag(pc0c1)+j*real(pc0c1));
pc90c1=double(ifftshift(ifft2(ifftshift(pc90.kSpace(:,:,1)))));
pc902=pc90c1; %(imag(pc90c1)+j*real(pc90c1));
pc180c1=double(ifftshift(ifft2(ifftshift(pc180.kSpace(:,:,1)))));
pc1802=pc180c1; %(imag(pc180c1)+j*real(pc180c1));
pc270c1=double(ifftshift(ifft2(ifftshift(pc270.kSpace(:,:,1)))));
pc2702=pc270c1; %(imag(pc270c1)+j*real(pc270c1));

GRE5c1=double(ifftshift(ifft2(ifftshift(GRE5.kSpace(:,:,1)))));
GRE10c1=double(ifftshift(ifft2(ifftshift(GRE10.kSpace(:,:,1)))));

%------------------------------------------------------------
%phase unwrap
%------------------------------------------------------------
scalingFactor = 20; %10; %20;
numImgs = 2;%6;
max_box_radius=10; %15;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

IM_mask=ones(size(GRE5c1));                     %Mask (if applicable)
promptMessage = sprintf('Is the phase correct? ');
titleBarCaption = 'Yes or No';
data{1} = GRE5c1;
data{2} = GRE10c1;
data{3}=pc0c1;
data{4}=pc90c1;
data{5}=pc180c1;
data{6}=pc270c1;

% correct phase offset
%theta = [0 1/4 1/2 3/4] *pi;
% theta = [0 0 0 0] *pi;
%      
%  data{3} = data{3} * exp(-1i*theta(1));
%  data{4} = data{4} * exp(-1i*theta(2));
%  data{5} = data{5} * exp(-1i*theta(3));
%  data{6} = data{6} * exp(-1i*theta(4));
 
 %GRE FieldMap
%GREFM=(angle(GRE5c1)-angle(GRE10c1))/(5*10^-3 * 2*pi) - 2.4; %no scalingFactor?
GREFM=(angle(data{1})-angle(data{2})); %no scalingFactor?


%M Field Map
% M = EllipticalModel2D(pc0c1, pc90c1, pc180c1, pc270c1);
% GSFM=angle(M)*2/(2*pi*.01);
u_M = EllipticalModel2D(data{3},data{4},data{5},data{6});
GSFM=angle(u_M)*2;
data{1} = GREFM;
data{2} = -GSFM;
   
if 1
    u_data={};
    for k=1:numImgs          
        IM = data{k};
        IM_mag=abs(IM);                             %Magnitude image
        IM_phase=IM; %angle(IM);                         %Phase image
        tryagain = true; %this really only makes sense if phaseUnwrapDebug is on . . .
        while tryagain
            close all
            residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
            branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
            [IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping-- this is where you can select start point of known phase
            %u_data{k} = IM_mag.*exp(1i*IM_unwrapped/scalingFactor); %rephase 10 is good for bottle
            u_data{k} = IM_unwrapped; %rephase 10 is good for bottle

            figure; imagesc(IM_phase),  axis square, axis off, title('Wrapped phase');
            figure; imagesc(angle(u_data{k})),  axis square, axis off, title('Unwrapped phase');
            figure; imagesc(IM_unwrapped), axis square, axis off, title('Golden Unwrapped phase'); %golden
            button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
            if strcmpi(button, 'Yes')
                tryagain = false;
            else
               tryagain = true;
            end
        end
    end
end
%------------
%
u_GSFM=u_data{2}*2/(6*10^-3* 2*pi);
u_GREFM=u_data{1}/(2*pi*5*10^-3);


%
% figure;
% subplot(1,2,1); imshow(GREFM,[]); title('GRE Field Map');
% subplot(1,2,2); imshow(GSFM,[]); title('Elliptical Signal Model Field Map');

figure;
u_CGREFM = u_GREFM; %imcrop(u_GREFM, [55 125 150 275]);
u_CGSFM = u_GSFM; %imcrop(u_GSFM, [55 125 150 275]);
% u_CGREFM = imcrop(u_GREFM, [5 170 250 230]);
% u_CGSFM = imcrop(u_GSFM, [5 170 250 230]);

maximum = max(max(max(abs(u_CGREFM))), max(max(abs(u_CGSFM))));

err = abs(u_CGREFM - (u_CGSFM+ (u_GREFM(270,100)-u_GSFM(270,100)))); %abs(u_CGREFM - (u_CGSFM-30.597)); (u_GREFM(273,102)-u_GSFM(273,102))
subplot(1,3,1); imagesc(u_CGREFM,[-500, 500]); title('GRE Field Map');
axis off;
subplot(1,3,2); imagesc((u_CGSFM+ (u_GREFM(270,100)-u_GSFM(270,100))),[-500, 500]); title('Elliptical Signal Model Field Map');
axis off;
subplot(1,3,3); imagesc(err, [0, 600]); title('Error');
axis off;


