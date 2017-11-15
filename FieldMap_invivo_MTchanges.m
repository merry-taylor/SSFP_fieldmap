%
close all;
clear all;
addpath('lib/phaseUnwrapping');

% % invivo- head ----------------------------------------------------------------

% pc0=load('Fieldmaps/meas_MID30_SSFPdjp_TR10_TE5_PC0_FA30_FID17087_Kspace');
% pc90=load('Fieldmaps/meas_MID31_SSFPdjp_TR10_TE5_PC90_FA30_FID17088_Kspace');
% pc180=load('Fieldmaps/meas_MID32_SSFPdjp_TR10_TE5_PC180_FA30_FID17089_Kspace');
% pc270=load('Fieldmaps/meas_MID33_SSFPdjp_TR10_TE5_PC270_FA30_FID17090_Kspace');
% GRE5=load('Fieldmaps/meas_MID35_gre_TE5_FID17092_Kspace');
% GRE10=load('Fieldmaps/meas_MID34_gre_TE10_FID17091_Kspace');

%phantom ---------------------------------------------------------------

pc0=load('Fieldmaps/meas_MID18_SSFPdjp_TR10_TE5_PC0_FA30_FID17075_Kspace');
pc90=load('Fieldmaps/meas_MID19_SSFPdjp_TR10_TE5_PC90_FA30_FID17076_Kspace');
pc180=load('Fieldmaps/meas_MID20_SSFPdjp_TR10_TE5_PC180_FA30_FID17077_Kspace');
pc270=load('Fieldmaps/meas_MID21_SSFPdjp_TR10_TE5_PC270_FA30_FID17078_Kspace');
GRE5=load('Fieldmaps/meas_MID23_gre_TE5_FID17080_Kspace');
GRE10=load('Fieldmaps/meas_MID22_gre_TE10_FID17079_Kspace');


pc0c1=double(ifftshift(ifft2(ifftshift(pc0.kSpace(:,:,1))))); %really showing averages
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
scalingFactor = 10; %20;
numImgs = 2;%6;
max_box_radius=15;                           %Maximum search box radius (pixels)
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
            [IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping
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
u_GSFM=u_data{2}/(10*10^-3 * 2*pi); % elliptical signal model field map- hertz %TR=10ms head and uniform phantom. knee, fat/water phantom TR=6ms
u_GREFM=u_data{1}/(2*pi*5*10^-3); % GRE field map - hertz


%
% figure;
% subplot(1,2,1); imshow(GREFM,[]); title('GRE Field Map');
% subplot(1,2,2); imshow(GSFM,[]); title('Elliptical Signal Model Field Map');

figure;
u_CGREFM = imcrop(u_GREFM, [55 125 150 275]);
u_CGSFM = imcrop(u_GSFM, [55 125 150 275]);
% u_CGREFM = imcrop(u_GREFM, [5 170 250 230]);
% u_CGSFM = imcrop(u_GSFM, [5 170 250 230]);

maximum = max(max(max(abs(u_CGREFM))), max(max(abs(u_CGSFM))));

err = abs(u_CGREFM - (u_CGSFM+(u_GREFM(223,132)-u_GSFM(223,132)))); %abs(u_CGREFM - (u_CGSFM-30.597));
subplot(1,3,1); imagesc(u_CGREFM,[-500, 500]); colorbar; title('GRE Field Map');
axis off;
subplot(1,3,2); imagesc((u_CGSFM+(u_GREFM(223,132)-u_GSFM(223,132))),[-500, 500]); colorbar; title('Elliptical Signal Model Field Map');
axis off;
subplot(1,3,3); imagesc(err, [0, 600]); colorbar; title('Error');
axis off;


