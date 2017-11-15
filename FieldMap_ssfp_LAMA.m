%2018 LAMA abstract
close all;
clear all;
addpath('lib/phaseUnwrapping');


%phantom ---------------------------------------------------------------
pc= load('ssfp_phantom_lama2.mat');
pc0=pc.im(:,:,1);
pc90=pc.im(:,:,2);
pc180=pc.im(:,:,3);
pc270=pc.im(:,:,4);



pc0c1=pc0; %double(ifftshift(ifft2(ifftshift(pc0(:,:,1))))); %really showing averages
pc02= pc0c1; %(imag(pc0c1)+j*real(pc0c1));
pc90c1=pc90;%double(ifftshift(ifft2(ifftshift(pc90(:,:,1)))));
pc902=pc90c1; %(imag(pc90c1)+j*real(pc90c1));
pc180c1=pc180;%double(ifftshift(ifft2(ifftshift(pc180(:,:,1)))));
pc1802=pc180c1; %(imag(pc180c1)+j*real(pc180c1));
pc270c1=pc270;%double(ifftshift(ifft2(ifftshift(pc270(:,:,1)))));
pc2702=pc270c1; %(imag(pc270c1)+j*real(pc270c1));


%------------------------------------------------------------
%phase unwrap
%------------------------------------------------------------
scalingFactor = 10; %20;
numImgs = 2;%6;
max_box_radius=15;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

IM_mask=ones(size(pc0c1));                     %Mask (if applicable)
promptMessage = sprintf('Is the phase correct? ');
titleBarCaption = 'Yes or No';
data{1} = 0;
data{2} = 0;
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
 



%M Field Map
% M = EllipticalModel2D(pc0c1, pc90c1, pc180c1, pc270c1);
% GSFM=angle(M)*2/(2*pi*.01);
u_M = EllipticalModel2D(data{3},data{4},data{5},data{6});
GSFM=angle(u_M)*2;
save(['GSfieldmap_phantom_LAMA2.mat'],'GSFM','-v7.3');
%data{1} = -GSFM;
%   
if 0
    u_data={};
    for k=1%1:numImgs          
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
% u_GSFM=u_data{1}/(10*10^-3 * 2*pi); % elliptical signal model field map- hertz %TR=10ms head and uniform phantom. knee, fat/water phantom TR=6ms
% 
% 
% 
% %
% % figure;
% % subplot(1,2,1); imshow(GREFM,[]); title('GRE Field Map');
% % subplot(1,2,2); imshow(GSFM,[]); title('Elliptical Signal Model Field Map');
% 
% figure;0
% u_CGSFM = imcrop(u_GSFM, [55 125 150 275]);
% % u_CGREFM = imcrop(u_GREFM, [5 170 250 230]);
% % u_CGSFM = imcrop(u_GSFM, [5 170 250 230]);
% 
% imagesc((u_CGSFM),[-500, 500]); colorbar; title('Elliptical Signal Model Field Map');
% axis off;



