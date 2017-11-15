function [ u_map ] = unwrapping_time( map )
%This function performs phase unwrapping 
%   map- input field map in radians
%   u_map- phase unwrapped field map in radians
IM_mask=ones(size(map));
scalingFactor = 10; %20;
numImgs = 1;%6;
max_box_radius=15;                           %Maximum search box radius (pixels)
threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

if 1
    u_data=[];          
        IM = map;
        IM_mag=abs(IM);                             %Magnitude image
        IM_phase=IM; %angle(IM);                         %Phase image
        tryagain = true; %this really only makes sense if phaseUnwrapDebug is on . . .
        while tryagain
            close all
            residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
            branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
            [IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping
            %u_data{k} = IM_mag.*exp(1i*IM_unwrapped/scalingFactor); %rephase 10 is good for bottle
            u_data = IM_unwrapped; %rephase 10 is good for bottle

            figure; imagesc(IM_phase),  axis square, axis off, title('Wrapped phase');
            figure; imagesc(angle(u_data)),  axis square, axis off, title('Unwrapped phase');
            figure; imagesc(IM_unwrapped), axis square, axis off, title('Golden Unwrapped phase'); %golden
%            button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
%             if strcmpi(button, 'Yes')
%                 tryagain = false;
%             else
%                tryagain = true;
%             end\
            u_map=u_data;
            tryagain=false;
        end
end


end

