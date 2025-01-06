% MRI phantom library. Syntax:
%
%            [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)
%
% Parameters:
%
%    ph_name  - character string giving the name of the
%               phantom (see the function text)
%
% Outputs:
%
%    R1Ph   - a cube of R1 values
%
%    R2Ph   - a cube of R2 values
%
%    PDPh   - a cube of PD values
%
%    dims   - row vector of three cube dimensions, m
%
%    npts   - row vector of three cube dimensions, points
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=phantoms.m>

function [R1Ph,R2Ph,PDPh,dims,npts]=phantoms(ph_name)

% Check consistency
grumble(ph_name);

% Get the path
P=mfilename('fullpath'); P=P(1:(end-9));

switch ph_name
    
    case 'boob'
        
        % Load the breast phantom
        load([P filesep 'mri_lab_boob.mat'],'boob','nvoxels');
        
        % Preallocate arrays
        R1Ph=zeros(nvoxels); R2Ph=zeros(nvoxels); PDPh=zeros(nvoxels);
        
        % Set number of points
        npts=size(R1Ph); 
        
        % Assume 1 mm resolution
        dims=npts.*[1e-3 1e-3 1e-3];
        
        % Assign water proton density (1.0 for pure water)
        PDPh(boob==-1.0)=0.00; % Air 
        PDPh(boob==-2.0)=0.40; % Skin
        PDPh(boob==-4.0)=0.50; % Muscle
        PDPh(boob==+1.1)=0.75; % Fibroconnective, high water
        PDPh(boob==+1.2)=0.70; % Fibroconnective, medium water
        PDPh(boob==+1.3)=0.65; % Fibroconnective, low water
        PDPh(boob==+2.0)=0.60; % Transitional
        PDPh(boob==+3.1)=0.55; % Fatty, low fat
        PDPh(boob==+3.2)=0.50; % Fatty, medium fat
        PDPh(boob==+3.3)=0.45; % Fatty, high fat
        
        % Assign longitudinal water relaxation rate (3 Tesla)
        R1Ph(boob==-1.0)=0.00; % Air 
        R1Ph(boob==-2.0)=0.82; % Skin
        R1Ph(boob==-4.0)=1.11; % Muscle
        R1Ph(boob==+1.1)=0.69; % Fibroconnective, high water
        R1Ph(boob==+1.2)=1.00; % Fibroconnective, medium water
        R1Ph(boob==+1.3)=1.50; % Fibroconnective, low water
        R1Ph(boob==+2.0)=0.90; % Transitional
        R1Ph(boob==+3.1)=2.73; % Fatty, low fat
        R1Ph(boob==+3.2)=2.67; % Fatty, medium fat
        R1Ph(boob==+3.3)=2.61; % Fatty, high fat

        % Assign transverse water relaxation rate (3 Tesla)
        R2Ph(boob==-1.0)= 0.00; % Air 
        R2Ph(boob==-2.0)= 6.49; % Skin
        R2Ph(boob==-4.0)=34.48; % Muscle
        R2Ph(boob==+1.1)=18.39; % Fibroconnective, high water
        R2Ph(boob==+1.2)=25.00; % Fibroconnective, medium water
        R2Ph(boob==+1.3)=35.00; % Fibroconnective, low water
        R2Ph(boob==+2.0)= 8.00; % Transitional
        R2Ph(boob==+3.1)=14.70; % Fatty, low fat
        R2Ph(boob==+3.2)=16.54; % Fatty, medium fat
        R2Ph(boob==+3.3)=18.39; % Fatty, high fat
    
    case 'brain-highres'
        
        % Load MRiLab phantom
        load([P filesep 'mri_lab_brain.mat'],'VObj');
        
        % Assign the variables
        R1Ph=1./VObj.T1; R2Ph=1./VObj.T2; PDPh=VObj.Rho;

        % Clean up infinities
        R1Ph(isinf(R1Ph))=0; R2Ph(isinf(R2Ph))=0;

        % Point counts
        npts=size(PDPh);

        % Physical dimensions
        dims=npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes];
    
    case 'brain-medres'
        
        % Load MRiLab phantom
        load([P filesep 'mri_lab_brain.mat'],'VObj');
        
       % Assign the variables
        R1Ph=1./VObj.T1; R2Ph=1./VObj.T2; PDPh=VObj.Rho;

        % Clean up infinities
        R1Ph(isinf(R1Ph))=0; R2Ph(isinf(R2Ph))=0;

        % Downsample
        R1Ph=R1Ph(1:2:end,1:2:end,1:2:end);
        R2Ph=R2Ph(1:2:end,1:2:end,1:2:end);
        PDPh=PDPh(1:2:end,1:2:end,1:2:end);

        % Point counts
        npts=size(PDPh);

        % Physical dimensions
        dims=2*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes];
        
    case 'brain-lowres'
        
        % Load MRiLab phantom
        load([P filesep 'mri_lab_brain.mat'],'VObj');

        % Assign the variables
        R1Ph=1./VObj.T1; R2Ph=1./VObj.T2; PDPh=VObj.Rho;

        % Clean up infinities
        R1Ph(isinf(R1Ph))=0; R2Ph(isinf(R2Ph))=0;
        
        % Downsample
        R1Ph=R1Ph(1:4:end,1:4:end,1:4:end);
        R2Ph=R2Ph(1:4:end,1:4:end,1:4:end);
        PDPh=PDPh(1:4:end,1:4:end,1:4:end);

        % Point counts
        npts=size(PDPh);

        % Physical dimensions
        dims=4*npts.*[VObj.YDimRes VObj.XDimRes VObj.ZDimRes];
        
    otherwise
        
        % Complain and bomb out
        error('unrecognised phantom name.');
        
end
        
end

% Consistency enforcement
function grumble(ph_name)
if ~ischar(ph_name)
    error('ph_name must be a character string.');
end
end

% Thank you for saving me from the Library. I can see 
% a lot of girls without boyfriends there, and all they
% can do is study. It is a terrible fate. But I have 
% you - the Library cannot have me now.
%
% One of IK's girlfriends

