% Exports magnetic resonance results with explicit physical axes. Syntax:
%
%      text=jcamp_signal(spin_system,axes,signal,units,names,info)
%
% Parameters:
%
%    spin_system - Spinach system structure, supplying its model field
%
%    axes,signal,units,names - see jcamp_grid; axes are in MATLAB array
%                  order, not acquisition-parameter order. Supports
%                  irregular times, frequency/field scans, mixed-domain
%                  arrays, and named quadrature or receiver components
%
%    info        - scalar structure with title, origin, owner, filename,
%                  detection ('CW' or 'PULSE'), method (JCAMP EMR core
%                  method), description (ASCII physical simulation
%                  parameters), and metadata (N-by-2 cell table).
%                  filename='' returns text without writing a file
%
% Outputs:
%
%    text - complete JCAMP character row, also written when requested
%
% All results are EMR SIMULATION data. Examples of core methods are FID,
% SPECTRUM, ENDOR, ESEEM, HYSCORE, ELDOR, KINETIC, and IMAGING. Describe
% sequence-specific names, e.g. DEER, in info.description or metadata.
% Units use JCAMP EMR keywords: SECOND, HERTZ, TESLA, etc. A custom scan
% must use physical coordinates, not indices. No interpolation, sorting,
% FFT, phase cycling, or normalisation occurs. Use jcamp_nmr for standard
% NMR outputs; this function covers explicitly sampled EMR data.
%
% talos@spindynamics.org

function text=jcamp_signal(spin_system,axes,signal,units,names,info)

% Check consistency
grumble(spin_system,info);

% Describe the simulated method and model field without guessing settings
metadata={'.DETECTION MODE',info.detection; '.METHOD',info.method;...
          '.SIMULATION SOURCE','Spinach'; '.SIMULATION PARAMETERS',info.description;...
          '$MODEL FIELD TESLA',spin_system.inter.magnet};
block=struct('title',info.title,'type','EMR SIMULATION',...
             'metadata',{[metadata; info.metadata]});

% Translate result arrays into the writer representation
blocks=jcamp_grid(axes,signal,units,names,block);
data=struct('title',info.title,'origin',info.origin,'owner',info.owner,'blocks',{blocks});
if ~isempty(info.filename), data.filename=info.filename; end
text=jcamp_export(data);

end

% Consistency enforcement
function grumble(spin_system,info)
if ~isstruct(spin_system)||~isscalar(spin_system)
    error('spin_system must be a scalar Spinach structure.');
end
if ~isstruct(info)||~isscalar(info)||...
   ~all(isfield(info,{'title','origin','owner','filename','detection','method','description','metadata'}))
    error('info must contain ownership, filename, detection, method, description, and metadata.');
end
if ~ischar(info.filename)||(~isempty(info.filename)&&~isrow(info.filename))||...
   ~iscell(info.metadata)||size(info.metadata,2)~=2
    error('info.filename must be a character row and metadata an N-by-2 cell table.');
end
end


