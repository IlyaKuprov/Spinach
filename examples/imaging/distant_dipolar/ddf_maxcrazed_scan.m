% Acquire the four water/acetone MAXCRAZED peaks for X, XY, and XZ.
% Inputs: npts is the [x y z] voxel count, max_step the largest LG4
% substep (s), delay the gap between the strong lobes (s), recovery
% the interval after the returning gradient before acquisition (s).
% Outputs: complex fid, spectrum, axes and peak frequencies (Hz), and
% the four peak magnitudes referenced to the strongest X peak.

function result=ddf_maxcrazed_scan(npts,max_step,delay,recovery)

% Check parameters before initialising the spin system
grumble(npts,max_step,delay,recovery);

% Water and acetone are separately diffusing proton resonances
gamma=spin('1H');
proton_freq=360e6;
sys.magnet=2*pi*proton_freq/gamma;
sys.isotopes={'1H','1H'};
sys.parallel={'processes',1};
sys.output='hush';
sys.disable={'hygiene'};
offset_hz=[1500 630];
carrier=mean(offset_hz);
inter.zeeman.scalar=num2cell(1e6*(offset_hz-carrier)/proton_freq);
inter.chem.parts={1,2};
inter.chem.concs=[1 1];

% Equal initial solvent volumes give the proton number densities
density=[997 785];
molar_mass=[18.01528 58.080]*1e-3;
number_density=0.5*6.02214076e23*density.*[2 6]./molar_mass;
hbar=1.054571817e-34;
kb=1.380649e-23;
temperature=298.15;
magn=number_density*(hbar*gamma)^2*sys.magnet/(4*kb*temperature);
diffusivity=[2 2]*1e-9;

% Build the two independent spin-half spaces
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
spin_system=basis(create(sys,inter),bas);
spin_system=assume(spin_system,'nmr');
dim=spin_system.bas.offsets(end);
drift=hamiltonian(spin_system);

% Use a bulk cell with one X pitch and two pitches along each strong axis
pitch=1e-3;
parameters.dims=pitch*[1 0.2 0.2];
parameters.groups={1,2};
parameters.frames=[1 2];

% Sample both dimensions on grids that put the four peaks on DFT bins
split=abs(diff(offset_hz));
nt1=8;
nt2=16;
dt1=1/(nt1*split);
dt2=2/(nt2*split);
result.fid=zeros(nt1,nt2,3);

% Build the chronological waveform for Fig. 1(c)
delta=1.5e-3;
pulse=10e-6;
g=2*pi/(gamma*delta*pitch);
strong=[0 0  0;...
        0 10 0;...
        0 0 10];
parameters.rf=zeros(8+nt2,2);
parameters.rf([1 4],:)=pi/(2*pulse);
parameters.acquire=[false(8,1);true(nt2,1)];

% Compare X, XY, and XZ with identical elapsed time before acquisition
for n=1:3

    % Collapse exactly uniform directions to their constant Fourier mode
    grid_size=npts;
    grid_size([false strong(n,2:3)==0])=1;
    parameters.npts=grid_size;
    nvox=prod(grid_size);
    spacing=parameters.dims./grid_size;
    symbol=ddf_periodic(grid_size,spacing,1/3);
    parameters.field=@(m)ddf_cartesian(m,grid_size,symbol);

    % Use separate diffusion blocks and a concentration-weighted receiver
    parameters.transport=sparse(dim*nvox,dim*nvox);
    source=sparse(dim,1);
    receiver=sparse(dim,1);
    parameters.scale=zeros(1,2);
    for q=1:2
        z=coil_state(spin_system,'Lz',q);
        source=source+z;
        parameters.scale(q)=magn(q)/real(z'*z);
        receiver=receiver+parameters.scale(q)*coil_state(spin_system,'L+',q);
        sites=(spin_system.bas.offsets(q)+1):spin_system.bas.offsets(q+1);
        proj=sparse(sites,sites,1,dim,dim);
        rate=ddf_diffusion(true(grid_size),spacing,diffusivity(q));
        parameters.transport=parameters.transport+1i*kron(rate,proj);
    end

    % Embed the spin operators only on the occupied spatial modes
    parameters.spc_dim=nvox;
    parameters.spn_dim=dim;
    parameters.rho0=kron(ones(nvox,1),source);
    parameters.coil=kron(ones(nvox,1)/nvox,receiver);
    H=kron(speye(nvox),drift);
    empty=sparse(nvox*dim,nvox*dim);
    G=cellfun(@inflate,g2fplanck(spin_system,parameters),'UniformOutput',false);

    % Add the balanced strong lobes after the 2X correlation lobe
    parameters.gradient=g*[zeros(2,3);1 0 0;0 0 0;2 0 0;strong(n,:);...
                           0 0 0;-strong(n,:);zeros(nt2,3)];
    for k=1:nt1

        % Increment t1 before the first correlation gradient
        parameters.duration=[pulse;(k-1)*dt1;delta;pulse;delta;delta;delay;delta;recovery;...
                             repmat(dt2,nt2-1,1)];
        parameters.steps=max(1,ceil(parameters.duration/max_step));
        result.fid(k,:,n)=transpose(ddf_sequence(spin_system,parameters,H,empty,empty,G,empty));

    end

end

% Use opposite indirect and direct Fourier signs for the selected DQ branch
result.f1=2*carrier+(-nt1/2:nt1/2-1)/(dt1*nt1);
result.f2=carrier+(-nt2/2:nt2/2-1)/(dt2*nt2);
result.spectrum=abs(fftshift(fftshift(fft(ifft(result.fid,[],1),[],2),1),2))/nt2;

% Read Table I frequencies and intensities in a common peak order
f1_peaks=[2*offset_hz(1);sum(offset_hz);sum(offset_hz);2*offset_hz(2)];
f2_peaks=offset_hz([1 1 2 2]).';
[~,row]=min(abs(result.f1-f1_peaks),[],2);
[~,col]=min(abs(result.f2-f2_peaks),[],2);
predicted=reshape(result.spectrum(sub2ind([nt1 nt2],row,col)+...
                                  (0:2)*nt1*nt2),4,3);
result.predicted=predicted/max(predicted(:,1));
result.f1_peaks=f1_peaks;
result.f2_peaks=f2_peaks;

end

% Consistency enforcement
function grumble(npts,max_step,delay,recovery)
if (~isnumeric(npts))||(~isequal(size(npts),[1 3]))||(~isreal(npts))||(any(~isfinite(npts(:))))||...
   (any(npts(:)~=fix(npts(:))))||(any(npts(:)<=0))
    error('npts must be a 1-by-3 row of positive integers.');
end
if (~isnumeric(max_step))||(~isscalar(max_step))||(~isreal(max_step))||(max_step<=0)||...
   (~isfinite(max_step))
    error('max_step must be a real positive finite numeric scalar.');
end
if (~isnumeric(delay))||(~isscalar(delay))||(~isreal(delay))||(delay<0)||(~isfinite(delay))
    error('delay must be a real nonnegative finite numeric scalar.');
end
if (~isnumeric(recovery))||(~isscalar(recovery))||(~isreal(recovery))||(recovery<0)||...
   (~isfinite(recovery))
    error('recovery must be a real nonnegative finite numeric scalar.');
end
end


