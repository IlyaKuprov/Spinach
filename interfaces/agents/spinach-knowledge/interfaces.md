# Spinach code index: interfaces

- Source root: `/home/kuprov/.openclaw/workspace/Spinach`
- Source commit: `f053e432a61d7144f3946d73d0a672e3ccfc3fc5`
- Source tree state: `clean`
- Path set: tracked MATLAB files from `git ls-files '*.m'`; untracked MATLAB files are excluded.
- Files indexed: **44** MATLAB files
- Generated: 2026-08-30T02:52:31

| File | Signature | Summary | LOC |
|---|---|---|---:|
| `interfaces/b2spinach.m` | `bdata=b2spinach(inpath)` | Imports time-domain NMR data recorded by Bruker instruments: reads the binary fid or ser file together with the acquisit | 464 |
| `interfaces/bootstrap.m` | `spin_system=bootstrap(volume)` | A minimal spin_system structure required to call many Spinach functions. Use this function if your system is not a spin  | 66 |
| `interfaces/castep/c2spinach.m` | `props=c2spinach(file_name)` | Reads the "new format" section of CASTEP .magres files. Syntax: props=c2spinach(file_name) Parameters: file_name -the na | 102 |
| `interfaces/comsol/comsol_conc.m` | `spin_system=comsol_conc(spin_system,file_name)` | Imports ASCII 2D concentration files produced by COMSOL. Syntax: spin_system=comsol_conc(spin_system,file_name) Paramete | 76 |
| `interfaces/comsol/comsol_import.m` | `mesh=comsol_import(comsol)` | COMSOL 2D mesh data import, cropping and preprocessing for Spinach. Syntax: mesh=comsol_import(comsol) Parameters: comso | 100 |
| `interfaces/comsol/comsol_mesh.m` | `mesh=comsol_mesh(file_name)` | Imports ASCII 2D mesh files produced by COMSOL. Syntax: mesh=comsol_mesh(file_name) Parameters: file_name -a character s | 146 |
| `interfaces/comsol/comsol_velo.m` | `mesh=comsol_velo(mesh,file_name)` | Imports ASCII 2D flow velocity files produced by COMSOL. Syntax: mesh=comsol_velo(mesh,file_name) Parameters: mesh -mesh | 74 |
| `interfaces/comsol/conc_plot.m` | `conc_plot(spin_system,conc,obs)` | 2D microfluidic concentration plotting function. Uses mesh tessellation information to plot concentrations as vertical b | 194 |
| `interfaces/comsol/mesh_crop.m` | `mesh=mesh_crop(mesh,ranges)` | 2D microfluidic mesh cropping. Updates the mesh object to remove anything outside the user-specified vertex coordi- nate | 109 |
| `interfaces/comsol/mesh_inact.m` | `mesh=mesh_inact(mesh,vertex_list)` | Marks 2D microfluidic mesh vertices as inactive in hydrodyna- mic and diffusive transport processes. Syntax: mesh=mesh_i | 63 |
| `interfaces/comsol/mesh_plot.m` | `mesh_plot(spin_system,qscale,nodelabels)` | 2D microfluidic mesh plotting function. Draws the mesh, its Vo- ronoi tessellation, and a quiver plot of velocities. Syn | 101 |
| `interfaces/comsol/mesh_preplot.m` | `mesh=mesh_preplot(mesh)` | Mesh preprocessing for drawing. Creates edge, triangle, and rectangle data structures needed for fast plotting later. Sy | 86 |
| `interfaces/comsol/mesh_vorn.m` | `mesh=mesh_vorn(mesh)` | Voronoi tessellation of a 2D COMSOL mesh. Syntax: mesh=mesh_vorn(mesh) Parameters: mesh -Spinach mesh object Outputs: me | 67 |
| `interfaces/g2spinach.m` | `[sys,inter]=g2spinach(props,particles,references,options)` | Makes Spinach data structures from parsed outputs of electronic structure theory packages, such as Gaussian and ORCA. Sy | 290 |
| `interfaces/gaussian/brokensymm.m` | `J=brokensymm(props_sing,props_trip)` | Exchange coupling estimation from a pair of DFT logs using Yamaguchi equation. The notation is: H=-2J*(Sa.Sb) Syntax: J= | 62 |
| `interfaces/gaussian/gparse.m` | `props=gparse(filename,options)` | A parser for Gaussian (03, 09, 16) calculation logs. Ex- tracts all potentially useful information. Syntax: props=gparse | 424 |
| `interfaces/gaussian/gslice.m` | `gslice()` | Slices a Gaussian geometry scan log into property calculation inputs at the energy minimum geometries. The function asks | 103 |
| `interfaces/gaussian/karplus_fit.m` | `[A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms)` | Fits a Karplus curve to a Gaussian dihedral angle scan. Syntax: [A,B,C,sA,sB,sC]=karplus_fit(dir_path,atoms) Parameters: | 124 |
| `interfaces/gissmo/gissmo2spinach.m` | `[sys,inter]=gissmo2spinach(filename,subsystem)` | Reads GISSMO files and forms Spinach data structures. Syntax: [sys,inter]=gissmo2spinach(file_name,subsystem) Parameters | 206 |
| `interfaces/hiper/spinach2hiper.m` | `spinach2hiper(file_name,amp,phi,off,dt)` | Exports phase-modulated optimal control waveforms into the format expected by Graham Smith's HiPER instrument. Syntax: s | 83 |
| `interfaces/jsonlab-1.5/jsonopt.m` | `val=jsonopt(key,default,varargin)` | No descriptive header found. | 35 |
| `interfaces/jsonlab-1.5/loadjson.m` | `data = loadjson(fname,varargin)` | No descriptive header found. | 503 |
| `interfaces/jsonlab-1.5/loadubjson.m` | `data = loadubjson(fname,varargin)` | No descriptive header found. | 457 |
| `interfaces/jsonlab-1.5/mergestruct.m` | `s=mergestruct(s1,s2)` | No descriptive header found. | 33 |
| `interfaces/jsonlab-1.5/savejson.m` | `json=savejson(rootname,obj,varargin)` | No descriptive header found. | 552 |
| `interfaces/jsonlab-1.5/saveubjson.m` | `json=saveubjson(rootname,obj,varargin)` | No descriptive header found. | 562 |
| `interfaces/jsonlab-1.5/struct2jdata.m` | `newdata=struct2jdata(data,varargin)` | No descriptive header found. | 96 |
| `interfaces/jsonlab-1.5/varargin2struct.m` | `opt=varargin2struct(varargin)` | No descriptive header found. | 40 |
| `interfaces/mestrenova/s2json.m` | `s2json(file_name,sys,inter,parameters,fid)` | Writes the parameters structure and the free induction decay into a JSON file that can be imported by MestreNova. Syntax | 70 |
| `interfaces/nmrpipe/fid2ascii.m` | `fid2ascii(file_name,fid)` | Writes out 1D, 2D and 3D free induction decays generated by Spinach as ASCII files in the following format [time1, time2 | 103 |
| `interfaces/nmrpipe/fid2pipe.m` | `fid2pipe(spin_system,file_root,fid,parameters,nmrpipe_root)` | Exports phase-sensitive 2D Spinach free induction decays into native NMRPipe time-domain files. Syntax: fid2pipe(spin_sy | 195 |
| `interfaces/orca/ocparse.m` | `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)` | ORCA cube file parser. Extracts the normalised probability density and the associated metric information from ORCA spin  | 95 |
| `interfaces/orca/oparse.m` | `props=oparse(file_name)` | A parser for ORCA text output logs, versions 2.6 to 6.1. Reads the geometry and every magnetic parameter that ORCA print | 479 |
| `interfaces/pdb_bmrb/nuclacid.m` | `[sys,inter]=nuclacid(pdb_file,shift_file,options)` | Nucleic acid data import function. Parses PDB and chemical shift data, runs a J-coupling guess using guess_j_nuc.m funct | 224 |
| `interfaces/pdb_bmrb/protein.m` | `[sys,inter,aux]=protein(pdb_file,bmrb_file,options)` | Protein data import function. Parses PDB and BMRB data, runs a J-coupl- ing guess, a CSA guess and outputs Spinach data  | 500 |
| `interfaces/pdb_bmrb/read_bmrb.m` | `[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)` | Reads BMRB file, extracts amino acid numbers, amino acid types, PDB atom identifiers and chemical shifts. Syntax: [aa_nu | 91 |
| `interfaces/pdb_bmrb/read_pdb_nuc.m` | `[res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)` | Reads the coordinates of all atoms from the user-specified PDB file and returns, for each atom, the residue number, the  | 95 |
| `interfaces/pdb_bmrb/read_pdb_pro.m` | `[aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,mod_id)` | Reads the a PDB file and returns amino acid numbers, amino acid types, PDB atom identifiers and Cartesian coordinates. S | 109 |
| `interfaces/retrieve_file.m` | `file_path=retrieve_file(file_url,file_name,dest_dir)` | Retrieves a file from an HTTPS link and stores it in a user- specified directory. Syntax: file_path=retrieve_file(file_u | 84 |
| `interfaces/spinjet/awg_interface.m` | `data=awg_interface(spin_system,awg_cmd,cmd_input)` | Interface to the Bruker SpinJet AWG, calling a library of Python scripts that inteface with the Bruker Xepr API. Syntax: | 141 |
| `interfaces/spinjet/py_run.m` | `arg_out=py_run(spin_system,pyscript,arg_in)` | Runs a python script from /interfaces/spinjet/Xepr_python/ folder of Bruker Xepr installation. Variable inputs and outpu | 126 |
| `interfaces/spinxml/parsexml.m` | `xml=parsexml(filename)` | Converts an XML file into a Matlab structure. Syntax: xml=parsexml(filename) This function is called by x2spinach() duri | 138 |
| `interfaces/spinxml/x2spinach.m` | `[sys,inter]=x2spinach(filename,shielding_refs)` | Reads SpinXML files and forms Spinach data structures. Syntax: [sys,inter]=x2spinach(file_name,shielding_refs) Parameter | 895 |
| `interfaces/v2spinach.m` | `vdata=v2spinach(inpath)` | Imports time-domain NMR data recorded by Varian and Agilent inst- ruments: reads the binary fid file and the procpar par | 329 |
