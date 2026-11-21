#include "../StdAfx.h"
#include "Gauss.h"

using namespace System;
using namespace System::IO;
using namespace System::Text::RegularExpressions;
using namespace System::Globalization;
using namespace System::Collections::Generic;

namespace {
    typedef cli::array<double> Values;
    typedef cli::array<String^> Lines;
    typedef SortedDictionary<int, Values^> AtomTable;
    typedef Dictionary<int, int> IsotopeTable;

    Lines^ words(String^ line)
    {
        return line->Split((cli::array<Char>^)nullptr, StringSplitOptions::RemoveEmptyEntries);
    }

    bool number(String^ text, double% value)
    {
        return Double::TryParse(text->Replace('D','E')->Replace('d','e'), NumberStyles::Float,
                                CultureInfo::InvariantCulture, value) && !Double::IsNaN(value) && !Double::IsInfinity(value);
    }

    Values^ diagonal(double value)
    {
        Values^ m=gcnew Values(9);
        m[0]=m[4]=m[8]=value;
        return m;
    }

    ::Matrix3x3^ matrix(Values^ values)
    {
        ::Matrix3x3^ m=gcnew ::Matrix3x3();
        for(int i=0;i<9;++i) m[i]=values[i];
        return m;
    }

    // Centre identifiers are external keys: ghost centres, dummies and lattice vectors are not spins.
    AtomTable^ orientation(Lines^ lines,int firstRow,int% n)
    {
        AtomTable^ atoms=gcnew AtomTable();
        for(int k=firstRow;k<lines->Length;++k) {
            Lines^ row=words(lines[k]); int id,z;
            if(row->Length==1 && row[0]->StartsWith("---")) return atoms;
            if((row->Length!=5 && row->Length!=6) || !Int32::TryParse(row[0],id) || !Int32::TryParse(row[1],z)) break;
            Values^ site=gcnew Values(4); site[0]=z;
            for(int a=0;a<3;++a) if(!number(row[row->Length-3+a],site[a+1])) return gcnew AtomTable();
            if(z>0) atoms[id]=site;
            n=k;
        }
        return gcnew AtomTable(); // Truncated or interleaved orientation: retain the last complete one.
    }

    Lines^ elementNames()
    {
        return (gcnew String("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"))->Split(' ');
    }

    // Some restart and reduced-output logs omit every orientation table. Read an
    // explicitly Cartesian archive or echoed input, never interpret internal coordinates as XYZ.
    AtomTable^ cartesian(Lines^ rows,bool archived,int start)
    {
        Lines^ elements=elementNames();
        AtomTable^ sites=gcnew AtomTable();
        for(int k=start;k<rows->Length;++k) {
            String^ line=rows[k]->Trim();
            if(!line->Length) break;
            // Gaussian input allows a symbol, an atomic number, or an ONIOM atom
            // specification with parenthesised PDB metadata and an optional link atom.
            Lines^ fields=words(Regex::Replace(line,"\\([^)]*\\)","")->Replace(',',' '));
            if(!fields->Length) break;
            int z;
            if(!Int32::TryParse(fields[0],z)) {
                Match^ label=Regex::Match(fields[0],"^([A-Za-z]+)(?=[-0-9]|$)");
                if(!label->Success) break;
                String^ symbol=label->Groups[1]->Value;
                z=System::Array::IndexOf(elements,symbol);
                if(z<0 && (symbol=="Bq" || symbol=="Tv" || symbol=="Gh")) z=0;
            }
            if(z<0 || z>=elements->Length) break;
            int end=fields->Length;
            if(!archived) for(int j=4;j<fields->Length && j<=5;++j) {
                if(fields[j]=="L" || fields[j]=="M" || fields[j]=="H") {end=j;break;}
            }
            if(end!=4 && !(end==5 && (fields[1]=="-1" || fields[1]=="0"))) return gcnew AtomTable();
            Values^ site=gcnew Values(4); site[0]=z;
            for(int j=1;j<=3;++j) if(!number(fields[end-4+j],site[j])) return gcnew AtomTable();
            if(z>0) sites[k-start+1]=site;
        }
        return sites;
    }

    AtomTable^ fallbackGeometry(Lines^ lines,int% multiplicity)
    {
        for(int n=lines->Length-1;n>=0;--n) {
            String^ line=lines[n]->Trim();
            if(line->StartsWith("1\\1\\") || line->StartsWith("1|1|")) {
                Char sep=line[1]; auto record=gcnew System::Text::StringBuilder(line);
                // Only the archive prefix through the atom list is needed (the Hessian may be huge).
                int boundaries=0; Char previous=' ';
                for(int k=n;k<lines->Length && boundaries<4;++k) {
                    String^ part=lines[k]->Trim();
                    if(k>n) {if(!part->Length) break; record->Append(part);}
                    for each(Char ch in part) {if(ch==sep && previous==sep) {++boundaries;previous=' ';} else previous=ch;}
                    if(part->EndsWith("@")) break;
                }
                Lines^ sections=record->ToString()->Split(gcnew Lines{gcnew String(sep,2)},StringSplitOptions::None);
                if(sections->Length<4) continue;
                Lines^ rows=sections[3]->Split(sep); if(rows->Length<2) continue;
                AtomTable^ found=cartesian(rows,true,1);
                if(found->Count) {
                    Lines^ state=rows[0]->Split(','); int mult;
                    if(!multiplicity && state->Length>1 && Int32::TryParse(state[1],mult)) multiplicity=mult;
                    return found;
                }
            } else if(line->StartsWith("Center") && line->EndsWith("Coordinates (Angstroms)") &&
                      n+2<lines->Length && Regex::IsMatch(lines[n+1],"^\\s*Number\\s+Number\\s+X\\s+Y\\s+Z\\s*$") &&
                      lines[n+2]->Trim()->StartsWith("---")) {
                // Standalone thermochemistry output prints this table without an orientation title.
                int end=n; AtomTable^ found=orientation(lines,n+3,end);
                if(found->Count) return found;
            } else if(line=="Input Z-matrix:" || line=="Redundant internal coordinates found in file." || line=="Z-Matrix found in chk file." || line=="Symbolic Z-matrix:") {
                int k=n+1; while(k<lines->Length && lines[k]->TrimStart()->StartsWith("Charge =")) ++k;
                AtomTable^ found=cartesian(lines,false,k);
                if(found->Count) return found;
            }
        }
        return gcnew AtomTable();
    }

    // Gaussian sometimes interchanges the axes of its reported rotation matrix. Fit the actual coordinates.
    Eigen::Matrix3d inertialRotation(AtomTable^ lab,AtomTable^ principal,Values^ axes)
    {
        if(!lab->Count || lab->Count!=principal->Count) throw gcnew FormatException("Gaussian principal-axis and molecular geometries do not match.");
        Eigen::MatrixXd a(lab->Count,3),b(lab->Count,3);
        int i=0;
        for each(int id in lab->Keys) {
            if(!principal->ContainsKey(id) || principal[id][0]!=lab[id][0]) throw gcnew FormatException("Unmatched Gaussian principal-axis centre.");
            for(int k=0;k<3;++k) {a(i,k)=principal[id][1+k]; b(i,k)=lab[id][1+k];} ++i;
        }
        a.rowwise()-=a.colwise().mean().eval(); b.rowwise()-=b.colwise().mean().eval();
        Eigen::Matrix3d cov=a.transpose()*b;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov,Eigen::ComputeFullU|Eigen::ComputeFullV);
        Eigen::Matrix3d sign=Eigen::Matrix3d::Identity(); sign(2,2)=(svd.matrixV()*svd.matrixU().transpose()).determinant()<0?-1:1;
        Eigen::Matrix3d r=svd.matrixV()*sign*svd.matrixU().transpose();
        // A linear or monatomic geometry cannot determine the transverse axes. Use the
        // printed principal basis when its direction reproduces the actual coordinates.
        if(svd.rank()<2 && axes!=nullptr) {
            Eigen::Matrix3d basis;
            for(int j=0;j<9;++j) basis(j/3,j%3)=axes[j];
            if((basis.transpose()*basis-Eigen::Matrix3d::Identity()).norm()<1e-4) {
                if((a*basis.transpose()-b).norm()<=1e-4*(1+b.norm())) r=basis;
                else if((a*basis-b).norm()<=1e-4*(1+b.norm())) r=basis.transpose();
                // Reversing the molecular axis can invert all coordinates without
                // changing a rank-two tensor. Keep the printed transverse basis.
                else if((a*basis.transpose()+b).norm()<=1e-4*(1+b.norm())) r=-basis;
                else if((a*basis+b).norm()<=1e-4*(1+b.norm())) r=-basis.transpose();
            }
        }
        if((a*r.transpose()-b).norm()>1e-3*(1+b.norm())) throw gcnew FormatException("Inconsistent Gaussian principal-axis coordinates.");
        return r;
    }

    // Keep the printed XYZ components by label, even when Gaussian prints columns instead of rows.
    Values^ components(Lines^ lines,int start,wchar_t first)
    {
        if(start+2>=lines->Length) return nullptr;
        Values^ m=gcnew Values(9); int seen=0;
        for(int k=start;k<start+3;++k) {
            MatchCollection^ fields=Regex::Matches(lines[k],"([XYZabc])([XYZabc])\\s*=\\s*([+\\-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[DdEe][+\\-]?[0-9]+)?|\\*+)");
            for each(Match^ f in fields) {
                int a=f->Groups[1]->Value[0]-first,b=f->Groups[2]->Value[0]-first;
                if(a<0 || a>2 || b<0 || b>2) return nullptr;
                int j=3*a+b;
                if((seen & (1<<j)) || !number(f->Groups[3]->Value,m[j])) return nullptr;
                seen |= 1<<j;
            }
        }
        return seen==511?m:nullptr;
    }

    void addTensor(Model^ model,Atom^ a,Atom^ b,InteractionKind kind,Unit^ unit,Values^ values)
    {
        Tensor^ t=gcnew Tensor(a,b,kind,matrix(values)); t->unit=unit;
        t->Frame=model->RefFrameCollection[1]; model->InteractionCollection->Add(t);
    }

    // Gaussian mass numbers and magnetic isotope labels need not agree (e.g. 12C geometry, 13C HFC).
    Isotope^ isotope(ArrayList^ library,int z,int mass,int magnetic)
    {
        Isotope^ first=nullptr; Isotope^ massMatch=nullptr;
        for each(Isotope^ iso in library) if(iso->AtomNumber==z) {
            if(first==nullptr) first=iso;
            if(iso->Mass==magnetic) return iso;
            if(iso->Mass==mass) massMatch=iso;
        }
        if(first==nullptr) throw gcnew FormatException("The GUI isotope table has no element with atomic number "+z+".");
        return massMatch!=nullptr?massMatch:first;
    }

    int massOf(IsotopeTable^ masses,int id)
    {
        int m; return masses->TryGetValue(id,m)?m:0;
    }

    bool sameSites(AtomTable^ a,AtomTable^ b)
    {
        if(a->Count!=b->Count) return false;
        for each(int id in a->Keys) if(!b->ContainsKey(id) || a[id][0]!=b[id][0]) return false;
        return true;
    }
}

Gauss::Gauss(String^ filename,Model^% model)
{
    SystemModel=model; FileName=filename; ReadingFile=gcnew StreamReader(filename); FormatName="Gauss";
}

bool Gauss::LoadFile()
{
    List<String^>^ buffer=gcnew List<String^>();
    try { while(!ReadingFile->EndOfStream) buffer->Add(ReadingFile->ReadLine()); }
    finally { ReadingFile->Close(); }
    Lines^ lines=buffer->ToArray();
    AtomTable^ geometry=gcnew AtomTable(); AtomTable^ principal=nullptr; AtomTable^ input=nullptr;
    AtomTable^ compositeDipolar=gcnew AtomTable(); Dictionary<int,double>^ compositeContact=gcnew Dictionary<int,double>();
    bool fullSystem=true, oniom=false; Lines^ elements=elementNames();
    IsotopeTable^ masses=gcnew IsotopeTable(); IsotopeTable^ hmass=gcnew IsotopeTable(); IsotopeTable^ qmass=gcnew IsotopeTable();
    AtomTable^ shielding=gcnew AtomTable(); AtomTable^ quadrupole=gcnew AtomTable();
    AtomTable^ dipolar=gcnew AtomTable(); Dictionary<int,double>^ contact=gcnew Dictionary<int,double>();
    typedef Tuple<int,int> Pair;
    Dictionary<Pair^,double>^ coupling=gcnew Dictionary<Pair^,double>();
    Values^ g=nullptr; Values^ chi=nullptr; Values^ axes=nullptr; int multiplicity=0; bool readCharge=false;
    for(int n=0;n<lines->Length;++n) {
        String^ line=lines[n]->Trim();
        if(line->Length==0) continue;
        // Restarts can omit the charge card and reuse the previous spin state.
        if(line->StartsWith("Entering Link 1 =")) {readCharge=false; fullSystem=true; oniom=false;}
        else if(line->StartsWith("ONIOM: generating point")) {oniom=true; fullSystem=line->EndsWith("on real system.");}
        else if(line=="ONIOM: calculating energy.") fullSystem=true;
        // Link atoms replace real sites in ONIOM subsystem calculations. Their magnetic
        // tensors must never overwrite a property's physical atom in the real system.
        else if(!fullSystem && line!="Isotropic Fermi Contact Couplings" && line!="ONIOM Magnetic shielding tensor (ppm):") continue;
        else if(!readCharge && line->StartsWith("Charge =")) {
            Match^ m=Regex::Match(line,"Multiplicity\\s*=\\s*(\\d+)");
            if(m->Success) {
                int next=Int32::Parse(m->Groups[1]->Value);
                if(multiplicity && next!=multiplicity) {contact->Clear(); dipolar->Clear(); hmass->Clear(); compositeDipolar->Clear(); compositeContact->Clear(); g=nullptr;}
                multiplicity=next; readCharge=true;
            }
        }
        else if(!multiplicity && line->Contains("alpha electrons") && line->Contains("beta electrons")) {
            Match^ m=Regex::Match(line,"^(\\d+)\\s+alpha electrons\\s+(\\d+)\\s+beta electrons");
            if(m->Success) multiplicity=1+Math::Abs(Int32::Parse(m->Groups[1]->Value)-Int32::Parse(m->Groups[2]->Value));
        }
        else if(line=="Input orientation:" || line=="Standard orientation:" || line=="Z-Matrix orientation:") {
            AtomTable^ next=orientation(lines,n+5,n);
            if(next->Count==0) continue;
            if(geometry->Count && !sameSites(geometry,next)) {
                shielding->Clear(); quadrupole->Clear(); qmass->Clear(); dipolar->Clear(); contact->Clear(); hmass->Clear(); masses->Clear();
                coupling->Clear(); compositeDipolar->Clear(); compositeContact->Clear(); input=nullptr; g=nullptr; chi=nullptr;
            }
            if(line=="Input orientation:") input=next;
            geometry=next; principal=nullptr; axes=nullptr;
        }
        else if(line=="Rotation matrix to Principal Axis frame:" && n+4<lines->Length) {
            Values^ basis=gcnew Values(9); bool valid=true;
            for(int i=0;i<3;++i) {
                Lines^ row=words(lines[n+2+i]);
                if(row->Length!=4) {valid=false; break;}
                for(int j=0;j<3;++j) if(!number(row[j+1],basis[i*3+j])) valid=false;
            }
            axes=valid?basis:nullptr;
        }
        else if(line=="Principal axis orientation:") principal=orientation(lines,n+5,n);
        else if(line->StartsWith("Atom ") && n+1<lines->Length && lines[n+1]->TrimStart()->StartsWith("IAtWgt=")) {
            Lines^ ids=words(line); Lines^ weights=words(lines[++n]->Split('=')[1]);
            for(int j=1;j<ids->Length && j<=weights->Length;++j) {int id,m; if(Int32::TryParse(ids[j],id) && Int32::TryParse(weights[j-1],m)) masses[id]=m;}
        }
        else if(line=="Magnetic shielding (ppm):" || line=="GIAO Magnetic shielding tensor (ppm):" || line=="ONIOM Magnetic shielding tensor (ppm):" || line=="SCF GIAO Magnetic shielding tensor (ppm):" || line=="MP2 GIAO Magnetic shielding tensor (ppm):") {
            shielding->Clear();
            for(int k=n+1;k<lines->Length && lines[k]->Contains("Isotropic =");) {
                int id; Lines^ fields=words(lines[k]);
                if(!Int32::TryParse(fields[0],id)) break;
                Values^ m=components(lines,k+1,'X');
                // Asterisks denote unprintable values, not zero: omit the affected tensor only.
                if(m!=nullptr) shielding[id]=m;
                k+=5;
                if(k<lines->Length && lines[k]->Trim()=="Eigenvectors:") k+=4;
                n=k-1;
            }
        }
        else if(line=="g tensor [g = g_e + g_RMC + g_DC + g_OZ/SOC]:" || line=="g tensor (ppm):") g=components(lines,n+1,'X');
        else if(line=="Magnetic susceptibility (cgs-ppm):" || line=="Magnetic susceptibility tensor (cgs-ppm):") {
            int k=n+1; if(k<lines->Length && lines[k]->Contains("Isotropic")) ++k;
            chi=components(lines,k,'X');
        }
        else if(line=="Total nuclear spin-spin coupling J (Hz):") {
            coupling->Clear(); List<int>^ columns=gcnew List<int>();
            for(int k=n+1;k<lines->Length;++k) {
                Lines^ row=words(lines[k]); int id;
                if(row->Length==0 || !Int32::TryParse(row[0],id)) break;
                if(row->Length==1 || !row[1]->Contains(".") && row[1]->IndexOfAny(gcnew cli::array<Char>{'D','E','*'})<0) {
                    columns->Clear(); for each(String^ field in row) {int col; if(!Int32::TryParse(field,col)) break; columns->Add(col);}
                } else {
                    for(int j=1;j<row->Length && j<=columns->Count;++j) {double v; if(columns[j-1]!=id && number(row[j],v)) coupling[Tuple::Create(id,columns[j-1])]=v;}
                }
                n=k;
            }
        }
        else if(line=="Isotropic Fermi Contact Couplings") {
            if(fullSystem) {contact->Clear(); hmass->Clear(); compositeDipolar->Clear(); compositeContact->Clear();}
            for(int k=n+2;k<lines->Length;++k) {
                Match^ m=Regex::Match(lines[k],"^\\s*(\\d+)\\s+([A-Za-z]+)\\((\\d+)\\)?\\s+\\S+\\s+\\S+\\s+(\\S+)");
                if(!m->Success) break;
                int id=Int32::Parse(m->Groups[1]->Value); double v;
                if(fullSystem || (geometry->ContainsKey(id) && geometry[id][0]<elements->Length && elements[(int)geometry[id][0]]==m->Groups[2]->Value)) {
                    hmass[id]=Int32::Parse(m->Groups[3]->Value);
                    if(fullSystem && number(m->Groups[4]->Value,v)) contact[id]=v;
                }
                n=k;
            }
        }
        else if(line=="Anisotropic Spin Dipole Couplings in Principal Axis System") {
            dipolar->Clear(); int k=n+1;
            while(k<lines->Length && (String::IsNullOrWhiteSpace(lines[k]) || lines[k]->TrimStart()->StartsWith("---") || lines[k]->TrimStart()->StartsWith("Atom"))) ++k;
            while(k+2<lines->Length && lines[k]->TrimStart()->StartsWith("Baa")) {
                Lines^ mid=words(lines[k+1]); int id;
                if(mid->Length<3 || !Int32::TryParse(mid[0],id)) break;
                Values^ m=gcnew Values(9); bool complete=true;
                for(int j=0;j<3;++j) {
                    Lines^ row=words(lines[k+j]); int start=System::Array::IndexOf(row,j==0?"Baa":(j==1?"Bbb":"Bcc")); double v; Values^ axis=gcnew Values(3);
                    if(start<0 || start+8!=row->Length) throw gcnew FormatException("Invalid Gaussian Anisotropic Couplings row.");
                    if(!number(row[start+3],v)) {if(!row[start+3]->Contains("*")) throw gcnew FormatException("Invalid Gaussian Anisotropic Couplings value."); complete=false;break;}
                    for(int a=0;a<3;++a) if(!number(row[start+5+a],axis[a])) complete=false;
                    for(int a=0;a<3;++a) for(int b=0;b<3;++b) m[a*3+b]+=v*axis[a]*axis[b];
                }
                if(complete) dipolar[id]=m;
                k+=3; while(k<lines->Length && String::IsNullOrWhiteSpace(lines[k])) ++k;
            }
            n=k-1;
        }
        else if(oniom && line=="Anisotropic Hyperfine tensors:") {
            // ONIOM's integrated six-component table is in the original input frame,
            // with rows XX, XY, YY, XZ, YZ, ZZ and Gaussian centre IDs as columns.
            compositeDipolar->Clear();
            for(int k=n+1;k+6<lines->Length;) {
                Lines^ fields=words(lines[k]); auto ids=gcnew cli::array<int>(fields->Length); bool valid=ids->Length>0;
                for(int col=0;col<ids->Length;++col) if(!Int32::TryParse(fields[col],ids[col])) valid=false;
                if(!valid) break;
                AtomTable^ panel=gcnew AtomTable(); for each(int id in ids) panel[id]=gcnew Values(9);
                for(int row=0;row<6 && valid;++row) {
                    Lines^ values=words(lines[k+row+1]); int index;
                    if(values->Length!=ids->Length+1 || !Int32::TryParse(values[0],index) || index!=row+1) {valid=false;break;}
                    const int lower[6]={0,3,4,6,7,8}; int a=lower[row]/3,b=lower[row]%3;
                    for(int col=0;col<ids->Length;++col) {double value; if(!number(values[col+1],value)) {panel->Remove(ids[col]);continue;}
                        if(panel->ContainsKey(ids[col])) panel[ids[col]][3*a+b]=panel[ids[col]][3*b+a]=value;}
                }
                if(!valid) break;
                for each(int id in panel->Keys) compositeDipolar[id]=panel[id];
                k+=7; n=k-1;
            }
        }
        else if(line=="ONIOM Fermi contact:") {
            compositeContact->Clear();
            for(int k=n+2;k<lines->Length;++k) {Lines^ f=words(lines[k]); int id; double v;
                if(f->Length!=3 || !Int32::TryParse(f[0],id)) break;
                if(number(f[2],v)) compositeContact[id]=v;
                n=k;
            }
        }
        else if(line=="Nuclear quadrupole coupling constants [Chi] (MHz):") {
            quadrupole->Clear(); qmass->Clear(); Eigen::Matrix3d r; bool fitted=false;
            for(int k=n+1;k+3<lines->Length;k+=4) {
                Match^ m=Regex::Match(lines[k],"^\\s*(\\d+)\\s+[A-Za-z]+\\(?(\\d+)\\)?\\s*$");
                if(!m->Success) break;
                int id=Int32::Parse(m->Groups[1]->Value); Values^ vals=components(lines,k+1,'a'); n=k+3;
                if(vals==nullptr || !geometry->ContainsKey(id)) continue;
                Eigen::Matrix3d t; for(int j=0;j<9;++j) t(j/3,j%3)=vals[j];
                if(t.isZero(0)) continue;
                if(!fitted) {
                    if(principal==nullptr) throw gcnew FormatException("Missing Gaussian principal-axis coordinates for quadrupolar tensor.");
                    r=inertialRotation(geometry,principal,axes); fitted=true;
                }
                t=(r*t*r.transpose()).eval(); t-=Eigen::Matrix3d::Identity()*(t.trace()/3);
                for(int j=0;j<9;++j) vals[j]=t(j/3,j%3);
                quadrupole[id]=vals; qmass[id]=Int32::Parse(m->Groups[2]->Value);
            }
        }
    }
    if(!geometry->Count) geometry=fallbackGeometry(lines,multiplicity);
    if(!geometry->Count && chi==nullptr && g==nullptr) return false;
    AtomsDictionary::CheckAtomCount(geometry->Count);
    ArrayList^ library=Isotopes::getAllIsotopes(); Dictionary<int,Atom^>^ atoms=gcnew Dictionary<int,Atom^>();
    SystemModel->AtomCollection->Refreshing=false;
    try {
        for each(int id in geometry->Keys) {
            Values^ row=geometry[id]; Atom^ a=gcnew Atom(); int mass=massOf(hmass,id); if(!mass) mass=massOf(qmass,id);
            a->isotope=isotope(library,(int)row[0],massOf(masses,id),mass);
            a->X=row[1]; a->Y=row[2]; a->Z=row[3]; atoms[id]=a; SystemModel->AtomCollection->Add(a);
            if(shielding->ContainsKey(id)) {addTensor(SystemModel,a,a,InteractionKind::CShielding,Units::ppm,shielding[id]); ToShift=true;}
            if(quadrupole->ContainsKey(id)) addTensor(SystemModel,a,a,InteractionKind::Quadrupolar,Units::MHz,quadrupole[id]);
        }
        for each(Pair^ pair in coupling->Keys) if(atoms->ContainsKey(pair->Item1) && atoms->ContainsKey(pair->Item2))
            addTensor(SystemModel,atoms[pair->Item1],atoms[pair->Item2],InteractionKind::Jcoupling,Units::Hz,diagonal(coupling[pair]));
        // Convert the integrated atomic-unit spin densities once, using the actual
        // nucleus's magnetic moment. Contact contains the 8*pi/3 Fermi prefactor.
        if(multiplicity>1 && (compositeDipolar->Count || compositeContact->Count)) {
            Eigen::Matrix3d rotation=Eigen::Matrix3d::Identity();
            if(compositeDipolar->Count && input!=nullptr && input!=geometry)
                rotation=inertialRotation(geometry,input,nullptr);
            const double hbar=1.054571817e-34,bohr=5.29177210903e-11;
            for each(int id in atoms->Keys) {
                double scale=1e-7*hbar*atoms[id]->isotope->Magnetogyric*1e4/(bohr*bohr*bohr);
                if(Double::IsNaN(scale) || Double::IsInfinity(scale)) continue;
                if(compositeContact->ContainsKey(id)) contact[id]=compositeContact[id]*scale*(8*acos(-1.0)/3)/(multiplicity-1);
                if(compositeDipolar->ContainsKey(id)) {
                    Eigen::Matrix3d t; for(int j=0;j<9;++j) t(j/3,j%3)=compositeDipolar[id][j];
                    t=(rotation*t*rotation.transpose()*scale).eval();
                    Values^ values=gcnew Values(9); for(int j=0;j<9;++j) values[j]=t(j/3,j%3);
                    dipolar[id]=values;
                }
            }
        }
        if(g!=nullptr || (multiplicity>1 && (contact->Count || dipolar->Count))) {
            Atom^ electron=gcnew Atom(); electron->isotope=Isotopes::FindIsotope("e",0,0); SystemModel->AtomCollection->Add(electron);
            addTensor(SystemModel,electron,electron,InteractionKind::GTensor,Units::Mi_b,g!=nullptr?g:diagonal(2.0023193043));
            for each(int id in atoms->Keys) {
                bool hasDipolar=multiplicity>1 && dipolar->ContainsKey(id);
                if(!contact->ContainsKey(id) && !hasDipolar) continue;
                Values^ h=hasDipolar?dipolar[id]:gcnew Values(9); double v;
                // Reduced output may first reveal the spin state in the closing archive.
                if(hasDipolar) for(int j=0;j<9;++j) h[j]/=(multiplicity-1);
                if(contact->TryGetValue(id,v)) {h[0]+=v; h[4]+=v; h[8]+=v;}
                addTensor(SystemModel,atoms[id],electron,InteractionKind::HFC,Units::Gauss,h);
            }
        }
        if(chi!=nullptr) {Atom^ carrier=gcnew Atom(); SystemModel->Susceptibility=gcnew Tensor(carrier,carrier,InteractionKind::Susceptibility,matrix(chi));}
    } finally {SystemModel->AtomCollection->Refreshing=true;}
    return true;
}
