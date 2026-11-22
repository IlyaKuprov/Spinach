#include "../StdAfx.h"
#include "Castep.h"

using namespace System::Collections::Generic;
using namespace System::Globalization;

typedef cli::array<String^> Fields;
typedef Dictionary<String^, Atom^> Sites;

namespace
{
    Fields^ fields(String^ line)
    {
        return line->Split((cli::array<Char>^)nullptr,StringSplitOptions::RemoveEmptyEntries);
    }

    double number(String^ token)
    {
        return Double::Parse(token->Replace(L'D',L'E')->Replace(L'd',L'e'),
                             NumberStyles::Float,CultureInfo::InvariantCulture);
    }

    String^ indexKey(String^ index)
    {
        return UInt32::Parse(index,NumberStyles::None,CultureInfo::InvariantCulture).ToString(CultureInfo::InvariantCulture);
    }

    // A glued alias can be ambiguous (e.g. O1 + 23 versus O12 + 3).
    Atom^ resolveSite(Fields^ row,int offset,int length,Sites^ sites)
    {
        String^ key=row[offset];
        if(length==2) {
            unsigned int index;
            if(!UInt32::TryParse(row[offset+1],NumberStyles::None,CultureInfo::InvariantCulture,index)) return nullptr;
            key+=" "+index.ToString(CultureInfo::InvariantCulture);
        }
        Atom^ atom;
        return sites->TryGetValue(key,atom) ? atom : nullptr;
    }

    Atom^ addSite(Model^ model,Sites^ sites,String^ species,String^ label,String^ index,Fields^ row,int start)
    {
        String^ id=indexKey(index);
        Atom^ atom=gcnew Atom();
        atom->isotope=Isotopes::FindIsotope(species->Split(':')[0],0,0);
        atom->Label=label+id;
        atom->X=number(row[start]);
        atom->Y=number(row[start+1]);
        atom->Z=number(row[start+2]);
        sites->Add(label+" "+id,atom);
        if(sites->ContainsKey(atom->Label)) sites[atom->Label]=nullptr;
        else sites->Add(atom->Label,atom);
        model->AtomCollection->Add(atom);
        return atom;
    }

    void store(Model^ model,Atom^ a,Atom^ b,InteractionKind kind,::Matrix3x3^ matrix,Unit^ unit)
    {
        Tensor^ tensor=gcnew Tensor(a,b,kind,matrix);
        tensor->unit=unit;
        tensor->Frame=model->RefFrameCollection[1];
        model->InteractionCollection->Add(tensor);
    }

    // castep2nqi: atomic EFG times the nuclear moment (barn), returned here in MHz.
    void quadrupole(Model^ model,Atom^ a,::Matrix3x3^ matrix,double moment)
    {
        double spin=a->isotope->Spin;
        if(spin<1) return;
        double scale=9.717362e21*1e-28*1.60217657e-19/6.62606957e-34*1e-6*moment/(2*spin*(2*spin-1));
        for(int k=0;k<9;k++) matrix[k]*=scale;
        store(model,a,a,InteractionKind::Quadrupolar,matrix,Units::MHz);
    }

    void checkUnits(Fields^ row)
    {
        if(row->Length!=3) throw gcnew FormatException("Invalid MAGRES units record.");
        String^ expected=nullptr;
        if(row[1]=="atom") expected="Angstrom";
        if(row[1]=="ms") expected="ppm";
        if(row[1]=="efg") expected="au";
        if(row[1]=="isc") expected="10^19.T^2.J^-1";
        if(expected!=nullptr && row[2]!=expected) throw gcnew FormatException("Unsupported "+row[1]+" units: "+row[2]);
    }
}

Castep::Castep(String^ filename,Model^% model)
{
    SystemModel=model;
    FileName=filename;
    ReadingFile=gcnew StreamReader(filename);
    FormatName="CASTEP";
}

bool Castep::LoadFile()
{
    SystemModel->AtomCollection->Refreshing=false;
    try {
        String^ first=ReadLine();
        while(first!=nullptr && (String::IsNullOrWhiteSpace(first) ||
              (first->TrimStart()->StartsWith("#") && !first->TrimStart()->StartsWith("#$")))) first=ReadLine();
        if(first==nullptr) return false;
        if(first->Trim()=="#$magres-abinitio-v1.0" || first->Trim()=="#$magres-abinitio") return NewFormat(nullptr);
        if(first->Trim()->StartsWith("#$magres-abinitio-v")) return false;
        first=first->Trim();
        if((first->StartsWith("<") && first->EndsWith(">")) || (first->StartsWith("[") && first->EndsWith("]")))
            return NewFormat(first->Substring(1,first->Length-2));
        return OldFormat(first);
    }
    catch(String^ message) {
        throw gcnew FormatException(message);
    }
    finally {
        ReadingFile->Close();
        SystemModel->AtomCollection->Refreshing=true;
    }
}

bool Castep::NewFormat(String^ block)
{
    Sites^ sites=gcnew Sites();
    List<Fields^>^ atomRows=gcnew List<Fields^>();
    List<Fields^>^ properties=gcnew List<Fields^>();
    String^ line;
    while((line=ReadLine())!=nullptr) {
        int hash=line->IndexOf('#');
        if(hash>=0) line=line->Substring(0,hash);
        line=line->Trim();
        if(line->Length==0) continue;
        if((line->StartsWith("[") && line->EndsWith("]")) || (line->StartsWith("<") && line->EndsWith(">"))) {
            String^ tag=line->Substring(1,line->Length-2);
            if(tag->StartsWith("/") && block==tag->Substring(1)) block=nullptr;
            else if(!tag->StartsWith("/") && block==nullptr) block=tag;
            else throw gcnew FormatException("Unbalanced MAGRES block: "+line);
            continue;
        }
        if(block!="atoms" && block!="magres") continue;
        Fields^ row=fields(line);
        if(row[0]=="units") checkUnits(row);
        else if(block=="atoms" && row[0]=="atom") {
            if(row->Length!=7) throw gcnew FormatException("Invalid MAGRES atom: "+line);
            atomRows->Add(row);
            AtomsDictionary::CheckAtomCount(atomRows->Count);
        } else if(block=="magres" && (row[0]=="ms" || row[0]=="efg" || row[0]=="isc")) properties->Add(row);
    }
    if(block!=nullptr) throw gcnew FormatException("Unclosed MAGRES block: "+block);
    for each(Fields^ row in atomRows) addSite(SystemModel,sites,row[1],row[2],row[3],row,4);
    if(sites->Count==0) return false;

    // Defer conversion of reduced K until all EFG records have selected their isotopes.
    typedef Tuple<int,int> Pair;
    Dictionary<Pair^,double>^ sums=gcnew Dictionary<Pair^,double>();
    Dictionary<Pair^,int>^ counts=gcnew Dictionary<Pair^,int>();
    for each(Fields^ row in properties) {
        bool pair=row[0]=="isc";
        int keys=row->Length-10;
        if(keys<1 || keys>(pair?4:2)) throw gcnew FormatException("Invalid MAGRES tensor width.");
        Atom^ a=nullptr; Atom^ b=nullptr;
        if(!pair) a=resolveSite(row,1,keys,sites);
        else for(int n=1;n<=2;n++) {
            int remaining=keys-n;
            if(remaining<1 || remaining>2) continue;
            Atom^ first=resolveSite(row,1,n,sites);
            Atom^ second=resolveSite(row,1+n,remaining,sites);
            if(first==nullptr || second==nullptr) continue;
            if(a!=nullptr) throw gcnew FormatException("Ambiguous MAGRES coupling sites.");
            a=first; b=second;
        }
        if(a==nullptr) throw gcnew FormatException("Unknown or ambiguous MAGRES site: "+String::Join(" ",row));
        ::Matrix3x3^ matrix=gcnew ::Matrix3x3();
        for(int k=0;k<9;k++) matrix[k]=number(row[keys+1+k]);
        if(row[0]=="ms") {
            store(SystemModel,a,a,InteractionKind::CShielding,matrix,Units::ppm);
            ToShift=true;
        } else if(!pair) {
            if(a->isotope->Spin<1) a->isotope=Isotopes::FindIsotopeWithHigherSpin(a->isotope);
            quadrupole(SystemModel,a,matrix,a->isotope->QuadrupoleMoment);
        } else if(a!=b) {
            Pair^ key=gcnew Pair(Math::Min(a->ID,b->ID),Math::Max(a->ID,b->ID));
            if(!sums->ContainsKey(key)) { sums->Add(key,0); counts->Add(key,0); }
            sums[key]+=(matrix[0]+matrix[4]+matrix[8])/3;
            counts[key]++;
        }
    }
    for each(Pair^ key in sums->Keys) {
        Atom^ a=SystemModel->AtomCollection[key->Item1];
        Atom^ b=SystemModel->AtomCollection[key->Item2];
        // J[Hz] = K[10^19 T^2/J] * 10^19 * h * gamma_A * gamma_B / (2*pi)^2.
        double hz=sums[key]/counts[key]*1e19*6.62607015e-34*
                  a->isotope->Magnetogyric*b->isotope->Magnetogyric/(4*PI*PI);
        ::Matrix3x3^ matrix=gcnew ::Matrix3x3();
        matrix[0]=hz; matrix[4]=hz; matrix[8]=hz;
        store(SystemModel,a,b,InteractionKind::Jcoupling,matrix,Units::Hz);
    }
    return true;
}

bool Castep::OldFormat(String^ line)
{
    Sites^ sites=gcnew Sites();
    Dictionary<String^,double>^ moments=gcnew Dictionary<String^,double>();
    Dictionary<String^,int>^ masses=gcnew Dictionary<String^,int>();
    Atom^ atom=nullptr;
    cli::array<double,2>^ lattice=nullptr;
    List<Fields^>^ geometry=nullptr;
    cli::array<double,2>^ geometryLattice=nullptr;
    while(line!=nullptr) {
        Fields^ row=fields(line);
        // The .castep report prints fractional coordinates; keep the last complete
        // table with its own real lattice, not a subsequent optimisation step's cell.
        if(line->Contains("Real Lattice(") && line->Contains("Reciprocal Lattice(")) {
            double scale;
            if(line->Contains("Real Lattice(A)")) scale=1.0;
            else if(line->Contains("Real Lattice(a0)")) scale=0.529177210903;
            else throw gcnew FormatException("Unsupported CASTEP lattice units: "+line);
            lattice=gcnew cli::array<double,2>(3,3);
            for(int i=0;i<3;i++) {
                line=ReadLine();
                if(line==nullptr || (row=fields(line))->Length<3) throw gcnew FormatException("Incomplete real lattice.");
                for(int j=0;j<3;j++) lattice[i,j]=number(row[j])*scale;
            }
        } else if(line->Contains("Fractional coordinates of atoms")) {
            if(lattice==nullptr) throw gcnew FormatException("No real lattice for fractional coordinates.");
            List<Fields^>^ frame=gcnew List<Fields^>();
            while((line=ReadLine())!=nullptr) {
                row=fields(line);
                if(row->Length==7 && row[0]=="x" && row[6]=="x") frame->Add(row);
                else if(frame->Count>0) break;
                else if(line->Contains("Number") || line->Contains("---") || String::IsNullOrWhiteSpace(line)) continue;
                else throw gcnew FormatException("Invalid CASTEP fractional coordinates.");
            }
            if(frame->Count>0 && line!=nullptr) {
                geometry=frame;
                geometryLattice=lattice;
            }
        } else if(row->Length>=7 && row[row->Length-4]=="Q" && row[row->Length-3]=="=" && row[row->Length-1]=="Barn") {
            moments[row[0]]=number(row[row->Length-2]);
            if(row[1]=="Isotope") masses[row[0]]=Int32::Parse(row[2],CultureInfo::InvariantCulture);
        } else if(row->Length==7 && row[2]=="Coordinates" && row[6]=="A") {
            if(!sites->TryGetValue(row[0]+" "+indexKey(row[1]),atom))
                atom=addSite(SystemModel,sites,row[0],row[0],row[1],row,3);
        } else if(atom!=nullptr && (line->Trim()=="TOTAL Shielding Tensor" || line->Trim()=="TOTAL tensor")) {
            bool shielding=line->Trim()=="TOTAL Shielding Tensor";
            ::Matrix3x3^ matrix=gcnew ::Matrix3x3();
            for(int i=0;i<3;i++) {
                do {line=ReadLine();} while(line!=nullptr && String::IsNullOrWhiteSpace(line));
                if(line==nullptr || (row=fields(line))->Length!=3) throw gcnew FormatException("Invalid legacy MAGRES tensor.");
                for(int j=0;j<3;j++) matrix[3*i+j]=number(row[j]);
            }
            if(shielding) {
                store(SystemModel,atom,atom,InteractionKind::CShielding,matrix,Units::ppm);
                ToShift=true;
            } else {
                String^ species=atom->isotope->Element;
                if(!moments->ContainsKey(species)) throw gcnew FormatException("No quadrupole moment for "+species);
                if(masses->ContainsKey(species)) atom->isotope=Isotopes::FindIsotope(species,0,masses[species]);
                else if(atom->isotope->Spin<1) atom->isotope=Isotopes::FindIsotopeWithHigherSpin(atom->isotope);
                quadrupole(SystemModel,atom,matrix,moments[species]);
            }
        }
        line=ReadLine();
    }
    if(geometry!=nullptr && sites->Count==0) {
        AtomsDictionary::CheckAtomCount(geometry->Count);
        for each(Fields^ row in geometry) {
            Atom^ site=gcnew Atom();
            site->isotope=Isotopes::FindIsotope(row[1]->Split(':')[0],0,0);
            site->Label=row[1]+row[2];
            double u=number(row[3]), v=number(row[4]), w=number(row[5]);
            site->X=u*geometryLattice[0,0]+v*geometryLattice[1,0]+w*geometryLattice[2,0];
            site->Y=u*geometryLattice[0,1]+v*geometryLattice[1,1]+w*geometryLattice[2,1];
            site->Z=u*geometryLattice[0,2]+v*geometryLattice[1,2]+w*geometryLattice[2,2];
            SystemModel->AtomCollection->Add(site);
        }
    }
    return SystemModel->AtomCollection->Count>0;
}
