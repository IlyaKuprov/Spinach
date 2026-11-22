#include "../StdAfx.h"
#include "Orca.h"

using namespace System::Collections::Generic;
using namespace System::Globalization;
using namespace System::Text::RegularExpressions;

namespace
{
    // Keep original ORCA atom indices, even when a dummy or capped site is
    // omitted from the physical spin system. Assign GUI IDs only at the end.
    ref class OrcaSite
    {
    public:
        String^ symbol;
        double x, y, z;
        int mass;
        Atom^ atom;
    };

    ref class NuclearProperties
    {
    public:
        String^ symbol;
        int mass, qMass;
        double spin, quadrupoleMoment;
        ::Matrix3x3^ hfc;
        ::Matrix3x3^ efg;
    };

    bool IsPhysicalSite(OrcaSite^ site)
    {
        return Regex::IsMatch(site->symbol,"^[A-Za-z]{1,2}$") && site->symbol!="X" && site->symbol!="DA";
    }

    bool Number(String^ text, double% value)
    {
        return Double::TryParse(text->Replace(L'D',L'E')->Replace(L'd',L'e'),
            NumberStyles::Float,CultureInfo::InvariantCulture,value) &&
            !Double::IsNaN(value) && !Double::IsInfinity(value);
    }

    cli::array<String^>^ Fields(String^ line)
    {
        return line->Split((cli::array<Char>^)nullptr,StringSplitOptions::RemoveEmptyEntries);
    }

    bool IsSectionRule(String^ line)
    {
        return line->Length>0 && line->Trim(L'-',L'=')->Length==0;
    }

    OrcaSite^ ParseSite(String^ line, bool atomicUnits, double scale)
    {
        auto f=Fields(line);
        int first=atomicUnits ? 5 : 1;
        if(f->Length!=first+3) return nullptr;
        OrcaSite^ site=gcnew OrcaSite();
        if(!Number(f[first],site->x) || !Number(f[first+1],site->y) || !Number(f[first+2],site->z)) return nullptr;
        site->symbol=f[atomicUnits ? 1 : 0];
        site->x*=scale; site->y*=scale; site->z*=scale;
        double mass;
        if(atomicUnits && Number(f[4],mass)) site->mass=Convert::ToInt32(mass);
        return site;
    }

    Isotope^ SiteIsotope(String^ symbol, int mass, bool exact)
    {
        // Atomic weights in the geometry table are not isotope declarations.
        // Use an exact table match when available, otherwise the GUI default.
        for each(Isotope^ candidate in Isotopes::getAllIsotopes())
            if(candidate->Element==symbol && candidate->Mass==mass) return candidate;
        if(exact) throw gcnew FormatException("ORCA requests "+mass+symbol+", absent from the GUI isotope table.");
        try { return Isotopes::FindIsotope(symbol,0,0); }
        catch(String^) { throw gcnew FormatException("ORCA element is not supported by the GUI: "+symbol); }
    }
}

Orca::Orca(String^ filename, Model^% model)
{
    SystemModel=model;
    FileName=filename;
    FormatName="Orca";
    ReadingFile=gcnew StreamReader(filename);
}

::Matrix3x3^ Orca::ReadMatrix()
{
    ::Matrix3x3^ matrix=gcnew ::Matrix3x3();
    for(int row=0; row<3; ++row)
    {
        String^ line;
        do {
            if(ReadingFile->EndOfStream) throw gcnew FormatException("Incomplete ORCA tensor at line "+LineCount);
            line=ReadLine()->Trim();
        } while(row==0 && (line=="" || IsSectionRule(line)));
        auto fields=Fields(line);
        if(fields->Length!=3) throw gcnew FormatException("Invalid ORCA tensor at line "+LineCount);
        for(int col=0; col<3; ++col)
        {
            double value;
            if(!Number(fields[col],value)) throw gcnew FormatException("Invalid ORCA tensor at line "+LineCount);
            matrix[3*row+col]=value;
        }
    }
    return matrix;
}

void Orca::AddTensor(Atom^ a, Atom^ b, InteractionKind kind, ::Matrix3x3^ matrix, Unit^ unit)
{
    Tensor^ tensor=gcnew Tensor(a,b,kind,matrix);
    tensor->unit=unit;
    tensor->Frame=SystemModel->RefFrameCollection[1];
    SystemModel->InteractionCollection->Add(tensor);
}

bool Orca::LoadFile()
{
    auto sites=gcnew List<OrcaSite^>();
    auto inputSites=gcnew List<OrcaSite^>();
    auto nuclear=gcnew SortedDictionary<int,NuclearProperties^>();
    auto shielding=gcnew SortedDictionary<int,::Matrix3x3^>();
    auto pairs=gcnew Dictionary<Tuple<int,int>^,double>();
    ::Matrix3x3^ g=nullptr;
    ::Matrix3x3^ zfs=nullptr;
    bool recognised=false, precedingAngstrom=false, inputBohr=false;
    enum class Section { None, G, ZFS, NMR, Hyperfine, J };
    Section section=Section::None;
    Regex^ nucleus=gcnew Regex("^Nucleus\\s+(\\d+)\\s*([a-zA-Z]+)");
    Regex^ massHeader=gcnew Regex("(?:ISTP|Isotope)\\s*=\\s*(\\d+)");
    Regex^ quadHeader=gcnew Regex("I=\\s*(\\S+)\\s+Q=\\s*(\\S+)\\s+barn");
    Regex^ pair=gcnew Regex("^NUCLEUS A\\s*=\\s*\\S+\\s+(\\d+)\\s+NUCLEUS B\\s*=\\s*\\S+\\s+(\\d+)");
    int atom=-1, other=-1;
    String^ pending=nullptr;

    try
    {
        while(pending!=nullptr || !ReadingFile->EndOfStream)
        {
            String^ line=pending==nullptr ? ReadLine()->Trim() : pending;
            pending=nullptr;
            if(line->StartsWith("Program Version "))
            {
                // Concatenated outputs can describe unrelated molecules.
                if(recognised)
                {
                    sites->Clear(); inputSites->Clear(); nuclear->Clear(); shielding->Clear(); pairs->Clear();
                    g=nullptr; zfs=nullptr; section=Section::None;
                    inputBohr=false; precedingAngstrom=false;
                }
                recognised=true;
                continue;
            }
            if(line->StartsWith("|") && line->Contains(">"))
            {
                String^ input=line->Substring(line->IndexOf('>')+1)->Trim();
                if(input->StartsWith("!") && Regex::IsMatch(input,"(?i)\\bbohrs?\\b")) inputBohr=true;
                if(Regex::IsMatch(input,"(?i)^\\*\\s*xyz\\s+"))
                {
                    inputSites->Clear();
                    while(!ReadingFile->EndOfStream)
                    {
                        line=ReadLine();
                        if(!line->Contains(">")) break;
                        input=line->Substring(line->IndexOf('>')+1)->Trim();
                        if(input=="*") break;
                        OrcaSite^ site=ParseSite(input,false,inputBohr ? 0.529177249 : 1.0);
                        if(site==nullptr) { inputSites->Clear(); break; }
                        inputSites->Add(site);
                    }
                }
                continue;
            }
            if(line=="CARTESIAN COORDINATES (ANGSTROEM)" || line=="CARTESIAN COORDINATES (A.U.)")
            {
                bool au=line->Contains("A.U.");
                auto frame=gcnew List<OrcaSite^>();
                while(!ReadingFile->EndOfStream)
                {
                    line=ReadLine()->Trim();
                    if(frame->Count==0 && (line=="" || IsSectionRule(line) || line->StartsWith("NO "))) continue;
                    OrcaSite^ site=ParseSite(line,au,au ? 0.529177249 : 1.0);
                    if(site==nullptr) { pending=line; break; }
                    frame->Add(site);
                }
                if(frame->Count)
                {
                    if(au && precedingAngstrom && frame->Count==sites->Count)
                    {
                        for(int i=0;i<sites->Count;i++) sites[i]->mass=frame[i]->mass;
                    }
                    else sites=frame;
                    precedingAngstrom=!au; // An A.U. table pairs with the preceding angstrom table.
                }
                section=Section::None;
                continue;
            }

            if(line=="ELECTRONIC G-MATRIX" || line=="ELECTRONIC G-MATRIX FROM EFFECTIVE HAMILTONIAN")
            { section=Section::G; g=nullptr; }
            else if(line->StartsWith("ZERO-FIELD-SPLITTING TENSOR"))
            { section=Section::ZFS; zfs=nullptr; }
            else if(line=="CHEMICAL SHIFTS" || line->StartsWith("CHEMICAL SHIELDINGS"))
            { section=Section::NMR; shielding->Clear(); atom=-1; }
            else if(line->StartsWith("ELECTRIC AND MAGNETIC HYPERFINE STRUCTURE"))
            { section=Section::Hyperfine; nuclear->Clear(); atom=-1; }
            else if(line->StartsWith("NMR SPIN-SPIN COUPLING CONSTANTS"))
            { section=Section::J; pairs->Clear(); atom=other=-1; }
            else if(line->StartsWith("MULLIKEN ATOMIC CHARGES") || line->StartsWith("LOEWDIN ATOMIC CHARGES") ||
                    line->StartsWith("ORCA TERMINATED") || line->Contains("****ORCA TERMINATED") ||
                    line->StartsWith("ORCA PROPERTY CALCULATIONS") || line->StartsWith("ORCA EULER ANGLE PROGRAM"))
                section=Section::None;
            else if(section==Section::G && (line=="The g-matrix:" || line=="g-matrix:"))
            { g=ReadMatrix(); section=Section::None; }
            else if(section==Section::ZFS && line->StartsWith("raw-matrix"))
            { zfs=ReadMatrix(); section=Section::None; }
            else if((section==Section::NMR || section==Section::Hyperfine) && line->StartsWith("Nucleus"))
            {
                Match^ m=nucleus->Match(line);
                atom=m->Success ? Int32::Parse(m->Groups[1]->Value) : -1;
                if(atom>=0 && section==Section::Hyperfine)
                {
                    NuclearProperties^ prop=gcnew NuclearProperties();
                    prop->symbol=m->Groups[2]->Value;
                    m=massHeader->Match(line);
                    if(m->Success) prop->mass=Int32::Parse(m->Groups[1]->Value);
                    nuclear[atom]=prop;
                }
            }
            else if(section==Section::NMR && atom>=0 &&
                    (line->StartsWith("Total shielding tensor") || line->StartsWith("Raw-matrix")))
            {
                auto matrix=ReadMatrix();
                // Early ORCA chemical-shift sections print dimensionless shielding.
                if(line->StartsWith("Raw-matrix"))
                    for(int i=0;i<9;i++) matrix[i]*=1e6;
                shielding[atom]=matrix;
            }
            else if(section==Section::Hyperfine && atom>=0)
            {
                auto prop=nuclear[atom];
                if(line->StartsWith("Q") && line->Contains("barn"))
                {
                    Match^ m=massHeader->Match(line);
                    if(m->Success) prop->qMass=Int32::Parse(m->Groups[1]->Value);
                    m=quadHeader->Match(line);
                    if(!m->Success || !Number(m->Groups[1]->Value,prop->spin) ||
                       !Number(m->Groups[2]->Value,prop->quadrupoleMoment))
                        throw gcnew FormatException("Invalid ORCA quadrupole isotope at line "+LineCount);
                }
                else if(line->StartsWith("Total HFC matrix") || line->StartsWith("Raw HFC matrix"))
                    prop->hfc=ReadMatrix();
                else if(line->StartsWith("Total EFG matrix") || line->StartsWith("Raw EFG matrix"))
                    prop->efg=ReadMatrix();
            }
            else if(section==Section::J && line->StartsWith("NUCLEUS A"))
            {
                Match^ m=pair->Match(line);
                atom=m->Success ? Int32::Parse(m->Groups[1]->Value) : -1;
                other=m->Success ? Int32::Parse(m->Groups[2]->Value) : -1;
            }
            else if(section==Section::J && atom>=0 && line->StartsWith("Total spin-spin coupling tensor"))
            {
                auto matrix=ReadMatrix();
                pairs[Tuple::Create(Math::Min(atom,other),Math::Max(atom,other))]=(matrix[0]+matrix[4]+matrix[8])/3;
            }
        }
    }
    finally { ReadingFile->Close(); }

    if(!recognised) return false;
    if(!sites->Count) sites=inputSites;
    if(!sites->Count) throw gcnew FormatException("ORCA output contains no Cartesian geometry.");
    int atomCount=0;
    for each(OrcaSite^ site in sites) if(IsPhysicalSite(site)) ++atomCount;
    AtomsDictionary::CheckAtomCount(atomCount);
    bool hasHfc=false;
    for(int i=0;i<sites->Count;i++)
    {
        OrcaSite^ site=sites[i];
        if(!IsPhysicalSite(site)) continue;
        NuclearProperties^ prop;
        int mass=site->mass;
        bool exact=false;
        if(nuclear->TryGetValue(i,prop))
        {
            if(!String::Equals(site->symbol,prop->symbol,StringComparison::OrdinalIgnoreCase))
                throw gcnew FormatException("ORCA nucleus disagrees with geometry at atom "+i);
            if(prop->hfc!=nullptr)
            {
                hasHfc=true;
                if(prop->mass) { mass=prop->mass; exact=true; }
            }
            else if(prop->efg!=nullptr && prop->qMass) { mass=prop->qMass; exact=true; }
        }
        auto a=gcnew Atom();
        String^ symbol=site->symbol->Substring(0,1)->ToUpperInvariant()+site->symbol->Substring(1)->ToLowerInvariant();
        a->isotope=SiteIsotope(symbol,mass,exact);
        a->X=site->x; a->Y=site->y; a->Z=site->z;
        SystemModel->AtomCollection->Add(a); site->atom=a;
    }
    Atom^ electron=nullptr;
    if(g!=nullptr || zfs!=nullptr || hasHfc)
    {
        electron=gcnew Atom();
        electron->isotope=Isotopes::FindIsotope("e",0,0);
        SystemModel->AtomCollection->Add(electron);
        if(g==nullptr && hasHfc)
        {
            g=gcnew ::Matrix3x3();
            for(int i=0;i<9;i+=4) g[i]=2.0023;
        }
        if(g!=nullptr) AddTensor(electron,electron,InteractionKind::GTensor,g,Units::Mi_b);
        if(zfs!=nullptr) AddTensor(electron,electron,InteractionKind::ZFS,zfs,Units::rev_cm);
    }
    for each(int i in shielding->Keys)
    {
        Atom^ a=sites[i]->atom;
        if(a!=nullptr) { AddTensor(a,a,InteractionKind::CShielding,shielding[i],Units::ppm); ToShift=true; }
    }
    for each(int i in nuclear->Keys)
    {
        Atom^ a=sites[i]->atom;
        if(a==nullptr) continue;
        NuclearProperties^ prop=nuclear[i];
        if(prop->hfc!=nullptr) AddTensor(a,electron,InteractionKind::HFC,prop->hfc,Units::MHz);
        // ORCA may print HFC and EFG for different isotopes at the same site.
        // A GUI atom describes one isotope. Never assign e.g. a 2H NQI to 1H.
        if(prop->efg!=nullptr && prop->spin>0.5 && prop->quadrupoleMoment!=0 &&
           a->isotope->Mass==prop->qMass)
        {
            const double efgToMHz=9.717362e21*1e-28*1.60217657e-19/6.62606957e-34/1e6;
            double scale=efgToMHz*prop->quadrupoleMoment/(2*prop->spin*(2*prop->spin-1));
            auto m=prop->efg;
            double iso=(m[0]+m[4]+m[8])/3;
            for(int k=0;k<9;k++) m[k]=scale*(m[k]-(k%4==0 ? iso : 0));
            AddTensor(a,a,InteractionKind::Quadrupolar,m,Units::MHz);
        }
    }
    for each(auto key in pairs->Keys)
    {
        Atom^ a=sites[key->Item1]->atom;
        Atom^ b=sites[key->Item2]->atom;
        if(a==nullptr || b==nullptr) continue;
        auto m=gcnew ::Matrix3x3();
        for(int i=0;i<9;i+=4) m[i]=pairs[key];
        AddTensor(a,b,InteractionKind::Jcoupling,m,Units::Hz);
    }
    return SystemModel->AtomCollection->Count>0;
}
