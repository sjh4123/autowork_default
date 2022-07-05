/*******************************************************************************
 PreSQ & AUTOWORK
 Copyright (c) Jui-Hsiang (Sean) Hung. All rights reserved.
 
 >>> SOURCE LICENSE >>>
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation (www.fsf.org); either version 2 of the
 License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 A copy of the GNU General Public License is available at
 http://www.fsf.org/licensing/licenses
 >>> END OF LICENSE >>>
 ******************************************************************************/

#include "moltemplatelmpdata.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

MoltemplateLmpData::MoltemplateLmpData(const StructureClass& sysVar):

/* bool */
is_restrictbynAtoms(true),
is_optimizeCharge(false),
is_mass_normalization(false),
is_correctAlkylDihedral(false),
is_HCHexisting(false),

/* int */
dop(1),
indexi(0),
indexii(0),
n_trials(sysVar.get_n_trial()),
n_total_atoms(0),
sequenceNum(24),
sequenceLen(30),
n_light(0),
n_heavy(0),

n_HCHbondType(0),
n_totalatomTypes(0),
atomType_CH3(0),
atomType_CH2(0),
atomType_CH(0),
atomType_H(0),

/* double */
FFmodify_CH3_charge(-0.222),
FFmodify_CH2_charge(-0.148),
FFmodify_CH_charge(-0.160),
FFmodify_H_charge(0.074),
FFmodify_HCH_charge(0.160),

moltemplateBoxSize(300.0),/** OPTIMIZED: will change accord. to actual packing **/
offset(4.0),   /** OPTIMIZED: will change accord. to constituent mers size **/
rotate(90.0),  /** OPTIMIZED: rotate this angle between consecutive mers **/
packingL(10.0),/** OPTIMIZED: will be equal to offset **/
offset_spacing(2.0),
packingL_spacing(5.0),

/* string */
path_cwd(return_SimulationFolderPath(sysVar)+"/lammps_inputs/start_data/moltemplate/"),
path_master(sysVar.get_Path()+"/"),
path_moltemplatesrc(path_master+"moltemplate/src/"),
path_oplsaaprm(path_master+"moltemplate/oplsaa.prm"),
path_monomerBank(path_master+"Monomer_bank/"),
path_sequenceBank(path_master+"sequence_bank/"),
path_scripts(path_master+"scripts/"),
merA(""),
merB(""),
copolymerType("random"), // "random" or "block"
tacticity("atactic"),    // "atactic" or "isotactic" or "syndiotactic"

/** STL containers **/
//------------------------------------------------------------------------------
merSet(),
chem_mono(),
merComp(),
mono_beads(),
sequenceSet(),
types_all(),
types_light(),
types_heavy(),
HCHatomIDs(),
FFmodify_alkylDihedral({-0.305938,2.697394,-0.896807,0.74567}),
n_types_all(),
n_types_light(),
n_types_heavy(),
HCHbondType()
{
    if (sysVar.get_is_moltempLmpData()) check_cwd(sysVar);
}





vector<vector<string>> MoltemplateLmpData::make_atacticSequenceSet()
{
    vector<vector<string>> sequenceSet;
    //AA_Monomer sequences(*this);
    //sequenceSet=sequences.determine_sequence(*this,chem_mono,mono_beads);
    sequenceNum=(int)sequenceSet.size();
    return sequenceSet;
}





int MoltemplateLmpData::n_monomerAtoms(const string& merltfile)
{
    int n_monomerAtoms=0;
    string in=path_monomerBank+merltfile;
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        string lineContent;
        string strData0;
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            if (strData0=="write(\"Data") {
                is_inside_block=true;
                getline(readFile,lineContent);//1st atom line
            } else if (strData0=="}") {
                is_inside_block=false;
            }
            if (is_inside_block) {
                ++n_monomerAtoms;
            }
        }
    } else {
        cout
        << "in MoltemplateLmpData::n_monomerAtoms():\n"
        << in <<" file cannot open.\n";
        exit(EXIT_FAILURE);
    } return n_monomerAtoms;
}





vector<vector<string>> MoltemplateLmpData::make_SequenceSet_adhoc()
{
    vector<vector<string>> sequenceSet;
    vector<string> merSet;
    int counter=0;
    int n_chainAtoms=0;
    int n_chains=0;
    int chain=0,mer=0;
    
    /** Loops through monomer types per chain **/
    for (indexi=0;indexi<(int)chem_mono.size();++indexi) {
        counter += mono_beads.at(indexi);
        if (sequenceLen>1) {
            n_chainAtoms += n_monomerAtoms(chem_mono.at(indexi)+"le.lt")*mono_beads.at(indexi);
        } else {
            n_chainAtoms += n_monomerAtoms(chem_mono.at(indexi)+".lt")*mono_beads.at(indexi);
        }
        //cout << n_monomerAtoms(chem_mono.at(indexi)+"le.lt") << "\n";
        //NOTE: approximate number of monomer atoms by le ltfile
    }
    /** sanity check **/
    if (counter!=sequenceLen) {
        cout
        << "in MoltemplateLmpData::make_SequenceSet_adhoc()\n"
        << "Number of mers per chain ("<<counter<<") "
        << "!= sequenceLen ("<<sequenceLen<<"). Please check.\n";
        exit(EXIT_FAILURE);
    }
    /** Determine number of chains **/
    if (n_total_atoms==0) {/** controlled by sequenceNum **/
        n_chains=sequenceNum;
    } else {/** controlled by n_total_atoms **/
        if (n_chainAtoms>n_total_atoms) {
            cout
            << "in MoltemplateLmpData::make_SequenceSet_adhoc()\n"
            << "n_chainAtoms("<<n_chainAtoms<<") > n_total_atoms("<<n_total_atoms<<")\n"
            << "The upper threshold of total number of atoms in system is reached.\n"
            << "Please either increase n_total_atoms or choose a shorter chain.\n";
            exit(EXIT_FAILURE);
        } n_chains=(int)((double)n_total_atoms/(double)n_chainAtoms);//floor
        //cout << n_chains << "\n";
    }
    /** sanity check **/
    if (n_chains<=0) {
        cout
        << "in MoltemplateLmpData::make_SequenceSet_adhoc()\n"
        << "n_chain = "<<n_chains<<", which is invalid. Please check.\n";
        exit(EXIT_FAILURE);
    }
    
    /** For generating small molecule **/
    //--------------------------------------------------------------------------
    if (sequenceLen==1)
    {
        if (chem_mono.size()>1) {
            cout
            << "sequenceLen = "<<sequenceLen<<", "
            << "merSet should only have one mer type! Please check.\n";
            exit(EXIT_FAILURE);
        }
        for (chain=0;chain<n_chains;++chain) {
            merSet.clear();
            for (mer=0;mer<sequenceLen;++mer) {
                merSet.push_back(chem_mono[0]+".lt");
            } sequenceSet.push_back(merSet);
        } return sequenceSet;
    }
    //--------------------------------------------------------------------------
    
    if (chem_mono.size()==1) { /** Homopolymer **/
        
        string chosenTac;
        if (randbool()) chosenTac="_T1.lt";
        else chosenTac=".lt";
        
        for (chain=0;chain<n_chains;++chain) {
            
            merSet.clear();
            
            if (tacticity=="atactic") /** random chirality **/
            {
                for (mer=0;mer<sequenceLen;++mer)
                {
                    if (mer==0) {//1st monomer: le ltfile
                        if (randbool()) {//true: T1
                            merSet.push_back(chem_mono[0]+"le_T1.lt");
                        } else {//false: normal
                            merSet.push_back(chem_mono[0]+"le.lt");
                        }
                    } else if (mer==(sequenceLen-1)) {//last monomer: re ltfile
                        if (randbool()) {//true: T1
                            merSet.push_back(chem_mono[0]+"re_T1.lt");
                        } else {//false: normal
                            merSet.push_back(chem_mono[0]+"re.lt");
                        }
                    } else {//intermediate monomers: i ltfiles
                        if (randbool()) {//true: T1
                            merSet.push_back(chem_mono[0]+"i_T1.lt");
                        } else {//false: normal
                            merSet.push_back(chem_mono[0]+"i.lt");
                        }
                    }
                }
            }
            else if (tacticity=="isotactic") /** monotonic chirality **/
            {
                for (mer=0;mer<sequenceLen;++mer)
                {
                    if (mer==0) {//1st monomer: le ltfile
                        merSet.push_back(chem_mono[0]+"le"+chosenTac);
                    } else if (mer==(sequenceLen-1)) {//last monomer: re ltfile
                        merSet.push_back(chem_mono[0]+"re"+chosenTac);
                    } else {//intermediate monomers: i ltfiles
                        merSet.push_back(chem_mono[0]+"i"+chosenTac);
                    }
                }
            }
            else if (tacticity=="syndiotactic") /** alternating chirality **/
            {
                string startTac,nextTac,currentTac;
                /** coin flip to decide start tacticity **/
                if (randbool()) {
                    startTac="_T1.lt";
                    nextTac =".lt";
                } else {
                    startTac=".lt";
                    nextTac ="_T1.lt";
                }
                for (mer=0;mer<sequenceLen;++mer)
                {
                    if (mer%2==0) {
                        currentTac=startTac;
                    } else {
                        currentTac=nextTac;
                    }
                    if (mer==0) {//1st monomer: le ltfile
                        merSet.push_back(chem_mono[0]+"le"+currentTac);
                    } else if (mer==(sequenceLen-1)) {//last monomer: re ltfile
                        merSet.push_back(chem_mono[0]+"re"+currentTac);
                    } else {//intermediate monomers: i ltfiles
                        merSet.push_back(chem_mono[0]+"i"+currentTac);
                    }
                }
            }
            else
            {
                cout
                << "in MoltemplateLmpData::make_SequenceSet_adhoc()\n"
                << "Please choose tacticity from {atactic,syndiotactic,isotactic}\n";
                exit(EXIT_FAILURE);
            }
            /** load in tacticity of current chain **/
            sequenceSet.push_back(merSet);
        }
        
    } else {
        
        if (copolymerType=="random") { /** Random Copolymer **/
            //TODO
        }
        
        if (copolymerType=="block") { /** Block Copolymer **/
            //TODO
        }
        
    } return sequenceSet;
}





void MoltemplateLmpData::make_LmpDataFilebyMoltemplate(StructureClass& sysVar)
{
    cout << "\nlmpdata prepared by Moltemplate.\n";system_wait(1);
    
    /** check number of molecules **/
    //if(sequenceNum>0) {
    //if ((int)sequenceSet.size()!=sequenceNum) {
    //    cout
    //    << "\nWarning: Number of Molecules != " << sequenceNum
    //    << "\n\n";
    //}
    //}
    cout << "\nNumber of Molecules = "<<sequenceSet.size()<<"\n\n";
    system_wait(1);
    
    /** set n_poly (now set n_poly in set_atomNumbers()) **/
    //sysVar.set_n_poly((int)sequenceSet.size());
    
    /** loop through all polymers and make corresponding polymer lt files **/
    for (indexi=0; indexi<(int)sequenceSet.size(); ++indexi) {
        
        /** check degrees of polymerization (DOP) of current polymer **/
        if (sequenceLen>0) {
            if (sequenceSet.at(indexi).size()!=sequenceLen) {
                cout
                << "\nWarning: At molecule#"<<indexi+1<<", "
                << "DOP="<<sequenceSet.at(indexi).size()<<" != "<<sequenceLen
                << "\n\n";
            }
        } else {
            cout
            << "MoltemplateLmpData::make_LmpDataFilebyMoltemplate():\n"
            << "sequenceLen("<<sequenceLen<<") should be >0. Please check.\n\n";
            exit(EXIT_FAILURE);
        }
        /** check monomer.lt's in the monomer bank **/
        for (indexii=0; indexii<(int)sequenceSet.at(indexi).size(); ++indexii)
        {
            if (check_monomerbank(sequenceSet.at(indexi).at(indexii))) {
                /** copy found monomer.lt from monomer bank to cwd **/
                string source=path_monomerBank+sequenceSet.at(indexi).at(indexii);
                copy_to_cwd(source);
            } else {
                /** if monomer.lt is NOT found, program will exit **/
                cout
                << "\nError: "
                << "At monomer#" <<indexii+1<<" ("<<sequenceSet.at(indexi).at(indexii)<<") "
                << "of molecule#"<<indexi+1 <<": "
                << "\nCan't find corresponding lt file in the monomer bank.\n\n"
                << "Automation terminated.\n\n";
                exit(EXIT_FAILURE);
            }
        }
        /** make oplsaa.lt (requiring oplsaa_subset.prm)
         NOTE: the oplsaa.lt is shared by the current polymer and all its
         constituent monomers; to make it more general, this function should
         generate an unique oplsaa.lt file for each polymer and its own
         monomers. This feature is NOT supported in the current version **/
        make_oplsaalt(sequenceSet.at(indexi));
        
        /** make poly.lt file **/
        if (sequenceLen>1) {
            make_polylt(indexi,sequenceSet.at(indexi));
        }
    }
    
    /** make system.lt file **/
    make_systemlt();
    
    if (sequenceLen>1)
    {
        /** correct alkyl dihedrals & charges in oplsaa.lt **/
        if (is_correctAlkylDihedral) {
            FFmodify_alkyl_dihedral_oplsaa();
            FFmodify_alkyl_charge_oplsaa();
            FFmodify_alkyl_findHCHbond_oplsaa();//find if H-CH exists
        }
    }
    
    /** invoke moltemplate to generate LAMMPS datafile **/
    invoke_moltemplate();
    
    /** evaluate proper starting box size **/
    evaluate_boxLen();
    
    if (sequenceLen>1)
    {
        /** correct hydrogens in system.data & system.in.charges **/
        if (is_correctAlkylDihedral) {
            if (is_HCHexisting) {
                FFmodify_alkyl_findHCHinBonds_systemdata();
                FFmodify_alkyl_writeHCHinMasses_systemdata();
                FFmodify_alkyl_writeHCHinAtoms_systemdata();
                FFmodify_alkyl_writeHCHincharges_systemcharges();
                FFmodify_alkyl_writeHCHinsettings_systemsettings();
            }
        }
        /** charge optimization (neutralization) **/
        if (is_optimizeCharge) optimize_charge();
    }
    
    if ((sequenceLen==1)||(!is_optimizeCharge))
    {
        getridof_full();
        getridof_ljcutcoullong();
    }
    
    /** renormalize atom mass **/
    if (is_mass_normalization) mass_normalization();
    
    /** move output files to due places **/
    mv_files();
    
    cout << "\nmoltemplate-enabled lmpdata made!\n";system_wait(1);
}





vector<vector<string>> MoltemplateLmpData::read_sequenceSetfromfile(const StructureClass& sysVar)
{
    vector<string> vs;
    vector<vector<string>> vvs;
    int countline=0,countnum=0;
    
    string input;
    input.append(sysVar.get_Path());
    input.append("/sequenceSet.txt");
    
    ifstream readFile(input.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        string monomer;
        while (getline(readFile,lineContent)) {
            ++countline;
            istringstream iss(lineContent);
            countnum=0; vs.clear();
            while (iss>>monomer) {
                ++countnum;
                vs.push_back(monomer);
            } vvs.push_back(vs);
        } readFile.close();
    } else {
        cout << "\nMoltemplateLmpData::read_sequenceSet: readFile can't open.\n";
    } return vvs;
}





bool MoltemplateLmpData::check_monomerbank(const string& monomer)
{
    string file;
    file.append(path_monomerBank);
    file.append(monomer);
    ifstream readFile(file.c_str());
    return readFile.good();
}





bool MoltemplateLmpData::check_is_heavyAtom(const int atomType)
{
    bool is_heavy=false;
    vector<int>::iterator itr;
    itr=find(types_heavy.begin(),types_heavy.end(),atomType);
    if (itr!=types_heavy.end()) {
        is_heavy=true;
    } else {
        is_heavy=false;
    }
    return is_heavy;
}





void MoltemplateLmpData::check_cwd(const StructureClass& sysVar)
{
    const string path_d   =sysVar.get_Path();
    const string simType_d=sysVar.get_simType();
    const string year_d   =sysVar.get_year();
    const string usicID_d =sysVar.get_usicID();
    const string usic_d   =sysVar.get_usic();
    string mkdir=
    "Automation="+path_d+";cd ${Automation};"
    "usic_dir=${Automation}/simulations/"+simType_d+"/"+year_d+"/"+usic_d+";"
    "simulations_dir=${usic_dir}/simulations;"
    "lammps_inputs=${simulations_dir}/lammps_inputs;"
    "moltemp_dir=${lammps_inputs}/start_data/moltemplate;"
    "test -e ${moltemp_dir};"
    "if [ $? -ne 0 ]; then "
    "   mkdir -p ${moltemp_dir};"
    "fi;";
    int ret=system(mkdir.c_str());
    if (ret!=0) {
        cout
        << "in MoltemplateLmpData::check_cwd():\n"
        << "system() returns "<<ret<<"\n\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::check_path(const string& path)
{
    string file;
    file.append(path);
    ifstream readFile(file.c_str());
    if (!readFile) {
        cout << path << " does NOT exist!\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::copy_to_cwd(const string& source)
{
    string bash="cp ";
    bash += source+" "+path_cwd;
    system(bash.c_str());
}





void MoltemplateLmpData::make_oplsaalt(const vector<string>& monomerSet)
{
    /** make oplsaa_subset.prm and put it in data file folder **/
    make_oplsaa_subset(monomerSet);
    
    /** invoke oplsaa_moltemplate.py to make oplsaa.lt **/
    string oplsaa_subset=path_cwd+"oplsaa_subset.prm";
    string oplsaa_py=path_moltemplatesrc+"oplsaa_moltemplate.py "+oplsaa_subset;
    string bash="cd "+path_cwd+"; "+oplsaa_py;
    system(bash.c_str());
}





void MoltemplateLmpData::make_polylt(const int polyindex,
                                     const vector<string>& monomerSet)
{
    string output;
    output.append(path_cwd+"/poly_"+to_string((long long int)polyindex+1)+".lt");
    
    ofstream writeFile(output.c_str());
    if (writeFile.is_open())
    {
        /** import oplsaa.lt **/
        writeFile << "import \"oplsaa.lt\"\n";
        
        /** import constituent monomer.lt's **/
        vector<string> used_monomers;
        vector<string>::iterator itr;
        for (indexii=0; indexii<(int)monomerSet.size(); ++indexii) {
            if (indexii>0) {
                itr=find(used_monomers.begin(),used_monomers.end(),monomerSet.at(indexii));
            } else {
                itr=used_monomers.end();
            }
            /** only import monomer if it is new (not repeating) **/
            if (itr==used_monomers.end()) {
                writeFile << "import \""+monomerSet.at(indexii)+"\"\n";
            } used_monomers.push_back(monomerSet.at(indexii));
        } writeFile << "\n";
        
        /** Define combined molecule (ex.polymer) **/
        writeFile << "poly_"<<polyindex+1<<" inherits OPLSAA {\n\n";
        writeFile << "    "<< "create_var {$mol}\n\n";
        
        vector<string> monomerSet_copy=monomerSet;
        double offset_cum=0;
        for (indexii=0; indexii<(int)monomerSet.size(); ++indexii)
        {
            /** erase .lt from name string **/
            monomerSet_copy.at(indexii).erase
            (monomerSet_copy.at(indexii).end()-3,monomerSet_copy.at(indexii).end());
            
            /** pack monomers along x-axis and rotate accordingly (1,0,0) **/
            writeFile
            << "    "
            << "monomer["<<indexii<<"] = new "<<monomerSet_copy.at(indexii);
            if (indexii>0) {
                writeFile
                << ".rot(" <<rotate*(indexii%2)<<",1,0,0)"
                << ".move("<<offset_cum<<",0,0)";
            } writeFile << "\n";
            
            /** evaluate offset distance based on C1-C2 of the pre-mer **/
            evaluate_offset(monomerSet_copy.at(indexii)+".lt");
            offset_cum+=offset;
        }
        
        /** add a list of bonds connecting propagating carbons **/
        writeFile << "\n    write('Data Bond List') {\n";
        for (indexii=0; indexii<(int)monomerSet.size()-1; ++indexii) { //NOTE: "DOP-1" bonds
            writeFile
            << "      "
            << "$bond:b"<<indexii+1<<"  "
            << "$atom:monomer["<<indexii<<"]/C2"<<"  "
            << "$atom:monomer["<<indexii+1<<"]/C1"<<"  "
            << "\n";
        } writeFile << "    }\n";
        
        /** end cap of poly.lt scope **/
        writeFile
        << "\n} # poly_"<<polyindex+1 << "\n";
        
        writeFile.close();
    } else {
        cout << "\nMoltemplateLmpData::make_polyltFile: writeFile can't open.\n";
    }
}





void MoltemplateLmpData::evaluate_offset(const string& merltfile)
{
    string in=path_monomerBank+merltfile;
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        vector<double> C1,C2;
        double dubVar=0;
        string lineContent,strVar;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strVar;
            if (strVar=="write(\"Data")
            {
                /** C1 coordinates **/
                getline(readFile,lineContent);
                iss.clear();iss.str(lineContent);
                for (int i=0;i<4;++i) iss >> strVar;
                for (int i=0;i<3;++i) {
                    iss >> dubVar;
                    C1.push_back(dubVar);
                }
                /** C2 coordinates **/
                getline(readFile,lineContent);
                iss.clear();iss.str(lineContent);
                for (int i=0;i<4;++i) iss >> strVar;
                for (int i=0;i<3;++i) {
                    iss >> dubVar;
                    C2.push_back(dubVar);
                }
                /** calculate C1-C2 distance **/
                double dx2=0,dy2=0,dz2=0,r2=0;
                dx2=pow(C2.at(0)-C1.at(0),2);
                dy2=pow(C2.at(1)-C1.at(1),2);
                dz2=pow(C2.at(2)-C1.at(2),2);
                r2=dx2+dy2+dz2;
                offset=pow(r2,0.5)+offset_spacing;
                //cout << merltfile <<"\n"<<"offset "<<offset<<"\n";
                return;
            }
        }
    } else {
        cout
        << "in MoltemplateLmpData::evaluate_offset():\n"
        << in <<" file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::evaluate_boxLen()
{
    string in=path_cwd+"system.data";
    string lineContent,strVar;
    double dubVar=0;
    double lmin=0,lmax=0;
    ifstream readFile(in.c_str());
    auto start_position=readFile.tellg();
    if (readFile.is_open())
    {
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> strVar;
            if (strVar=="Atoms") {
                getline(readFile,lineContent);//blank
                while (getline(readFile,lineContent))
                {
                    iss.clear();iss.str(lineContent);
                    if (lineContent==""||strVar=="Bonds") break;
                    for (int i=0;i<4;++i) iss >> strVar;
                    for (int i=0;i<3;++i) {
                        iss >> dubVar;
                        if(dubVar<lmin)lmin=dubVar;
                        if(dubVar>lmax)lmax=dubVar;
                    } moltemplateBoxSize=(lmax-lmin)+10.0;
                }
            }
        }
        cout << "Starting Box Size = "<< moltemplateBoxSize << "\n\n";
        system_wait(2);
    } else {
        cout
        << "in MoltemplateLmpData::evaluate_boxLen():\n"
        << in <<" file cannot open.\n";
        exit(EXIT_FAILURE);
    }
    readFile.clear();
    readFile.seekg(start_position);
    string out=path_cwd+"system.tmp";
    ofstream writeFile(out.c_str());
    while (getline(readFile,lineContent)) {
        istringstream iss(lineContent);
        for (int i=0;i<2;++i) iss >> strVar;
        iss >> strVar;
        if (strVar=="xlo"&&lineContent!="") {
            writeFile
            << "     "
            << lmin-5.0 << " "
            << lmax+5.0 << " "
            << "xlo xhi \n";
        } else if (strVar=="ylo"&&lineContent!="") {
            writeFile
            << "     "
            << lmin-5.0 << " "
            << lmax+5.0 << " "
            << "ylo yhi \n";
        } else if (strVar=="zlo"&&lineContent!="") {
            writeFile
            << "     "
            << lmin-5.0 << " "
            << lmax+5.0 << " "
            << "zlo zhi \n";
        } else {
            writeFile << lineContent << "\n";
        }
    } readFile.close(); writeFile.close();
    string mv=
    "cd "+path_cwd+";"+"rm "+in+";"+"mv "+out+" "+in+";"
    "cd "+path_cwd+";";
    system(mv.c_str());
}





void MoltemplateLmpData::make_systemlt()
{
    string output;
    output.append(path_cwd+"/system.lt");
    
    ofstream writeFile(output.c_str());
    if (writeFile.is_open())
    {
        int n_poly=(int)sequenceSet.size();
        if (sequenceLen>1) {
            for (indexi=0; indexi<n_poly; ++indexi) {
                writeFile
                << "import \"poly_"+to_string((long long int)indexi+1)+".lt\"\n";
            } writeFile << "\n";
        } else {
            if (chem_mono.size()>1) {
                cout
                << "sequenceLen = "<<sequenceLen<<", "
                << "merSet should only have one mer type! Please check.\n";
                exit(EXIT_FAILURE);
            }
            /** import constituent monomer.lt's **/
            vector<string> used_monomers;
            vector<string>::iterator itr;
            for (indexi=0;indexi<(int)sequenceSet.size();++indexi) {
                itr=find(used_monomers.begin(),used_monomers.end(),sequenceSet.at(indexi).at(0));
                /** only import monomer if it is new (not repeating) **/
                if (itr==used_monomers.end()) {
                    writeFile << "import \""+sequenceSet.at(indexi).at(0)+"\"\n";
                } used_monomers.push_back(sequenceSet.at(indexi).at(0));
            } writeFile << "\n";
        }
        
        /** Pack molecules in square spiral shape **/
        //----------------------------------------------------------------------
        packingL=offset+packingL_spacing;
        int counter=0;
        int n=0,bndl=0,bndh=0;
        int n_now=0,n_pre=0;
        int signy=1,signz=-1;
        int timey=0,timez=0;
        double valy=0,valz=0;
        double offset_x=0;
        if (sequenceLen>1) {
            offset_x=-50;
            writeFile << "polymer_1 = new poly_1.move("<<offset_x<<",0,0)\n";
        } else {
            //packingL=10;
            offset_x=-5;
            writeFile << "molecule_1 = new "+chem_mono.at(0)+".move("<<offset_x<<",0,0)\n";
        }
        for (indexi=1; indexi<n_poly; ++indexi) {
            n=0;
            while (true) {
                ++n;
                bndl=(n-1)*n;
                bndh=n*(n+1);
                if (bndl<indexi && indexi<=bndh) break;
            } n_now=n;
            if (n_now!=n_pre) {
                counter=0;
                signy*=-1;
                signz*=-1;
            }
            if (counter<n_now) {
                timey=1; valy+=packingL*(double)signy*(double)timey;
                timez=0;
            } else {
                timey=0;
                timez=1; valz+=packingL*(double)signz*(double)timez;
            }
            if (sequenceLen>1) {
                writeFile
                << "polymer_"<<indexi+1<<" = new "
                << "poly_"<<indexi+1<<".move("<<offset_x<<","<<valy<<","<<valz<<")"<< "\n";
            } else {
                writeFile
                << "molecule_"<<indexi+1<<" = new "
                << chem_mono.at(0)+".move("<<offset_x<<","<<valy<<","<<valz<<")"<< "\n";
            } n_pre=n_now; ++counter;
        } writeFile << "\n";
        //----------------------------------------------------------------------
        
        double hbox=moltemplateBoxSize*0.5;
        double fbox=moltemplateBoxSize;
        
        if (true) {
            writeFile
            << "write_once(\"Data Boundary\") {"     << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  xlo xhi" << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  ylo yhi" << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  zlo zhi" << "\n"
            << "}"                                   << "\n\n";
        } else {
            writeFile
            << "write_once(\"Data Boundary\") {"     << "\n"
            << "   0.0  "<<fbox<<"  xlo xhi"         << "\n"
            << "   0.0  "<<fbox<<"  ylo yhi"         << "\n"
            << "   0.0  "<<fbox<<"  zlo zhi"         << "\n"
            << "}"                                   << "\n\n";
        } writeFile.close();
        
    } else {
        cout << "\nMoltemplateLmpData::make_systemltFile: writeFile can't open.\n";
    }
}





void MoltemplateLmpData::FFmodify_alkyl_dihedral_oplsaa()
{
    string in=path_cwd+"oplsaa.lt";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        int    intData=0;
        double dubData=0;
        string lineContent;
        string strData0,strData1;
        vector<int> atomTypes;
        vector<double> old_dihedral;
        vector<double> vecData;
        
        string out=path_cwd+"oplsaa_tmp.lt";
        ofstream writeFile(out.c_str());
        
        /** correct alkyl dihedral coeffs **/
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            if (!is_inside_block) {
                writeFile << lineContent << "\n";
            } iss >> strData0;
            if (strData0=="write_once(\"In") {
                getline(readFile,lineContent);
                iss.clear(); iss.str(lineContent);
                iss >> strData0;
                if (strData0=="dihedral_coeff") {
                    is_inside_block=true;
                } else {
                    writeFile << lineContent << "\n";
                }
            } else if (strData0=="}") {
                if (is_inside_block) {
                    writeFile << lineContent << "\n";
                } is_inside_block=false;
            }
            if (is_inside_block)
            {
                iss >> strData0;
                strData0.erase(strData0.begin(),strData0.begin()+10);//@dihedral:
                for (indexi=0;indexi<strData0.size();++indexi) {
                    if (strData0.at(indexi)!='-') {
                        strData1.push_back(strData0.at(indexi));
                    } else {
                        strData1.push_back(' ');
                    }
                }
                istringstream iss1(strData1);
                strData1.clear(); atomTypes.clear();//NOTE
                while (iss1>>strData0) {
                    intData=atoi(strData0.c_str());
                    atomTypes.push_back(intData);
                }
                int counter=0;
                for (indexi=0;indexi<(int)atomTypes.size();++indexi) {
                    //NOTE: allows repeated, ex. 80-80-80-81
                    if (atomTypes.at(indexi)==80) ++counter;//CH3
                    if (atomTypes.at(indexi)==81) ++counter;//CH2
                    if (atomTypes.at(indexi)==82) ++counter;//CH
                }
                if (counter==(int)atomTypes.size()) {
                    iss >> strData0;//opls
                    old_dihedral.clear();
                    while (iss>>dubData) old_dihedral.push_back(dubData);
                    writeFile << "    " << "dihedral_coeff @dihedral:";
                    for (indexi=0;indexi<(int)atomTypes.size();++indexi) {
                        if (indexi!=(int)atomTypes.size()-1) {
                            writeFile << atomTypes.at(indexi) << "-";
                        } else {
                            writeFile << atomTypes.at(indexi) << " ";
                        }
                    } writeFile << "opls ";
                    /** replace old with new dihedral coeffs **/
                    for (indexi=0;indexi<(int)FFmodify_alkylDihedral.size();++indexi) {
                        writeFile << FFmodify_alkylDihedral.at(indexi) << " ";
                    } writeFile << "\n";
                } else {
                    writeFile << lineContent << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./oplsaa.lt;"
        "mv ./oplsaa_tmp.lt ./oplsaa.lt"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_dihedral_oplsaa():\n"
        << "oplsaa.lt file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_charge_oplsaa()
{
    string in=path_cwd+"oplsaa.lt";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        int    atomType=0;
        //double dubData=0;
        string lineContent;
        string strData0,strData1;
        
        string out=path_cwd+"oplsaa_tmp.lt";
        ofstream writeFile(out.c_str());
        
        /** correct alkyl/hydrogen charges **/
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            if (!is_inside_block) {
                writeFile << lineContent << "\n";
            } iss >> strData0;
            if (strData0=="write_once(\"In") {
                getline(readFile,lineContent);
                iss.clear(); iss.str(lineContent);
                iss >> strData0;
                if (strData0=="set") {
                    is_inside_block=true;
                } else {
                    writeFile << lineContent << "\n";
                }
            } else if (strData0=="}") {
                if (is_inside_block) {
                    writeFile << lineContent << "\n";
                    n_totalatomTypes=atomType;
                    //cout << n_totalatomTypes << "\n";
                } is_inside_block=false;
            }
            if (is_inside_block)
            {
                ++atomType;
                iss >> strData0;//type
                iss >> strData0;
                strData0.erase(strData0.begin(),strData0.begin()+6);//@atom:
                bool   is_alkyl=false;
                double new_charge=0;
                if (atoi(strData0.c_str())==80) //CH3
                {
                    is_alkyl=true;
                    new_charge=FFmodify_CH3_charge;
                    atomType_CH3=atomType;
                    //cout << "atomType_CH3 " << atomType_CH3 << "\n";
                }
                if (atoi(strData0.c_str())==81) //CH2
                {
                    is_alkyl=true;
                    new_charge=FFmodify_CH2_charge;
                    atomType_CH2=atomType;
                    //cout << "atomType_CH2 " << atomType_CH2 << "\n";
                }
                if (atoi(strData0.c_str())==82) //CH
                {
                    is_alkyl=true;
                    new_charge=FFmodify_CH_charge;
                    atomType_CH=atomType;
                    //cout << "atomType_CH " << atomType_CH << "\n";
                }
                if (atoi(strData0.c_str())==85) //H
                {
                    is_alkyl=true;
                    new_charge=FFmodify_H_charge;
                    atomType_H=atomType;
                    //cout << "atomType_H " << atomType_H << "\n";
                }
                if (is_alkyl) {
                    writeFile<<"    "<<"set type @atom:"<<strData0;
                    writeFile<<" charge "<<new_charge<<"\n";
                } else {
                    writeFile << lineContent << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./oplsaa.lt;"
        "mv ./oplsaa_tmp.lt ./oplsaa.lt"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_charge_oplsaa():\n"
        << "oplsaa.lt file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_findHCHbond_oplsaa()
{
    string in=path_cwd+"oplsaa.lt";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        int    intData=0;
        int    bondType=0;
        //double dubData=0;
        string lineContent;
        string strData0,strData1;
        vector<int> atomTypes;
        
        /** find H-CH bond **/
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            if (strData0=="write_once(\"In") {
                getline(readFile,lineContent);
                iss.clear(); iss.str(lineContent);
                iss >> strData0;
                if (strData0=="bond_coeff") is_inside_block=true;
            } else if (strData0=="}") {
                is_inside_block=false;
            }
            if (is_inside_block)
            {
                ++bondType;
                iss >> strData0;
                strData0.erase(strData0.begin(),strData0.begin()+6);//@bond:
                for (indexi=0;indexi<strData0.size();++indexi) {
                    if (strData0.at(indexi)!='-') {
                        strData1.push_back(strData0.at(indexi));
                    } else {
                        strData1.push_back(' ');
                    }
                }
                istringstream iss1(strData1);
                strData1.clear(); atomTypes.clear();//NOTE
                while (iss1>>strData0) {
                    intData=atoi(strData0.c_str());
                    atomTypes.push_back(intData);
                }
                bool is_HCHbond=false;
                //NOTE: not allows repeated, only 82-85 or 85-82 --> H-CH bond
                if(atomTypes.at(0)==82)if(atomTypes.at(1)==85)is_HCHbond=true;
                if(atomTypes.at(0)==85)if(atomTypes.at(1)==82)is_HCHbond=true;
                
                if (is_HCHbond) {
                    is_HCHexisting=true;
                    ++n_HCHbondType;
                    HCHbondType.push_back({n_HCHbondType,bondType});
                }
                if (n_HCHbondType>1) {
                    cout
                    << "in MoltemplateLmpData::FFmodify_alkyl_findHCHbond_oplsaa():\n"
                    << "n_HCHbondType>1, please check "<<in<<"\n";
                    exit(EXIT_FAILURE);
                }
            }
        } readFile.close();
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_findHCHbond_oplsaa():\n"
        << "oplsaa.lt file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_findHCHinBonds_systemdata()
{
    string in=path_cwd+"system.data";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        string lineContent;
        string strData0,strData1;
        
        HCHatomIDs.clear();
        
        /** find H-CH bond in Bonds section of datafile **/
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            if (strData0=="Bonds") {
                is_inside_block=true;
                getline(readFile,lineContent);//blank
                getline(readFile,lineContent);
            } else if (lineContent=="") {//blank
                is_inside_block=false;
            }
            if (is_inside_block)
            {
                iss.clear(); iss.str(lineContent);
                int bondID=0,bondType=0,atom1=0,atom2=0;
                iss >> bondID;
                iss >> bondType;
                iss >> atom1;
                iss >> atom2;
                if (bondType==HCHbondType[0][1]) {
                    HCHatomIDs.push_back(atom1);
                    HCHatomIDs.push_back(atom2);
                }
            }
        } readFile.close();
        
        string cp=
        "cd "+path_cwd+";"
        "cp ./system.data ./system_org.data"; system(cp.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_findHCHinBonds_systemdata():\n"
        << "system.data file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_writeHCHinMasses_systemdata()
{
    string in=path_cwd+"system.data";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        int    counter=0;
        int    atomType=0;
        string lineContent;
        string strData0,strData1,strData2;
        string atominfo,Hatominfo;
        
        string out=path_cwd+"system_tmp.data";
        ofstream writeFile(out.c_str());
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            iss >> strData1;
            iss >> strData2;
            if (!is_inside_block) {
                if (lineContent=="") {//blank
                    writeFile << "\n";
                } else {
                    if (strData1=="atom" && strData2=="types") {
                        writeFile
                        <<"     "<<n_totalatomTypes+1<<"  "<<strData1<<" "<<strData2
                        <<"\n";
                    } else {
                        writeFile << lineContent << "\n";
                    }
                }
            }
            //NOTE: if lineContent is void, strData0 will retain previously stored value
            //cout << "strData0        " << strData0 << "\n";
            //cout << "lineContent     " << lineContent << "\n";
            //cout << "lineContent_pre " << lineContent_pre << "\n\n";
            
            if (strData0=="Masses") {
                is_inside_block=true;
                getline(readFile,lineContent);//blank
                writeFile << "\n";
                getline(readFile,lineContent);
            } else if (lineContent=="") {//blank
                if (is_inside_block) {
                    writeFile << "\n";
                } is_inside_block=false;
            }
            
            /** write H of HCH in Masses section of datafile **/
            if (is_inside_block)
            {
                ++counter;
                iss.clear(); iss.str(lineContent);
                iss >> atomType;
                atominfo.clear();
                if (atomType==atomType_H) {
                    for (indexi=0;indexi<6;++indexi) {
                        iss>>strData0;
                        Hatominfo+=strData0;
                        Hatominfo+=" ";
                    }
                } else {
                    for (indexi=0;indexi<6;++indexi) {
                        iss>>strData0;
                        atominfo+=strData0;
                        atominfo+=" ";
                    }
                }
                /** existing atomTypes **/
                //--------------------------------------------------------------
                if (atomType==atomType_CH3) {
                    writeFile
                    <<"    "<<atomType_CH3<<" "
                    <<atominfo<<"charge="<<FFmodify_CH3_charge<<"\n";
                } else if (atomType==atomType_CH2) {
                    writeFile
                    <<"    "<<atomType_CH2<<" "
                    <<atominfo<<"charge="<<FFmodify_CH2_charge<<"\n";
                } else if (atomType==atomType_CH) {
                    writeFile
                    <<"    "<<atomType_CH<<" "
                    <<atominfo<<"charge="<<FFmodify_CH_charge<<"\n";
                } else if (atomType==atomType_H) {
                    writeFile
                    <<"    "<<atomType_H<<" "
                    <<Hatominfo<<"charge="<<FFmodify_H_charge<<"\n";
                } else {
                    writeFile << lineContent << "\n";
                }
                /** newly added H-CH type **/
                //--------------------------------------------------------------
                if (counter==n_totalatomTypes) {
                    writeFile
                    <<"    "<<n_totalatomTypes+1<<" "
                    <<Hatominfo<<"charge="<<FFmodify_HCH_charge<<"\n";
                }
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./system.data;"
        "mv ./system_tmp.data ./system.data"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_writeHCHinMasses_systemdata():\n"
        << "system.data file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}




void MoltemplateLmpData::FFmodify_alkyl_writeHCHinAtoms_systemdata()
{
    string in=path_cwd+"system.data";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        bool   is_inside_block=false;
        string lineContent;
        string strData0;
        
        string out=path_cwd+"system_tmp.data";
        ofstream writeFile(out.c_str());
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            if (!is_inside_block) {
                if (lineContent=="") {//blank
                    writeFile << "\n";
                } else {
                    writeFile << lineContent << "\n";
                }
            }
            if (strData0=="Atoms") {
                is_inside_block=true;
                getline(readFile,lineContent);//blank
                writeFile << "\n";
                getline(readFile,lineContent);
            } else if (lineContent=="") {//blank
                if (is_inside_block) {
                    writeFile << "\n";
                } is_inside_block=false;
            }
            
            /** find H of HCH in Atoms section of datafile **/
            if (is_inside_block)
            {
                //full: atom-ID molecule-ID atom-type q x y z
                int atomID=0,molID=0,atomType=0; string q,x,y,z;
                iss.clear(); iss.str(lineContent);
                iss >> atomID;
                iss >> molID;
                iss >> atomType;
                iss >> q; iss >> x; iss >> y; iss >> z;
                vector<int>::iterator itr;
                itr=find(HCHatomIDs.begin(),HCHatomIDs.end(),atomID);
                if (itr!=HCHatomIDs.end()) {
                    if (atomType==atomType_H) {
                        writeFile
                        << atomID << " "
                        << molID  << " "
                        << n_totalatomTypes+1 << " "
                        << q << " " << x << " " << y << " " << z << "\n";
                    } else {
                        writeFile << lineContent << "\n";
                    }
                } else {
                    writeFile << lineContent << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./system.data;"
        "mv ./system_tmp.data ./system.data"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_writeHCHinAtoms_systemdata():\n"
        << "system.data file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_writeHCHincharges_systemcharges()
{
    string in=path_cwd+"system.in.charges";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        int    counter=0;
        string lineContent;
        string strData0;
        
        string out=path_cwd+"system_tmp.in.charges";
        ofstream writeFile(out.c_str());
        
        while (getline(readFile,lineContent))
        {
            ++counter;
            
            /** existing atomTypes **/
            //------------------------------------------------------------------
            writeFile << lineContent << "\n";
            
            /** newly added H-CH type **/
            //------------------------------------------------------------------
            if (counter==n_totalatomTypes) {
                writeFile
                << "    set type "<<n_totalatomTypes+1<<" charge "
                << FFmodify_HCH_charge<<"\n";
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./system.in.charges;"
        "mv ./system_tmp.in.charges ./system.in.charges"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_writeHCHincharges_systemcharges():\n"
        << "system.in.charges file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::FFmodify_alkyl_writeHCHinsettings_systemsettings()
{
    string in=path_cwd+"system.in.settings";
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        int    counter=0;
        string lineContent;
        string strData0;
        string pairCoeffs;
        
        string out=path_cwd+"system_tmp.in.settings";
        ofstream writeFile(out.c_str());
        
        while (getline(readFile,lineContent))
        {
            ++counter;
            
            /** existing atomTypes **/
            //------------------------------------------------------------------
            istringstream iss(lineContent);
            if (counter==atomType_H) {
                writeFile << lineContent << "\n";
                iss >> strData0;//pair_coeff
                iss >> strData0;//n
                iss >> strData0;//n
                while (iss>>strData0) {
                    pairCoeffs+=strData0;
                    pairCoeffs+=" ";
                }
            } else {
                writeFile << lineContent << "\n";
            }
            /** newly added H-CH type **/
            //------------------------------------------------------------------
            if (counter==n_totalatomTypes) {
                writeFile
                << "    pair_coeff "<<n_totalatomTypes+1<<" "<<n_totalatomTypes+1<<" "
                << pairCoeffs << "\n";
            }
        } readFile.close(); writeFile.close();
        
        string mv=
        "cd "+path_cwd+";rm ./system.in.settings;"
        "mv ./system_tmp.in.settings ./system.in.settings"; system(mv.c_str());
        
    } else {
        cout
        << "in MoltemplateLmpData::FFmodify_alkyl_writeHCHinsettings_systemsettings():\n"
        << "system.in.settings file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::invoke_moltemplate()
{
    /** NOTE: system.lt is in cwd **/
    string bash =
    "cd "+path_cwd+"; "+path_moltemplatesrc+"moltemplate.sh ./system.lt";
    system(bash.c_str());
}





void MoltemplateLmpData::mv_files()
{
    string datafile,incharge,insetting;
    if ((sequenceLen>1)&&(is_optimizeCharge)) {
        datafile="system_updated.data";
        incharge="system_updated.in.charges";
        insetting="system_updated.in.settings";
    } else {
        datafile="system.data";
        incharge="system.in.charges";
        insetting="system.in.settings";
    }
    for (indexii=0; indexii<n_trials; ++indexii) {
        string data=
        "cd "+path_cwd+";"
        "cp "+datafile+" ../; cd ../;"
        "mv "+datafile+" "
        "input_00"+to_string((long long int)indexii)+".data";
        system(data.c_str());
    }
    string init="cd "+path_cwd+";"
    "cp "+incharge+" "+insetting+" system.in system.in.init ../../generation/;"
    "cd ../../generation/;";
    if ((sequenceLen>1)&&(is_optimizeCharge)) {
        init+=
        "mv system_updated.in.charges  system.in.charges;"
        "mv system_updated.in.settings system.in.settings";
    } system(init.c_str());
    
    string output="cd "+path_cwd+";";
    if ((sequenceLen>1)&&(is_optimizeCharge)) {
        output+="mkdir output_updated; mv *_updated.* output_updated/;";
    } output+="mkdir output; mv system.in* system*data output_ttree output/";
    system(output.c_str());
    
    string input="cd "+path_cwd+"; mkdir input; mv *.lt *.prm *.txt input/";
    system(input.c_str());
    
    //string scriptsource=path_master+"scripts";
    string simfolder=path_cwd+"../../../../simulations/";
    string anafolder=path_cwd+"../../../../analysis/";
    string cpfiles=
    "cp "+path_scripts+"generation.inp "   +simfolder+"lammps_inputs/generation/;"
    "cp "+path_scripts+"quench.inp "       +simfolder+"lammps_inputs/quench/;"
    "cp "+path_scripts+"equilibration.inp "+simfolder+"lammps_inputs/equilibration/;"
    "cp "+path_scripts+"production.inp "   +simfolder+"lammps_inputs/production/;"
    "cp "+path_scripts+"amdat.inp "        +anafolder+"AMDAT_inputs/;";
    system(cpfiles.c_str());
}





int MoltemplateLmpData::show_positinindex(const int a, const std::vector<int>& va)
{
    int positionindex=0;
    for (indexi=0; indexi<(int)va.size(); ++indexi) {
        if (va.at(indexi)==a) {
            positionindex=indexi;
            break;
        }
    } return positionindex;
}





void MoltemplateLmpData::initialize_n_typesVec()
{
    n_types_all.clear();
    n_types_heavy.clear();
    n_types_light.clear();
    
    vector<int> tmpvi;
    for (indexi=0; indexi<(int)types_all.size(); ++indexi) {
        tmpvi.clear();
        tmpvi.push_back(types_all.at(indexi)); // type_all
        tmpvi.push_back(0); // intialize with 0
        n_types_all.push_back(tmpvi);
        n_types_heavy.push_back(tmpvi);
        n_types_light.push_back(tmpvi);
    }
}
void MoltemplateLmpData::finalize_n_typesVec(const int n_poly)
{
    int n_all=0,n_heavy=0,n_light=0;
    for (indexi=0; indexi<(int)types_all.size(); ++indexi) {
        n_all += n_types_all.at(indexi).at(1);
        n_heavy += n_types_heavy.at(indexi).at(1);
        n_light += n_types_light.at(indexi).at(1);
    }
    //cout << "n_poly  " << n_poly  << "\n";
    //cout << "n_all   " << n_all   << "\n";
    //cout << "n_heavy " << n_heavy << "\n";
    //cout << "n_light " << n_light << "\n";
    
    /** to get atom_types and their respective numbers in each molecule **/
    if(false)
    {
        for (indexi=0; indexi<(int)types_all.size(); ++indexi) {
            n_types_all.at(indexi).at(1) /= n_poly;
            n_types_heavy.at(indexi).at(1) /= n_poly;
            n_types_light.at(indexi).at(1) /= n_poly;
        }
        /* check the number to see if it's consistent (should be integer values,
         meaning they are completely divided) */
        n_all=0,n_heavy=0,n_light=0;
        for (indexi=0; indexi<(int)types_all.size(); ++indexi) {
            n_all += n_types_all.at(indexi).at(1);
            n_heavy += n_types_heavy.at(indexi).at(1);
            n_light += n_types_light.at(indexi).at(1);
        }
    }
}





void MoltemplateLmpData::add_n_types_all(const int atomType)
{
    vector<int>::iterator itr;
    itr=find(types_all.begin(),types_all.end(),atomType);
    if(itr!=types_all.end()) {
        indexi=show_positinindex(atomType,types_all);
        n_types_all.at(indexi).at(1) += 1;
    }
}
void MoltemplateLmpData::add_n_types_heavy(const int atomType)
{
    vector<int>::iterator itr;
    itr=find(types_all.begin(),types_all.end(),atomType);
    if(itr!=types_all.end()) {
        indexi=show_positinindex(atomType,types_all);
        n_types_heavy.at(indexi).at(1) += 1;
    }
}
void MoltemplateLmpData::add_n_types_light(const int atomType)
{
    vector<int>::iterator itr;
    itr=find(types_all.begin(),types_all.end(),atomType);
    if(itr!=types_all.end()) {
        indexi=show_positinindex(atomType,types_all);
        n_types_light.at(indexi).at(1) += 1;
    }
}





void MoltemplateLmpData::set_atomTypes(StructureClass& sysVar)
{
    /** threshold mass for light atoms (a.m.u.) **/
    double light=3.0;
    
    types_light.clear(); types_heavy.clear();
    string data=path_cwd;
    string datafile;
    if (sequenceLen>1) datafile="output_updated/system_updated.data";
    else datafile="output/system.data";
    
    if (is_mass_normalization) {
        string test=path_cwd+"output/system_massOrg.data";
        ifstream readFile(test.c_str());
        if(readFile.is_open()) {
            data+="output/system_massOrg.data";
        } else {
            data+="output/system.data";
        }
    } else data+=datafile;
    
    ifstream readFile(data.c_str());
    if(readFile.is_open())
    {
        string lineContent;
        string strvar;
        double dubvar=0;
        bool is_inblock=false;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> strvar;
            if (strvar=="Masses") {
                is_inblock=true;
                continue;
            } else if (strvar=="Atoms") {
                break;
            }
            if(is_inblock)
            {
                iss >> dubvar;
                types_all.push_back(atoi(strvar.c_str()));
                if (dubvar<=light) {
                    types_light.push_back(atoi(strvar.c_str()));
                } else {
                    types_heavy.push_back(atoi(strvar.c_str()));
                }
            }
        } readFile.close();
        check_lasttwo(types_all); sysVar.set_types_all(types_all);
        check_lasttwo(types_light); sysVar.set_types_light(types_light);
        check_lasttwo(types_heavy); sysVar.set_types_heavy(types_heavy);
    } else {
        if (is_path_exist(path_cwd)) {
            cout
            << "in MoltemplateLmpData::set_atomTypes():\n"
            << data <<" can NOT open. Please check.\n";            
            exit(EXIT_FAILURE);
        }
    }
}





void MoltemplateLmpData::set_atomNumbers(StructureClass& sysVar)
{
    n_light=0; n_heavy=0;
    int chainLen=sysVar.get_chainLen();
    
    string data=path_cwd;
    string datafile;
    if (sequenceLen>1) datafile="output_updated/system_updated.data";
    else datafile="output/system.data";
    
    if (is_mass_normalization) {
        string test=path_cwd+"output/system_massOrg.data";
        ifstream readFile(test.c_str());
        if(readFile.is_open()) {
            data+="output/system_massOrg.data";
        } else {
            data+="output/system.data";
        }
    } else data+=datafile;
    
    /** create the data structure used to store atom_types and thier respective
     ** numbers (in each molecule) **/
    initialize_n_typesVec();
    
    ifstream readFile(data.c_str());
    if(readFile.is_open())
    {
        string lineContent;
        string strvar,str_now,str_pre;
        int intvar=0,n_poly=0;
        bool is_inblock=false;
        while (getline(readFile,lineContent)) {
            istringstream iss(lineContent);
            iss >> strvar; str_now=strvar;
            if (strvar=="Atoms") {
                is_inblock=true;
                continue;
            } else if (strvar=="Bonds") {
                break;
            }
            if(is_inblock)
            {
                if (str_now==str_pre) break;
                iss >> intvar; // molID
                n_poly = intvar;
                iss >> intvar; // atomType
                add_n_types_all(intvar);
                if (check_is_heavyAtom(intvar)) {
                    ++n_heavy;
                    add_n_types_heavy(intvar);
                } else {
                    ++n_light;
                    add_n_types_light(intvar);
                }
            } str_pre=str_now;
        } readFile.close();
        
        if (sysVar.get_is_moltempLmpData()) {
            if (n_poly!=(int)sequenceSet.size()) {
                cout
                << "in MoltemplateLmpData::set_atomNumbers(), "
                << "molid != n_poly" << "\n";
                exit(EXIT_FAILURE);
            }
        }
        //cout << "n_total " << n_heavy+n_light << "\n";
        //cout << "n_heavy " << n_heavy << "\n";
        //cout << "n_light " << n_light << "\n";
        
        /** set number of polymers (i.e. molecules) */
        sysVar.set_n_poly(n_poly);
        
        /** set number of heavy and light atoms **/
        n_heavy /= n_poly;
        n_light /= n_poly;
        sysVar.set_n_light(n_light);
        sysVar.set_n_heavy(n_heavy);
        
        /** According to the design of monomer structures in the monomer bank,
         ** right-end and left-end monomers have one more H than intermediate
         ** monomers; n_typeSet stores the number of heavy and light atoms in
         ** polymers constructed from using monomers in the monomer bank for
         ** later use in wrapper function for AMDAT for multibody analysis **/
        sysVar.set_n_typeSet
        ({
            {n_heavy/chainLen,((n_light-2)/chainLen)+1}, // le monomer
            {n_heavy/chainLen,(n_light-2)/chainLen},     // im monomers
            {n_heavy/chainLen,((n_light-2)/chainLen)+1}  // re monomer
        });
        
        /** n_typeVecs have the struture that shows the number of atoms for each
         ** corresponding type of atom;
         ** i.e. per 1-D element: (type,number of type per molecule) **/
        finalize_n_typesVec(n_poly);
        sysVar.set_n_types_all(n_types_all);
        sysVar.set_n_types_heavy(n_types_heavy);
        sysVar.set_n_types_light(n_types_light);
        
    } else {
        if (is_path_exist(path_cwd)) {
            cout
            << "in MoltemplateLmpData::set_atomNumbers():\n"
            << data <<" can NOT open. Please check.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void MoltemplateLmpData::check_lasttwo(std::vector<int>& vi)
{
    size_t svi=vi.size();
    if (vi.at(svi-1)==vi.at(svi-2)) vi.pop_back();
}
void MoltemplateLmpData::check_lasttwo(std::vector<double>& vd)
{
    size_t svd=vd.size();
    if (vd.at(svd-1)==vd.at(svd-2)) vd.pop_back();
}
void MoltemplateLmpData::check_lasttwo(std::vector<std::vector<int>>& vvi)
{
    size_t svvi=vvi.size();
    if (vvi.at(svvi-1).at(0)==vvi.at(svvi-2).at(0)) vvi.pop_back();
}
void MoltemplateLmpData::check_lasttwo(std ::vector<std::vector<double>>& vvd)
{
    size_t svvd=vvd.size();
    /* pop back if the last two types are repeated */
    if (vvd.at(svvd-1).at(0)==vvd.at(svvd-2).at(0)) vvd.pop_back();
}




void MoltemplateLmpData::make_oplsaa_subset(const vector<string>& monomerSet)
{
    /** Primary contributor: Venkatesh Meenakshisundaram
     ** Integration into master code: Jui-Hsiang Hung (SJH) **/
    
    /** Monomer set to be used for creating oplsaa_subset **/
    const vector<string> mono_inputs = monomerSet;
    //const vector<string> mono_inputs = {"M001i.lt","M025i.lt"}; // test line
    
    /** path to oplsaa_subset.prm file **/
    const string opls_subset_file = path_cwd+"oplsaa_subset.prm";
    
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////
    /** vector to store all atom types including the repeats **/
    vector<string> atom_keys;
    for(int vecii=0; vecii<mono_inputs.size(); vecii++)
    {
        /** path to monomer.lt in monomer bank **/
        const string mono = path_monomerBank+mono_inputs[vecii];
        
        bool read_switch= false;
        ifstream read_mono(mono.c_str());
        if(read_mono.is_open())
        {
            string mono_line;
            while(getline(read_mono,mono_line))
            {
                if(mono_line.size()!=0)
                {
                    /** Modified: SJH **/
                    string key;
                    istringstream iss(mono_line);
                    iss >> key;
                    //if(mono_line == "	write(\"Data Atoms\") {")
                    if(key=="write(\"Data")
                    {
                        read_switch = true;
                        continue;
                    }
                    else if(key=="}")
                    {
                        read_switch = false;
                        break;
                    }
                    /** Determine atom types, element names and the raw_charges as
                     given in the opls table **/
                    if(read_switch)
                    {
                        vector<string> stringvector;
                        string stringelements;
                        istringstream checkstring (mono_line);
                        while (checkstring >> stringelements)
                        {
                            stringvector.push_back(stringelements);
                            /** storing every entity of the line (checkstring) in
                             a vector; The vector is overwritten for every line **/
                        }
                        string load_line;
                        load_line.clear();
                        bool load_switch=false;
                        for(int readii=0; readii<stringvector.at(2).size(); readii++)
                        {
                            /** Modified: SJH **/
                            if(stringvector.at(2).at(readii)==':')
                            {
                                load_switch = true;
                                continue;
                            }
                            if(load_switch)
                            {
                                load_line += stringvector.at(2).at(readii);
                            }
                        }
                        atom_keys.push_back(load_line);
                    }
                }
            }
            read_mono.close();
        }
        else
        {
            cout
            << "Monomer ("<< mono_inputs[vecii] << ") does NOT exist. \n"
            << "Please check the following path to the file\n" << mono << "\n";
            exit(EXIT_FAILURE);
        }
    }
    /** Cleaning up the stored data. Remove duplicate atoms types **/
    for(int ii=0;ii<atom_keys.size();ii++)
    {
        for(int iii=ii+1;iii<atom_keys.size();iii++)
        {
            if(atom_keys[ii]==atom_keys[iii])
            {
                atom_keys.erase(atom_keys.begin()+iii);
                iii--;
            }
        }
    }
    /** Convert the vectors string to vector int in order to sort the atom_types
     in ascending order **/
    vector<int>atom_types;
    for(int ii=0; ii<atom_keys.size();ii++)
    {
        atom_types.push_back(atoi(atom_keys[ii].c_str()));
    }
    sort(atom_types.begin(),atom_types.end());
    
    /** Read the master opls file and store the ones that match the atom_types
     into new subset file **/
    ofstream write_subset_file;
    write_subset_file.open(opls_subset_file.c_str());
    
    ifstream read_prm(path_oplsaaprm.c_str());
    if(read_prm.is_open())
    {
        string prm_line;
        bool check_switch=false;
        while(getline(read_prm,prm_line))
        {
            if(prm_line.size()!=0)
            {
                if(prm_line == "      ##  Atom Type Definitions  ##")
                {
                    check_switch = true;
                    write_subset_file << prm_line << "\n";
                    /** Modified: SJH **/
                    getline(read_prm,prm_line);
                    write_subset_file << prm_line << "\n";
                    getline(read_prm,prm_line);
                    write_subset_file << prm_line << "\n";
                    continue;
                }
                /** Modified: SJH **/
                else if (prm_line=="      ################################")
                {
                    check_switch = false;
                    write_subset_file << prm_line << "\n";
                    continue;
                }
                else if(check_switch)
                {
                    vector<string> stringvector;
                    string stringelements;
                    istringstream checkstring(prm_line);
                    while (checkstring >> stringelements)
                    {
                        //cout << stringelements << " ";
                        stringvector.push_back(stringelements);
                        /** storing every entity of the line (checkstring) in a
                         vector; The vector is overwritten for every line **/
                    } //cout << "\n";
                    for(int checkii=0; checkii<atom_types.size(); checkii++)
                    {
                        if(atom_types.at(checkii)==atoi(stringvector.at(1).c_str()))
                        {
                            write_subset_file << prm_line << "\n";
                            break;
                        }
                        /** Modified: SJH **/
                        //else if (stringvector[0][0]=='#')
                        //{
                        //    write_subset_file << prm_line << "\n";
                        //    break;
                        //}
                    }
                }
                else
                {
                    write_subset_file << prm_line << "\n";
                }
            }
            else
            {
                write_subset_file << prm_line << "\n";
            }
        }
        read_prm.close();
    }
    else
    {
        cout
        << "The oplsaa.prm source file don't exist. \n"
        << "Please check the following path to the file\n"
        << path_oplsaaprm << "\n";
        exit(EXIT_FAILURE);
    }
    write_subset_file.close();
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
}





void MoltemplateLmpData::getridof_full()
{
    string in=path_cwd+"system.data";
    string out=path_cwd+"tmp.data";
    ofstream writeFile(out.c_str());
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        string lineContent,str;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> str;
            if (lineContent=="") {
                writeFile << "\n";
            } else if (str=="Atoms") {
                writeFile << "Atoms\n";
            } else {
                writeFile << lineContent << "\n";
            }
        } readFile.close(); writeFile.close();
    } else {
        cout
        << "in MoltemplateLmpData::getridof_full():\n"
        << in+" cannot open. Please check.\n";
        exit(EXIT_FAILURE);
    }
    string mv="rm "+in+";mv "+out+" "+in;
    system(mv.c_str());
}





void MoltemplateLmpData::getridof_ljcutcoullong()
{
    string in=path_cwd+"system.in.settings";
    string out=path_cwd+"tmp.data";
    ofstream writeFile(out.c_str());
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        string lineContent,str;
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> str;
            if (lineContent=="") {
                writeFile << "\n";
            } else if (str=="pair_coeff") {
                writeFile << "    pair_coeff ";
                while (iss>>str) {
                    if (str=="lj/cut/coul/long") {
                        continue;
                    } else {
                        writeFile << str << " ";
                    }
                } writeFile << "\n";
            } else {
                writeFile << lineContent << "\n";
            }
        } readFile.close(); writeFile.close();
    } else {
        cout
        << "in MoltemplateLmpData::getridof_ljcutcoullong():\n"
        << in+" cannot open. Please check.\n";
        exit(EXIT_FAILURE);
    }
    string mv="rm "+in+";mv "+out+" "+in;
    system(mv.c_str());
}





void MoltemplateLmpData::mass_normalization()
{
    string in=path_cwd;
    if ((sequenceLen>1)&&(is_optimizeCharge)) in+="system_updated.data";
    else in+="system.data";
    
    vector<int> typeVec;
    vector<double> massVec;
    vector<vector<int>> typeNum;
    double mass_normalized=0;
    
    /** -- read Masses section for atomTypes and masses -- **/
    
    ifstream readFile(in.c_str());
    auto start_position=readFile.tellg();
    if (readFile.is_open())
    {
        string lineContent;
        string strData0;
        bool   is_inblock=false;
        int    type=0;
        double mass=0;
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strData0;
            if (strData0=="Masses") {
                is_inblock=true;
                getline(readFile,lineContent);//blank
                getline(readFile,lineContent);
            } else if (lineContent=="") {//blank
                if (is_inblock) break;
            }
            if(is_inblock)
            {
                iss.clear(); iss.str(lineContent);
                iss >> type;
                iss >> mass;
                typeVec.push_back(type);
                typeNum.push_back({type,0});
                massVec.push_back(mass);
            }
        }
        if (typeVec.size()!=massVec.size()) {
            cout <<
            "in MoltemplateLmpData::mass_normalization(), "
            "typeVec.size()!=massVec.size() \n"
            "Please check. \n";
            exit(EXIT_FAILURE);
        }
        //for (indexi=0; indexi<typeVec.size();++indexi) {
        //    cout
        //    << typeVec.at(indexi) << " "
        //    << massVec.at(indexi) << "\n";
        //}
    } else {
        cout
        << "in MoltemplateLmpData::mass_normalization():\n"
        << in+" cannot open. Please check.\n";
        exit(EXIT_FAILURE);
    }
    
    string lineContent;
    string strData0;
    int    atomID=0,molID=0,type=0;
    bool   is_inblock=false;
    double totalmass=0;
    
    /** -- read Atoms section to add up total mass -- **/
    
    readFile.clear();
    readFile.seekg(start_position);
    while (getline(readFile,lineContent))
    {
        istringstream iss(lineContent);
        iss >> strData0;
        if (strData0=="Atoms") {
            is_inblock=true;
            getline(readFile,lineContent);//blank
            getline(readFile,lineContent);
        } else if (lineContent=="") {//blank
            if (is_inblock) break;
        }
        if(is_inblock)
        {
            iss.clear(); iss.str(lineContent);
            iss >> atomID;
            iss >> molID;
            iss >> type;
            typeNum.at(type-1).at(1) += 1;
            totalmass += massVec.at(show_positinindex(type,typeVec));
        }
    } mass_normalized=totalmass/(double)atomID;
    
    if (molID!=(int)sequenceSet.size()) {
        cout << "molID  " << molID << "\n";
        cout << "n_poly " << sequenceSet.size() << "\n";
        cout
        << "in MoltemplateLmpData::mass_normalization(),\n"
        << "molid != n_poly" << "\n";
        exit(EXIT_FAILURE);
    }
    
    double mass_verify=0;
    for (indexi=0; indexi<(int)typeNum.size(); ++indexi) {
        cout
        << "atomType"<<typeNum.at(indexi).at(0)<<" "
        << massVec.at(indexi)<<" "
        << typeNum.at(indexi).at(1)<<"\n";
        mass_verify += massVec.at(indexi)*typeNum.at(indexi).at(1);
    }
    cout << "totalmass        " << totalmass       << "\n";
    cout << "totalmass_verify " << mass_verify     << "\n";
    cout << "n_total          " << atomID          << "\n";
    cout << "mass_normalized  " << mass_normalized << "\n";
    
    /** -- replace original masses by normalized mass -- **/
    
    string out=path_cwd+"system_tmp.data";
    
    ofstream writeFile(out.c_str());
    readFile.clear();
    readFile.seekg(start_position);
    is_inblock=false;
    strData0.clear();
    int    atomType=0;
    double atomMass=0;
    while (getline(readFile,lineContent))
    {
        istringstream iss(lineContent);
        iss >> strData0;
        if (!is_inblock) {
            if (lineContent=="") {//blank
                writeFile << "\n";
            } else {
                writeFile << lineContent << "\n";
            }
        }
        if (strData0=="Masses") {
            is_inblock=true;
            getline(readFile,lineContent);//blank
            writeFile << "\n";
            getline(readFile,lineContent);
        } else if (lineContent=="") {//blank
            if (is_inblock) {
                writeFile << "\n";
            } is_inblock=false;
        }
        if(is_inblock)
        {
            iss.clear(); iss.str(lineContent);
            iss >> atomType;
            iss >> atomMass;
            writeFile << "    ";
            writeFile << atomType << " ";
            writeFile << mass_normalized << " ";//replace mass by normalized value
            while (iss>>strData0) {
                writeFile << strData0 << " ";
            } writeFile << "\n";
        }
    } readFile.close(); writeFile.close();
    
    /** store a copy of original mass & mv new mass to system_updated.data **/
    string mv=
    "mv "+in+" "+path_cwd+"system_massOrg.data;"
    "mv "+out+" "+in; system(mv.c_str());
}





void MoltemplateLmpData::optimize_charge()
{
    /** Primary contributor: Venkatesh Meenakshisundaram
     ** Integration into master code: Jui-Hsiang Hung (SJH) **/
    
    /** Modified: SJH **/
    const string lammps_file   = path_cwd+"system.data";
    const string pair_file     = path_cwd+"system.in.settings";
    const string o_lammps_file = path_cwd+"system_updated.data";
    const string o_charge_file = path_cwd+"system_updated.in.charges";
    const string o_pair_file   = path_cwd+"system_updated.in.settings";
    
    //--------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////
    /**READING LAMMPS data file to get all associated data of the system**/
    //Identify the types of atoms in the system and the element they are associated with
    vector<int> atom_types;		//number of atom types
    //vector<int> num_atoms;	//number of atoms in each type
    vector<string> atom_names;	//The elements associated with each atom type
    vector<double> type_charges;//The charges associated with each atom type
    bool mass_switch = false;
    
    ifstream lmp_mass_obj(lammps_file.c_str());
    if(lmp_mass_obj.is_open())
    {
        string lmp_line;
        while(getline(lmp_mass_obj,lmp_line))
        {
            if(lmp_line.size()!=0)
            {
                if(lmp_line == "Masses")
                {
                    mass_switch = true;
                    continue;
                }
                else if((lmp_line == "Atoms") || (lmp_line == "Atoms # full"))
                {
                    mass_switch = false;
                    break;
                }
                /*Determine atom types, element names and the raw_charges as
                 given in the opls table*/
                if(mass_switch)
                {
                    vector<string> stringvector;
                    string stringelements;
                    istringstream checkstring (lmp_line);
                    while (checkstring >> stringelements)
                    {
                        stringvector.push_back(stringelements);
                        /*storing every entity of the line (checkstring) in a
                         vector; The vector is overwritten for every line*/
                    }
                    //Load the atom types
                    atom_types.push_back(atoi(stringvector[0].c_str()));
                    
                    //Load the element name associated with each atom type
                    char ele_char[stringvector[2].size()];
                    strcpy(ele_char,stringvector[2].c_str());
                    string temp_ele;
                    if(islower(ele_char[2])) {
                        temp_ele = ele_char[1] + ele_char[2];
                    } else {
                        temp_ele = ele_char[1];
                    }
                    atom_names.push_back(temp_ele);
                    
                    //Load the charges associated with each type of atoms
                    char chr_char[stringvector[stringvector.size()-1].size()];
                    strcpy(chr_char,stringvector[stringvector.size()-1].c_str());
                    bool key = false;
                    string temp_chr;
                    temp_chr.clear();
                    for(int ii=0; ii<stringvector[stringvector.size()-1].size();ii++)
                    {
                        if(chr_char[ii]=='=') {
                            key = true;
                            continue;
                        }
                        if(key==true) {
                            temp_chr+=chr_char[ii];
                        }
                    }
                    type_charges.push_back(atof(temp_chr.c_str()));
                    //num_atoms.resize(atom_types.size());
                }
            }
        }
        lmp_mass_obj.close();
    }
    else
    {
        /** Modified: SJH **/
        cout
        << "The lammps data file don't exist. Please check the following path"
        " to the file\n" << lammps_file << "\n";
        exit(EXIT_FAILURE);
    }
    
    //Identify all the atoms in the system and its atom type
    vector<int> all_atom_id;
    vector<int> all_atom_type;
    bool atom_switch = false;
    
    ifstream lmp_atom_obj(lammps_file.c_str());
    if(lmp_atom_obj.is_open())
    {
        string lmp_line;
        while(getline(lmp_atom_obj,lmp_line))
        {
            if(lmp_line.size()!=0)
            {
                if((lmp_line == "Atoms") || (lmp_line == "Atoms # full")) {
                    atom_switch = true;
                    continue;
                } else if(lmp_line == "Bonds") {
                    atom_switch = false;
                    break;
                }
                /*Determine atom types, element names and the raw_charges as
                 given in the opls table*/
                if(atom_switch)
                {
                    vector<string> stringvector;
                    string stringelements;
                    istringstream checkstring (lmp_line);
                    while (checkstring >> stringelements)
                    {
                        stringvector.push_back(stringelements);
                        /*storing every entity of the line (checkstring) in a
                         vector; The vector is overwritten for every line*/
                    }
                    //Load the atom id
                    all_atom_id.push_back(atoi(stringvector[0].c_str()));
                    //Load the atom type associated with each atom ID
                    all_atom_type.push_back(atoi(stringvector[2].c_str()));
                }
            }
        }
        lmp_atom_obj.close();
    }
    else
    {
        /** Modified: SJH **/
        cout << "The lammps data file don't exist. Please check the following "
        "path to the file\n" << lammps_file << "\n";
        exit(EXIT_FAILURE);
    }
    
    /*Collect all the bond infomation from lammps file
     (Atom connections-> Reconstruct the molecules)*/
    vector<int> central_atom;
    vector<vector<int>> connect_atom;
    bool bond_switch = false;
    ifstream lmp_bond_obj(lammps_file.c_str());
    if(lmp_bond_obj.is_open())
    {
        string lmp_line;
        vector<int> temp_connect;
        while(getline(lmp_bond_obj,lmp_line))
        {
            if(lmp_line.size()!=0)
            {
                if(lmp_line == "Bonds") {
                    bond_switch = true;
                    continue;
                } else if(lmp_line == "Angles") {
                    bond_switch = false;
                    break;
                }
                /*Determine atom types, element names and the raw_charges as
                 given in the opls table*/
                if(bond_switch)
                {
                    vector<string> stringvector;
                    string stringelements;
                    istringstream checkstring (lmp_line);
                    while (checkstring >> stringelements)
                    {
                        stringvector.push_back(stringelements);
                        /*storing every entity of the line (checkstring) in a
                         vector; The vector is overwritten for every line*/
                    }
                    
                    if(central_atom.size()!=0)
                    {
                        int temp_center = atoi(stringvector[2].c_str());
                        if(central_atom[central_atom.size()-1] == temp_center) {
                            temp_connect.push_back(atoi(stringvector[3].c_str()));
                        } else {
                            connect_atom.push_back(temp_connect);
                            temp_connect.clear();
                            central_atom.push_back(atoi(stringvector[2].c_str()));
                            temp_connect.push_back(atoi(stringvector[3].c_str()));
                        }
                    } else {
                        central_atom.push_back(atoi(stringvector[2].c_str()));
                        temp_connect.push_back(atoi(stringvector[3].c_str()));
                    }
                }
            }
        }
        connect_atom.push_back(temp_connect);
        lmp_bond_obj.close();
    }
    else
    {
        /** Modified: SJH **/
        cout << "The lammps data file don't exist. Please check the following "
        "path to the file\n" << lammps_file << "\n";
        exit(EXIT_FAILURE);
    }
    //--------------------------------------------------------------------------
    
    /**Cleaning up the stored data**/
    //Combine all atom ID together
    for(int centerii=0; centerii<central_atom.size(); centerii++)
    {
        for(int bii=centerii+1; bii<central_atom.size();bii++)
        {
            if(central_atom[centerii]==central_atom[bii])
            {
                for(int inii=0; inii<connect_atom[bii].size(); inii++)
                {
                    connect_atom[centerii].push_back(connect_atom[bii][inii]);
                }
                connect_atom.erase(connect_atom.begin()+bii);
                central_atom.erase(central_atom.begin()+bii);
                bii--;
            }
        }
    }
    /*If an atom is shared by two carbons, remove the atom from the vector of
     latter carbon and keep it only the former carbon. NOTE: This ensures that
     the charges are not double counted*/
    for(int centerii=0; centerii<central_atom.size(); centerii++)
    {
        for(int bii=0; bii<connect_atom[centerii].size();bii++)
        {
            for(int connectii=centerii+1; connectii<central_atom.size();connectii++)
            {
                for(int cii=0; cii<connect_atom[connectii].size();cii++)
                {
                    if(connect_atom[centerii][bii]==connect_atom[connectii][cii])
                    {
                        connect_atom[connectii].erase(connect_atom[connectii].begin()+cii);
                        cii--;
                    }
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    
    /**If an Oxygen or a Nitrogen is the central atom, club the atoms to the
     connected oxygen or nitrogen as one single super atom**/
    for(int centerii=0; centerii<central_atom.size(); centerii++)
    {
        //Determine the type of the central atom
        int temp_type=0;
        for(int aaii=0; aaii<all_atom_id.size(); aaii++)
        {
            if(all_atom_id[aaii]==central_atom[centerii])
            {
                temp_type = all_atom_type[aaii];
                break;
            }
        }
        string temp_atom;
        for(int atii=0; atii<atom_types.size(); atii++)
        {
            if(atom_types[atii]==temp_type)
            {
                temp_atom = atom_names[atii];
                break;
            }
        }
        if(temp_atom == "O")
        {
            for(int cenii=0;cenii<central_atom.size();cenii++)
            {
                for(int conii=0;conii<connect_atom[cenii].size();conii++)
                {
                    if(central_atom[centerii]==connect_atom[cenii][conii])
                    {
                        cout << "An oxygen atom was found to be a central atom."
                        " This atom and other atoms connected to this oxygen "
                        "will be added to Carbon atom to which this oxygen is connected\n";
                        for(int transii=0;transii<connect_atom[centerii].size();transii++)
                        {
                            connect_atom[cenii].push_back(connect_atom[centerii][transii]);
                        }
                        central_atom.erase(central_atom.begin()+centerii);
                        connect_atom.erase(connect_atom.begin()+centerii);
                        centerii--;
                    }
                }
            }
        }
        else if(temp_atom == "N")
        {
            for(int cenii=0;cenii<central_atom.size();cenii++)
            {
                for(int conii=0;conii<connect_atom[cenii].size();conii++)
                {
                    if(central_atom[centerii]==connect_atom[cenii][conii])
                    {
                        cout << "A nitrogen atom was found to be a central atom. "
                        "This atom and other atoms connected to this nitrogen "
                        "will be added to Carbon atom to which this nitrogen is connected\n";
                        for(int transii=0;transii<connect_atom[centerii].size();transii++)
                        {
                            connect_atom[cenii].push_back(connect_atom[centerii][transii]);
                        }
                        central_atom.erase(central_atom.begin()+centerii);
                        connect_atom.erase(connect_atom.begin()+centerii);
                        centerii--;
                    }
                }
            }
        }
        else if(temp_atom!="C")
        {
            /** Modified: SJH **/
            cout << "This program is designed to work on molecules that have "
            "Carbon, Oxygen, and Nitrogen as center atoms. Please check the "
            "following atom_ID in the lammps input file:\n";
            cout << "Atom ID: " << central_atom[centerii]
            << "\tAtom type name: " << temp_atom << "\n";
            exit(EXIT_FAILURE);
        }
    }
    //--------------------------------------------------------------------------
    
    /**Reconstructing the molecules to identify the charge build up in the
     system and adjust the charges to make the system neutral**/
    vector<int> new_atom_types;
    vector<double> new_atom_charges;
    vector<int> duplication_atom_t_info;
    vector<double> duplication_atom_charge;
    
    vector<int> charged_center_atoms;
    vector<int> num_connect_atoms;
    vector<double> net_charge_store;
    
    for(int centerii=0; centerii<central_atom.size(); centerii++)
    {
        //Determine the type of the central atom
        int temp_type=0;
        for(int aaii=0; aaii<all_atom_id.size(); aaii++)
        {
            if(all_atom_id[aaii]==central_atom[centerii])
            {
                temp_type = all_atom_type[aaii];
                break;
            }
        }
        string temp_atom;
        double center_atom_charge =0;
        for(int atii=0; atii<atom_types.size(); atii++)
        {
            if(atom_types[atii]==temp_type)
            {
                temp_atom = atom_names[atii];
                center_atom_charge = type_charges[atii];
                break;
            }
        }
        if(temp_atom != "C")
        {
            /** Modified: SJH **/
            cout << "This program is designed to work on molecules that have "
            "carbon as center atom. Please check the following atom_ID in the "
            "lammps input file:\n";
            cout << "Atom ID: " << central_atom[centerii] << "\tAtom type name: "
            << temp_atom << "\n";
            exit(EXIT_FAILURE);
        }
        int connect_counter=0;
        double connect_charges=0;
        for(int connectii=0; connectii<connect_atom[centerii].size(); connectii++)
        {
            //Determine the type of the central atom
            int temp_connect_type=0;
            for(int aaii=0; aaii<all_atom_id.size(); aaii++)
            {
                if(all_atom_id[aaii]==connect_atom[centerii][connectii])
                {
                    temp_connect_type = all_atom_type[aaii];
                    break;
                }
            }
            for(int atii=0; atii<atom_types.size(); atii++)
            {
                if(atom_types[atii]==temp_connect_type)
                {
                    if(atom_names[atii]!= "C")
                    {
                        connect_charges += type_charges[atii];
                        connect_counter++;
                    }
                    break;
                }
            }
        }
        double net_central_charge= center_atom_charge+connect_charges;
        //double net_central_charge= abs(center_atom_charge)-abs(connect_charges);
        if((round(net_central_charge*1000)/1000)!=0)
        {
            charged_center_atoms.push_back(central_atom[centerii]);
            num_connect_atoms.push_back(connect_counter);
            net_charge_store.push_back(net_central_charge);
            //cout << "There is a charge on carbon " << central_atom[centerii]
            //<< ". The charge will be adjusted to make that node neutral.\n";
        }
    }
    //--------------------------------------------------------------------------
    
    /**Isolating the nodes that have charges**/
    vector<int> temp_center_iso;
    vector<int> temp_connect_iso;
    vector<double> temp_charge_iso;
    /*Pairwise charge neurtalization (check upto 4 consecutive charges to see if
     they cancel out, if not store the charge for neutralization)*/
    for(int chargeii=0;chargeii<charged_center_atoms.size();chargeii++)
    {
        if(((chargeii+1)<charged_center_atoms.size())&&
           (((round((net_charge_store[chargeii]+
                     net_charge_store[chargeii+1])*1000))/1000)==0))
        {
            chargeii++;
        }
        else if(((chargeii+2)<charged_center_atoms.size())&&
                (((round((net_charge_store[chargeii]+
                          net_charge_store[chargeii+1]+
                          net_charge_store[chargeii+2])*1000))/1000)==0))
        {
            chargeii+=2;
        }
        else if(((chargeii+3)<charged_center_atoms.size())&&
                (((round((net_charge_store[chargeii]+
                          net_charge_store[chargeii+1]+
                          net_charge_store[chargeii+2]+
                          net_charge_store[chargeii+3])*1000))/1000)==0))
        {
            chargeii+=3;
        }
        else
        {
            temp_center_iso.push_back(charged_center_atoms[chargeii]);
            temp_connect_iso.push_back(num_connect_atoms[chargeii]);
            temp_charge_iso.push_back(net_charge_store[chargeii]);
            //cout << "Charge found on node surrounding carbon "
            //<< charged_center_atoms[chargeii] << ". This will be neutralized. "
            //<< net_charge_store[chargeii] << "\n";
        }
    }
    vector<vector<int>> center_isolation;
    vector<vector<int>> connect_isolation;
    vector<vector<double>> charge_isolation;
    vector<int> center_loader;
    vector<int> connect_loader;
    vector<double> charge_loader;
    for(int ii=0; ii<temp_center_iso.size();ii++)
    {
        if((ii!=temp_center_iso.size()-1) && ((temp_center_iso[ii+1]-temp_center_iso[ii])<3))
        {
            center_loader.push_back(temp_center_iso[ii]);
            connect_loader.push_back(temp_connect_iso[ii]);
            charge_loader.push_back(temp_charge_iso[ii]);
        }
        else
        {
            center_loader.push_back(temp_center_iso[ii]);
            connect_loader.push_back(temp_connect_iso[ii]);
            charge_loader.push_back(temp_charge_iso[ii]);
            
            center_isolation.push_back(center_loader);
            connect_isolation.push_back(connect_loader);
            charge_isolation.push_back(charge_loader);
            
            center_loader.clear();
            connect_loader.clear();
            charge_loader.clear();
        }
    }
    vector<int> altered_atom_ID;
    vector<double> adjust_charge;
    vector<int> associated_type;
    for(int ii=0; ii<center_isolation.size();ii++)
    {
        double sum_up_charge=0;
        int max_connected = *max_element(connect_isolation[ii].begin(),
                                         connect_isolation[ii].end());
        bool switcher=true;
        for(int jj=0; jj<center_isolation[ii].size();jj++)
        {
            if((max_connected==connect_isolation[ii][jj]) && (switcher==true))
            {
                cout
                << "Charge found on node surrounding carbon "
                << center_isolation[ii][jj] << ". This will be neutralized.\n";
                switcher=false;
                altered_atom_ID.push_back(center_isolation[ii][jj]);
                for(int aidii=0; aidii<all_atom_id.size();aidii++)
                {
                    if(center_isolation[ii][jj]==all_atom_id[aidii])
                    {
                        associated_type.push_back(all_atom_type[aidii]);
                    }
                }
            }
            sum_up_charge+=charge_isolation[ii][jj];
        }
        adjust_charge.push_back(sum_up_charge);
    }
    //--------------------------------------------------------------------------
    
    /**Create new atom types depending on carbons that needs to be neutralized**/
    /*NOTE:The original atom type is not deleted, an additional type will be
     added as a duplicate. Also the type of the atom that is altered is updated.*/
    for(int ii=0; ii<altered_atom_ID.size(); ii++)
    {
        //Check if an already updated atom type is being updated again;
        bool type_creation=false;
        int type_repeat=0;
        for(int checkii=0; checkii<duplication_atom_t_info.size(); checkii++)
        {
            //Check the atom type
            if(associated_type[ii]==duplication_atom_t_info[checkii])
            {
                //Check the net central charge around that node
                if(adjust_charge[ii]==duplication_atom_charge[checkii])
                {
                    type_creation=true;
                    type_repeat= checkii;
                    break;
                }
            }
        }
        int new_atom_t =0;
        double charge_adjusted =0;
        if(!type_creation)
        {
            for(int oii=0; oii< atom_types.size();oii++)
            {
                if(associated_type[ii]==atom_types[oii])
                {
                    charge_adjusted = type_charges[oii] + (-1*adjust_charge[ii]);
                }
            }
            if(new_atom_types.size()==0)
            {
                //Create a new atom type
                new_atom_t = atom_types[atom_types.size()-1]+1;
            }
            else
            {
                //Create a new atom type
                new_atom_t = new_atom_types[new_atom_types.size()-1]+1;
            }
            //Store the created new atom type
            new_atom_types.push_back(new_atom_t);
            //Store the charge associated with the new atom type
            new_atom_charges.push_back(charge_adjusted);
            //Store the original atom type information that was updated
            duplication_atom_t_info.push_back(associated_type[ii]);
            //Store the orginal charge associated with that node that was updated
            duplication_atom_charge.push_back(adjust_charge[ii]);
        }
        else
        {
            new_atom_t = new_atom_types[type_repeat];
        }
        
        for(int updateii=0; updateii<all_atom_id.size(); updateii++)
        {
            if(altered_atom_ID[ii]==all_atom_id[updateii])
            {
                //Update the atom type to new atom type;
                all_atom_type[updateii] = new_atom_t;
                break;
            }
        }
    }
    //--------------------------------------------------------------------------
    
    /**Write new lammps data file**/
    ofstream write_lmp_file;
    write_lmp_file.open(o_lammps_file.c_str());
    bool mass_section= false;
    bool atom_section= false;
    int aatii=0;
    vector<string> mass_type_charge;
    mass_type_charge.resize(duplication_atom_t_info.size());
    ifstream read_lmp_file(lammps_file.c_str());
    if(read_lmp_file.is_open())
    {
        string lmp_line;
        while(getline(read_lmp_file,lmp_line))
        {
            if(lmp_line.size()!=0)
            {
                if(lmp_line== "Masses") {
                    mass_section=true;
                    atom_section=false;
                    write_lmp_file << lmp_line << "\n";
                    continue;
                } else if((lmp_line == "Atoms") || (lmp_line == "Atoms # full")) {
                    mass_section=false;
                    atom_section=true;
                    /** Modified: SJH **/
                    //write_lmp_file << lmp_line << "\n";
                    write_lmp_file << "Atoms" << "\n";
                    continue;
                } else if((lmp_line == "Bonds")) {
                    mass_section=false;
                    atom_section=false;
                    write_lmp_file << lmp_line << "\n";
                    continue;
                }
                vector<string> stringvector;
                string stringelements;
                istringstream checkstring (lmp_line);
                while (checkstring >> stringelements)
                {
                    /*storing every entity of the line (checkstring) in a vector;
                     The vector is overwritten for every line*/
                    stringvector.push_back(stringelements);
                }
                if(mass_section)
                {
                    for(int dupii=0; dupii<duplication_atom_t_info.size(); dupii++)
                    {
                        string temp_store;
                        temp_store.clear();
                        if(duplication_atom_t_info[dupii]==atoi(stringvector[0].c_str()))
                        {
                            for(int storeii=1; storeii<stringvector.size();storeii++)
                            {
                                if(storeii == stringvector.size()-1)
                                {
                                    for(int strii=0; strii<stringvector[stringvector.size()-1].size(); strii++)
                                    {
                                        if(stringvector[stringvector.size()-1][strii]=='=')
                                        {
                                            temp_store+=stringvector[stringvector.size()-1][strii];
                                            break;
                                        }
                                        else
                                        {
                                            temp_store+=stringvector[stringvector.size()-1][strii];
                                        }
                                    }
                                }
                                else
                                {
                                    temp_store+=(stringvector[storeii] + " ");
                                }
                            }
                            mass_type_charge[dupii] = temp_store;
                        }
                    }
                    if(atom_types[atom_types.size()-1] == atoi(stringvector[0].c_str()))
                    {
                        write_lmp_file << lmp_line << "\n";
                        for(int writeii=0; writeii<new_atom_types.size(); writeii++)
                        {
                            write_lmp_file
                            << "    "
                            << new_atom_types[writeii] << " "
                            << mass_type_charge[writeii]
                            << new_atom_charges[writeii] << "\n";
                        }
                    } else {
                        write_lmp_file << lmp_line << "\n";
                    }
                }
                else if(atom_section)
                {
                    write_lmp_file
                    << stringvector[0] << " "
                    << stringvector[1] << " "
                    << all_atom_type[aatii] << " "
                    << stringvector[3] << " "
                    << stringvector[4] << " "
                    << stringvector[5] << " "
                    << stringvector[6] << "\n";
                    aatii++;
                }
                else if(stringvector.size()==3)
                {
                    if((stringvector[1]=="atom")&&(stringvector[2]=="types"))
                    {
                        write_lmp_file
                        << "     "
                        << (atom_types.size()+new_atom_types.size())
                        << "  atom types\n";
                    } else {
                        write_lmp_file << lmp_line << "\n";
                    }
                }
                else
                {
                    write_lmp_file << lmp_line << "\n";
                }
            }
            else
            {
                write_lmp_file << lmp_line << "\n";
            }
        }
        read_lmp_file.close();
    }
    else
    {
        /** Modified: SJH **/
        cout << "The lammps data file don't exist. Please check the following "
        "path to the file\n" << lammps_file << "\n";
        exit(EXIT_FAILURE);
    }
    write_lmp_file.close();	//Close the newly written lammps data input
    //--------------------------------------------------------------------------
    
    /**Write the charges output file in lammps accepted format**/
    ofstream write_charge_file;
    write_charge_file.open(o_charge_file.c_str());
    for (int ii=0; ii<atom_types.size();ii++)
    {
        write_charge_file
        << "    set type "
        << atom_types[ii] << " charge "
        << ((round(type_charges[ii]*1000))/1000) << "\n";
    }
    for (int ii=0; ii<new_atom_types.size();ii++)
    {
        write_charge_file
        << "    set type "
        << new_atom_types[ii] << " charge "
        << ((round(new_atom_charges[ii]*1000))/1000) << "\n";
    }
    write_charge_file.close();	//Close the newly written charge data file
    //--------------------------------------------------------------------------
    
    /**Write the lammps setting file in lammps accepted format**/
    vector<string> setting_duplication;
    setting_duplication.resize(duplication_atom_t_info.size());
    ofstream write_pair_file;
    write_pair_file.open(o_pair_file.c_str());
    ifstream read_pair_file(pair_file.c_str());
    if(read_pair_file.is_open())
    {
        string pair_line;
        while(getline(read_pair_file,pair_line))
        {
            /** Modified: SJH **/
            vector<string> stringvector;
            if(pair_line.size()!=0)
            {
                //vector<string> stringvector;
                string stringelements;
                istringstream checkstring (pair_line);
                while (checkstring >> stringelements)
                {
                    /*storing every entity of the line (checkstring) in a vector;
                     The vector is overwritten for every line*/
                    stringvector.push_back(stringelements);
                }
                for(int dupii=0; dupii<duplication_atom_t_info.size(); dupii++)
                {
                    string temp_store;
                    temp_store.clear();
                    if((duplication_atom_t_info[dupii]==atoi(stringvector[2].c_str()))&&
                       (stringvector[0] == "pair_coeff"))
                    {
                        for(int storeii=3; storeii<stringvector.size(); storeii++)
                        {
                            /** Modified: SJH **/
                            if(stringvector.at(storeii)=="lj/cut/coul/long") continue;
                            temp_store+= (stringvector[storeii] + " ");
                        }
                        setting_duplication[dupii] = temp_store;
                    }
                }
                if((atom_types[atom_types.size()-1]==atoi(stringvector[2].c_str()))&&
                   (stringvector[0] == "pair_coeff"))
                {
                    /** Modified: SJH **/
                    for (indexi=0; indexi<(int)stringvector.size(); ++indexi) {
                        if (stringvector.at(indexi)=="lj/cut/coul/long") continue;
                        write_pair_file << stringvector.at(indexi) << " ";
                    } write_pair_file << "\n";
                    //write_pair_file << pair_line << "\n";
                    
                    for(int writeii=0; writeii<new_atom_types.size(); writeii++)
                    {
                        /** Modified: SJH **/
                        write_pair_file
                        //<< "    pair_coeff "
                        << "pair_coeff "
                        << new_atom_types[writeii] << " "
                        << new_atom_types[writeii] << " "
                        << setting_duplication[writeii] << "\n";
                    }
                }
                else
                {
                    /** Modified: SJH **/
                    for (indexi=0; indexi<(int)stringvector.size(); ++indexi) {
                        if (stringvector.at(indexi)=="lj/cut/coul/long") continue;
                        write_pair_file << stringvector.at(indexi) << " ";
                    } write_pair_file << "\n";
                    //write_pair_file << pair_line << "\n";
                }
            }
            else
            {
                /** Modified: SJH **/
                for (indexi=0; indexi<(int)stringvector.size(); ++indexi) {
                    if (stringvector.at(indexi)=="lj/cut/coul/long") continue;
                    write_pair_file << stringvector.at(indexi) << " ";
                } write_pair_file << "\n";
                //write_pair_file << pair_line << "\n";
            }
        }
        read_pair_file.close();
    }
    else
    {
        /** Modified: SJH **/
        cout << "The lammps settings file don't exist. "
        "Please check the following path to the file\n" << pair_file << "\n";
        exit(EXIT_FAILURE);
    }
    write_pair_file.close();
    ////////////////////////////////////////////////////////////////////////////
    //--------------------------------------------------------------------------
}





/** public setters **/
//----------------------------------------------------------------------
/* bool */
void MoltemplateLmpData::set_is_restrictbynAtoms(const bool b)
{is_restrictbynAtoms=b;}
void MoltemplateLmpData::set_is_optimizeCharge(const bool b)
{is_optimizeCharge=b;}
void MoltemplateLmpData::set_is_mass_normalization(const bool b)
{is_mass_normalization=b;}
void MoltemplateLmpData::set_is_correctAlkylDihedral(const bool b)
{is_correctAlkylDihedral=b;}
/* int */
void MoltemplateLmpData::set_dop(const int i){dop=i;}
void MoltemplateLmpData::set_sequenceNum(const int i){sequenceNum=i;}
void MoltemplateLmpData::set_sequenceLen(const int i){sequenceLen=i;}
void MoltemplateLmpData::set_n_total_atoms(const int i){n_total_atoms=i;}
/* double */
void MoltemplateLmpData::set_moltemplateBoxSize(const double d){moltemplateBoxSize=d;}
/* string */
void MoltemplateLmpData::set_merA(const string& str){merA=str;}
void MoltemplateLmpData::set_merB(const string& str){merB=str;}
void MoltemplateLmpData::set_copolymerType(const string& str){copolymerType=str;}
void MoltemplateLmpData::set_tacticity(const string& str){tacticity=str;}
/* vector<string> */
void MoltemplateLmpData::set_merSet(const vector<string>& vstr){merSet=vstr;}
void MoltemplateLmpData::set_chem_mono(const vector<string>& vstr){chem_mono=vstr;}
/* vector<int> */
void MoltemplateLmpData::set_merComp(const vector<int>& vi){merComp=vi;}
void MoltemplateLmpData::set_mono_beads(const vector<int>& vi){mono_beads=vi;}
/* vector<vector<string>> */
void MoltemplateLmpData::set_sequenceSet(const vector<vector<string>>& vvs)
{
    string out=path_cwd+"sequenceSet.txt";
    ofstream writeFile(out.c_str());
    sequenceSet=vvs;
    for (indexi=0; indexi<(int)vvs.size(); ++indexi) {
        //writeFile<<indexi+1<<" ";
        for (indexii=0; indexii<(int)vvs.at(indexi).size(); ++indexii) {
            //writeFile <<"#"<<indexii+1<<"("<<vvs.at(indexi).at(indexii)<<") ";
            writeFile << vvs.at(indexi).at(indexii) << " ";
        } writeFile << "\n";
    } writeFile.close();
}




