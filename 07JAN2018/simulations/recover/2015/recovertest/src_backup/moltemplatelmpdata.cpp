//
//  moltemplatelmpdata.cpp
//  cppWork
//
//  Created by SJH on 12/18/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#include "moltemplatelmpdata.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

MoltemplateLmpData::MoltemplateLmpData(const StructureClass& sysVar):
/* string */
path_master(sysVar.get_Path()+"/"),
/* int */
sequenceNum(0),
sequenceLen(0),
/* double */
moltemplateBoxSize(300.0),
offset(4.0), /**offset distance between each consecutive monomer**/
rotate(180.0) /**rotating angle between each consecutive monomer**/
{
    /** Path to moltemplate/src/ (NOTE: need "/" at the end of path) **/
    path_moltemplatesrc = path_master+"moltemplate/src/";
    
    /** Path to the monomer bank (NOTE: need "/" at the end of path) **/
    path_monomerBank = path_master+"Monomer_bank/";
    
    /** Path to the oplsaa.prm file (full oplsaa parameter table) **/
    path_oplsaaprm = path_master+"moltemplate/oplsaa.prm";
    
    
    /*** NOTE: <DON'T CHANGE!> start_data/moltemplate/ is the cwd ***/
    path_cwd=return_SimulationFolderPath(sysVar)+"/lammps_inputs/start_data/moltemplate/";
    n_trials=sysVar.get_n_trial();
}





void MoltemplateLmpData::make_LmpDataFilebyMoltemplate(const StructureClass& sysVar)
{
    /** check number of molecules **/
    if(sequenceNum>0) {
        if ((int)sequenceSet.size()!=sequenceNum) {
            cout
            << "\nWarning: Number of Molecules != " << sequenceNum
            << "\n\n";
        }
    }
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
        }
        /** check monomer.lt's in the monomer bank **/
        for (indexii=0; indexii<(int)sequenceSet.at(indexi).size(); ++indexii) {
            
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
        make_polylt(indexi,sequenceSet.at(indexi));
    }
    
    /** make system.lt file **/
    make_systemlt();
    
    /** invoke moltemplate to generate LAMMPS datafile **/
    invoke_moltemplate();
    
    /** charge optimization **/
    optimize_charge();
    
    /** move output files to due places **/
    mv_files();
}





vector<vector<string>> MoltemplateLmpData::read_sequenceSetfromfile(const StructureClass& sysVar)
{
    vector<string> vs;
    vector<vector<string>> vvs;
    int countline=0,countnum=0;
    
    string input;
    input.append(sysVar.get_Path());
    input.append("/test_sequences.txt");
    //input.append("/test_sequences2.txt");
    
    ifstream readfile(input.c_str());
    if (readfile.is_open())
    {
        string lineContent;
        string monomer;
        while (getline(readfile,lineContent)) {
            ++countline;
            istringstream iss(lineContent);
            countnum=0; vs.clear();
            while (iss>>monomer) {
                ++countnum;
                vs.push_back(monomer);
            } vvs.push_back(vs);
        } readfile.close();
    } else {
        cout << "\nMoltemplateLmpData::read_sequenceSet: readfile can't open.\n";
    } return vvs;
}





bool MoltemplateLmpData::check_monomerbank(const string& monomer)
{
    string file;
    file.append(path_monomerBank);
    file.append(monomer);
    ifstream readfile(file.c_str());
    return readfile.good();
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





void MoltemplateLmpData::check_path(const string& path)
{
    string file;
    file.append(path);
    ifstream readfile(file.c_str());
    if (!readfile) {
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
    
    ofstream writefile(output.c_str());
    if (writefile.is_open())
    {
        /** import oplsaa.lt **/
        writefile << "import \"oplsaa.lt\"\n";
        
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
                writefile << "import \""+monomerSet.at(indexii)+"\"\n";
            } used_monomers.push_back(monomerSet.at(indexii));
        } writefile << "\n";
        
        /** Define combined molecule (ex.polymer) **/
        writefile << "poly_"<<polyindex+1<<" inherits OPLSAA {\n\n";
        writefile << "    "<< "create_var {$mol}\n\n";
        
        vector<string> monomerSet_copy=monomerSet;
        for (indexii=0; indexii<(int)monomerSet.size(); ++indexii) {
            
            /** erase .lt from name string **/
            monomerSet_copy.at(indexii).erase
            (monomerSet_copy.at(indexii).end()-3,monomerSet_copy.at(indexii).end());
            
            /** pack monomers along x-axis and rotate accordingly (1,0,0) **/
            writefile
            << "    "
            << "monomer["<<indexii<<"] = new "<<monomerSet_copy.at(indexii);
            if (indexii>0) {
                writefile
                << ".rot(" <<rotate<<",1,0,0)"
                << ".move("<<offset*indexii<<",0,0)";
            } writefile << "\n";
        }
        
        /** add a list of bonds connecting propagating carbons **/
        writefile << "\n    write('Data Bond List') {\n";
        for (indexii=0; indexii<(int)monomerSet.size()-1; ++indexii) { //NOTE: "DOP-1" bonds
            writefile
            << "      "
            << "$bond:b"<<indexii+1<<"  "
            << "$atom:monomer["<<indexii<<"]/C2"<<"  "
            << "$atom:monomer["<<indexii+1<<"]/C1"<<"  "
            << "\n";
        } writefile << "    }\n";
        
        /** end cap of poly.lt scope **/
        writefile
        << "\n} # poly_"<<polyindex+1 << "\n";
        
        writefile.close();
    } else {
        cout << "\nMoltemplateLmpData::make_polyltFile: writefile can't open.\n";
    }
}





void MoltemplateLmpData::make_systemlt()
{
    string output;
    output.append(path_cwd+"/system.lt");
    
    ofstream writefile(output.c_str());
    if (writefile.is_open())
    {
        int n_poly=(int)sequenceSet.size();
        for (indexi=0; indexi<n_poly; ++indexi) {
            writefile
            << "import \"poly_"+to_string((long long int)indexi+1)+".lt\"\n";
        } writefile << "\n";
        
        double dely=25.0;
        double delz=20.0;
        double sign=1.0;
        int county=0,countz=0;
        for (indexi=0; indexi<n_poly; ++indexi) {
            writefile
            << "polymer_"<<indexi+1<<" = new "
            << "poly_"<<indexi+1<<".move(-50,"
            << sign*dely*(double)county<<","
            << sign*delz*(double)countz<<")"<< "\n";
            ++county;
            if (county%4==0) {
                sign *= -1.0;
                county=0;
                ++countz;
            }
        } writefile << "\n";
        
        double hbox=moltemplateBoxSize*0.5;
        double fbox=moltemplateBoxSize;
        
        if (true) {
            writefile
            << "write_once(\"Data Boundary\") {"     << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  xlo xhi" << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  ylo yhi" << "\n"
            << "   -"<<hbox<<"  "<<hbox<<"  zlo zhi" << "\n"
            << "}"                                   << "\n\n";
        } else {
            writefile
            << "write_once(\"Data Boundary\") {"     << "\n"
            << "   0.0  "<<fbox<<"  xlo xhi"         << "\n"
            << "   0.0  "<<fbox<<"  ylo yhi"         << "\n"
            << "   0.0  "<<fbox<<"  zlo zhi"         << "\n"
            << "}"                                   << "\n\n";
        }
        
        writefile.close();
    } else {
        cout << "\nMoltemplateLmpData::make_systemltFile: writefile can't open.\n";
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
    for (indexii=0; indexii<n_trials; ++indexii) {
        string data=
        "cd "+path_cwd+"; cp system_updated.data ../; "+
        "cd ../; "
        "mv system_updated.data input_00"+to_string((long long int)indexii)+".data";
        system(data.c_str());
    }
    
    string init=
    "cd "+path_cwd+";"
    "cp system_updated.in.charges system_updated.in.settings system.in "
    "system.in.init ../../generation/;"
    "cd ../../generation/;"
    "mv system_updated.in.charges  system.in.charges;"
    "mv system_updated.in.settings system.in.settings";
    system(init.c_str());
    
    string output=
    "cd "+path_cwd+";"
    "mkdir output_updated; mv *_updated.* output_updated/;"
    "mkdir output; mv system.data system.in system.in.charges system.in.init "
    "system.in.settings output/";
    system(output.c_str());
    
    string input=
    "cd "+path_cwd+";"
    "mkdir input; mv *.lt *.prm input/";
    system(input.c_str());
    
    
    string scriptsource=path_master+"scripts";
    string simfolder=path_cwd+"../../../../simulations/";
    string anafolder=path_cwd+"../../../../analysis/";
    string cpfiles=
    "cp "+scriptsource+"/generation.inp "   +simfolder+"lammps_inputs/generation/;"
    "cp "+scriptsource+"/quench.inp "       +simfolder+"lammps_inputs/quench/;"
    "cp "+scriptsource+"/equilibration.inp "+simfolder+"lammps_inputs/equilibration/;"
    "cp "+scriptsource+"/production.inp "   +simfolder+"lammps_inputs/production/;"
    "cp "+scriptsource+"/amdat.inp "        +anafolder+"AMDAT_inputs/;";
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
    /** thresholding mass (a.m.u.) **/
    double light=3.0;
    
    types_light.clear(); types_heavy.clear();
    string data=path_cwd+"output_updated/system_updated.data";
    ifstream readfile(data.c_str());
    if(readfile.is_open())
    {
        string lineContent;
        string strvar;
        double dubvar=0;
        bool is_inblock=false;
        while (getline(readfile,lineContent)) {
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
        } readfile.close();
        check_lasttwo(types_all); sysVar.set_types_all(types_all);
        check_lasttwo(types_light); sysVar.set_types_light(types_light);
        check_lasttwo(types_heavy); sysVar.set_types_heavy(types_heavy);
    } else {
        cout <<
        "in MoltemplateLmpData::set_atomTypes(), readfile can NOT open. "
        "Please check.\n";
        //exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::set_atomNumbers(StructureClass& sysVar)
{
    string data=path_cwd+"output_updated/system_updated.data";
    ifstream readfile(data.c_str());
    n_light=0; n_heavy=0;
    
    /** create the data structure used to store atom_types and thier respective
     numbers (in each molecule) **/
    initialize_n_typesVec();
    
    if(readfile.is_open())
    {
        string lineContent;
        string strvar,str_now,str_pre;
        int intvar=0,n_poly=0;
        bool is_inblock=false;
        while (getline(readfile,lineContent)) {
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
        } readfile.close();
        
        if (n_poly!=(int)sequenceSet.size()) {
            cout
            << "in MoltemplateLmpData::set_atomNumbers, "
            << "molid != n_poly" << "\n";
            exit(EXIT_FAILURE);
        }
        //cout << "n_total " << n_heavy+n_light << "\n";
        //cout << "n_heavy " << n_heavy << "\n";
        //cout << "n_light " << n_light << "\n";
        
        sysVar.set_n_poly((int)sequenceSet.size());
        n_heavy /= n_poly;
        n_light /= n_poly;
        sysVar.set_n_light(n_light);
        sysVar.set_n_heavy(n_heavy);
        
        /** n_typeVecs have the struture that shows the number of atoms for each
         corresponding type of atom;
         i.e. per 1-D element: (type,number of type per molecule) **/
        finalize_n_typesVec(n_poly);
        sysVar.set_n_types_all(n_types_all);
        sysVar.set_n_types_heavy(n_types_heavy);
        sysVar.set_n_types_light(n_types_light);
        
    } else {
        cout <<
        "in MoltemplateLmpData::set_atomNumbers(), readfile can NOT open. "
        "Please check.\n";
        //exit(EXIT_FAILURE);
    }
}





void MoltemplateLmpData::check_lasttwo(std::vector<int>& vi)
{
    size_t svi=vi.size();
    if (vi.at(svi-1)==vi.at(svi-2)) vi.pop_back();
}





void MoltemplateLmpData::make_oplsaa_subset(const vector<string>& monomerSet)
{
    /** Primary contributor: Venkatesh Meenakshisundaram
     Integration into master code: Jui-Hsiang Hung (SJH) **/
    
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
                    if(mono_line == "	write(\"Data Atoms\") {")
                    {
                        read_switch = true;
                        continue;
                    }
                    else if(mono_line =="	}")
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
                        for(int readii=0; readii<stringvector[2].size(); readii++)
                        {
                            if(stringvector[2][readii]==':')
                            {
                                load_switch = true;
                                continue;
                            }
                            if(load_switch)
                            {
                                load_line += stringvector[2][readii];
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





void MoltemplateLmpData::optimize_charge()
{
    /** Primary contributor: Venkatesh Meenakshisundaram
     Integration into master code: Jui-Hsiang Hung (SJH) **/
    
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





/*==( public setters )==*/
void MoltemplateLmpData::set_sequenceNum(const int i){sequenceNum=i;}
void MoltemplateLmpData::set_sequenceLen(const int i){sequenceLen=i;}
/* double */
void MoltemplateLmpData::set_moltemplateBoxSize(const double d){moltemplateBoxSize=d;}
/* vector<vector<string>> */
void MoltemplateLmpData::set_sequenceSet(const std::vector<std::vector<std::string>>& vvs)
{
    bool is_cout=false;
    sequenceSet=vvs;
    for (indexi=0; indexi<(int)vvs.size(); ++indexi) {
        if (is_cout) cout<<indexi+1<<" ";
        for (indexii=0; indexii<(int)vvs.at(indexi).size(); ++indexii) {
            if (is_cout) cout<<"#"<<indexii+1<<"("<<vvs.at(indexi).at(indexii)<<") ";
        } if (is_cout) cout << "\n";
    }
}


