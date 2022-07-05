/*******************************************************************************
 PreSQ & AUTOWORK
 Copyright (c) Jui-Hsiang "Sean" Hung (SJH). All rights reserved.
 
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

#include "userinterface.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

UserInterface::UserInterface(const vector<string>& args):
//NOTE: Class() --> 'value-initialize' every member in class
var(VarClass()),
sysVar(StructureClass()),
moldata(MoltemplateLmpData()),
ds(FitDielectric()),
ws(WorkScripts()),
aa(AmdatAnalysis()),
fd(FitData()),
tt(TheoryTest()),

/* bool */
is_inpVar_successful(false),
is_inforloop(false),
is_action(false),
/* int */
current_line_number(0),
indi(0),indii(0),
argi(0),argii(0),//NOTE:argi is reserved for use in parse_commandLineArgs()
n_fors(0),n_ifs(0),
/* string */
current_readFileName(),
lineContent(),
command(),
/* STL */
strArgs(args),
veclineContents()
{
    sysVar.set_startdatetime(return_datetime());
    ws.set_numCores({1,1,1,1,1});
    ws.set_run_cores({1,1,1,1,1});
}





void UserInterface::show_version()
{
    cout
    << "\n\n"
    << "************************************************************"
    << "\n"
    << "*             -- PreSQ & AUTOWORK ("<<VERSION<<") --           *"
    << "\n"
    << "*  Author: Jui-Hsiang (Sean) Hung, Simmons Research Group  *"
    << "\n"
    << "*          http://blogs.uakron.edu/simmonsgroup/           *"
    << "\n"
    << "************************************************************"
    << "\n\n\n";
    system_wait(1);
}





void UserInterface::show_help_info()
{
    cout
    << "Current supported options:\n"
    << "NOTE: what's inside [ ] is the argument\n\n"
    << "-i or -inp [path]: path to the input file\n\n"
    << "-v or -var [name value1 value2 ...]: variable name and its value(s) to "
    << "be replaced in input file\n\n"
    << "-s or -screen [path]: path to the screen file\n\n"
    << "-version: check software version\n\n";
}





void UserInterface::check_var_format(const string& var)
{
    char i=var[0];
    if (i=='0'||i=='1'|i=='2'|i=='3'|i=='4'|i=='5'|i=='6'|i=='7'|i=='8'|i=='9')
    {
        cout<<"Error: for variable ("+var+"), it should start with a letter.\n";
        exit(EXIT_FAILURE);
    }
}





double UserInterface::check_n_args_comkeyword()
{
    int offset=1,n_args=0;
    if ((argi+offset)>strArgs.size()-1||strArgs.at(argi+offset).at(0)=='-') {
        cout<< "Error: No arg found for commandline keyword ("<<strArgs.at(argi)<<").\n";
        exit(EXIT_FAILURE);
    }
    do { ++n_args; ++offset;
        if ((argi+offset)>strArgs.size()-1) break;
    } while (strArgs.at(argi+offset).at(0)!='-');
    
    return n_args;
}





void UserInterface::show_inplineContent()
{
    cout<<"line "<<current_line_number<<" "<<lineContent<<"\n";
}





void UserInterface::check_readFile_n_args(const std::string& lineContent,const int n_args)
{
    int n_args_read=0;
    string strVar,command;
    istringstream iss(lineContent);
    while (iss>>strVar) {
        if(n_args_read==0)command=strVar;
        ++n_args_read;
    }
    if (n_args_read!=n_args) {
        cout
        <<"in "<<current_readFileName<<":\n"
        <<"at/around line "<<current_line_number<<":\n"
        <<"current command ("+command+") takes "<<n_args-1<<" arguments;\n"
        <<n_args_read-1<<" arguments are found. Please check.\n\n";
        exit(EXIT_FAILURE);
    }
}





void UserInterface::check_readFile_min_args(const std::string& lineContent,const int min)
{
    int n_args_read=0;
    string strVar,command;
    istringstream iss(lineContent);
    while (iss>>strVar) {
        if(n_args_read==0)command=strVar;
        ++n_args_read;
    }
    if (n_args_read<min) {
        cout
        <<"in "<<current_readFileName<<":\n"
        <<"at/around line "<<current_line_number<<":\n"
        <<"current command ("+command+") takes at least "<<min-1<<" arguments;\n"
        <<n_args_read-1<<" arguments are found. Please check.\n\n";
        exit(EXIT_FAILURE);
    }
}





int UserInterface::get_n_args(const std::string& lineContent)
{
    int n_args_read=0;
    string strVar,command;
    istringstream iss(lineContent);
    while (iss>>strVar) {
        if(n_args_read==0)command=strVar;
        ++n_args_read;
    } return n_args_read;
}





std::string UserInterface::get_lineContent_at_line(const string& in,const int line)
{
    int linecount=0;
    string lineContent;
    string str;
    ifstream readFile(in.c_str());
    if (readFile.is_open()) {
        while (getline(readFile,lineContent)) {
            ++linecount;
            istringstream iss(lineContent);
            if (linecount==line) {
                str=get_effectiveStr(lineContent); break;
            }
        } readFile.close();
    } else {
        cout
        << "in UserInterface::get_lineContent_at_line():\n"
        << in<<" cannot open.\n";
        exit(EXIT_FAILURE);
    } return str;
}





string UserInterface::get_partialstr(const string& str,
                                     const int beg,
                                     const int end)
{
    string ret;
    for (int i=beg;i<=end;++i) {
        ret+=str.at(i);
    } return ret;
}





string UserInterface::get_line_comp(const string& lineContent,
                                    const int pos,
                                    const bool check)
{
    int n_args_read=0;
    string strVar,command;
    vector<string> strVec;
    istringstream iss(lineContent);
    while (iss>>strVar) {
        if(n_args_read==0)command=strVar;
        strVec.push_back(strVar);
        ++n_args_read;
    }
    if (pos<=n_args_read-1&&pos>=0) {
        return str_from_parse_rule00(strVec.at(pos));
    } else {
        if (check) {
            cout
            << "in UserInterface::get_line_comp():\n"
            << "at/around line "<<current_line_number<<": "
            << "allowable component position is between [0,"<<n_args_read-1<<"]\n\n";
            exit(EXIT_FAILURE);
        } else {
            return "";
        }
    }
}





string UserInterface::get_effectiveStr(const string& inpstr)
{
    // tokenize a composite string into smaller pieces of strings;
    // ignore content after #
    string str;
    vector<string> tstrings;
    istringstream iss(inpstr);
    while (iss>>str) {
        if (str.at(0)!='#') {
            tstrings.push_back(str);
        } else {
            break;
        }
    } str.clear();
    if (tstrings.size()==0) {
        str="";
    } else {
        for (size_t i=0;i<tstrings.size();++i) {
            str+=tstrings.at(i);
            str+=" ";//dlm:space
        }
    } return str;
}





void UserInterface::parse_commandLineArgs()
{
    if (strArgs.size()>1)
    {
        /** read vars to be replaced in input file (if any) **/
        for (argi=0;argi<(int)strArgs.size();++argi)
        {
            string longKeyword;
            if (strArgs.at(argi).size()>2) {
                for (int s=1;s<(int)strArgs.at(argi).size();++s) {
                    longKeyword+=strArgs.at(argi).at(s);
                }
            }
            if ((strArgs.at(argi).size()==2&&strArgs.at(argi).at(1)=='v')||
                longKeyword=="var")
            {
                /** store vars to be replaced in input file **/
                //NOTE:
                //variables use in commandline are stored in the same container
                //as those defined in input script. If same variable name is
                //used in both command line and input script, the one defined in
                //input script will be used.
                if (check_n_args_comkeyword()>=2) {
                    string var_name=strArgs.at(argi+1);
                    int n_vals=check_n_args_comkeyword()-1;
                    vector<string> var_vals;
                    for (int i=1;i<=n_vals;++i) {
                        if (strArgs.at(argi+1+i).at(0)=='!') {//NOTE: igonore first '!'
                            string tmp=strArgs.at(argi+1+i);
                            tmp.erase(tmp.begin(),tmp.begin()+1);
                            strArgs.at(argi+1+i)=tmp;
                        } var_vals.push_back(strArgs.at(argi+1+i));
                    } var.add_new_var(var_name,var_vals);
                } else {
                    check_var_format(strArgs.at(argi+1));
                    cout<<"no value(s) for variable ("<<strArgs.at(argi+1)<<").\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        /** read input file **/
        for (argi=0;argi<(int)strArgs.size();++argi)
        {
            if (strArgs.at(argi).at(0)=='-')
            {
                string longKeyword;
                if (strArgs.at(argi).size()>2) {
                    for (int s=1;s<(int)strArgs.at(argi).size();++s) {
                        longKeyword+=strArgs.at(argi).at(s);
                    }
                }
                if ((strArgs.at(argi).size()==2&&strArgs.at(argi).at(1)=='i')||
                    longKeyword=="inp")
                {
                    if (check_n_args_comkeyword()!=1) {
                        cout<<"There should be one input file. Please check.\n\n";
                        exit(EXIT_FAILURE);
                    }
                    /** call input file parser **/
                    parse_inputFile(strArgs.at(argi+1));
                    is_inpVar_successful=true;
                }
                else if ((strArgs.at(argi).size()==2&&strArgs.at(argi).at(1)=='s')||
                         longKeyword=="screen")
                {
                    /** screen file collects both cpp and bash output **/
                    //TODO
                }
                else if ((strArgs.at(argi).size()==2&&strArgs.at(argi).at(1)=='h')||
                         longKeyword=="help")
                {
                    /** help shows supported command-line options and other usefuls **/
                    show_help_info();
                    exit(EXIT_SUCCESS);
                }
                else if (longKeyword=="version")
                {
                    cout << VERSION << "\n\n";
                    exit(EXIT_SUCCESS);
                }
                else
                {
                    cout
                    << "Option ("<<strArgs.at(argi)<<") not found. Please check. "
                    << "Use -h for help.\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    else
    {
        is_inpVar_successful=false;
    }
}





void UserInterface::parse_inputFile(const string& in)
{
    current_readFileName=in;
    
    ifstream readFile(in.c_str());
    //auto start_position=readFile.tellg();
    
    if (readFile.is_open())
    {
        /** I. STORE LINES FOR LATER USE **/
        veclineContents.clear();
        while (getline(readFile,lineContent)) {
            veclineContents.push_back(lineContent);
        } readFile.close();
        //readFile.clear();
        //readFile.seekg(start_position);
        
        string strVar,keyword;
        int n_total_lines=(int)veclineContents.size();
        
        if (false)
        {
            /** II. SYNTAX CHECK; NO ACTION **/
            is_action=false;
            current_line_number=0;
            for (current_line_number=1;current_line_number<=n_total_lines;
                 ++current_line_number)
            {
                lineContent=veclineContents.at(current_line_number-1);
                lineContent=get_effectiveStr(lineContent);
                //show_inplineContent();
                if (lineContent!="") {
                    command=get_line_comp(lineContent,0);
                    command_parser(command);
                }
            }
        }
        
        /** III. EXECUTE COMMANDS; TAKE ACTION **/
        is_action=true;
        for (current_line_number=1;current_line_number<=n_total_lines;
             ++current_line_number)
        {
            lineContent=veclineContents.at(current_line_number-1);
            lineContent=get_effectiveStr(lineContent);
            if (lineContent!="") {
                command=get_line_comp(lineContent,0);
                command_parser(command);
            }
        }
        
    } else {
        cout << "\ninput file to AUTOWORK ("+in+") cannot open.\n\n";
        exit(EXIT_FAILURE);
    }
}





void UserInterface::command_parser(const string& command)
{
    //show_inplineContent();
    
    string strVar,keyword,arg;
    
    if (command=="exit") {
        //if (is_action) {
        cout
        <<"\nat line "<<current_line_number<<": "
        <<"AUTOWORK stopped at user's request.\n";
        exit(EXIT_SUCCESS);
        //}
    }
    
    //TODO: COMMAND LISTS (1.necessary,2.prompts)
    
    /** STATIC: SYS_CONFIG **/
    //**************************************************************
    if (command=="typeB") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_typeB(atoi(strVar.c_str()));
    }
    else if (command=="typeS") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_typeS(atoi(strVar.c_str()));
    }
    else if (command=="userName") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_userName(strVar);
    }
    else if (command=="masterFolder") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_masterFolder(strVar);
    }
    else if (command=="simType") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_simType(strVar);
    }
    else if (command=="year") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_year(strVar);
    }
    else if (command=="usicID") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_usicID(strVar);
    }
    else if (command=="nameString") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        //sysVar.set_is_use_named_nameString(true);
        sysVar.set_nameString(strVar);
    }
    else if (command=="n_trial") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_n_trial(atoi(strVar.c_str()));
    }
    else if (command=="systemUnit") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_systemUnit(strVar);
    }
    else if (command=="startTemp") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_startTemp(atof(strVar.c_str()));
    }
    else if (command=="crossTemp") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_crossTemp(atof(strVar.c_str()));
    }
    else if (command=="hTres") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_hTres(atof(strVar.c_str()));
    }
    else if (command=="lTres") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_lTres(atof(strVar.c_str()));
    }
    else if (command=="n_kuhnsegs") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_n_kuhnsegs(atof(strVar.c_str()));
    }
    else if (command=="n_equ_blocks") {
        check_readFile_min_args(lineContent,2);
        double n_equ_blocks=0;
        double n_kuhnsegs=0;
        double n_relaxation=0;
        if (get_n_args(lineContent)==2) {
            strVar=get_line_comp(lineContent,1);
            n_equ_blocks=atof(strVar.c_str());
            n_kuhnsegs=sysVar.get_n_kuhnsegs();
            n_relaxation=n_equ_blocks*pow(n_kuhnsegs,2);
            sysVar.set_n_equ_blocks(n_equ_blocks);
            sysVar.set_n_relaxation(n_relaxation);
        } else {
            //strVar=get_line_comp(lineContent,1);
            //n_equ_blocks=atof(strVar.c_str());
            //n_kuhnsegs=sysVar.get_n_kuhnsegs();
            //n_relaxation=n_equ_blocks*pow(n_kuhnsegs,2);
            //sysVar.set_n_equ_blocks(n_equ_blocks);
            //sysVar.set_n_relaxation(n_relaxation);
            
            /** Set n_equ_blocks for all regimes **/
            double n_equ=0; vector<double> n_equs;
            double n_rel=0; vector<double> n_relx;
            n_kuhnsegs=sysVar.get_n_kuhnsegs();
            for (int argi=1;argi<get_n_args(lineContent);++argi) {
                strVar=get_line_comp(lineContent,argi);
                n_equ=atof(strVar.c_str());
                n_equs.push_back(n_equ);                
                n_rel=n_equ*pow(n_kuhnsegs,2);
                n_relx.push_back(n_rel);
            }
            sysVar.set_n_equ_blocks_plus(n_equs);
            sysVar.set_n_relaxation_plus(n_relx);
            
            /* error check for number of elements consistency */
            if (sysVar.get_equilibration_times().size()>0)
            {
                size_t n_equs=sysVar.get_n_equ_blocks_plus().size();
                size_t n_equ_times=sysVar.get_equilibration_times().size();
                if (n_equs!=n_equ_times) {
                    cout
                    << "Number of elements in n_equ_blocks ("<<n_equs<<") is inconsistent with "
                    << "that in logtau_targets ("<<n_equ_times<<").\n\n"; exit(EXIT_FAILURE);
                }
            }
        }
    }
    else if (command=="n_prd_blocks") {
        check_readFile_min_args(lineContent,2);
        if (get_n_args(lineContent)==2) {
            strVar=get_line_comp(lineContent,1);
            sysVar.set_n_prd_blocks(atof(strVar.c_str()));
        } else {
            strVar=get_line_comp(lineContent,1);
            sysVar.set_n_prd_blocks(atof(strVar.c_str()));
            
            /** Set n_equ_blocks for all regimes **/
            double n_prd=0;vector<double> n_prds;
            for (int argi=1;argi<get_n_args(lineContent);++argi) {
                strVar=get_line_comp(lineContent,argi);
                n_prd=atof(strVar.c_str());
                n_prds.push_back(n_prd);
            } sysVar.set_n_prd_blocks_plus(n_prds);
            
            /* error check for number of elements consistency */
            if (sysVar.get_equilibration_times().size()>0)
            {
                size_t n_prds=sysVar.get_n_prd_blocks_plus().size();
                size_t n_equ_times=sysVar.get_equilibration_times().size();
                if (n_prds!=n_equ_times) {
                    cout
                    << "Number of elements in n_prd_blocks ("<<n_prds<<") is inconsistent with "
                    << "that in logtau_targets ("<<n_equ_times<<").\n\n"; exit(EXIT_FAILURE);
                }
            }
        }
    }
    else if (command=="n_digits") {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_n_digits(atoi(strVar.c_str()));
    }
    else if (command=="logtau_targets")
    {
        double nequ=sysVar.get_n_equ_blocks();
        double nrel=sysVar.get_n_relaxation();
        if (nequ<=0) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"n_equ_blocks has to be set before logtau_targets "
            <<"with a positive value.\n\n";
            exit(EXIT_FAILURE);
        }
        
        /** Set relaxation times for all regimes **/
        double tau=0;vector<double> equ;
        if (sysVar.get_n_equ_blocks_plus().size()>0)
        {
            for (int argi=1;argi<get_n_args(lineContent);++argi) {
                tau=pow(10,atof(get_line_comp(lineContent,argi).c_str()));
                nrel=sysVar.get_n_relaxation_plus().at(argi-1);
                equ.push_back(nrel*tau);
            }
        } else {
            for (int argi=1;argi<get_n_args(lineContent);++argi) {
                tau=pow(10,atof(get_line_comp(lineContent,argi).c_str()));
                equ.push_back(nrel*tau);
            }
        } sysVar.set_equilibration_times(equ);
        
        /* error check for number of elements consistency */
        if (sysVar.get_n_equ_blocks_plus().size()>0||
            sysVar.get_n_prd_blocks_plus().size()>0)
        {
            size_t n_equs=sysVar.get_n_equ_blocks_plus().size();
            size_t n_prds=sysVar.get_n_prd_blocks_plus().size();
            size_t n_equ_times=sysVar.get_equilibration_times().size();
            if (n_equs!=n_equ_times) {
                cout
                << "Number of elements in n_equ_blocks ("<<n_equs<<") is inconsistent with "
                << "that in logtau_targets ("<<n_equ_times<<").\n\n"; exit(EXIT_FAILURE);
            }
            if (n_prds!=n_equ_times) {
                cout
                << "Number of elements in n_prd_blocks ("<<n_prds<<") is inconsistent with "
                << "that in logtau_targets ("<<n_equ_times<<").\n\n"; exit(EXIT_FAILURE);
            }
        }
    }
    else if (command=="n_regime_temps")
    {
        int nregimes=(int)sysVar.get_equilibration_times().size();
        if (nregimes<=0) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"logtau_targets has to be set before n_regime_temps.\n\n";
            exit(EXIT_FAILURE);
        }
        check_readFile_n_args(lineContent,nregimes+1);
        vector<int> ntemps;
        for (int argi=1;argi<get_n_args(lineContent);++argi) {
            ntemps.push_back(atoi(get_line_comp(lineContent,argi).c_str()));
        }
        if (ntemps.size()!=nregimes) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"size of n_regime_temps ("<<ntemps.size()<<") and "
            <<"size of logtau_targets ("<<nregimes<<") don't match.\n\n";
            exit(EXIT_FAILURE);
        } sysVar.set_n_regime_temps(ntemps);
    }
    else if (command=="ts_regime")
    {
        int nregimes=(int)sysVar.get_equilibration_times().size();
        if (nregimes<=0) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"logtau_targets has to be set before ts_regime.\n\n";
            exit(EXIT_FAILURE);
        }
        check_readFile_n_args(lineContent,nregimes+1);
        vector<double> tsregime;
        for (int argi=1;argi<get_n_args(lineContent);++argi) {
            tsregime.push_back(atof(get_line_comp(lineContent,argi).c_str()));
        }
        if (tsregime.size()!=nregimes) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"size of ts_regime("<<tsregime.size()<<") and "
            <<"size of logtau_targets("<<nregimes<<") don't match.\n\n";
            exit(EXIT_FAILURE);
        } sysVar.set_ts_regime(tsregime);
    }
    else if (command=="regime_beg")
    {
        check_readFile_n_args(lineContent,2);
        int nregimes=(int)sysVar.get_equilibration_times().size();
        int regime=atoi(get_line_comp(lineContent,1).c_str());
        if (regime<0||regime>nregimes-1) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"regime_beg should be within [0,"<<nregimes-1<<"]\n\n";
            exit(EXIT_FAILURE);
        } sysVar.set_regime_beg(regime);
    }
    else if (command=="regime_end")
    {
        check_readFile_n_args(lineContent,2);
        int nregimes=(int)sysVar.get_equilibration_times().size();
        int regime=atoi(get_line_comp(lineContent,1).c_str());
        if (regime<0||regime>nregimes-1) {
            cout
            <<"in "<<current_readFileName<<":\n"
            <<"at/around line "<<current_line_number<<":\n"
            <<"regime_end should be within [0,"<<nregimes-1<<"]\n\n";
            exit(EXIT_FAILURE);
        } sysVar.set_regime_end(regime);
    }
    else if (command=="computeCluster")
    {
        check_readFile_n_args(lineContent,2);
        strVar=get_line_comp(lineContent,1);
        sysVar.set_computeCluster(strVar);
    }
    //**************************************************************
    /** END STATIC: SYS_CONFIG **/
    
    
    /** STATIC KEYWORDS **/
    //**************************************************************
    
    //**************************************************************
    /** END STATIC KEYWORDS **/
    
    
    /** ACTION: MAKE **/
    //**************************************************************
    else if (command=="make")
    {
        UserInterface inistate=*this;
        keyword=get_line_comp(lineContent,1);
        
        if (is_action)
        {
            if (keyword=="filesystem")
            {
                check_readFile_n_args(lineContent,2);
                sysVar.set_is_makeFileSystem(1);
                sysVar.set_is_return(1);
            }
            else if (keyword=="lmpdata")
            {
                string arg=get_line_comp(lineContent,2);
                
                if (arg=="default")
                {
                    check_readFile_n_args(lineContent,3);
                    sysVar.set_is_defaultLmpData(1);
                    sysVar.set_is_return(1);
                }
                else if (arg=="moltemplate")
                {
                    check_readFile_min_args(lineContent,3);
                    sysVar.set_is_moltempLmpData(1);
                    sysVar.set_is_return(1);
                    
                    if (get_n_args(lineContent)>3)
                    {
                        int offset=2;//offset to position at "moltemplate"
                        string argkeyword;
                        for (int argi=offset+1;argi<get_n_args(lineContent);++argi)
                        {
                            argkeyword=get_line_comp(lineContent,argi);
                            if (argkeyword=="dop") {
                                int dop=atoi(get_line_comp(lineContent,argi+1).c_str());
                                moldata.set_dop(dop);
                                argi+=1;
                            }
                            if (argkeyword=="merSet") {//merSet <n> <{S(n)}>
                                int n_merSet=atoi(get_line_comp(lineContent,argi+1).c_str());
                                vector<string> tmpstr;
                                vector<int> tmpvi;
                                for (int i=0;i<n_merSet*2;i+=2) {
                                    string mer=get_line_comp(lineContent,argi+2+i);
                                    int comp=atoi(get_line_comp(lineContent,argi+2+i+1).c_str());
                                    tmpstr.push_back(mer);
                                    tmpvi.push_back(comp);
                                }
                                moldata.set_merSet(tmpstr);
                                moldata.set_merComp(tmpvi);
                                argi=argi+1+n_merSet*2;
                            }
                            if (argkeyword=="copolymerType") {
                                string str=get_line_comp(lineContent,argi+1);
                                moldata.set_copolymerType(str);
                                argi+=1;
                            }
                            if (argkeyword=="tacticity") {
                                string str=get_line_comp(lineContent,argi+1);
                                moldata.set_tacticity(str);
                                argi+=1;
                            }
                            if (argkeyword=="n_total_atoms") {
                                int tmpint=atoi(get_line_comp(lineContent,argi+1).c_str());
                                moldata.set_is_restrictbynAtoms(true);//NOTE
                                moldata.set_n_total_atoms(tmpint);
                                argi+=1;
                            }
                            if (argkeyword=="sequenceNum") {
                                int tmpint=atoi(get_line_comp(lineContent,argi+1).c_str());
                                moldata.set_is_restrictbynAtoms(false);//NOTE
                                moldata.set_n_total_atoms(0);//NOTE
                                moldata.set_sequenceNum(tmpint);
                                argi+=1;
                            }
                        }
                    }
                }
                else {
                    cout
                    <<"in "<<current_readFileName<<":\n"
                    <<"at/around line "<<current_line_number<<":\n"
                    <<"arg ("+arg+") for "+keyword+"not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="scripts")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="default") {
                    sysVar.set_is_LammpsPhase(1);
                    sysVar.set_is_return(1);
                }
                else if (arg=="prepared") {
                    sysVar.set_is_use_prepScripts(1);
                    sysVar.set_is_return(1);
                }
                else {
                    cout
                    <<"in "<<current_readFileName<<":\n"
                    <<"at/around line "<<current_line_number<<":\n"
                    <<"arg ("+arg+") for "+keyword+"not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                cout
                << "in "<<current_readFileName<<":\n"
                << "at/around line "<<current_line_number<<": "
                << "keyword ("+keyword+") not found.\n";
                exit(EXIT_FAILURE);
            }
            submit_jobs(this);
            *this=inistate;
        }
    }
    //**************************************************************
    /** END ACTION: MAKE **/
    
    
    /** ACTION: RUN **/
    //**************************************************************
    else if (command=="run")
    {
        UserInterface inistate=*this;
        check_readFile_min_args(lineContent,2);
        
        if (is_action)
        {
            int offset=0;
            vector<string> foundList;
            
            for (int argi=offset+1;argi<get_n_args(lineContent);++argi)
            {
                keyword=get_line_comp(lineContent,argi);
                foundList.push_back(keyword);
                
                if (keyword=="fitDielectric")
                {
                    if (argi!=1) {
                        cout << "fitDielectric can not be used with other modules.\n";
                        exit(EXIT_FAILURE);
                    }
                    sysVar.set_is_fit_dielectric(1);
                    
                    if (get_n_args(lineContent)>2)
                    {
                        int offset=1;//offset to position at "dielectric_analysis"
                        string argkeyword;
                        for (int argii=offset+1;argii<get_n_args(lineContent);++argii)
                        {
                            argkeyword=get_line_comp(lineContent,argii);
                            
                            if (argkeyword=="form") {
                                string str=get_line_comp(lineContent,argii+1);
                                if (str!="NPIC"&&str!="NIST"&&str!="manual") {
                                    cout
                                    << "in "<<current_readFileName<<":\n"
                                    << "at/around line "<<current_line_number<<":\n"
                                    << "form should be NPIC, NIST, or manual.\n";
                                    exit(EXIT_FAILURE);
                                } argii+=1;
                                
                                ds.set_form(str);
                                cout << "str "<<str<<"\n";
                                ds.set_sfreq_range({0,inf});
                                
                                if (str=="manual") {
                                    ds.set_is_manualFitSingleDS(1);
                                } else {
                                    ds.set_is_manualFitSingleDS(0);
                                }
                            }
                            else if (argkeyword=="sfreq_range") {
                                string lo=get_line_comp(lineContent,argii+1);
                                string hi=get_line_comp(lineContent,argii+2);
                                ds.set_sfreq_range({atof(lo.c_str()),atof(hi.c_str())});
                                argii+=2;
                            }
                            else if (argkeyword=="dc_free") {
                                string str=get_line_comp(lineContent,argii+1);
                                if (str=="yes") {
                                    ds.set_is_use_dc_free(1);
                                } else if (str=="no") {
                                    ds.set_is_use_dc_free(0);
                                } else {
                                    cout
                                    << "in "<<current_readFileName<<":\n"
                                    << "at/around line "<<current_line_number<<": "
                                    << "arg for "+argkeyword+" should be 'yes' or 'no'.\n";
                                    exit(EXIT_FAILURE);
                                } argii+=1;
                            }
                            else if (argkeyword=="func") {
                                string str=get_line_comp(lineContent,argii+1);
                                ds.set_func(str);
                                argii+=1;
                            }
                            else if (argkeyword=="model") {
                                string str=get_line_comp(lineContent,argii+1);
                                ds.set_model(str);
                                argii+=1;
                            }
                            else {
                                cout
                                << "in "<<current_readFileName<<":\n"
                                << "at/around line "<<current_line_number<<": "
                                << "arg ("+argkeyword+") for fitDielectric not found.\n";
                                exit(EXIT_FAILURE);
                            }
                        }
                    } break;//NOTE
                }
                else if (keyword=="simulation")
                {
                    sysVar.set_is_Simulations(1);
                }
                else if (keyword=="aging")
                {
                    sysVar.set_is_aging(1);
                }
                else if (keyword=="analysis")
                {
                    sysVar.set_is_AMDAT(1);
                }
                else if (keyword=="fitdata")
                {
                    sysVar.set_is_fitData(1);
                }
                else if (keyword=="cutoff_traject")
                {
                    sysVar.set_is_cutoff_traject(1);
                }
                else if (keyword=="alglibfit")
                {
                    sysVar.set_is_alglibfit(1);
                }
                else if (keyword=="theorytest")
                {
                    sysVar.set_is_theorytest(1);
                }
                else if (keyword=="sub")
                {
                    string arg=get_line_comp(lineContent,argi+1);
                    if (arg=="yes") {
                        sysVar.set_is_directSub(1);
                    } else if (arg=="no") {
                        sysVar.set_is_directSub(0);
                    } else {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                        exit(EXIT_FAILURE);
                    } argi+=1;
                }
                else if (keyword=="hold")
                {
                    string arg=get_line_comp(lineContent,argi+1);
                    if (arg=="yes") {
                        sysVar.set_is_watch_hold(1);
                    } else if (arg=="no") {
                        sysVar.set_is_watch_hold(0);
                    } else {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                        exit(EXIT_FAILURE);
                    } argi+=1;
                }
                else if (keyword=="amdatinp")
                {
                    string arg=get_line_comp(lineContent,argi+1);
                    if (arg=="yes") {
                        sysVar.set_is_amdatinp(1);
                    } else if (arg=="no") {
                        sysVar.set_is_amdatinp(0);
                    } else {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                        exit(EXIT_FAILURE);
                    } argi+=1;
                }
                else if (keyword=="sim_restart")
                {
                    int arg=atoi(get_line_comp(lineContent,argi+1).c_str());
                    if (arg<0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be >=0.\n";
                        exit(EXIT_FAILURE);
                    } else {
                        sysVar.set_sim_restart(arg);
                        cout<<"sim_restart "<<arg<<"\n";
                    } argi+=1;
                    
                    while (get_line_comp(lineContent,argi+1,false)=="-p")
                    {
                        string arg=get_line_comp(lineContent,argi+2);
                        int restart=atoi(get_line_comp(lineContent,argi+3).c_str());
                        if (arg=="gen") {
                            sysVar.set_gen_restart(restart);
                            cout<<"gen_restart "<<restart<<"\n";
                        }
                        else if (arg=="qch") {
                            sysVar.set_qch_restart(restart);
                            cout<<"qch_restart "<<restart<<"\n";
                        }
                        else if (arg=="equ") {
                            sysVar.set_equ_restart(restart);
                            cout<<"equ_restart "<<restart<<"\n";
                        }
                        else if (arg=="res") {
                            sysVar.set_res_restart(restart);
                            cout<<"res_restart "<<restart<<"\n";
                        }
                        else if (arg=="prd") {
                            sysVar.set_prd_restart(restart);
                            cout<<"prd_restart "<<restart<<"\n";
                        }
                        else {
                            cout
                            << "in "<<current_readFileName<<":\n"
                            << "at/around line "<<current_line_number<<": "
                            << "arg ("+arg+") not found.\n"
                            << "Available args are: gen,qch,equ,res,prd\n";
                            exit(EXIT_FAILURE);
                        } argi+=3;
                    }
                }
                else if (keyword=="ana_restart")
                {
                    int arg=atoi(get_line_comp(lineContent,argi+1).c_str());
                    if (arg<0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be >=0.\n";
                        exit(EXIT_FAILURE);
                    } else {
                        sysVar.set_ana_restart(arg);
                        cout<<"ana_restart "<<arg<<"\n";
                    } argi+=1;
                }
                else if (keyword=="fd_res")
                {
                    int arg=atoi(get_line_comp(lineContent,argi+1).c_str());
                    if (arg<0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be >=0.\n";
                        exit(EXIT_FAILURE);
                    } else {
                        sysVar.set_fd_res(arg);
                        cout<<"fd_res "<<arg<<"\n";
                    } argi+=1;
                }
                else if (keyword=="tt_res")
                {
                    int arg=atoi(get_line_comp(lineContent,argi+1).c_str());
                    if (arg<0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "arg for "+keyword+" should be >=0.\n";
                        exit(EXIT_FAILURE);
                    } else {
                        sysVar.set_tt_res(arg);
                        cout<<"tt_res "<<arg<<"\n";
                    } argi+=1;
                }
                else
                {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "keyword ("+keyword+") not found.\n";
                    exit(EXIT_FAILURE);
                }
            } check_run(foundList);
            cout
            << "\n"
            << "run regime "<<sysVar.get_regime_beg()<<" to "
            << "regime "<<sysVar.get_regime_end()<<"...\n";
            //system_countdown(3);
            submit_jobs(this);
            *this=inistate;
        }
    }
    //**************************************************************
    /** END ACTION: RUN **/
    
    
    /** ACTION: INTERACTION KEYWORDS **/
    //**************************************************************
    else if (command=="cp")
    {
        check_readFile_min_args(lineContent,3);
        if (is_action)
        {
            sysVar.set_is_copy_from_to(1);
            int offset=0;
            /** source includes 2nd arg to second to last arg **/
            string arg,sourceall;
            for (int argi=offset+1;argi<get_n_args(lineContent)-1;++argi) {
                arg=get_line_comp(lineContent,argi);
                sourceall+=arg;
                sourceall+=" ";
            } string source=sourceall;
            /** destination is the final arg **/
            string destin=get_line_comp(lineContent,get_n_args(lineContent)-1);
            sysVar.set_copy_source(source);
            sysVar.set_copy_destin(destin);
        }
    }
    else if (command=="define")
    {
        check_readFile_min_args(lineContent,2);
        if (is_action)
        {
            string str,var_name=get_line_comp(lineContent,1);
            vector<string> var_vals;
            for (int argi=2;argi<get_n_args(lineContent);++argi) {
                string vals=get_line_comp(lineContent,argi);
                istringstream iss(vals);
                while (iss>>str) var_vals.push_back(str);
            } var.add_new_var(var_name,var_vals);
        }
    }
    else if (command=="next")
    {
        check_readFile_n_args(lineContent,2);
        if (is_action) var.update_val_pos(get_line_comp(lineContent,1));
    }
    else if (command=="print")
    {
        check_readFile_min_args(lineContent,2);
        if (is_action) print_command();
    }
    else if (command=="show")
    {
        check_readFile_min_args(lineContent,2);
        keyword=get_line_comp(lineContent,1);
        if (is_action)
        {
            if (keyword=="all_var_names")
            {
                check_readFile_n_args(lineContent,2);
                var.show_all_var_names();
            }
            else if (keyword=="all_var_vals")
            {
                check_readFile_n_args(lineContent,2);
                var.show_all_var_vals();
            }
            else if (keyword=="var_vals")
            {
                check_readFile_n_args(lineContent,3);
                var.show_var_vals(get_line_comp(lineContent,2));
            }
            else if (keyword=="all_tinfo")
            {
                check_readFile_n_args(lineContent,2);
                sysVar.show_inptinfo();
            }
            else if (keyword=="regime_tinfo")
            {
                check_readFile_n_args(lineContent,3);
                int nregimes=(int)sysVar.get_equilibration_times().size();
                int regime=atoi(get_line_comp(lineContent,2).c_str());
                if (regime<0||regime>nregimes-1) {
                    cout
                    <<"in "<<current_readFileName<<":\n"
                    <<"at/around line "<<current_line_number<<":\n"
                    <<"regime index should be within [0,"<<nregimes-1<<"]\n\n";
                    exit(EXIT_FAILURE);
                } sysVar.show_regime_tinfo(regime);
            }
            else
            {
                cout
                << "in "<<current_readFileName<<":\n"
                << "at/around line "<<current_line_number<<": "
                << "keyword ("+keyword+") not found.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    else if (command=="set")
    {
        check_readFile_min_args(lineContent,2);
        keyword=get_line_comp(lineContent,1);
        
        if (is_action)
        {
            /** tinfo setter **/
            //------------------------------------------------------------------
            if (keyword=="tinfo")
            {
                check_readFile_min_args(lineContent,3);
                arg=get_line_comp(lineContent,2);
                if (arg=="initialize")
                {
                    sysVar.clear_inptinfo();
                    cout << "initializing inptinfo...\n";
                    system_wait(1);
                    
                    /* get number of regimes to be initialized */
                    
                    // 1. number of Ts by default set to n_regimes
                    //----------------------------------------------------------
                    int n_regimes=(int)sysVar.get_equilibration_times().size();
                    int n_temps=0,n_pre=0,n_now=0;
                    sysVar.set_max_initialized_regime(n_regimes-1);//0-based
                    
                    // 2. number of Ts by user specification
                    //----------------------------------------------------------
                    if (get_n_args(lineContent)==4) {
                        int tmp=atoi(get_line_comp(lineContent,3).c_str());
                        if (tmp>=0) {
                            n_regimes=tmp+1;//NOTE:+1 for number of regimes
                            sysVar.set_max_initialized_regime(tmp);
                        }
                    }
                    if (n_regimes<0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "initialize tinfo failed: n_regimes<0\n";
                        exit(EXIT_FAILURE);
                    } else {
                        n_temps=sysVar.get_n_regime_temps().at(0);
                        for (int i=1;i<n_regimes;++i) {
                            n_pre=sysVar.get_n_regime_temps().at(i-1);
                            n_now=sysVar.get_n_regime_temps().at(i);
                            if (n_now>n_pre) {
                                n_temps=n_now;
                            } n_pre=n_now;
                        }
                    }
                    /* set temperature vector 1d: temperatures */
                    //----------------------------------------------------------
                    vector<double> vd;
                    if (n_temps==0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "initialize tinfo failed: n_temps=0\n";
                        exit(EXIT_FAILURE);
                    } else {
                        for (int i=0;i<n_temps;++i) vd.push_back(0);
                    }
                    /* set temperature vector 2d: trials */
                    //----------------------------------------------------------
                    vector<vector<double>> vvd;
                    int n_trial=sysVar.get_n_trial();
                    if (n_trial<=0) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "initialize tinfo failed: n_trial<=0\n";
                        exit(EXIT_FAILURE);
                    } else {
                        for (int i=0;i<n_trial;++i) vvd.push_back(vd);
                    }
                    /* set temperature vector 3d: regimes */
                    //----------------------------------------------------------
                    vector<vector<vector<double>>> tmp;
                    for (int i=0;i<n_regimes;++i) {
                        tmp.push_back(vvd);
                    }
                    /* assigned 3d initialized vector to sysVar */
                    //----------------------------------------------------------
                    sysVar.set_inptinfo(tmp);
                    cout << "inptinfo successfully initialized!\n";
                    system_wait(1);
                }
                else if (arg=="regime")
                {
                    check_readFile_n_args(lineContent,4);
                    int regime=atoi(get_line_comp(lineContent,3).c_str());
                    if (regime>sysVar.get_max_initialized_regime()) {
                        cout
                        << "in "<<current_readFileName<<":\n"
                        << "at/around line "<<current_line_number<<": "
                        << "Max initializable regime is set to be regime_"
                        << sysVar.get_max_initialized_regime()<<"; "
                        << "current set tinfo regime has exceeded the limit.\n";
                        exit(EXIT_FAILURE);
                    }
                    string str;
                    vector<vector<double>> tmp;
                    for (int i=0;i<sysVar.get_n_trial();++i) {
                        ++current_line_number;
                        str=get_lineContent_at_line(current_readFileName,current_line_number);
                        str=str_from_parse_rule02(str);
                        if (str=="") {
                            cout
                            << "in "<<current_readFileName<<":\n"
                            << "at/around line "<<current_line_number<<": "
                            << "set tinfo failed: missing data!\n";
                            exit(EXIT_FAILURE);
                        }
                        vector<double> vd;
                        istringstream iss(str);
                        while (iss>>str) vd.push_back(atof(str.c_str()));
                        tmp.push_back(vd);
                    } sysVar.set_regime_tinfo(regime,tmp);
                }
                else
                {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            //------------------------------------------------------------------
            /** end tinfo setter **/
            
            
            /** sysinfo setters **/
            //------------------------------------------------------------------
            else if (keyword=="is_GPU")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    sysVar.set_is_GPU(1);
                } else if (arg=="no") {
                    sysVar.set_is_GPU(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="which_GPU")
            {
                if (sysVar.get_n_trial()!=1) {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "to use this keyword ("+keyword+"), n_trial has to be 1."
                    << "\n"; exit(EXIT_FAILURE);
                }
                check_readFile_n_args(lineContent,3);
                int arg=atoi(get_line_comp(lineContent,2).c_str());
                sysVar.set_which_GPU(arg);
            }
            else if (keyword=="is_makeNewAnanlysisFolder")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    sysVar.set_is_makeNewAnanlysisFolder(1);
                } else if (arg=="no") {
                    sysVar.set_is_makeNewAnanlysisFolder(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_backupfitdata")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    sysVar.set_is_backupfitdata(1);
                } else if (arg=="no") {
                    sysVar.set_is_backupfitdata(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            //------------------------------------------------------------------
            /** sysinfo setters **/
            
            
            /** simulation setters **/
            //------------------------------------------------------------------
            else if (keyword=="n_steps")
            {
                check_readFile_n_args(lineContent,4);
                arg=get_line_comp(lineContent,2);
                if (arg=="gen")
                {
                    ws.set_steps_gen(atof(get_line_comp(lineContent,3).c_str()));
                }
                else
                {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="quench_rate")
            {
                check_readFile_n_args(lineContent,3);
                ws.set_quenchRate(atof(get_line_comp(lineContent,2).c_str()));
            }
            else if (keyword=="compute_node")
            {
                check_readFile_min_args(lineContent,3);
                string str;
                vector<int> node0,node1,node2;
                vector<string> var_vals;
                int n_args=get_n_args(lineContent);
                for (int argi=2;argi<n_args;++argi) {
                    str=get_line_comp(lineContent,argi);
                    var_vals.push_back(str);
                    int dashind=-1;
                    string which,node;
                    for (int i=0;i<(int)str.size();++i) {
                        if (str.at(i)=='-') {
                            dashind=i; break;
                        }
                    }
                    if (dashind<0) {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"compute node should be specified by syantax: "
                        <<"(compute.id)-(node.id), ex. 0-10,1-11,2-0 \n\n";
                        exit(EXIT_FAILURE);
                    } else {
                        for(int i=0;i<dashind;++i)which+=str.at(i);
                        for(int i=dashind+1;i<(int)str.size();++i)node+=str.at(i);
                        if (atoi(which.c_str())==0) {//Node-0
                            ws.set_is_node0(1);
                            sysVar.set_is_node0_info(1);
                            node0.push_back(atoi(node.c_str()));
                        } else if (atoi(which.c_str())==1) {//Node-1
                            ws.set_is_node1(1);
                            sysVar.set_is_node1_info(1);
                            node1.push_back(atoi(node.c_str()));
                        } else if (atoi(which.c_str())==2) {//Node-2
                            ws.set_is_node2(1);
                            sysVar.set_is_node2_info(1);
                            node2.push_back(atoi(node.c_str()));
                        }
                    }
                }
                /* node2 cannot be mixed in use with node0 or node1 */
                if ((ws.get_is_node0()||ws.get_is_node1())&&(ws.get_is_node2())) {
                    cout
                    << "Node-2 can NOT be used together with Node-0 or Node-1. "
                    << "Please re-specify the compute nodes.\n\n";
                    exit(EXIT_FAILURE);
                }
                ws.set_node0(node0);
                ws.set_node1(node1);
                ws.set_node2(node2);
                sysVar.set_node0_info(node0);
                sysVar.set_node1_info(node1);
                sysVar.set_node2_info(node2);
                var.add_new_var(command,var_vals);
            }
            else if (keyword=="n_cores_run")
            {
                check_readFile_n_args(lineContent,4);
                string arg=get_line_comp(lineContent,2);
                int val=atoi(get_line_comp(lineContent,3).c_str());
                vector<int> run_cores=ws.get_run_cores();
                vector<int> numCores=ws.get_numCores();
                if (arg=="gen") {
                    numCores.at(0) =val; run_cores.at(0)=val;
                } else if (arg=="qch") {
                    numCores.at(1) =val; run_cores.at(1)=val;
                } else if (arg=="equ") {
                    numCores.at(2) =val; run_cores.at(2)=val;
                } else if (arg=="res") {
                    numCores.at(3) =val; run_cores.at(3)=val;
                } else if (arg=="prd") {
                    numCores.at(4) =val; run_cores.at(4)=val;
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg ("+arg+") not found.\n"
                    << "Available args are: gen,qch,equ,res,prd\n";
                    exit(EXIT_FAILURE);
                }
                ws.set_run_cores(run_cores);
                ws.set_numCores(numCores);
            }
            else if (keyword=="n_cores_requested")
            {
                check_readFile_n_args(lineContent,4);
                string arg=get_line_comp(lineContent,2);
                int val=atoi(get_line_comp(lineContent,3).c_str());
                vector<int> numCores=ws.get_numCores();
                if (arg=="gen") {
                    if (val>=ws.get_run_cores().at(0)) {
                        numCores.at(0)=val;
                    } else {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"n_cores_requested ("<<val<<") cannot be < "
                        <<"n_cores_run ("<<ws.get_run_cores().at(0)<<").\n";
                        exit(EXIT_FAILURE);
                    }
                } else if (arg=="qch") {
                    if (val>=ws.get_run_cores().at(1)) {
                        numCores.at(1)=val;
                    } else {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"n_cores_requested ("<<val<<") cannot be < "
                        <<"n_cores_run ("<<ws.get_run_cores().at(1)<<").\n";
                        exit(EXIT_FAILURE);
                    }
                } else if (arg=="equ") {
                    if (val>=ws.get_run_cores().at(2)) {
                        numCores.at(2)=val;
                    } else {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"n_cores_requested ("<<val<<") cannot be < "
                        <<"n_cores_run ("<<ws.get_run_cores().at(2)<<").\n";
                        exit(EXIT_FAILURE);
                    }
                } else if (arg=="res") {
                    if (val>=ws.get_run_cores().at(3)) {
                        numCores.at(3)=val;
                    } else {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"n_cores_requested ("<<val<<") cannot be < "
                        <<"n_cores_run ("<<ws.get_run_cores().at(3)<<").\n";
                        exit(EXIT_FAILURE);
                    }
                } else if (arg=="prd") {
                    if (val>=ws.get_run_cores().at(4)) {
                        numCores.at(4)=val;
                    } else {
                        cout
                        <<"in "<<current_readFileName<<":\n"
                        <<"at/around line "<<current_line_number<<":\n"
                        <<"n_cores_requested ("<<val<<") cannot be < "
                        <<"n_cores_run ("<<ws.get_run_cores().at(4)<<").\n";
                        exit(EXIT_FAILURE);
                    }
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg ("+arg+") not found.\n"
                    << "Available args are: gen,qch,equ,res,prd\n";
                    exit(EXIT_FAILURE);
                } ws.set_numCores(numCores);
            }
            //------------------------------------------------------------------
            /** end simulation setters **/
            
            
            /** analysis setters **/
            //------------------------------------------------------------------
            else if (keyword=="analysispart")
            {
                check_readFile_n_args(lineContent,3);
                aa.set_analysispart(get_line_comp(lineContent,2));
            }
            else if (keyword=="relaxation_target")
            {
                check_readFile_n_args(lineContent,3);
                aa.set_relaxation_target(get_line_comp(lineContent,2));
            }
            else if (keyword=="c_amdat")
            {
                check_readFile_n_args(lineContent,3);
                aa.set_amdat_numCores(atoi(get_line_comp(lineContent,2).c_str()));
            }
            else if (keyword=="is_mbodies")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    aa.set_is_mbodies(1);
                } else if (arg=="no") {
                    aa.set_is_mbodies(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_monoStruct")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    aa.set_is_monoStruct(1);
                } else if (arg=="no") {
                    aa.set_is_monoStruct(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_strings")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    aa.set_is_strings(1);
                } else if (arg=="no") {
                    aa.set_is_strings(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_changeAtomType")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    aa.set_is_changeAtomType(1);
                } else if (arg=="no") {
                    aa.set_is_changeAtomType(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_use_voroNeighbors")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    aa.set_is_use_voroNeighbors(1);
                } else if (arg=="no") {
                    aa.set_is_use_voroNeighbors(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="logtaubeta")
            {
                check_readFile_n_args(lineContent,3);
                aa.set_logtaubeta(atof(get_line_comp(lineContent,2).c_str()));
            }
            //------------------------------------------------------------------
            /** end analysis setters **/
            
            
            /** fd setters **/
            //------------------------------------------------------------------
            else if (keyword=="is_use_FG")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_use_FG(1);
                } else if (arg=="no") {
                    fd.set_is_use_FG(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_Fs_by_spline")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_Fs_by_spline(1);
                } else if (arg=="no") {
                    fd.set_is_fit_Fs_by_spline(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_full_alpha")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_full_alpha(1);
                } else if (arg=="no") {
                    fd.set_is_fit_full_alpha(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_use_gammafunc")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_use_gammafunc(1);
                } else if (arg=="no") {
                    fd.set_is_use_gammafunc(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_sExp")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_sExp(1);
                } else if (arg=="no") {
                    fd.set_is_fit_sExp(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_lnFs")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_lnFs(1);
                } else if (arg=="no") {
                    fd.set_is_fit_lnFs(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_1stmoment_tau")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_1stmoment_tau(1);
                } else if (arg=="no") {
                    fd.set_is_1stmoment_tau(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_Arrhenius")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_Arrhenius(1);
                } else if (arg=="no") {
                    fd.set_is_fit_Arrhenius(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_fit_tauFit")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_fit_tauFit(1);
                } else if (arg=="no") {
                    fd.set_is_fit_tauFit(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_applytauFitcut")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_applytauFitcut(1);
                } else if (arg=="no") {
                    fd.set_is_applytauFitcut(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_use_KWWassist")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_use_KWWassist(1);
                } else if (arg=="no") {
                    fd.set_is_use_KWWassist(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_find_DWF")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_find_DWF(1);
                } else if (arg=="no") {
                    fd.set_is_find_DWF(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_find_NGP")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_find_NGP(1);
                } else if (arg=="no") {
                    fd.set_is_find_NGP(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_calc_thermoData")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_calc_thermoData(1);
                } else if (arg=="no") {
                    fd.set_is_calc_thermoData(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_normalModes")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_normalModes(1);
                } else if (arg=="no") {
                    fd.set_is_normalModes(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_qspectrum")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    fd.set_is_qspectrum(1);
                } else if (arg=="no") {
                    fd.set_is_qspectrum(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="definitionofTA")
            {
                check_readFile_n_args(lineContent,3);
                fd.set_definitionofTA(get_line_comp(lineContent,2));
            }
            else if (keyword=="fcorr_model")
            {
                check_readFile_n_args(lineContent,3);
                fd.set_fcorr_model(get_line_comp(lineContent,2));
            }
            else if (keyword=="extrp_model")
            {
                check_readFile_n_args(lineContent,3);
                fd.set_extrp_model(get_line_comp(lineContent,2));
            }
            else if (keyword=="presq_model")
            {
                check_readFile_n_args(lineContent,3);
                fd.set_presq_model(get_line_comp(lineContent,2));
            }
            else if (keyword=="xtauA")
            {
                check_readFile_n_args(lineContent,3);
                tt.set_is_reltaur(1);//NOTE:set values for var in TheoryTest
                fd.set_xtauA(atof(get_line_comp(lineContent,2).c_str()));
            }
            else if (keyword=="xTA")
            {
                check_readFile_n_args(lineContent,3);
                fd.set_xTA(atof(get_line_comp(lineContent,2).c_str()));
            }
            else if (keyword=="waveindices")
            {
                check_readFile_n_args(lineContent,5);
                int beg=atoi(get_line_comp(lineContent,2).c_str());
                int end=atoi(get_line_comp(lineContent,3).c_str());
                int set=atoi(get_line_comp(lineContent,4).c_str());
                fd.set_waveindices({beg,end,set});
            }
            else if (keyword=="logtselect")
            {
                check_readFile_min_args(lineContent,3);
                vector<double> logts;
                for (int argi=2;argi<get_n_args(lineContent);++argi) {
                    logts.push_back(atof(get_line_comp(lineContent,argi).c_str()));
                } fd.set_logtselect(logts);
            }
            //------------------------------------------------------------------
            /** end fd setters **/
            
            
            /** tt setters **/
            //------------------------------------------------------------------
            else if (keyword=="is_AGtest")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    tt.set_is_AGtest(1);
                } else if (arg=="no") {
                    tt.set_is_AGtest(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_GLMtest")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    tt.set_is_GLMtest(1);
                } else if (arg=="no") {
                    tt.set_is_GLMtest(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_Leporinitest")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    tt.set_is_Leporinitest(1);
                } else if (arg=="no") {
                    tt.set_is_Leporinitest(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (keyword=="is_HWtest")
            {
                check_readFile_n_args(lineContent,3);
                string arg=get_line_comp(lineContent,2);
                if (arg=="yes") {
                    tt.set_is_HWtest(1);
                } else if (arg=="no") {
                    tt.set_is_HWtest(0);
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "arg for "+keyword+" should be 'yes' or 'no'.\n";
                    exit(EXIT_FAILURE);
                }
            }
            //------------------------------------------------------------------
            /** end tt setters **/
            
            else
            {
                cout
                << "in "<<current_readFileName<<":\n"
                << "at/around line "<<current_line_number<<": "
                << "keyword ("+keyword+") not found.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    //**************************************************************
    /** END ACTION: INTERACTION KEYWORDS **/
    
    
    /** IF STATEMENT **/
    //**************************************************************
    else if (command=="if")
    {
        check_readFile_min_args(lineContent,2);
        if (is_action)
        {
            bool is_true=false;
            int n_args=get_n_args(lineContent);
            string token=get_line_comp(lineContent,2);
            if (n_args==2)
            {
                if (token=="true"||atoi(token.c_str())>0) {
                    is_true=true;
                } else if ((token=="false"||atoi(token.c_str())==0)) {
                    is_true=false;
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "token ("+token+") not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else if (n_args==4)
            {
                double lhs=atof(get_line_comp(lineContent,1).c_str());
                double rhs=atof(get_line_comp(lineContent,3).c_str());
                cout << lhs<<" "<<token<<" "<<rhs<<"\n";
                if (token==">") {
                    if (lhs>rhs) {
                        is_true=true;
                    } else {
                        is_true=false;
                    }
                } else if (token=="<") {
                    if (lhs<rhs) {
                        is_true=true;
                    } else {
                        is_true=false;
                    }
                } else if (token=="==") {
                    if (lhs==rhs) {
                        is_true=true;
                    } else {
                        is_true=false;
                    }
                } else {
                    cout
                    << "in "<<current_readFileName<<":\n"
                    << "at/around line "<<current_line_number<<": "
                    << "token ("+token+") not found.\n";
                    exit(EXIT_FAILURE);
                }
            }
            else
            {
                cout
                << "in "<<current_readFileName<<":\n"
                << "at/around line "<<current_line_number<<": "
                << "n_args ("<<n_args<<") not correct.\n";
                exit(EXIT_FAILURE);
            }
            
            ++n_ifs;
            
            int begline=current_line_number+1;//line after "if"
            int endline=0;//line before "endif"
            find_if_content(begline,endline);
            
            if (is_true)
            {
                /** EXECUTE LOOP CONTENT **/
                execute_if_content(begline,endline);
                
                /** CONTINUE AFTER "END" LINE **/
                current_line_number=endline+2;
            }
            else
            {
                /** CONTINUE AFTER "END" LINE **/
                current_line_number=endline+2;
            }
        }
    }
    else if (command=="endif")
    {
        check_readFile_n_args(lineContent,1);
        --n_ifs;
    }
    //**************************************************************
    /** END IF STATEMENT **/
    
    
    /** FOR LOOP STRUCTURE **/
    //**************************************************************
    // NOTE: Currently, for loop does not support nested structure;
    // i.e. only one level of for loop is supported
    else if (command=="for")
    {
        check_readFile_min_args(lineContent,4);
        if (is_action)
        {
            int n_args=get_n_args(lineContent);
            int begind=0,endind=0,incind=1;
            string forVar=get_line_comp(lineContent,1);;
            begind=atoi(get_line_comp(lineContent,2).c_str());
            endind=atoi(get_line_comp(lineContent,3).c_str());
            if (n_args==4) {
                incind=1;
                if (begind>endind) {
                    cout
                    <<"in "<<current_readFileName<<":\n"
                    <<"at/around line "<<current_line_number<<":\n"
                    <<"end index should >= begin index.\n";
                    exit(EXIT_FAILURE);
                }
            } else if (n_args==5) {
                incind=atoi(get_line_comp(lineContent,4).c_str());
            }
            
            ++n_fors;
            
            vector<string> var_vals;
            for (int i=begind;i<endind;i+=incind) {
                var_vals.push_back(to_string((long long int)i));
            } var.add_new_var(forVar,var_vals);
            
            /** FIND LOOP CONTENT: IN-BETWEEN "FOR" and "END" **/
            int begline=current_line_number+1;//line after "for"
            int endline=0;//line before "end"
            find_loop_content(begline,endline);
            
            /** EXECUTE LOOP CONTENT **/
            for (int i=begind;i<endind;i+=incind) {
                execute_loop_content(begline,endline,forVar);
            }
            
            /** CONTINUE AFTER "END" LINE **/
            current_line_number=endline+2;
        }
    }
    else if (command=="endfor")
    {
        check_readFile_n_args(lineContent,1);
        --n_fors;
    }
    //**************************************************************
    /** END FOR LOOP STRUCTURE **/
    
    
    /** COMMAND NOT FOUND **/
    //**************************************************************
    else
    {
        cout
        << "in "<<current_readFileName<<":\n"
        << "at/around line "<<current_line_number<<": "
        << "command ("+command+") not found.\n";
        exit(EXIT_FAILURE);
    }
    //**************************************************************
    /** END COMMAND_PARSER() **/ return;
}





void UserInterface::find_if_content(const int begline,
                                    int& endline)
{
    n_ifs=1;
    bool def=is_action;
    is_action=false;
    int lineid=begline-1;//0-based
    do {
        lineContent=veclineContents.at(lineid);
        lineContent=get_effectiveStr(lineContent);
        if (lineContent!="") {
            command=get_line_comp(lineContent,0);
            command_parser(command);
        } ++lineid;
    } while (n_ifs>0);
    endline=lineid-1;//NOTE: the line before "end"
    is_action=def;
}





void UserInterface::execute_if_content(const int begline,
                                       const int endline)
{
    bool def=is_action;
    is_action=true;
    for (int i=begline;i<=endline;++i) {
        lineContent=veclineContents.at(i-1);
        lineContent=get_effectiveStr(lineContent);
        if (lineContent!="") {
            command=get_line_comp(lineContent,0);
            command_parser(command);
        }
    } is_action=def;
}





void UserInterface::find_loop_content(const int begline,
                                      int& endline)
{
    n_fors=1;
    bool def=is_action;
    is_action=false;
    int lineid=begline-1;//0-based
    do {
        lineContent=veclineContents.at(lineid);
        lineContent=get_effectiveStr(lineContent);
        if (lineContent!="") {
            command=get_line_comp(lineContent,0);
            command_parser(command);
        } ++lineid;
    } while (n_fors>0);
    endline=lineid-1;//NOTE: the line before "end"
    is_action=def;
}





void UserInterface::execute_loop_content(const int begline,
                                         const int endline,
                                         const string& forVar)
{
    bool def=is_action;
    is_action=true;
    for (int i=begline;i<=endline;++i) {
        lineContent=veclineContents.at(i-1);
        lineContent=get_effectiveStr(lineContent);
        if (lineContent!="") {
            command=get_line_comp(lineContent,0);
            command_parser(command);
        }
    } is_action=def; var.update_val_pos(forVar);
}





void UserInterface::check_run(const vector<string>& foundList)
{
    auto itr1=find(foundList.begin(),foundList.end(),"simulation");
    auto itr2=find(foundList.begin(),foundList.end(),"analysis");
    auto itr3=find(foundList.begin(),foundList.end(),"fitdata");
    bool is_missing_fitdata=
    itr1!=foundList.end()&&itr2!=foundList.end()&&itr3==foundList.end()&&
    sysVar.get_is_directSub()&&sysVar.get_is_watch_hold()&&
    (abs(sysVar.get_regime_end()-sysVar.get_regime_beg())>0);
    if (is_missing_fitdata) {
        cout
        << "in "<<current_readFileName<<":\n"
        << "at/around line "<<current_line_number<<": "
        << "fitdata is missing.\n";
        exit(EXIT_FAILURE);
    } else {
        return;
    }
}





string UserInterface::evaluate_var(const std::string& var_name)
{
    bool is_found=false;
    string value=var.evaluate_var(var_name,is_found);
    if (is_found) {
        return value;
    } else {
        cout
        << "in "<<current_readFileName<<":\n"
        << "at (or around) line "<<current_line_number<<": "
        << "variable ("+var_name+") not found! \n";
        exit(EXIT_FAILURE);
    }
}





string UserInterface::str_from_parse_rule00(const std::string& str)
{
    // rule:
    // will evaluate var value
    
    bool evlflag=false;
    int varbeg=0,varend=0;
    string ret,var;
    for (int indi=0;indi<str.size();++indi) {
        if (str.at(indi)=='$'&&str.at(indi+1)=='{') {
            varbeg=indi+2;
            for (int indii=varbeg;indii<str.size();++indii) {
                if (str.at(indii)=='}') {
                    varend=indii-1;
                    evlflag=!evlflag;
                    var.clear();
                    var=evaluate_var(get_partialstr(str,varbeg,varend));
                    indi=indii;
                    break;
                }
            }
        }
        if (evlflag) {
            ret+=var;
            evlflag=!evlflag;
        } else {
            ret+=str.at(indi);
        }
    } return ret;
}





string UserInterface::str_from_parse_rule01(const std::string& str)
{
    // rule:
    // strings to parse are in ""
    // will evaluate var value
    
    bool strflag=false;
    bool evlflag=false;
    int strbeg=0,strend=0;
    int varbeg=0,varend=0;
    string ret,var;
    for (int indi=0;indi<str.size();++indi) {
        if (str.at(indi)=='"') {
            strbeg=indi+1;
            for (int indii=strbeg;indii<str.size();++indii) {
                if (str.at(indii)=='"') {
                    strend=indii-1;
                    strflag=!strflag;
                    break;
                }
            } continue;
        }
        if (str.at(indi)=='$'&&str.at(indi+1)=='{') {
            varbeg=indi+2;
            for (int indii=varbeg;indii<str.size();++indii) {
                if (str.at(indii)=='}') {
                    varend=indii-1;
                    evlflag=!evlflag;
                    var.clear();
                    var=evaluate_var(get_partialstr(str,varbeg,varend));
                    indi=indii;
                    break;
                }
            }
        }
        if (indi>=strbeg&&indi<=strend)
        {
            if (strflag||evlflag) {
                if (evlflag) {
                    ret+=var;
                    evlflag=!evlflag;
                } else {
                    ret+=str.at(indi);
                }
            }
        }
    } return ret;
}





string UserInterface::str_from_parse_rule02(const std::string& str)
{
    // rule:
    // strings to be parse are in {}
    // content elements are deliminated by commas(,)
    
    bool strflag=false;
    int strbeg=0,strend=0;
    string ret,var;
    for (int indi=0;indi<str.size();++indi) {
        if (str.at(indi)=='{') {
            strbeg=indi+1;
            for (int indii=strbeg;indii<str.size();++indii) {
                if (str.at(indii)=='}') {
                    strend=indii-1;
                    strflag=!strflag;
                    break;
                }
            } continue;
        }
        if (indi>=strbeg&&indi<=strend)
        {
            if (strflag) {
                if (str.at(indi)==',') {
                    ret+=" ";
                } else {
                    ret+=str.at(indi);
                }
            }
        }
    } return ret;
}





void UserInterface::print_command()
{
    string str;
    int n_args=get_n_args(lineContent);
    for (int argi=1;argi<n_args;++argi) {
        str+=get_line_comp(lineContent,argi);
        str+=" ";
    } //cout << str<<"\n"; exit(0);
    cout << str_from_parse_rule01(str) << "\n";
}





UserInterface * UserInterface::get_inpVar()
{
    if (is_inpVar_successful) {
        return this;
    } else {
        return nullptr;
    }
}




