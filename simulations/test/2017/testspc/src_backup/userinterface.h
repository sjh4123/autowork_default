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

#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "varclass.h"
#include "structureclass.h"
#include "moltemplatelmpdata.h"

#include "fitdielectric.h"
#include "workscripts.h"
#include "amdatanalysis.h"
#include "fitdata.h"

//#include "lmpscripts.h"
#include "theorytest.h"

#define VERSION "24AUG2017"

namespace autoWork
{
    class UserInterface
    {
        VarClass var;
        StructureClass sysVar;
        MoltemplateLmpData moldata;
        FitDielectric ds;
        WorkScripts ws;
        AmdatAnalysis aa;
        FitData fd;
        TheoryTest tt;
        
        bool is_inpVar_successful;
        bool is_inforloop;
        bool is_action;
        
        int current_line_number;
        int indi,indii;
        int argi,argii;//NOTE:argi is reserved for use in parse_commandLineArgs()
        int n_fors,n_ifs;
        
        std::string current_readFileName;
        std::string lineContent;
        std::string command;
        std::vector<std::string> strArgs;
        std::vector<std::string> veclineContents;
        
        void show_help_info();
        void check_var_format(const std::string&);
        double check_n_args_comkeyword();
        void show_inplineContent();
        void parse_inputFile(const std::string&);
        void command_parser(const std::string&);
        void check_readFile_n_args(const std::string&,const int);
        void check_readFile_min_args(const std::string&,const int);
        void print_command();
        void find_if_content(const int,int&);
        void find_loop_content(const int,int&);
        void execute_if_content(const int,const int);
        void execute_loop_content(const int,const int,const std::string&);
        void check_run(const std::vector<std::string>&);
        int  get_n_args(const std::string&);
        
        std::string evaluate_var(const std::string&);
        std::string str_from_parse_rule00(const std::string&);
        std::string str_from_parse_rule01(const std::string&);
        std::string str_from_parse_rule02(const std::string&);
        
        std::string get_lineContent_at_line(const std::string&,const int);
        std::string get_partialstr(const std::string&,const int,const int);
        std::string get_line_comp(const std::string&,const int,const bool check=true);
        std::string get_effectiveStr(const std::string&);
        
    public:
        
        UserInterface()=delete;
        UserInterface(const std::vector<std::string>&);
        UserInterface(const UserInterface&)=default;
        UserInterface& operator= (const UserInterface&)=default;
        ~UserInterface()=default;
        
        void show_version();
        void parse_commandLineArgs();
        
        const StructureClass& get_sysVar() const {return sysVar;}
        const MoltemplateLmpData& get_moldata() const {return moldata;}
        const FitDielectric& get_ds() const {return ds;}
        const WorkScripts& get_ws() const {return ws;}
        const AmdatAnalysis& get_aa() const {return aa;}
        const FitData& get_fd() const {return fd;}
        const TheoryTest& get_tt() const {return tt;}
        
        UserInterface * get_inpVar();
    };
}
#endif /* USERINTERFACE_H */




