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

#ifndef VARCLASS_H
#define VARCLASS_H

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

namespace autoWork
{
    class VarClass
    {
        bool is_val_end;
        std::string strVar;
        std::string var_name;
        std::vector<int> val_pos;
        std::vector<std::string> vecStrVar;
        std::vector<std::string> var_names;
        std::vector<std::vector<std::string>> var_vals;
        
    public:
        
        VarClass();
        VarClass(const VarClass&)=default;
        VarClass& operator= (const VarClass&)=default;
        ~VarClass()=default;
        
        void add_new_var(const std::string&,const std::vector<std::string>&);
        void add_a_new_var_name(const std::string&);
        void add_a_new_var_vals(const std::vector<std::string>&);
        void add_a_new_val_pos(const int);
        void update_val_pos(const std::string&);
        void show_all_var_names();
        void show_all_var_vals();
        void show_var_name_by_index(const int);
        void show_var_vals_by_index(const int);
        void show_var_vals(const std::string&);
        int  get_var_name_index(const std::string&);
        std::string evaluate_var(const std::string&,bool&);
        
        /** public setters **/
        void set_strVar(const std::string&);
        void set_var_name(const std::string&);
        void set_val_pos(const std::vector<int>&);
        void set_vecStrVar(const std::vector<std::string>&);
        void set_var_names(const std::vector<std::string>&);
        void set_var_vals(const std::vector<std::vector<std::string>>&);
        
        /** public getters **/
        const std::string get_strVar() const {return strVar;}
        const std::string get_var_name() const {return var_name;}
        /* STL */
        const std::vector<int>& get_val_pos() const {return val_pos;}
        const std::vector<std::string>& get_vecStrVar() const
        {return vecStrVar;}
        const std::vector<std::string>& get_var_names() const
        {return var_names;}
        const std::vector<std::vector<std::string>>& get_var_vals() const
        {return var_vals;}
    };
}
#endif /* VARCLASS_H */




