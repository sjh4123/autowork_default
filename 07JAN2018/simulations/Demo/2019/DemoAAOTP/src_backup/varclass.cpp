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

#include "varclass.h"

using namespace std;
using namespace autoWork;

VarClass::VarClass():
/* bool */
is_val_end(false),
/* int */
/* double */
/* string */
strVar(""),
var_name(""),
/* STL */
val_pos(),
vecStrVar(),
var_names(),
var_vals()
{
    //
}



void VarClass::add_new_var(const string& var_name,
                           const vector<std::string>& var_vals)
{
    // var_names, var_vals, and val_pos eash is a container itself but together
    // form a peudo data container, where elements from these three containers
    // have 1-to-1 mapping relationship;
    // NOTE: var_name can't be repeated; if a var_name already exists,
    // then var_namw & var_vals will be replaced, and val_pos be set to zero.
    
    //cout << var_name;
    auto itr=find(var_names.begin(),var_names.end(),var_name);
    if (itr==var_names.end()) {
        //cout << ": NEW var \n";
        add_a_new_var_name(var_name);
        add_a_new_var_vals(var_vals);
        add_a_new_val_pos(0);
    } else {
        //cout << ": OLD VAR \n";
        string varname=*itr;
        int varid=get_var_name_index(varname);
        vector<vector<string>> tmpvals=get_var_vals();
        vector<vector<string>> oldvals=tmpvals;
        tmpvals.at(varid)=var_vals;
        set_var_vals(tmpvals);
        vector<int> tmppos=get_val_pos();
        vector<int> oldpos=tmppos;
        tmppos.at(varid)=0;
        set_val_pos(tmppos);
        if (true) {
            cout << "variable ("+varname+") value has been changed from (";
            for (size_t i=0;i<oldvals.at(varid).size();++i) {
                cout << oldvals.at(varid).at(i);
                if (i!=oldvals.at(varid).size()-1) cout<<" ";
            } cout << ") to (";
            for (size_t i=0;i<var_vals.size();++i) {
                cout << var_vals.at(i);
                if (i!=var_vals.size()-1) cout<<" ";
            } cout << ")\n";
        }
    }
}



void VarClass::add_a_new_var_name(const string& str)
{
    var_names.push_back(str);
}



void VarClass::add_a_new_var_vals(const vector<string>& vals)
{
    var_vals.push_back(vals);
}



void VarClass::add_a_new_val_pos(const int new_pos)
{
    val_pos.push_back(new_pos);
}



void VarClass::update_val_pos(const string& var_name)
{
    int varnameid=get_var_name_index(var_name);
    vector<int> tmp=val_pos;
    tmp.at(varnameid)+=1;
    int valpos=tmp.at(varnameid);
    if (valpos>(var_vals.at(varnameid).size()-1)) {
        is_val_end=true;
    } else {
        is_val_end=false;
    }
    if (!is_val_end) {
        val_pos=tmp;
    } else {
        tmp.at(varnameid)=(int)var_vals.at(varnameid).size()-1;
        val_pos=tmp;
    }
}

void VarClass::show_all_var_names()
{
    for (int i=0;i<(int)var_names.size();++i) {
        cout << var_names.at(i) << "\n";
    }
}



void VarClass::show_all_var_vals()
{
    for (size_t i=0;i<var_vals.size();++i) {
        for (size_t ii=0;ii<var_vals.at(i).size();++ii) {
            cout << var_vals.at(i).at(ii) << " ";
        } cout << "\n";
    }
}



void VarClass::show_var_name_by_index(const int ind)
{
    cout << var_names.at(ind) << "\n";
}


void VarClass::show_var_vals_by_index(const int ind)
{
    for (size_t i=0;i<var_vals.at(ind).size();++i) {
        cout << var_vals.at(ind).at(i) << " ";
    } cout << "\n";
}



void VarClass::show_var_vals(const std::string& var_name)
{
    cout << var_name << " ";
    int varind=get_var_name_index(var_name);
    for (int ii=0;ii<(int)var_vals.at(varind).size();++ii) {
        cout << var_vals.at(varind).at(ii) << " ";
    } cout << "\n";
}



int VarClass::get_var_name_index(const string& name)
{
    int varind=-1;
    for (int i=0;i<(int)var_names.size();++i) {
        if (var_names.at(i)==name) {
            varind=i; break;
        }
    }
    if (varind>=0) {
        return varind;
    } else {
        cout
        << "in VarClass::get_var_name_index():\n"
        << "variable ("+name+") not found.\n";
        exit(EXIT_FAILURE);
    }
}



string VarClass::evaluate_var(const std::string& var_name,bool& is_found)
{
    string val;
    int varind=0,valpos=0;
    auto itr=find(var_names.begin(),var_names.end(),var_name);
    if (itr==get_var_names().end()) {
        is_found=false;
    } else {
        is_found=true;
        varind=get_var_name_index(*itr);
        valpos=val_pos.at(varind);
        val=var_vals.at(varind).at(valpos);
    } return val;
}



/** public setters **/
void VarClass::set_strVar(const string& str){strVar=str;}
void VarClass::set_var_name(const string& str){var_name=str;}
void VarClass::set_vecStrVar(const vector<string>& vs){vecStrVar=vs;}
void VarClass::set_var_names(const vector<string>& vs){var_names=vs;}
void VarClass::set_val_pos(const std::vector<int>& vi){val_pos=vi;}
void VarClass::set_var_vals(const vector<vector<string>>& vvs){var_vals=vvs;}




