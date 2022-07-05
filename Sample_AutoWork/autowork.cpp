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

#include "userinterface.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

int main(int argc, char** argv)
{
    //cout << "argv[0] = " << argv[0] << "\n";
    //cout << "argc =    " << argc << "\n";
    
    vector<string> strArgs;
    for (int i=0;i<argc;++i) strArgs.push_back(argv[i]);
    if (false) {
        cout << strArgs.size() << "\n";
        for (int i=0;i<(int)strArgs.size();++i) {
            cout << strArgs.at(i) << " ";
        } cout << "\n";
    }
    
    /** if use redirection to a file with specific file name **/
    //ifstream in("in.txt");
    //auto cinbuf_default = std::cin.rdbuf(in.rdbuf());
    //streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    //cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!
    //ofstream out("out.txt");
    //auto coutbuf_default = std::cout.rdbuf(out.rdbuf());
    //streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    //std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!    
    
    static UserInterface ui(strArgs);
    ui.show_version();
    ui.parse_commandLineArgs();
    
    //std::cin.rdbuf(cinbuf_default);   //reset to standard input again
    //std::cout.rdbuf(coutbuf_default); //reset to standard output again
    
    string forMasterWatch;
    if (ui.get_inpVar()==nullptr) {
        forMasterWatch=submit_jobs();
    } //cout << forMasterWatch << "\n\n";
    system("date");
    cout<<"AUTOWORK finished successfully.\n\n";
    return 0;
}




