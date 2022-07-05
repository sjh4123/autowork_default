//
//  Created by Sean Jh H. on 6/3/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "functions.h"

/*******************************************************************************
 NOTE:
 1) If you are using GNU g++ compiler make an alias for g++ like:
    g++ -std=c++0x -std=gnu++0x
 2) The packmol executable should be first generated on your machine;
    Also, make sure the packmol exe has sufficient permission (ex. chmod to 755)
 *******************************************************************************/

using namespace std;

int main(int argc, const char **argv)
{
    //cout << argc << "\n";
    
    /* command line arguments */
    vector<string> variable_names;
    vector<string> variable_content;
    int n_variables=0;
    
    for(int argii=0; argii<argc; ++argii)
    {
        if(argv[argii][0] == '-') // flags
        {
            if(argv[argii][1] == 'v') // variable flag
            {
                variable_names.push_back(argv[argii+1]);
                variable_content.push_back(argv[argii+2]);
                ++n_variables;
            }
        }
    }
    
    const string forMasterWatch = autoWork::submit_jobs();
    //cout << forMasterWatch << "\n";
    
    return 0;
}