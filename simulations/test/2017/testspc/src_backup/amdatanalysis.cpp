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

#include "amdatanalysis.h"
#include "functions.h"

using namespace std;
using namespace autoWork;
using namespace alglib;

AmdatAnalysis::AmdatAnalysis(const StructureClass& sysVar,
                             const WorkScripts& ws):
/* bool */
is_check_amdat_version(false),//AMDAT version consistency check
is_changeAtomOrder(false),
is_changeAtomType(false),
is_recoverAtomType(false),
is_changeFrameSteps(false),
is_deleteZerothFrames(false),
is_NPT(true),
is_fixed_q(true),
is_fixed_maxL(true),
is_auto_qvalue(false),
is_mbodies(false),
is_keep_new_trajectory(false),
is_use_displacement_list(true),
is_use_find_fast(false),
is_write_pdb_file(false),
is_write_list_trajectory(false),//NOTE
is_use_strFac_frame(false),
is_rouseModes(false),
is_use_voroNeighbors(false),
is_stringmb(true),
is_streamlined_strings(true),
is_monoStruct(true),//NOTE
is_mb_all_molecule(true),
is_binning(false),
/** typical AMDAT analyses **/
is_strFac(false),
is_rdf(false),
is_msd(false),
is_ngp(false),
is_isfs(false),
/** more advanced AMDAT analyses **/
is_baf(false),
is_composition(false),
is_u2dist(false),
is_stiffness_dist(false),
is_isf(false),
is_strings(false),
is_peak_frame(true),
/* int */
indexi(0),
indexii(0),
amdat_numCores(1),
amdat_run_cores(1),
amdat_priority(0),
fullblock(0),
rdf_nbins(100),
vhs_nbins(75),
u2dist_nbins(100),
stiffness_nbins(100),
prd_blocksize(ws.get_prd_blocksize()),
waveindex(19),
n_poly(sysVar.get_n_poly()),
chainLen(sysVar.get_chainLen()),
n_moments(3),
strFac_frame(-2),//n in arg means (n-2)th frame
n_species_types(1),
n_rouseModes(1),
/* double */
prd_exp_base(ws.get_prd_exp_base()),
n_prd_blocks(ws.get_n_prd_blocks()),
DWF_time(0),
maxLenScale(55),
vhs_rangebinned(5.0),
max_u2(1.5),
max_stiffness(1.0),
threshold_percentile_lo(45.0),
threshold_percentile_hi(55.0),
strings_threshold(0.7),//NOTE: this is the a of a*sigma(TA)=threshold
threshold_percentile(93.5),//NOTE
logtaubeta(0),
/* string */
amdat_svn("137"),//NOTE
amdat_exe("/home/"+sysVar.get_userName()+"/AMDAT"+amdat_svn+"/AMDAT"),
threshold_keyword("greater"),
symmetry("symmetric"),
geometry("xyz"),
relaxation_target("isfs"),
centertype("centroid"),
species("KGPolymer"),
segmode("mode12"),
analysispart("all"),
analysispart_bin("0_0_0"),
//KGPolymer,KG12Rouse,binaryLJ,Cu4Ag6,OTP1mb,OTP3mbs,UAPS30,MartiniPS30,SiO2
/* STL */
/** speciesName-speciesType now only supports 1-to-1 species-type relation **/
speciesName({"all"}),//"all" and other species if any ("all" should be the last)
speciesType({"1"}),//types used for sigma matrix (same order as in speciesName)
n_typeSet(sysVar.get_n_typeSet()),
mbodyList({"mbList1"}),//"mbList1","mbList2",...
mbodyType({"mbType1"}),//"mbType1","mbType2",...
sigmaMatrix({{1}})
{
    if (sysVar.get_is_AMDAT()) {
        cout << "\nAMDAT in use is "<<amdat_exe<<"\n";
    }
    if (sysVar.get_systemUnit()=="real") {
        logtaubeta  =  3;
        waveindex   = 19;
        maxLenScale = 55;
    } else if (sysVar.get_systemUnit()=="lj") {
        logtaubeta  =  0;
        waveindex   = 26;
        maxLenScale = 12.43848;
    } else if (sysVar.get_systemUnit()=="metal") {
        logtaubeta  =  0;
        waveindex   = 20;
        maxLenScale = 25;
    } DWF_time=pow(10,logtaubeta);
}





void AmdatAnalysis::node_allocation(const StructureClass& sysVar,
                                    const WorkScripts& ws,
                                    const int counter,
                                    std::ofstream& writefile,
                                    const int n_trl)
{
    int cores_d=amdat_numCores;
    
    /* at least one node needs to be used */
    if (!(ws.get_is_node0()||ws.get_is_node1()||ws.get_is_node2())) {
        cout
        << "\n"
        << "All compute nodes are switched off! Please switch on at least one."
        << "\n\n";
        exit(EXIT_FAILURE);
    }
    /* node2 can NOT be mixed in use with node0 or node1 */
    if ((ws.get_is_node0()||ws.get_is_node1())&&(ws.get_is_node2())) {
        cout
        << "Node-2 can NOT be used together with Node-0 or Node-1. "
        << "Please re-specify the compute nodes.\n\n";
        exit(EXIT_FAILURE);
    }
    
    int numNode0=0;//number of nodes used on node0
    int numNode1=0;//number of nodes used on node1
    int numNode2=0;//number of nodes used on node2
    
    if (ws.get_is_node0())numNode0=(int)(ws.get_node0().size());
    if (ws.get_is_node1())numNode1=(int)(ws.get_node1().size());
    if (ws.get_is_node2())numNode2=(int)(ws.get_node2().size());
    
    // sumNodes is the total number of nodes available for job allocation
    int sumNode=numNode0+numNode1+numNode2;
    int noml=counter/(ws.get_n_cores_node()/cores_d);
    int tmpNode=noml%sumNode;
    int tmpNode1=0;
    
    /** Node-0 and Node-1 assignment **/
    if (ws.get_is_node0()&&ws.get_is_node1()) { /* using both node0 and node1 */
        if (tmpNode<numNode0) {
            writefile << "#$ -q all.q@compute-0-" << ws.get_node0().at(tmpNode) << ".local\n";
        } else {
            tmpNode1=(tmpNode-numNode0);
            writefile << "#$ -q all.q@compute-1-" << ws.get_node1().at(tmpNode1) << ".local\n";
        }
    } else if (ws.get_is_node0()) { /* only using node0 */
        writefile << "#$ -q all.q@compute-0-" << ws.get_node0().at(tmpNode) << ".local\n";
    } else if (ws.get_is_node1()) { /* only using node1 */
        writefile << "#$ -q all.q@compute-1-" << ws.get_node1().at(tmpNode) << ".local\n";
    }
    /** Node-2 assignment **/
    if (ws.get_is_node2()) { /* only using node2 */
        writefile << "#$ -q all.q@compute-2-" << ws.get_node2().at(tmpNode) << ".local\n";
    } writefile << "\n";
    
}





const string AmdatAnalysis::get_filenameString()
{
    return "${usic}_${trial}_${namestring}_T${Temp}";
}





void AmdatAnalysis::systemBlock_definition(ofstream& ai)
{
    if (species=="KGPolymer") {
        
        n_species_types=1;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "polymer 400"                                  << "\n"
        << "1"                                            << "\n"
        << "20"                                           << "\n\n";
        ai
        << "create_list all" << "\n"
        << "all"             << "\n\n";
        
    }
    else if (species=="KG12Rouse") {
        
        is_rouseModes=true;
        is_mb_all_molecule=false;
        n_species_types=1;
        
        int chainlength=12;
        
        //vector<int> modes={2};
        vector<int> modes={96,48,24,12,6,4,3,2,1};
        n_rouseModes=(int)modes.size();
        mbodyList.clear();
        for (int i=0;i<modes.size();++i) {
            mbodyList.push_back("mode"+to_string((long long int)modes.at(i)));
        }
        
        if (false) //Overlapping Rouse Modes (More States)
        {
            bool is_endtoend=true;
            
            ai
            << "system_np"                                    << "\n"
            << "custom"                                       << "\n"
            << "${path_prdcustom}"                            << "\n"
            << "exponential ${n_prd_blocks} ${blocksize} "
            << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
            << "polymer 800"                                  << "\n"
            << "1"                                            << "\n"
            << "12"                                           << "\n\n";
            ai
            << "create_list all" << "\n"
            << "all"             << "\n\n";
            
            ai << "## D.O.P:       "<<chainlength<<"\n";
            ai << "## Rouse Modes: ";
            for (size_t i=0;i<modes.size();++i) {
                ai << modes.at(i);
                if (i!=modes.size()-1) ai << ",";
            } ai << "\n\n";
            
            
            if (is_endtoend)
            {
                for (size_t i=0;i<modes.size();++i)
                {
                    int intralength=chainlength/modes.at(i);
                    int n_mbodies=chainlength-intralength+1;
                    int startind=0;
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai
                        << "create_multibodies "
                        << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" "
                        << "mbType1"<<" "
                        << centertype << " "
                        << "species_atomlist polymer ";
                        for (int iii=0;iii<intralength;++iii) {
                            ai << "1 "<<startind+iii<<" ";
                        } ++startind; ai << "\n";
                    } ai << "\n";
                    
                    ai << "combine_multibody_lists mode"<<modes.at(i)<<" ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n";
                    ai << "combine_trajectories mode"<<modes.at(i)<<" ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n\n";
                } ai << "\n\n\n\n";
                
                for (size_t i=0;i<modes.size();++i)
                {
                    if (modes.at(i)==chainlength) continue;
                    
                    int intralength=chainlength/modes.at(i);
                    int n_mbodies=chainlength-intralength+1;
                    int startind=0;
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai
                        << "create_multibodies "
                        << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" "
                        << "mbType1"<<" "
                        << centertype << " "
                        << "species_atomlist polymer ";
                        ai
                        << "1 "<<startind<<" "
                        << "1 "<<startind+intralength-1<<" ";
                        ++startind; ai << "\n";
                    } ai << "\n";
                    
                    ai << "combine_multibody_lists mode"<<modes.at(i)<<"_endtoend ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n";
                    ai << "combine_trajectories mode"<<modes.at(i)<<"_endtoend ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n\n";
                } ai << "\n\n\n\n";
            }
            else
            {
                for (size_t i=0;i<modes.size();++i)
                {
                    int intralength=chainlength/modes.at(i);
                    int n_mbodies=chainlength-intralength+1;
                    int startind=0;
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai
                        << "create_multibodies "
                        << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" "
                        << "mbType1"<<" "
                        << centertype << " "
                        << "species_atomlist polymer ";
                        for (int iii=0;iii<intralength;++iii) {
                            ai << "1 "<<startind+iii<<" ";
                        } ++startind; ai << "\n";
                    } ai << "\n";
                    
                    ai << "combine_multibody_lists mode"<<modes.at(i)<<" ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n";
                    ai << "combine_trajectories mode"<<modes.at(i)<<" ";
                    for (int ii=1;ii<=n_mbodies;++ii) {
                        ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                    } ai << "\n\n";
                } ai << "\n\n\n\n";
            }
        }
        else if (false) //Non-Overlapping Rouse Modes (Less States)
        {
            ai
            << "system_np"                                    << "\n"
            << "custom"                                       << "\n"
            << "${path_prdcustom}"                            << "\n"
            << "exponential ${n_prd_blocks} ${blocksize} "
            << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
            << "polymer 800"                                  << "\n"
            << "1"                                            << "\n"
            << "12"                                           << "\n\n";
            ai
            << "create_list all" << "\n"
            << "all"             << "\n\n";
            
            ai << "## D.O.P:       "<<chainlength<<"\n";
            ai << "## Rouse Modes: ";
            for (size_t i=0;i<modes.size();++i) {
                ai << modes.at(i);
                if (i!=modes.size()-1) ai << ",";
            } ai << "\n\n";
            
            for (size_t i=0;i<modes.size();++i)
            {
                int intralength=chainlength/modes.at(i);
                int n_mbodies=modes.at(i);
                int beadind=0;
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai
                    << "create_multibodies "
                    << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" "
                    << "mbType1"<<" "
                    << centertype << " "
                    << "species_atomlist polymer ";
                    for (int iii=0;iii<intralength;++iii) {
                        ai << "1 "<<beadind<<" ";++beadind;
                    } ai << "\n";
                } ai << "\n";
                
                ai << "combine_multibody_lists mode"<<modes.at(i)<<" ";
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                } ai << "\n";
                ai << "combine_trajectories mode"<<modes.at(i)<<" ";
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                } ai << "\n\n";
            } ai << "\n\n\n\n";
        }
        else if (false) //random sampled groups of mers in the system
        {
            ai
            << "system_np"                                    << "\n"
            << "custom"                                       << "\n"
            << "${path_prdcustom}"                            << "\n"
            << "exponential ${n_prd_blocks} ${blocksize} "
            << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
            << "polymer 1"                                    << "\n"
            << "1"                                            << "\n"
            << "9600"                                         << "\n\n";
            ai
            << "create_list all" << "\n"
            << "all"             << "\n\n";
            
            ai << "## n: ";
            for (size_t i=0;i<modes.size();++i) {
                ai << modes.at(i);
                if (i!=modes.size()-1) ai << ",";
            } ai << "\n\n";
            
            for (size_t i=0;i<modes.size();++i)
            {
                int n_mers=9600;//number of beads in the system
                int intralength=modes.at(i);//for group size (n)
                int n_mbodies=100;//samping size
                //int n_mbodies=n_mers/intralength;//samping size
                
                int mer=0;
                vector<int> group;
                vector<vector<int>> group_all;
                vector<int>::iterator itr;
                
                std::random_device rd;
                std::mt19937 mt(rd());
                std::uniform_int_distribution<int> dist(0,n_mers-1);
                
                /** NOTE: no repeat of mer index within groups of mers **/
                for (int ii=0;ii<n_mbodies;++ii) {
                    group.clear();
                    group.push_back(dist(mt));
                    if (intralength>1) {
                        for (int iii=1;iii<intralength;++iii) {
                            do {
                                mer=dist(mt);
                                itr=find(group.begin(),group.end(),mer);
                            } while(itr!=group.end());
                            group.push_back(mer);
                        }
                        /** sort group (mer index in ascending order) **/
                        sort(group.begin(),group.end());
                    } group_all.push_back(group);
                }
                /** NOTE: mers cannot be fully identical between groups **/
                for (int i=0;i<group_all.size()-1;++i) {
                    for (int j=i+1;j<group_all.size();++j) {
                        int count1=0;
                        for (int k=0;k<group_all.at(i).size();++k) {
                            if (group_all.at(i).at(k)==group_all.at(j).at(k)) {
                                ++count1;
                            }
                        }
                        //if found repeat at jth group
                        if (count1==(int)group_all.at(i).size()) {
                            bool is_retry=false;
                            do {
                                is_retry=false;
                                group.clear();
                                group.push_back(dist(mt));
                                if (intralength>1) {
                                    for (int k=1;k<intralength;++k) {
                                        do {
                                            mer=dist(mt);
                                            itr=find(group.begin(),group.end(),mer);
                                        } while(itr!=group.end());
                                        group.push_back(mer);
                                    } sort(group.begin(),group.end());
                                }
                                for (int s=0;s<=i;++s) {//NOTE: [0,i]
                                    int count2=0;
                                    for (int l=0;l<group.size();++l) {
                                        if (group.at(l)==group_all.at(s).at(l)) {
                                            ++count2;
                                        }
                                    }
                                    if (count2==(int)group.size()) {
                                        is_retry=true;
                                        break;
                                    }
                                }
                            } while (is_retry);
                            group_all.at(j)=group;
                        }
                    }
                } sort(group_all.begin(),group_all.end());//sort ascending
                
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai
                    << "create_multibodies "
                    << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" "
                    << "mbType1"<<" "
                    << centertype << " "
                    << "species_atomlist polymer ";
                    for (int iii=0;iii<intralength;++iii) {
                        ai << "1 "<<group_all.at(ii-1).at(iii)<<" ";
                    } ai << "\n";
                } ai << "\n";
                
                ai << "combine_multibody_lists mode"<<modes.at(i)<<" ";
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                } ai << "\n";
                ai << "combine_trajectories mode"<<modes.at(i)<<" ";
                for (int ii=1;ii<=n_mbodies;++ii) {
                    ai << "mode"<<modes.at(i)<<"_mb"<<ii<<"_list"<<" ";
                } ai << "\n\n";
            } ai << "\n\n\n\n";
        }
        else
        {
            cout
            << "At AmdatAnalysis::systemBlock_definition():\n"
            << "species="<<species<<", no def!\n";exit(EXIT_FAILURE);
        }
    }
    else if (species=="binaryLJ") {
        
        n_species_types=2;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "type1 6400 type2 1600"                        << "\n"
        << "1 2"                                          << "\n"
        << "1 0"                                          << "\n"
        << "0 1"                                          << "\n\n";
        ai
        << "create_list all"        << "\n"
        << "all"                    << "\n\n";
        if (speciesName.size()==3) {
            for (indexi=0;indexi<speciesName.size()-1;++indexi) {
                ai
                << "create_list  "
                << speciesName.at(indexi)+"_list"<< "\n"
                << "type_species "
                << speciesName.at(indexi)<<" "
                << speciesType.at(indexi)<<"\n";
            }
        } else {
            cout << species << " should have 3 species!\n";
            exit(EXIT_FAILURE);
        } ai << "\n";
        
    }
    else if (species=="Cu4Ag6") {
        
        n_species_types=2;
        
        ai
        << "system_nv"                                    << "\n"// NOTE: system_nv
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "Cu 1600 Ag 2400"                              << "\n"
        << "1 2"                                          << "\n"
        << "1 0"                                          << "\n"
        << "0 1"                                          << "\n\n";
        ai
        << "create_list all"      << "\n"
        << "all"                  << "\n\n";
        if (speciesName.size()==3) {
            for (indexi=0;indexi<speciesName.size()-1;++indexi) {
                ai
                << "create_list  "
                << speciesName.at(indexi)+"_list"<< "\n"
                << "type_species "
                << speciesName.at(indexi)<<" "
                << speciesType.at(indexi)<<"\n";
            }
        } else {
            cout << species << " should have 3 species!\n";
            exit(EXIT_FAILURE);
        } ai << "\n";
    }
    else if (species=="OTP1mb") {
        
        is_mb_all_molecule=true;
        n_species_types=1;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "OTP 512"                                      << "\n"
        << "1 2"                                          << "\n"
        << "18 14"                                        << "\n\n";
        ai
        << "create_list all" << "\n"
        << "all"             << "\n\n";
        
    }
    else if (species=="OTP3mbs") {
        
        is_mb_all_molecule=false;
        n_species_types=1;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "OTP 512"                                      << "\n"
        << "1 2"                                          << "\n"
        << "18 14"                                        << "\n\n";
        ai
        << "create_list all" << "\n"
        << "all"             << "\n\n";
        
        int current_C=0;
        int current_H=0;
        //int n_total_C=18;
        
        /* 1st ring (6C+5H) */
        ai
        << "create_multibodies "
        << "Ring1 mbType1 "<<centertype<<" "
        << "species_atomlist OTP ";
        for (int i=0;i<6;++i) {
            ai << "1 "<< i <<" ";
        } current_C += 6;
        for (int i=0;i<5;++i) {
            ai << "2 "<< i <<" ";
        } current_H += 5;
        /* 2nd ring (6C+4H) */
        ai
        << "\n"
        << "create_multibodies "
        << "Ring2 mbType1 "<<centertype<<" "
        << "species_atomlist OTP ";
        for (int i=0;i<6;++i) {
            ai << "1 "<< i+current_C <<" ";
        } current_C += 6;
        for (int i=0;i<4;++i) {
            ai << "2 "<<i+current_H<<" ";
        } current_H += 4;
        /* 3rd ring (6C+5H) */
        ai
        << "\n"
        << "create_multibodies "
        << "Ring3 mbType1 "<<centertype<<" "
        << "species_atomlist OTP ";
        for (int i=0;i<6;++i) {
            ai << "1 "<< i+current_C <<" ";
        }
        for (int i=0;i<5;++i) {
            ai << "2 "<<i+current_H<<" ";
        } ai << "\n\n";
        
        ai << "combine_multibody_lists mbList1 Ring1 Ring2 Ring3 \n";
        ai << "combine_trajectories mbList1 Ring1 Ring2 Ring3 \n\n";
    }
    else if (species=="UAPS30") {
        
        is_mb_all_molecule=false;
        n_species_types=1;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "UAPS30 24"                                    << "\n"
        << "1 2"                                          << "\n"
        << "60 180"                                       << "\n\n";
        ai
        << "create_list all" << "\n"
        << "all"             << "\n\n";
        
        int dop=30;
        int current_hC=0;
        int current_bC=0;
        
        for (int i=1;i<=dop;++i) {
            ai
            << "create_multibodies "
            << "mer"<<i<<" "<<"mbType1"<<" "<<centertype<<" "
            << "species_atomlist UAPS30 ";
            /* 2 hydroC + 6 benzoC */
            for (int ii=0;ii<2;++ii) {
                ai << "1 "<< ii+current_hC <<" ";
            } current_hC+=2;
            for (int ii=0;ii<6;++ii) {
                ai << "2 "<< ii+current_bC <<" ";
            } current_bC+=6;
            ai << "\n";
        } ai << "\n";
        ai << "combine_multibody_lists mbList1 ";
        for (int i=1;i<=dop;++i) ai << "mer"<<i << " ";
        ai << "\n";
        ai << "combine_trajectories mbList1 ";
        for (int i=1;i<=dop;++i) ai << "mer"<<i << " ";
        ai << "\n\n";
    }
    else if (species=="MartiniPS30") {
        
        is_mb_all_molecule=false;
        n_species_types=1;
        
        ai
        << "system_np"                                    << "\n"
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "polymer 32"                                   << "\n"
        << "1 2"                                          << "\n"
        << "30 90"                                        << "\n\n";
        ai
        << "create_list all \n"
        << "all"
        << "\n\n";
        
        int dop=30;
        int current_hC=0;
        int current_bC=0;
        
        for (int i=1;i<=dop;++i) {
            ai
            << "create_multibodies "
            << "mer"<<i<<" "<<"mbType1"<<" "<<centertype<<" "
            << "species_atomlist polymer ";
            /* 2 hydroC + 6 benzoC */
            for (int ii=0;ii<1;++ii) {
                ai << "1 "<< ii+current_hC <<" ";
            } current_hC+=1;
            for (int ii=0;ii<3;++ii) {
                ai << "2 "<< ii+current_bC <<" ";
            } current_bC+=3;
            ai << "\n";
        } ai << "\n";
        ai << "combine_multibody_lists mbList1 ";
        for (int i=1;i<=dop;++i) ai << "mer"<<i << " ";
        ai << "\n";
        ai << "combine_trajectories mbList1 ";
        for (int i=1;i<=dop;++i) ai << "mer"<<i << " ";
        ai << "\n\n";
        
    }
    else if (species=="SiO2") {
        
        n_species_types=2;
        
        ai
        << "system_nv"                                    << "\n" // NOTE: system_nv
        << "custom"                                       << "\n"
        << "${path_prdcustom}"                            << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                   << "\n"
        << "Si 512 O 1024"                                << "\n"
        << "1 2"                                          << "\n"
        << "1 0"                                        << "\n"
        << "0 1"                                       << "\n\n";
        ai
        << "create_list all" << "\n"
        << "all"             << "\n\n";
        if (speciesName.size()==3) {
            for (indexi=0;indexi<speciesName.size()-1;++indexi) {
                ai
                << "create_list  "
                << speciesName.at(indexi)+"_list"<< "\n"
                << "type_species "
                << speciesName.at(indexi)<<" "
                << speciesType.at(indexi)<<"\n";
            }
        } else {
            cout << species << " should have 3 species!\n";
            exit(EXIT_FAILURE);
        } ai << "\n";
    }
    else {
        cout
        << "in AmdatAnalysis::systemBlock_definition():\n"
        << "species specified ("<<species<<") NOT found!\n";
        exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::systemBlock_multibody_allmolecule(const StructureClass& sysVar,
                                                      ofstream& ai)
{
    test_AnalysisFolder(sysVar,*this);
    
    if (is_monoStruct) /** use struct in monomer bank **/
    {
        ai
        << "system_np"                                 << "\n"
        << "custom"                                    << "\n"
        << "${path_prdcustom}"                         << "\n"
        << "exponential ${n_prd_blocks} ${blocksize} "
        << "${prd_exp_base} 0 0 ${ts} "                << "\n";
        for (int polyindex=1;polyindex<=n_poly;++polyindex) {//NOTE: 1-based index
            ai
            << "leMer_"<<polyindex<<" 1 "
            << "imMer_"<<polyindex<<" "<<chainLen-2<<" "
            << "reMer_"<<polyindex<<" 1 ";
        } ai << "\n";
        
        ai << "1 2" << "\n";
        for (int polyindex=0;polyindex<n_poly;++polyindex) {
            ai
            << n_typeSet.at(0).at(0) << " " << n_typeSet.at(0).at(1) << "\n"
            << n_typeSet.at(1).at(0) << " " << n_typeSet.at(1).at(1) << "\n"
            << n_typeSet.at(2).at(0) << " " << n_typeSet.at(2).at(1) << "\n";
        } ai << "\n";
        
        ai
        << "create_list all \n"
        << "all"
        << "\n\n";
        
    } else { /** user-specified atom order **/
        
        systemBlock_definition(ai);
    }
    
    /** NOTE:
     ** In multibody, there can be multiple multibody lists with only one
     ** trajectory type ; multibodies can be combined to be in a same list and
     ** share a same trajectory **/
    
    //if (mbodyList.size()>1) {
    if (0) {
        cout
        << "in AmdatAnalysis::systemBlock_multibody_allmolecule():\n"
        << "mbodyList >1 \n"
        << "Curently the program wrapper only supports mbodyList = 1.\n";
        exit(EXIT_FAILURE);
    }
    //if (mbodyType.size()>1) {
    if (0) {
        cout
        << "in AmdatAnalysis::systemBlock_multibody_allmolecule():\n"
        << "mbodyType >1 \n"
        << "Curently the program wrapper only supports mbodyType = 1.\n";
        exit(EXIT_FAILURE);
    }
    
    if (is_mb_all_molecule)
    {
        ai
        << "create_multibodies "
        << mbodyList.at(0)<<" "
        << mbodyType.at(0)<<" "
        << centertype<<" "
        << "all_molecule"
        << "\n\n";
    }
    
    if (is_write_list_trajectory)
    {
        ai
        << "write_list_trajectory all "
        << "./statistics/mbTraj/all_"
        << get_filenameString()
        << ".xyz"
        << "\n";
        ai
        << "write_list_trajectory "
        << mbodyList.at(0)<<" "
        << "./statistics/mbTraj/"<<mbodyList.at(0)<<"_"
        << get_filenameString()
        << ".xyz"
        << "\n\n";
        
    } else {
        
        ai
        << "##write_list_trajectory all "
        << "./statistics/mbTraj/all_"
        << get_filenameString()
        << ".xyz"
        << "\n";
        ai
        << "##write_list_trajectory "
        << mbodyList.at(0)<<" "
        << "./statistics/mbTraj/"<<mbodyList.at(0)<<"_"
        << get_filenameString()
        << ".xyz"
        << "\n\n";
    }
}





void AmdatAnalysis::make_sigmamatrix(const StructureClass& sysVar)
{
    test_AnalysisFolder(sysVar,*this);
    
    vector<string> targetspecies;
    if (is_mbodies) {
        targetspecies=mbodyType;
    } else {
        targetspecies=speciesType;
    }
    /* sanity check of sigma matrix size */
    if (sigmaMatrix.size()!=targetspecies.size()) {
        cout << "sigmaMatrix size should equal number of types.\n";
        exit(EXIT_FAILURE);
    }
    double sigma=0;
    for (indexi=0;indexi<(int)sigmaMatrix.size();++indexi) {
        if (sigmaMatrix.at(indexi).size()!=sigmaMatrix.size()) {
            cout << "sigmaMatrix should be a square matrix.\n";
            exit(EXIT_FAILURE);
        }
        for (indexii=0;indexii<(int)sigmaMatrix.at(indexi).size();++indexii) {
            sigma=sigmaMatrix.at(indexi).at(indexii);
            if (sigma<0) {
                cout << "sigma value for string analysis is wrong.\n";
                exit(EXIT_FAILURE);
            }
        }
    }
    string out;
    out.append(return_AnalysisFolderPath(sysVar));
    out.append("/statistics/strings/sigma_matrix.txt");
    ofstream writeFile(out.c_str());
    /** NOTE:
     ** The dimension of the sigma matrix is based on trajectory types of
     ** interest **/
    for (indexi=0;indexi<(int)sigmaMatrix.size();++indexi) {
        writeFile << targetspecies.at(indexi) << " ";
        for (indexii=0;indexii<(int)sigmaMatrix.at(indexi).size();++indexii) {
            writeFile << sigmaMatrix.at(indexi).at(indexii) << " ";
        } writeFile << "\n";
    }
}





void AmdatAnalysis::make_amdatInputFile(const StructureClass& sysVar,
                                        const WorkScripts& ws)
{
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_inputs");
    
    if (is_strings) {
        o.append("/amdat_strings.inp");
    } else {
        o.append("/amdat.inp");
    } ofstream ai(o.c_str());
    
    string path;
    if (is_keep_new_trajectory) {
        path.append(return_SimulationFolderPath(sysVar));
        path.append("/production/trajectory");
        path.append("/new_trajectory_");//NOTE: "new_"
        path.append(sysVar.get_usic());
        path.append("_${trial}");
        path.append("_${namestring}");
        path.append("_T${Temp}");
        path.append(".path.custom");
    } else {
        path.append(return_SimulationFolderPath(sysVar));
        path.append("/production/trajectory");
        path.append("/trajectory_");
        path.append(sysVar.get_usic());
        path.append("_${trial}");
        path.append("_${namestring}");
        path.append("_T${Temp}");
        path.append(".prd.custom");
    }
    
    if (is_mbodies) {
        systemBlock_multibody_allmolecule(sysVar,ai);
        if (is_strings) make_sigmamatrix(sysVar);
    } else {
        systemBlock_definition(ai);
        if (is_strings) make_sigmamatrix(sysVar);
    }
    
    /** Typical AMDAT Analyses **/
    //--------------------------------------------------------------------------
    for (int rouse=0;rouse<n_rouseModes;++rouse)
    {
        write_amdat_typical(ai,sysVar,rouse);
    }
    
    /** String Analyses **/
    //--------------------------------------------------------------------------
    if (is_strings) write_amdat_stringInput(ai,sysVar);
    
    ai.close();
}





void AmdatAnalysis::write_amdat_typical(ofstream& ai,
                                        const StructureClass& sysVar,
                                        const int rouse)
{
    /** Mean Square Displacement **/
    //--------------------------------------------------------------------------
    if (is_msd) {
        
        if (is_mbodies) {
            
            if (rouse==0) {
                ai
                << "msd "
                << "./statistics/msd/msd_"
                << get_filenameString()
                << ".all.dat"
                << "\n";
                ai
                << "list all"
                << "\n\n";
            }
            ai
            << "msd "
            << "./statistics/msd/msd_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n\n";
            
        } else {
            
            for (indexi=0;indexi<(int)speciesName.size();++indexi) {
                ai
                << "msd "
                << "./statistics/msd/msd_"
                << get_filenameString()
                << "."+speciesName.at(indexi)+".dat"
                << "\n";
                if (speciesName.at(indexi)=="all") {
                    ai
                    << "list "+speciesName.at(indexi)
                    << "\n";
                } else {
                    ai
                    << "list "+speciesName.at(indexi)+"_list"
                    << "\n";
                }
            }
            
        } ai << "\n";
    }
    /** Calculates Non-Gaussian Parameter from MSD **/
    //--------------------------------------------------------------------------
    if (is_ngp) {
        
        if (is_mbodies) {
            
            if (rouse==0) {
                ai
                << "ngp "
                << "./statistics/ngp/ngp_"
                << get_filenameString()
                << ".all.dat"
                << "\n";
                ai
                << "list all"
                << "\n\n";
            }
            ai
            << "ngp "
            << "./statistics/ngp/ngp_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n\n";
            
        } else {
            
            for (indexi=0;indexi<(int)speciesName.size();++indexi) {
                
                ai
                << "ngp "
                << "./statistics/ngp/ngp_"
                << get_filenameString()
                << "."+speciesName.at(indexi)+".dat"
                << "\n";
                if (speciesName.at(indexi)=="all") {
                    ai
                    << "list "+speciesName.at(indexi)
                    << "\n";
                } else {
                    ai
                    << "list "+speciesName.at(indexi)+"_list"
                    << "\n";
                }
            }
            
        } ai << "\n";
    }
    /** Structure Factor **/
    //--------------------------------------------------------------------------
    if (is_strFac) {
        
        if (is_mbodies) {
            
            if (rouse==0) {
                ai
                << "structure_factor "
                << "./statistics/strFac/strFac_"
                << get_filenameString()
                << ".all.dat" << " "
                << symmetry   << " "
                << geometry   << " ";
                if (is_use_strFac_frame) {
                    ai << "0" << " ";
                    ai << strFac_frame;
                } else {
                    ai << "0" << " ";
                    ai << fullblock;
                } ai << "\n";
                ai
                << "list all"
                << "\n\n";
            }
            ai
            << "structure_factor "
            << "./statistics/strFac/strFac_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << symmetry   << " "
            << geometry   << " ";
            if (is_use_strFac_frame) {
                ai << "0" << " ";
                ai << strFac_frame;
            } else {
                ai << "0" << " ";
                ai << fullblock;
            } ai << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n\n";
            
        } else {
            
            for (indexi=0;indexi<(int)speciesName.size();++indexi) {
                
                if (speciesName.at(indexi)=="all")
                {
                    if (n_species_types>1)
                    {
                        /** asymmetric **/
                        for (int i=0;i<n_species_types-1;++i)
                        {
                            if (speciesName.at(i+1)=="all") {
                                cout
                                << "in AmdatAnalysis::make_amdatInputFile():\n"
                                << "strFac asymmetric pair is wrong!\n\n";
                                exit(EXIT_FAILURE);
                            }
                            ai
                            << "structure_factor "
                            << "./statistics/strFac/strFac_"
                            << get_filenameString()
                            << "."+speciesName.at(i)+speciesName.at(i+1)+".dat"<<" "
                            << "asymmetric "
                            << geometry   << " ";
                            if (is_use_strFac_frame) {
                                ai << "0" << " ";
                                ai << strFac_frame;
                            } else {
                                ai << "0" << " ";
                                ai << fullblock;
                            } ai << "\n";
                            ai
                            << "list "+speciesName.at(i)+"_list \n"
                            << "list "+speciesName.at(i+1)+"_list \n";
                        }
                        ai
                        << "structure_factor "
                        << "./statistics/strFac/strFac_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat" << " "
                        << symmetry   << " "
                        << geometry   << " ";
                        if (is_use_strFac_frame) {
                            ai << "0" << " ";
                            ai << strFac_frame;
                        } else {
                            ai << "0" << " ";
                            ai << fullblock;
                        } ai << "\n";
                        ai
                        << "list all"
                        << "\n\n";
                    } else {
                        ai
                        << "structure_factor "
                        << "./statistics/strFac/strFac_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat" << " "
                        << symmetry   << " "
                        << geometry   << " ";
                        if (is_use_strFac_frame) {
                            ai << "0" << " ";
                            ai << strFac_frame;
                        } else {
                            ai << "0" << " ";
                            ai << fullblock;
                        } ai << "\n";
                        ai
                        << "list all" //NOTE: "all"
                        << "\n\n";
                    }
                }
                else {
                    ai
                    << "structure_factor "
                    << "./statistics/strFac/strFac_"
                    << get_filenameString()
                    << "."+speciesName.at(indexi)+".dat" << " "
                    << symmetry   << " "
                    << geometry   << " ";
                    if (is_use_strFac_frame) {
                        ai << "0" << " ";
                        ai << strFac_frame;
                    } else {
                        ai << "0" << " ";
                        ai << fullblock;
                    } ai << "\n";
                    ai
                    << "list "+speciesName.at(indexi)+"_list"
                    << "\n\n";
                }
            }
            
        } ai << "\n";
    }
    /** Radial Distribution Function **/
    //--------------------------------------------------------------------------
    if (is_rdf) {
        
        if (is_mbodies) {
            
            if (rouse==0) {
                ai
                << "\n"
                << "rdf "
                << "./statistics/rdf/rdf_"
                << get_filenameString()
                << ".all.dat" << " "
                << symmetry   << " "
                << rdf_nbins  << " "
                << fullblock  << " "
                //<< maxLenScale/2.0 << " "
                << "0 "
                << "\\srdf"
                << "_${usic}"
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << "_all"
                << "\n";
                ai
                << "list all"
                << "\n";
                ai
                << "structure_factor_from_rdf "
                << "./statistics/rdf/strFac_rdf_"
                << get_filenameString()
                << ".all.dat" << " "
                << rdf_nbins << " "
                << "rdf"
                << "_${usic}"
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << "_all"
                << "\n\n";
            }
            ai
            << "rdf "
            << "./statistics/rdf/rdf_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << symmetry   << " "
            << rdf_nbins  << " "
            << fullblock  << " "
            //<< maxLenScale/2.0 << " "
            << "0 "
            << "\\srdf_"
            << get_filenameString()
            << "_"<<mbodyList.at(rouse)
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n";
            ai
            << "structure_factor_from_rdf "
            << "./statistics/rdf/strFac_rdf_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"<<" "
            << rdf_nbins << " "
            << "rdf"
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}"
            << "_"<<mbodyList.at(rouse)
            << "\n";
            
        } else {
            
            for (indexi=0;indexi<(int)speciesName.size();++indexi) {
                
                if (speciesName.at(indexi)=="all")
                {
                    if (n_species_types>1)
                    {
                        /** asymmetric **/
                        for (int i=0;i<n_species_types-1;++i)
                        {
                            if (speciesName.at(i+1)=="all") {
                                cout
                                << "in AmdatAnalysis::make_amdatInputFile():\n"
                                << "rdf asymmetric pair is wrong!\n\n";
                                exit(EXIT_FAILURE);
                            }
                            ai
                            << "rdf "
                            << "./statistics/rdf/rdf_"
                            << get_filenameString()
                            << "."+speciesName.at(i)+speciesName.at(i+1)+".dat"<<" "
                            << "asymmetric "
                            << rdf_nbins  << " "
                            << fullblock  << " "
                            << "0 "
                            << "\\srdf"
                            << "_${usic}"
                            << "_${trial}"
                            << "_${namestring}"
                            << "_T${Temp}"
                            << "_"<<speciesName.at(i)+speciesName.at(i+1)
                            << "\n";
                            ai
                            << "list "+speciesName.at(i)+"_list \n"
                            << "list "+speciesName.at(i+1)+"_list \n";
                            ai
                            << "structure_factor_from_rdf "
                            << "./statistics/rdf/strFac_rdf_"
                            << get_filenameString()
                            << "."+speciesName.at(i)+speciesName.at(i+1)+".dat"<<" "
                            << rdf_nbins  << " "
                            << "rdf"
                            << "_${usic}"
                            << "_${trial}"
                            << "_${namestring}"
                            << "_T${Temp}"
                            << "_"<<speciesName.at(i)+speciesName.at(i+1)<<"\n"
                            << "\n";
                        }
                        ai
                        << "rdf "
                        << "./statistics/rdf/rdf_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat"<<" "
                        << symmetry   << " "
                        << rdf_nbins  << " "
                        << fullblock  << " "
                        << "0 "
                        << "\\srdf"
                        << "_${usic}"
                        << "_${trial}"
                        << "_${namestring}"
                        << "_T${Temp}"
                        << "_"<<speciesName.at(indexi)
                        << "\n";
                        ai
                        << "list all"
                        << "\n";
                        ai
                        << "structure_factor_from_rdf "
                        << "./statistics/rdf/strFac_rdf_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat"<<" "
                        << rdf_nbins << " "
                        << "rdf"
                        << "_${usic}"
                        << "_${trial}"
                        << "_${namestring}"
                        << "_T${Temp}"
                        << "_"<<speciesName.at(indexi)<< "\n"
                        << "\n";
                    } else {
                        ai
                        << "rdf "
                        << "./statistics/rdf/rdf_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat" << " "
                        << symmetry   << " "
                        << rdf_nbins  << " "
                        << fullblock  << " "
                        << "0 "
                        << "\\srdf"
                        << "_${usic}"
                        << "_${trial}"
                        << "_${namestring}"
                        << "_T${Temp}"
                        << "_"<<speciesName.at(indexi)<< "\n";
                        ai
                        << "list all" //NOTE: "all"
                        << "\n";
                        ai
                        << "structure_factor_from_rdf "
                        << "./statistics/rdf/strFac_rdf_"
                        << get_filenameString()
                        << "."+speciesName.at(indexi)+".dat"<<" "
                        << rdf_nbins << " "
                        << "rdf"
                        << "_${usic}"
                        << "_${trial}"
                        << "_${namestring}"
                        << "_T${Temp}"
                        << "_"<<speciesName.at(indexi)<< "\n"
                        << "\n";
                    }
                    
                } else {
                    ai
                    << "rdf "
                    << "./statistics/rdf/rdf_"
                    << get_filenameString()
                    << "."+speciesName.at(indexi)+".dat" << " "
                    << symmetry   << " "
                    << rdf_nbins  << " "
                    << fullblock  << " "
                    << "0 "
                    << "\\srdf"
                    << "_${usic}"
                    << "_${trial}"
                    << "_${namestring}"
                    << "_T${Temp}"
                    << "_"<<speciesName.at(indexi)<<"\n";
                    ai
                    << "list "+speciesName.at(indexi)+"_list"
                    << "\n";
                    ai
                    << "structure_factor_from_rdf "
                    << "./statistics/rdf/strFac_rdf_"
                    << get_filenameString()
                    << "."+speciesName.at(indexi)+".dat"<<" "
                    << rdf_nbins << " "
                    << "rdf"
                    << "_${usic}"
                    << "_${trial}"
                    << "_${namestring}"
                    << "_T${Temp}"
                    << "_"<<speciesName.at(indexi)<<"\n"
                    << "\n";
                }
            }
        } ai << "\n";
    }
    /** Calculates the composition and number density **/
    //--------------------------------------------------------------------------
    if (is_composition) {
        
        if (is_mbodies) {
            
            ai
            << "composition "
            << "./statistics/composition/comp_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n";
            
        } else {
            
            ai
            << "composition "
            << "./statistics/composition/comp_"
            << get_filenameString()
            << ".all.dat"
            << "\n"
            << "list all"
            << "\n";
            
        } ai << "\n";
    }
    /** Calculates distribution of square displacements at a specified time **/
    //--------------------------------------------------------------------------
    if (is_u2dist) {
        
        if (is_mbodies) {
            
            ai
            << "u2dist "
            << "./statistics/u2dist/u2dist_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << u2dist_nbins
            << " "
            << max_u2
            << " "
            << "${dwf_frame}"
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n";
            
        } else {
            
            ai
            << "u2dist "
            << "./statistics/u2dist/u2dist_"
            << get_filenameString()
            << ".all.dat"
            << " "
            << u2dist_nbins
            << " "
            << max_u2
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list all"
            << "\n";
            
        } ai << "\n";
    }
    /** Calculates distribution of inverse Debye-Waller factor values 1/u2 **/
    //--------------------------------------------------------------------------
    if (is_stiffness_dist) {
        
        if (is_mbodies) {
            
            ai
            << "stiffness_dist "
            << "./statistics/stiffness/stiffness_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << stiffness_nbins
            << " "
            << max_stiffness
            << " "
            << "${dwf_frame}"
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n";
            
        } else {
            
            ai
            << "stiffness_dist "
            << "./statistics/stiffness/stiffness_"
            << get_filenameString()
            << ".all.dat"
            << " "
            << stiffness_nbins
            << " "
            << max_stiffness
            << " "
            << "${dwf_frame}"
            << "\n"
            << "list all"
            << "\n";
            
        } ai << "\n";
    }
    /** Bond Autocorrelation Function **/
    //--------------------------------------------------------------------------
    if (is_baf) {
        
        if (is_mbodies) {
            
            ai
            << "baf "
            << "./statistics/baf/baf_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << mbodyList.at(rouse)<<"_endtoend"
            << "\n\n";
            
        } else {
            
            cout
            << "in AmdatAnalysis::write_amdat_typical():\n"
            << "baf requires defining multibodies.\n";
            exit(EXIT_FAILURE);
        }
        
    }
    /** Self Part of Intermediate Scattering Function **/
    //--------------------------------------------------------------------------
    if (is_isfs) {
        
        if (is_mbodies) {
            
            if (rouse==0) {
                ai
                << "isfs "
                << "./statistics/isfs/isfs_"
                << get_filenameString()
                << ".all.dat" << " "
                << "${waveindex}"   << " "
                << "${waveindex}"   << " "
                << geometry   << " "
                << "${maxLenScale}" << " "
                << fullblock
                << "\n"
                << "list all"
                << "\n\n";
            }
            ai
            << "isfs "
            << "./statistics/isfs/isfs_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << "${waveindex}"   << " "
            << "${waveindex}"   << " "
            << geometry   << " "
            << "${maxLenScale}" << " "
            << fullblock
            << "\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n\n";
            
        } else {
            
            for (indexi=0;indexi<(int)speciesName.size();++indexi) {
                ai
                << "isfs "
                << "./statistics/isfs/isfs_"
                << get_filenameString()
                << "."+speciesName.at(indexi)+".dat" << " "
                << "${waveindex}"   << " "
                << "${waveindex}"   << " "
                << geometry   << " "
                << "${maxLenScale}" << " "
                << fullblock
                << "\n";
                if (speciesName.at(indexi)=="all") {
                    ai
                    << "list "+speciesName.at(indexi)
                    << "\n";
                } else {
                    ai
                    << "list "+speciesName.at(indexi)+"_list"
                    << "\n";
                }
            } ai << "\n";
        }
    }
    /** Full Intermediate Scattering Function **/
    //--------------------------------------------------------------------------
    if (is_isf) {
        
        if (is_mbodies) {
            
            ai
            << "isf_list "
            << "./statistics/isf/isf_"
            << get_filenameString()
            << "."<<mbodyList.at(rouse)<<".dat"
            << " "
            << mbodyList.at(rouse)
            << " "
            << geometry
            << " "
            << "${waveindex}"
            << " "
            << "${waveindex}"
            << " "
            << "\n";
            
        } else {
            
            ai
            << "isf_list "
            << "./statistics/isf/isf_"
            << get_filenameString()
            << ".all.dat"
            << " "
            << "all"
            << " "
            << geometry
            << " "
            << "${waveindex}"
            << " "
            << "${waveindex}"
            << " "
            << "\n";
            
        } ai << "\n";
    }
}





void AmdatAnalysis::write_amdat_loopstring(ofstream& ai,
                                           const StructureClass& sysVar,
                                           const int rouse)
{
    if (is_use_displacement_list)
    {
        if (is_mbodies)
        {
            // is_use_displacement_list = true
            // is_mbodies = true
            
            ai
            /* create analysis file (WITH .dat appended) */
            << "displacement_list "
            << "./statistics/strings/displacement_stats_"
            << get_filenameString()<< "_frame${frame}";
            if (is_mbodies) ai << "."<<mbodyList.at(rouse)<<".dat";//NOTE: .dat
            ai << " ";
            /* create list obj. (WITHOUT .dat appended) */
            ai
            << "displacement_list_"
            << get_filenameString()<< "_frame${frame}";
            if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
            ai << " ";
            ai << "${frame}\n";
            ai
            << "list "<<mbodyList.at(rouse)
            << "\n";
            
        } else {
            
            // is_use_displacement_list = true
            // is_mbodies = false
            
            ai
            /* create analysis file (WITH .dat appended) */
            << "displacement_list "
            << "./statistics/strings/displacement_stats_"
            << get_filenameString()
            << ".all.dat"
            << " "
            /* create list obj. (WITHOUT .dat appended) */
            << "displacement_list_"
            << get_filenameString()<< "_frame${frame}"
            << " "
            << "${frame}";
            ai << "\n";
            ai
            << "list all"
            << "\n";
            
        } ai << "\n";
        
        if (is_write_pdb_file) {
            ai
            << "value_list "
            << "write_pdb"
            << " "
            /* use list object */
            << "displacement_list_"<< get_filenameString()
            << "_frame${frame}";
            if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
            ai << " ";
            ai << "${frame}\n";
            ai
            << " "
            /* create pdb filename stem */
            << "./statistics/strings/displacement_pdb_"
            << get_filenameString()<< "_frame${frame}";
            if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
            ai
            << " "
            << "0 0"
            << "\n\n";
        }
        
        if (threshold_keyword=="between")
        {
            if (is_mbodies)
            {
                // is_use_displacement_list = true
                // threshold_keyword == "between"
                // is_mbodies = true
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString()<< "_frame${frame}";
                if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi << "_"
                << get_filenameString()<< "_frame${frame}";
                if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile_lo
                << " "
                << threshold_percentile_hi
                << "\n\n";
                
            } else {
                
                // is_use_displacement_list = true
                // threshold_keyword == "between"
                // is_mbodies = false
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString()<< "_frame${frame}";
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi
                << "_${usic}"
                << "_${trial}"
                << "_${namestring}"
                << "_T${Temp}"
                << "_frame${frame}";
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile_lo
                << " "
                << threshold_percentile_hi
                << "\n\n";
            }
            
        } /* End if (threshold_keyword=="between") */
        else
        {
            if (is_mbodies)
            {
                // is_use_displacement_list = true
                // threshold_keyword=="greater" or "less"
                // is_mbodies = true
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString()<< "_frame${frame}";
                if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << get_filenameString()<< "_frame${frame}";
                if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile
                << "\n\n";
                
            } else {
                
                // threshold_keyword=="greater" or "less"
                // is_mbodies = false
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString()<<"_frame${frame}";
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << get_filenameString()<<"_frame${frame}";
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile
                << "\n\n";
            }
        }
    } /* End if (threshold_keyword) */
    
    if (is_use_find_fast) {
        
        if (is_mbodies)
        {
            cout
            << "in AmdatAnalysis::write_amdat_stringInput():\n"
            << "mbodies not implemented for using find_fast approach.\n";
            exit(EXIT_FAILURE);
        }
        ai
        << "vhs "
        << "./statistics/strings/vhs_"
        << get_filenameString()<< "_frame${frame}";
        ai
        << ".all.dat"
        << " "
        << vhs_rangebinned
        << " "
        << vhs_nbins
        << "\n"
        << "list all"
        << "\n";
        ai
        << "compare_gaussian ./statistics/strings/gaussian_"
        << get_filenameString()<<"_frame${frame}";
        ai
        << ".all.dat"
        << " "
        << "${frame}"
        << "\n";
        ai
        << "list all"
        << "\n";
        ai
        << "find_fast "
        << "fast_"
        << get_filenameString()
        << " "
        << "./statistics/strings/fast_"
        << get_filenameString()<<"_frame${frame}"
        << ".all.dat"
        << "\n"
        << "list all"
        << "\n";
    }
    
    // analysis based on voronoi neighbors
    //--------------------------------------------------------------------------
    if (is_use_voroNeighbors)
    {
        ai
        << "neighbor_decorrelation_function "
        << "./statistics/voro/ndf_"<< get_filenameString()<<"_frame${frame}";
        if (is_mbodies) {
            ai<<"."<<mbodyList.at(rouse)<<".dat";
        } else {
            ai<<"."<<analysispart<<".dat";
        } ai << " vneigh\n";
        displacement_value_list(ai,rouse);
        
        ai
        << "persistent_neighbors "
        << "${frame} vneigh p_neigh_frame${frame} per_neigh_frame${frame} "
        << centertype<<"\n";
        displacement_value_list(ai,rouse);
        
        ai
        << "size_statistics "
        << "./statistics/strings/strings_"
        << get_filenameString()<<"_frame${frame}";
        if (is_mbodies) {
            ai<<"."<<mbodyList.at(rouse)<<".dat";
        } else {
            ai<<"."<<analysispart<<".dat";
        } ai << " p_neigh_frame${frame} "<<n_moments<<"\n\n";
    }
    //--------------------------------------------------------------------------
    
    if (is_use_displacement_list) { // based on "value_list"
        
        if (is_mbodies) { //for multibodies
            
            if (!is_use_voroNeighbors)
            {
                if (false)
                {
                    // MSD
                    //----------------------------------------------------------
                    ai
                    << "msd "
                    << "./statistics/msd/msd_"
                    << get_filenameString()<<"_frame${frame}"
                    << "."+analysispart+".dat"
                    << "\n";
                    displacement_value_list(ai,rouse);
                    
                    // ISFS
                    //----------------------------------------------------------
                    ai
                    << "isfs "
                    << "./statistics/isfs/isfs_"
                    << get_filenameString()<<"_frame${frame}"
                    << "."+analysispart+".dat" << " "
                    << "${waveindex}"   << " "
                    << "${waveindex}"   << " "
                    << geometry   << " "
                    << "${maxLenScale}" << " "
                    << fullblock
                    << "\n";
                    displacement_value_list(ai,rouse);
                }
                
                // STRINGS
                //--------------------------------------------------------------
                ai
                << "streamlined_strings "
                << "./statistics/strings/strings_"
                << get_filenameString()<<"_frame${frame}"
                << "."+analysispart+".dat";
                //if (is_mbodies) ai <<"."<<mbodyList.at(rouse)<<".dat";//NOTE: .dat
                //else ai << "."+analysispart+".dat";
                ai
                << " "
                << "${frame}"
                << " "
                << strings_threshold
                << " "
                << "./statistics/strings/sigma_matrix.txt"
                << " "
                << n_moments
                << "\n";
                displacement_value_list(ai,rouse);
                
                
                // Create multibodies of strings
                //--------------------------------------------------------------
                ai
                << "string_multibodies null ${frame}"
                << " "
                << strings_threshold
                << " "
                << "./statistics/strings/sigma_matrix.txt"
                << " "
                << "stringmbs_frame${frame}"
                << " "
                << "string_frame${frame}"
                << " "
                << centertype
                << "\n";
                displacement_value_list(ai,rouse);
                
                
                // Calc Rg of strings
                //--------------------------------------------------------------
                ai
                << "gyration_radius "
                << "./statistics/strings/Rg_strings_"
                << get_filenameString()<<"_frame${frame}"
                << "."+analysispart+".dat"
                << " "
                << "stringmbs_frame${frame}"
                << "\n\n";
            }
            ai
            << "delete_trajectory_list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile << "_"
            << get_filenameString()
            << "_frame${frame}";
            if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
            ai << "\n";
            ai
            << "delete_valuelist "
            << "displacement_list_"
            << get_filenameString()
            << "_frame${frame}";
            if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
            ai << "\n\n";
            
        } else { //for non-multibodies
            
            if (!is_use_voroNeighbors)
            {
                if (false)
                {
                    // MSD
                    //----------------------------------------------------------
                    ai
                    << "msd "
                    << "./statistics/msd/msd_"
                    << get_filenameString()<<"_frame${frame}"
                    << "."+analysispart+".dat"
                    << "\n";
                    displacement_value_list(ai,rouse);
                    
                    // ISFS
                    //----------------------------------------------------------
                    ai
                    << "isfs "
                    << "./statistics/isfs/isfs_"
                    << get_filenameString()<<"_frame${frame}"
                    << "."+analysispart+".dat" << " "
                    << "${waveindex}"   << " "
                    << "${waveindex}"   << " "
                    << geometry   << " "
                    << "${maxLenScale}" << " "
                    << fullblock
                    << "\n";
                    displacement_value_list(ai,rouse);
                }
                
                // STRINGS
                //--------------------------------------------------------------
                ai
                << "streamlined_strings "
                << "./statistics/strings/strings_"
                << get_filenameString()
                <<"_frame${frame}"
                << "."+analysispart+".dat";
                ai << " "
                << "${frame}"
                << " "
                << strings_threshold
                << " "
                << "./statistics/strings/sigma_matrix.txt"
                << " "
                << n_moments
                << "\n";
                displacement_value_list(ai,rouse);
                
                
                // Create multibodies of strings
                //--------------------------------------------------------------
                ai
                << "string_multibodies null ${frame}"
                << " "
                << strings_threshold
                << " "
                << "./statistics/strings/sigma_matrix.txt"
                << " "
                << "stringmbs_frame${frame}"
                << " "
                << "string_frame${frame}"
                << " "
                << centertype
                << "\n";
                displacement_value_list(ai,rouse);
                
                
                // Calc Rg of strings
                //--------------------------------------------------------------
                ai
                << "gyration_radius "
                << "./statistics/strings/Rg_strings_"
                << get_filenameString()<<"_frame${frame}"
                << "."+analysispart+".dat"
                << " "
                << "stringmbs_frame${frame}"
                << "\n\n";
            }
            ai
            << "delete_trajectory_list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile << "_"
            << get_filenameString()
            << "_frame${frame}"
            << "\n";
            ai
            << "delete_valuelist "
            << "displacement_list_"<< get_filenameString()
            <<"_frame${frame}"
            << "\n\n";
        }
    }
    
    if (is_use_find_fast) { // based on "find_fast"
        
        // is_use_displacement_list = false
        // is_use_find_fast = true
        
        ai
        << "streamlined_strings "
        << "./statistics/strings/strings_"
        << get_filenameString()<<"_frame${frame}"
        << "."+analysispart+".dat";
        ai << " "
        << "${frame}"
        << " "
        << strings_threshold
        << " "
        << "./statistics/strings/sigma_matrix.txt"
        << " "
        << n_moments
        << "\n";
        
        ai
        << "list "
        << "fast_"<<get_filenameString()<<"_frame${frame}"
        << "\n\n";
        
        ai
        << "delete_trajectory_list "
        << "fast_"<<get_filenameString()<<"_frame${frame}"
        << "\n\n";
    }
}





void AmdatAnalysis::write_amdat_stringInput(ofstream& ai,
                                            const StructureClass& sysVar,
                                            const int frame_index)
{
    if (is_use_displacement_list)
    {
        if (is_mbodies)
        {
            // is_use_displacement_list = true
            // is_mbodies = true
            
            ai
            /* create analysis file (WITH .dat appended) */
            << "displacement_list "
            << "./statistics/strings/displacement_stats_"
            << get_filenameString();
            if (frame_index>0) ai << "_frame" << frame_index;
            if (is_mbodies) ai << "."<<mbodyList.at(0)<<".dat";//NOTE: .dat
            ai << " ";
            /* create list obj. (WITHOUT .dat appended) */
            ai
            << "displacement_list_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            if (is_mbodies) ai <<"_"<<mbodyList.at(0);
            ai
            << " ";
            if (frame_index==0) {
                ai << "${frame}";
            } else {
                ai << frame_index;
            } ai << "\n";
            
            ai
            << "list "<<mbodyList.at(0)
            << "\n";
            
        } else {
            
            // is_use_displacement_list = true
            // is_mbodies = false
            
            ai
            /* create analysis file (WITH .dat appended) */
            << "displacement_list "
            << "./statistics/strings/displacement_stats_"
            << get_filenameString()
            << ".all.dat"
            << " "
            /* create list obj. (WITHOUT .dat appended) */
            << "displacement_list_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai
            << " ";
            if (frame_index==0) {
                ai << "${frame}";
            } else {
                ai << frame_index;
            } ai << "\n";
            
            ai
            << "list all"
            << "\n";
            
        } ai << "\n";
        
        if (is_write_pdb_file) {
            ai
            << "value_list "
            << "write_pdb"
            << " "
            /* use list object */
            << "displacement_list_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai
            << " "
            /* create pdb filename stem */
            << "./statistics/strings/displacement_pdb_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai
            << " "
            << "0 0"
            << "\n\n";
        }
        
        if (threshold_keyword=="between")
        {
            if (is_mbodies)
            {
                // is_use_displacement_list = true
                // threshold_keyword == "between"
                // is_mbodies = true
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi << "_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile_lo
                << " "
                << threshold_percentile_hi
                << "\n\n";
                
            } else {
                
                // is_use_displacement_list = true
                // threshold_keyword == "between"
                // is_mbodies = false
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile_lo << "_"
                << threshold_percentile_hi << "_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile_lo
                << " "
                << threshold_percentile_hi
                << "\n\n";
            }
            
        } /* End if (threshold_keyword=="between") */
        else
        {
            if (is_mbodies)
            {
                // is_use_displacement_list = true
                // threshold_keyword=="greater" or "less"
                // is_mbodies = true
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile
                << "\n\n";
                
            } else {
                
                // threshold_keyword=="greater" or "less"
                // is_mbodies = false
                
                ai
                << "value_list "
                << "threshold_percentile"
                << " "
                << "displacement_list_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                ai
                << " "
                << "list_"
                << threshold_keyword << "_"
                << threshold_percentile << "_"
                << get_filenameString();
                if (frame_index>0) ai <<"_frame"<<frame_index;
                ai
                << " "
                << threshold_keyword
                << " "
                << threshold_percentile
                << "\n\n";
            }
        }
    } /* End if (threshold_keyword) */
    
    if (is_use_find_fast) {
        
        if (is_mbodies)
        {
            cout
            << "in AmdatAnalysis::write_amdat_stringInput():\n"
            << "mbodies not implemented for using find_fast approach.\n";
            exit(EXIT_FAILURE);
        }
        
        ai
        << "vhs "
        << "./statistics/strings/vhs_"
        << get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        ai
        << ".all.dat"
        << " "
        << vhs_rangebinned
        << " "
        << vhs_nbins
        << "\n"
        << "list all"
        << "\n";
        
        ai
        << "compare_gaussian ./statistics/strings/gaussian_"
        << get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        ai
        << ".all.dat"
        << " ";
        if (frame_index==0) {
            ai << "${frame}";
        } else {
            ai << frame_index;
        } ai << "\n";
        ai
        << "list all"
        << "\n";
        
        ai
        << "find_fast "
        << "fast_"
        << get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        ai
        << " "
        << "./statistics/strings/fast_"
        << get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        ai
        << ".all.dat"
        << "\n"
        << "list all"
        << "\n";
    }
    
    if (is_use_displacement_list) { // based on "value_list"
        
        if (is_mbodies)
        {
            // is_use_displacement_list = true
            // is_mbodies = true
            
            find_strings_syntax(ai,frame_index);
            
            if (!is_stringmb) //NOTE: is_stringmb=false
            {
                if (threshold_keyword=="between") { // "between"
                    ai
                    << "list "
                    << "list_"
                    << threshold_keyword << "_"
                    << threshold_percentile_lo << "_"
                    << threshold_percentile_hi << "_"
                    << get_filenameString();
                    if (frame_index>0) ai <<"_frame"<<frame_index;
                    if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                    ai
                    << "\n\n";
                } else { // "greater" or "less"
                    ai
                    << "list "
                    << "list_"
                    << threshold_keyword << "_"
                    << threshold_percentile << "_"
                    << get_filenameString();
                    if (frame_index>0) ai <<"_frame"<<frame_index;
                    if (is_mbodies) ai <<"_"<<mbodyList.at(0);
                    ai
                    << "\n\n";
                }
            }
            ai
            << "delete_valuelist "
            << "displacement_list_"<< get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            if (is_mbodies) ai <<"_"<<mbodyList.at(0);
            ai << "\n\n";
            
        } else {
            
            // is_use_displacement_list = true
            // is_mbodies = false
            
            find_strings_syntax(ai,frame_index);
            
            if (!is_stringmb) //NOTE: is_stringmb=false
            {
                if (threshold_keyword=="between") { // "between"
                    ai
                    << "list "
                    << "list_"
                    << threshold_keyword << "_"
                    << threshold_percentile_lo << "_"
                    << threshold_percentile_hi
                    << "_${usic}"
                    << "_${trial}"
                    << "_${namestring}"
                    << "_T${Temp}";
                    if (frame_index>0) ai <<"_frame"<<frame_index;
                    ai << "\n\n";
                } else { // "greater" or "less"
                    ai
                    << "list "
                    << "list_"
                    << threshold_keyword << "_"
                    << threshold_percentile << "_"
                    << get_filenameString();
                    if (frame_index>0) ai <<"_frame"<<frame_index;
                    ai << "\n\n";
                }
            }
            /** delete used displacement list **/
            ai
            << "delete_valuelist "
            << "displacement_list_"<< get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai << "\n\n";
        }
    }
    
    if (is_use_find_fast) { // based on "find_fast"
        
        // is_use_displacement_list = false
        // is_use_find_fast = true
        
        find_strings_syntax(ai,frame_index);
        
        ai
        << "list "
        << "fast_"<< get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        ai << "\n\n";
    }
}





void AmdatAnalysis::find_strings_syntax(std::ofstream& ai,
                                        const int frame_index)
{
    //NOTE: indexi is changed accordingly with outer loop calling this function
    
    if (is_stringmb) {
        stringmb_syntax(ai,frame_index);
    } else {
        ai
        << "strings "
        << "./statistics/strings/strings_"<< get_filenameString();
        if (frame_index>0) ai <<"_frame"<<frame_index;
        if (is_mbodies) ai <<"."<<mbodyList.at(0)<<".dat";//NOTE: .dat
        else ai << "."+analysispart+".dat";
        ai << " ";
        if (frame_index==0) {
            ai << "${frame}";
        } else {
            ai << frame_index;
        }
        ai
        << " "
        << strings_threshold
        << " "
        << "./statistics/strings/sigma_matrix.txt"
        << " "
        << "stringTrajList_"<< get_filenameString();
        if (frame_index>0) ai << "_frame" << frame_index;
        ai << "\n";
    }
}





void AmdatAnalysis::stringmb_syntax(std::ofstream& ai,
                                    const int frame_index)
{
    //NOTE: indexi is changed accordingly with outer loop calling this function
    
    ai
    << "string_multibodies "
    << "null ";//NOTE: a null argument
    if (frame_index==0) {
        ai << "${frame}";
    } else {
        ai << frame_index;
    }
    ai
    << " "
    << strings_threshold
    << " "
    << "./statistics/strings/sigma_matrix.txt"
    << " "
    << "stringmb_"<< get_filenameString();
    if (frame_index>0) ai <<"_frame"<<frame_index;
    if (is_mbodies) ai <<"_"<<mbodyList.at(0);
    ai
    << " "
    << "stringmbType"
    << " "
    << centertype
    << "\n";
    if (is_mbodies) {
        if (threshold_keyword=="between") { // "between"
            ai
            << "list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile_lo << "_"
            << threshold_percentile_hi << "_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            if (is_mbodies) ai <<"_"<<mbodyList.at(0);
            ai << "\n\n";
        } else { // "greater" or "less"
            ai
            << "list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile
            << "_${usic}"
            << "_${trial}"
            << "_${namestring}"
            << "_T${Temp}";
            if (frame_index>0) ai <<"_frame"<<frame_index;
            if (is_mbodies) ai <<"_"<<mbodyList.at(0);
            ai << "\n\n";
        }
    } else {
        if (threshold_keyword=="between") { // "between"
            ai
            << "list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile_lo << "_"
            << threshold_percentile_hi << "_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai << "\n\n";
        } else { // "greater" or "less"
            ai
            << "list "
            << "list_"
            << threshold_keyword << "_"
            << threshold_percentile << "_"
            << get_filenameString();
            if (frame_index>0) ai <<"_frame"<<frame_index;
            ai << "\n\n";
        }
    }
    /** find string distribution by size stats **/
    ai
    << "size_statistics "
    << "./statistics/strings/strings_"<< get_filenameString();
    if (frame_index>0) ai <<"_frame"<<frame_index;
    if (is_mbodies) ai <<"."<<mbodyList.at(0)<<".dat";//NOTE: .dat
    else ai << "."+analysispart+".dat";
    ai << " ";
    ai
    << "stringmb_"<< get_filenameString();
    if (frame_index>0) ai <<"_frame"<<frame_index;
    if (is_mbodies) ai <<"_"<<mbodyList.at(0);
    ai << " ";
    ai << n_moments;
    ai << "\n\n";
}





void AmdatAnalysis::displacement_value_list(ofstream& ai,const int rouse)
{
    if (threshold_keyword=="between") { // "between"
        ai
        << "list "
        << "list_"
        << threshold_keyword << "_"
        << threshold_percentile_lo << "_"
        << threshold_percentile_hi << "_"
        << get_filenameString()
        <<"_frame${frame}";
        if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
        ai << "\n\n";
    } else { // "greater" or "less"
        ai
        << "list "
        << "list_"
        << threshold_keyword << "_"
        << threshold_percentile << "_"
        << get_filenameString()
        <<"_frame${frame}";
        if (is_mbodies) ai <<"_"<<mbodyList.at(rouse);
        ai << "\n\n";
    }
}





std::vector<int> AmdatAnalysis::make_amdatInputFile_strings_loop(const StructureClass& sysVar,
                                                                 const WorkScripts& ws,
                                                                 const int n_trl,
                                                                 const int n_sys,
                                                                 const double Temp_d,
                                                                 const int peakframe,
                                                                 const int blockframe)
{
    string output;
    output.append(return_AnalysisFolderPath(sysVar));
    output.append("/AMDAT_inputs");
    output.append("/amdat_strings_");
    output.append(sysVar.get_usic());
    output.append("_00"+to_string((long long int)n_trl));
    output.append("_"+sysVar.get_nameString(n_sys));
    output.append("_T"+to_string((long long int)Temp_d));
    output.append(".inp");
    ofstream ai(output.c_str());
    
    bool is_use_all_inblockFrames=true;
    int  beg_frame=0;
    int  end_frame=0;
    
    if (is_use_all_inblockFrames)
    {
        beg_frame=1;//NOTE: DON'T use zero!
        end_frame=blockframe;
        //end_frame=blockframe-1;//NOTE: to avoid final frame issues
    }
    else
    {
        int pre_frames=20; // # of frames before NGP peak
        int pos_frames=20;  // # of frames including and after NGP peak
        beg_frame=peakframe-pre_frames;
        end_frame=peakframe+(pos_frames-1);
        if (beg_frame<1) { // should at least start from 1st frame
            beg_frame=1;
        }
        if (end_frame>=(blockframe-2)) { // avoid last frames
            end_frame=blockframe-2;
        }
    }
    
    /* System Block Definition */
    if (is_mbodies) {
        systemBlock_multibody_allmolecule(sysVar,ai);
        if (is_strings) make_sigmamatrix(sysVar);
    } else {
        systemBlock_definition(ai);
        if (is_strings) make_sigmamatrix(sysVar);
    }
    
    ai
    << "## peakframe  " << peakframe  << "\n"
    << "## blockframe " << blockframe << "\n"
    << "## beg_frame  " << beg_frame  << "\n"
    << "## end_frame  " << end_frame  << "\n";
    ai
    << "\n\n";
    
    if (true) // since amdat r128
    {
        //cout << "n_rouseModes "<<n_rouseModes<<"\n";
        for (int rouse=0;rouse<n_rouseModes;++rouse) {
            if (is_use_voroNeighbors) {
                if (is_mbodies) {
                    ai
                    << "create_voronoi_neighborlist vneigh -1\n"
                    << "list "<<mbodyList.at(rouse)
                    << "\n\n";
                    ai
                    << "value_statistics "
                    << "./statistics/voro/vneigh_"<< get_filenameString();
                    if(is_mbodies)ai<<"."<<mbodyList.at(rouse)<<".dat";
                    ai << " vneigh "<<n_moments<<"\n\n";
                    ai
                    << "neighbor_decorrelation_function "
                    << "./statistics/voro/ndf_"<< get_filenameString();
                    if(is_mbodies)ai<<"."<<mbodyList.at(rouse)<<".dat";
                    ai << " vneigh\n";
                    ai << "list "<<mbodyList.at(rouse)<<"\n\n";
                } else {
                    ai
                    << "create_voronoi_neighborlist vneigh -1\n"
                    << "list "+analysispart+"\n\n";
                    ai
                    << "value_statistics "
                    << "./statistics/voro/vneigh_"<< get_filenameString()
                    << "."+analysispart+".dat"
                    << " vneigh "<<n_moments<<"\n\n";
                    ai
                    << "neighbor_decorrelation_function "
                    << "./statistics/voro/ndf_"<< get_filenameString()
                    << "."+analysispart+".dat"
                    << " vneigh\n"
                    << "list "+analysispart+"\n\n";
                }
            }
            ai << "for frame "<< beg_frame << " " << end_frame+1 << "\n\n";
            write_amdat_loopstring(ai,sysVar,rouse);
            ai << "end\n\n";
        }
        
    } else {
        
        for (int frame_index=beg_frame; frame_index<=end_frame; ++frame_index)
        {
            /** OLD wrpper for amdat **/
            write_amdat_stringInput(ai,sysVar,frame_index); ai << "\n";
        }
    } ai.close();
    
    return {beg_frame,end_frame};
}





void AmdatAnalysis::make_amdatSubFile(const StructureClass& sysVar,
                                      const WorkScripts& ws,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d,
                                      const int frame)
{
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/AMDAT_submission_files");
    o.append("/AMDAT_");
    o.append(sysVar.get_usic());
    o.append("_00"+to_string((long long int)n_trl));
    o.append("_"+sysVar.get_nameString(n_sys));
    o.append("_T"+to_string((long long int)Temp_d));
    o.append(".qsub");
    //cout << o << "\n";
    
    ofstream as0(o.c_str());
    if(as0.is_open())	{
        
        as0
        << "#!/bin/bash"                                                << "\n"
        << "#$ -V"                                                      << "\n"
        << "#$ -cwd"                                                    << "\n"
        << "#$ -j y"                                                    << "\n"
        << "#$ -pe orte " << amdat_numCores                             << "\n"
        << "#$ -p " << amdat_priority                                   << "\n";
        
        streamsize ss=as0.precision();
        as0 << fixed << setprecision(0);
        
        as0
        << "#$ -hold_jid prd_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T"	<< Temp_d
        << "\n";
        
        as0
        << "#$ -N amdat_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << "\n";
        
        as0
        << "#$ -o ./AMDAT_submission_files/cluster_out/out_"
        << sysVar.get_usic()
        << "_00" << n_trl
        << "_" << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".amdat.o"
        << "\n";
        
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /** Compute node allocation **/
        //----------------------------------------------------------------------
        static int counter=0;
        node_allocation(sysVar,ws,counter,as0,n_trl);
        ++counter;
        
        /** AMDAT executable **/
        //----------------------------------------------------------------------
        as0 << "mpirun -np " << amdat_run_cores << " " << amdat_exe;
        
        as0 << fixed << setprecision(0);
        
        /** AMDAT inputfile path **/
        //----------------------------------------------------------------------
        if (is_strings) {
            if (is_peak_frame) {
                as0
                << " -i ./AMDAT_inputs/amdat_strings.inp";
            } else {
                as0
                << " -i ./AMDAT_inputs/amdat_strings_"
                << sysVar.get_usic()
                << "_00" << n_trl << "_"
                << sysVar.get_nameString(n_sys)
                << "_T" << Temp_d
                << ".inp";
            }
        } else {
            as0
            << " -i ./AMDAT_inputs/amdat.inp";
        }
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /** AMDAT temporary file **/
        //----------------------------------------------------------------------
        as0
        << fixed << setprecision(0)
        << " -t ./AMDAT_"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".tmp";
        
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        int n_poly_d  = sysVar.get_n_poly();
        int backLen_d = 0;
        int sideLen_d = 0;
        int numSide   = 0;
        
        if (sysVar.get_is_defaultLmpData()) {
            n_poly_d  = sysVar.get_n_polyVec()[n_sys];
            backLen_d = sysVar.get_backLenVec()[n_sys];
            sideLen_d = sysVar.get_sideLenVec()[n_sys];
            numSide=sideLen_d*backLen_d;
        }
        
        /** defaultly set in StructureClass **/
        //if (is_fixed_maxL) {maxLenScale=sysVar.get_maxLenScale()[n_sys];}
        //if (is_fixed_q) {waveindex=sysVar.get_waveindex()[n_sys];}
        
        string path_prdcustom;
        string backup="Regime_"+to_string(sysVar.get_current_regime());
        
        if (is_keep_new_trajectory) {
            if (sysVar.get_is_aging()) {
                path_prdcustom=
                return_SimulationFolderPath(sysVar)+"/production/trajectory/"+
                backup+"/new_trajectory_"+sysVar.get_usic()+
                "_00"+to_string((long long int)n_trl)+
                "_"+sysVar.get_nameString(n_sys)+
                "_T"+to_string((long long int)Temp_d)+
                ".prd.custom";
            } else {
                path_prdcustom=
                return_SimulationFolderPath(sysVar)+"/production/trajectory/"+
                "new_trajectory_"+sysVar.get_usic()+
                "_00"+to_string((long long int)n_trl)+
                "_"+sysVar.get_nameString(n_sys)+
                "_T"+to_string((long long int)Temp_d)+
                ".prd.custom";
            }
        } else {
            if (sysVar.get_is_aging()) {
                path_prdcustom=
                return_SimulationFolderPath(sysVar)+"/production/trajectory/"+
                backup+"/trajectory_"+sysVar.get_usic()+
                "_00"+to_string((long long int)n_trl)+
                "_"+sysVar.get_nameString(n_sys)+
                "_T"+to_string((long long int)Temp_d)+
                ".prd.custom";
            } else {
                path_prdcustom=
                return_SimulationFolderPath(sysVar)+"/production/trajectory/"+
                "trajectory_"+sysVar.get_usic()+
                "_00"+to_string((long long int)n_trl)+
                "_"+sysVar.get_nameString(n_sys)+
                "_T"+to_string((long long int)Temp_d)+
                ".prd.custom";
            }
        }
        
        
        /** input varaibles **/
        //======================================================================
        as0
        << " -c usic "           << sysVar.get_usic()
        << " -c trial 00"        << n_trl
        << " -c namestring "     << sysVar.get_nameString(n_sys)
        << " -c Regime "         << sysVar.get_current_regime()
        << " -c prd_exp_base "   << prd_exp_base
        << " -c n_prd_blocks "   << ws.get_n_prd_blocks()
        << " -c blocksize "      << ws.get_prd_blocksize()
        << " -c n_poly "         << n_poly_d
        << " -c maxLenScale "    << maxLenScale
        << " -c waveindex "      << waveindex
        << " -c path_prdcustom " << path_prdcustom;
        
        if (sysVar.get_is_defaultLmpData()) {
            as0
            << " -c backLen "       << backLen_d
            << " -c numSide "       << numSide;
        }
        if (sysVar.get_n_heavy()>0 ||
            sysVar.get_n_light()>0) {
            as0
            << " -c n_heavy " << sysVar.get_n_heavy()
            << " -c n_light " << sysVar.get_n_light();
        }
        if (is_strings) {
            if (is_peak_frame) {
                as0 << " -c frame " << frame;
            }
        }
        as0
        << " -c ts "            << ws.get_timestep_size()
        << fixed << setprecision(0)
        << " -c Temp "          << Temp_d;
        //======================================================================
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        /* to find the first waveindex from strFac file */
        //if (is_isfs||is_isf) {
        //    as0
        //    << " -c waveindex ";
        //    if (get_is_fixed_q()) {
        //        as0 << sysVar.get_waveindex()[n_sys];}
        //    else {
        //        as0 << find_waveindex(sysVar,n_trl,n_sys,Temp_d);}
        //}
        
        /** screen file path **/
        as0 << fixed << setprecision(0);
        //----------------------------------------------------------------------
        as0
        << " > ./screen/"
        << sysVar.get_usic()
        << "_00" << n_trl << "_"
        << sysVar.get_nameString(n_sys)
        << "_T" << Temp_d
        << ".amdat.screen"
        << "\n";
        
        as0.precision(ss);
        as0 << resetiosflags(ios::fixed|ios::showpoint);
        
        as0.close();
    }
    else cout << "amdatSubFile: 'amdat.sh' cannot open." << "\n";
    
}





int AmdatAnalysis::find_waveindex(const StructureClass& sysVar,
                                  const double Temp_d)
{
    /** find wavenumber used overall structure factor information **/
    string analysispart="all";
    
    /** number of discarded points due to NPT noise **/
    int n_discard=5;
    
    //FitData fitobj;
    
    for (int n_trl=0; n_trl<sysVar.get_n_trial(); ++n_trl) {
        for (int n_sys=sysVar.get_n_sys_beg(); n_sys<=sysVar.get_n_sys_end(); ++n_sys) {
            
            string in;
            in.append(return_AnalysisFolderPath(sysVar));
            in.append("/statistics/strFac/strFac_");
            in.append(sysVar.get_usic());
            in.append("_00"+to_string((long long int)n_trl));
            in.append("_"+sysVar.get_nameString(n_sys));
            in.append("_T"+to_string((long long int)Temp_d));
            in.append("."+analysispart+".dat");
            
            string out;
            out.append(return_AnalysisFolderPath(sysVar));
            out.append("/statistics/strFac/strFac_");
            out.append(sysVar.get_usic());
            out.append("_00"+to_string((long long int)n_trl));
            out.append("_"+sysVar.get_nameString(n_sys));
            out.append("_T"+to_string((long long int)Temp_d));
            out.append("."+analysispart+".interpolant");
            
            vector<vector<double>> interpolantdata;
            
            ifstream readfile(out.c_str());
            if (readfile.is_open()) {
                string lineContent;
                double dubVar=0;
                double wavetmp=0;
                double strFactmp=0;
                /* trash lines */
                getline(readfile,lineContent);
                getline(readfile,lineContent);
                /* run over data points to discard */
                for (indexi=0; indexi<n_discard; ++indexi) {
                    getline(readfile,lineContent);
                }
                /** Format of strFac file (amdat r89)
                 ** 1st column: wavenumber, q;
                 ** 7th column: structure factor, S(q) **/
                while (getline(readfile,lineContent)) {
                    istringstream iss(lineContent);
                    iss >> wavetmp; // 1st
                    iss >> dubVar;
                    iss >> dubVar;
                    iss >> dubVar;
                    iss >> dubVar;
                    iss >> dubVar;
                    iss >> strFactmp; // 7th
                    interpolantdata.push_back({wavetmp,strFactmp});
                } readfile.close();
                //fitobj.fitdata_processing_1d(interpolantdata);
                //real_1d_array x_in=fitobj.fit_xData.c_str();
                //real_1d_array y_in=fitobj.fit_yData.c_str();
                //alglib::spline1dinterpolant s;
                //alglib::spline1dbuildcubic(x_in,y_in,s);
            } else {
                cout
                << "in AmdatAnalysis::find_waveindex():\n"
                << in << " cannot open.\n";
                exit(EXIT_FAILURE);
            }
            
        }
    } return waveindex;
}





int AmdatAnalysis::find_peak_ngp_frame(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d)
{
    string analysispart="all";
    
    string in;
    in.append(return_AnalysisFolderPath(sysVar));
    in.append("/statistics/ngp/ngp_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append("."+analysispart+".dat");
    
    double time_threshold=0; // in-block time cutoff
    
    const int current_regime=sysVar.get_current_regime();
    double previous_teq,current_teq,teq;
    if (current_regime>0) {
        previous_teq = sysVar.get_equilibration_times().at(current_regime-1);
        current_teq  = sysVar.get_equilibration_times().at(current_regime);
        teq=max(previous_teq,current_teq);
    } else {
        teq=sysVar.get_equilibration_times().at(current_regime);
    } time_threshold=teq/sysVar.get_n_equ_blocks();
    
    int frame=0;
    vector<double> data_decreasing;
    ifstream readFile(in.c_str());
    if (readFile.is_open())
    {
        string lineContent;
        vector<double> time;
        vector<double> ngp;
        double timetmp=0,ngptmp=0;
        
        /* This is the first trash line in file */
        getline(readFile,lineContent);
        
        while (getline(readFile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> timetmp; time.push_back(timetmp); // 1st column: time
            iss >> ngptmp;  ngp.push_back(ngptmp);   // 2nd column: ngp
            
            if (ngp.size()>=3) {
                int    midindex=(int)ngp.size()-2;//0-based
                double mid=ngp.at(midindex);
                double forwdiff=ngp.at(midindex)-ngp.at(midindex+1);
                double backdiff=ngp.at(midindex)-ngp.at(midindex-1);
                
                if ((mid>time_threshold)&&(forwdiff>0)&&(backdiff>0)) {
                    frame = midindex;//0-based
                    break;
                }
            }
        } readFile.close();
    } else {
        cout
        << "in AmdatAnalysis::find_peak_ngp_frame():\n"
        << in << " file cannot open." << "\n";
    } return frame;
}





void AmdatAnalysis::make_amdatSubScript(const StructureClass& sysVar,
                                        const std::vector<std::vector<double>> tinfo) const
{
    const int n_trial  = sysVar.get_n_trial();
    const int n_system = sysVar.get_n_system();
    
    string o;
    o.append(return_AnalysisFolderPath(sysVar));
    o.append("/submission_scripts");
    o.append("/amdatSub_");
    o.append(sysVar.get_usic());
    o.append(".sh");
    ofstream as1(o.c_str());
    
    as1 << "#!/bin/bash" << "\n\n";
    
    /* Trial */
    //--------------------------------------------------------------------------
    as1 << "trial=(";
    for (int i=0; i<n_trial; ++i) {
        as1 << i;
        if (i!=(n_trial-1)) as1 << " ";
    } as1 << ")\n";
    
    /* namestring */
    //--------------------------------------------------------------------------
    as1 << "namestring=(";
    for (int i=sysVar.get_n_sys_beg(); i<=sysVar.get_n_sys_end(); ++i) {
        as1 << sysVar.get_nameString(i);
        if (i!=sysVar.get_n_sys_end()) as1 << " ";
    } as1 << ")\n";
    
    /* Temperature */
    //--------------------------------------------------------------------------
    streamsize ss=as1.precision();
    as1 << fixed << setprecision(0);
    as1 << "Temp=(";
    for (int n_trl=0; n_trl<n_trial; ++n_trl) {
        for (int n_sys=0; n_sys<n_system; ++n_sys) {
            int index=n_trl*n_system+n_sys;
            int n_Temp=(int)tinfo[index].size();
            as1 << "\"";
            for (int T=0; T<n_Temp; ++T) {
                as1
                << tinfo[index][T];
                if (T!=(n_Temp-1)) as1 << " ";
            }
            if (n_trl!=(n_trial-1)) as1 << "\"" << " ";
        }
    } as1 << "\")\n\n";
    as1.precision(ss);
    as1 << resetiosflags(ios::fixed|ios::showpoint);
    
    
    /* make 'simulations' folder the current working directory */
    //--------------------------------------------------------------------------
    as1
    << "cd "
    << return_AnalysisFolderPath(sysVar)                            << "\n";
    
    as1 << "\n";
    
    as1
    << "for i in ${trial[@]}; do"                                   << "\n"
    << "  for ii in ${namestring}; do"                              << "\n"
    << "    for iii in ${Temp[$i]}; do"                             << "\n"
    << "      qsub_File=./AMDAT_submission_files/AMDAT_"
    << sysVar.get_usic()
    << "_00${i}"
    << "_${ii}"
    << "_T${iii}"
    << ".qsub"                                                      << "\n"
    
    << "      test -e ${qsub_File}"                                 << "\n"
    << "      if [ \"$?\" -eq \"0\" ]; then"                        << "\n"
    << "        qsub ${qsub_File}"                                  << "\n"
    << "      else"                                                 << "\n"
    << "        continue"                                           << "\n"
    << "      fi"                                                   << "\n"
    << "    done"                                                   << "\n"
    << "  done"                                                     << "\n"
    << "done"                                                       << "\n";
    //--------------------------------------------------------------------------
    as1.close();
    
    
    /** Job submission **/
    //--------------------------------------------------------------------------
    if (sysVar.get_is_directSub()) {
        string sub;
        sub.append("bash ");
        sub.append(return_AnalysisFolderPath(sysVar));
        sub.append("/submission_scripts");
        sub.append("/amdatSub_");
        sub.append(sysVar.get_usic());
        sub.append(".sh");
        int ret=system(sub.c_str());
        if (ret!=0) {
            cout
            << "in AmdatAnalysis::make_amdatSubScript():\n"
            << "job submissoin failed.\n";
            exit(EXIT_FAILURE);
        }
    }
}





void AmdatAnalysis::change_atom_order(const StructureClass& sysVar,
                                      const int n_trl,
                                      const int n_sys,
                                      const double Temp_d)
{
    string in;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/trajectory/");
    in.append("trajectory_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.custom");
    ifstream readFile(in.c_str());
    
    string out;
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/trajectory/");
    out.append("newOrder_");//NOTE
    out.append("trajectory_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.custom");
    ofstream writeFile(out.c_str());
    
    int n_poly=32;
    int backlen=30;
    int back_per_monomer=1;
    int side_per_monomer=3;
    int monomersize=back_per_monomer+side_per_monomer;
    
    int n_columns_of_data=7; // # of columns of data to be stored
    
    if (readFile.is_open())
    {
        bool   is_lineForAtomNum=false;
        int    countStepkeyword=0;
        int    n_atoms=0;
        int    line_index=0;
        string lineContent;
        string strData0,strData1;
        
        vector<double> vecData;
        vector<vector<double>> trajectoryinfo;
        
        while (getline(readFile,lineContent))
        {
            ++line_index;
            istringstream iss(lineContent);
            
            // write out "type x y z ix iy iz" in custom file
            //------------------------------------------------------------------
            if ((countStepkeyword!=0)&&(countStepkeyword<=n_atoms)) {
                
                for (int i=0; i<n_columns_of_data; ++i) {
                    vecData.push_back(0);
                    iss >> vecData[i];
                }
                trajectoryinfo.push_back(vecData);
                
                /* NOTE: clean container content after each use */
                vecData.clear();
                
                // when all trajectories are read into storing vector
                // write out the intended format to file
                if (countStepkeyword==n_atoms) {
                    
                    for (int i=0; i<n_poly; ++i) {
                        
                        int offset_index=i*(backlen*monomersize);
                        
                        for (int ii=offset_index; ii<(offset_index+backlen); ++ii) {
                            
                            /* backbone bead */
                            int back_index=ii;
                            
                            for (int j=0; j<back_per_monomer; ++j) {
                                for (int iii=0; iii<trajectoryinfo[back_index+j].size(); ++iii) {
                                    writeFile
                                    << trajectoryinfo[back_index+j][iii] << " ";
                                } writeFile << "\n";
                            }
                            
                            /* side group beads */
                            int sum_sides=side_per_monomer*(back_index-offset_index);
                            int side_index=(offset_index+backlen)+sum_sides;
                            
                            for (int j=0; j<side_per_monomer; ++j) {
                                for (int iii=0; iii<trajectoryinfo[side_index+j].size(); ++iii) {
                                    writeFile
                                    << trajectoryinfo[side_index+j][iii] << " ";
                                } writeFile << "\n";
                            }
                        }
                    }
                    /* NOTE: clean container content after each use */
                    trajectoryinfo.clear();
                }
            } else {
                writeFile << lineContent << "\n";
            }
            //------------------------------------------------------------------
            
            iss >> strData0;
            
            if (is_lineForAtomNum) {
                n_atoms=stoi(strData0); // use # of atoms defined in file
                //cout << n_atoms << "\n";
                is_lineForAtomNum=false;
            }
            
            if (strData0=="ITEM:") {
                iss >> strData1;
                if (strData1=="NUMBER") {
                    is_lineForAtomNum=true;
                    //cout << line_index << "\n";
                    
                } else if (strData1=="ATOMS") {
                    // NOTE:
                    // if the read-in data type (ex. double) is not the same as
                    // the receiving variable type (ex. string), then the
                    // receiving varibale will not store the read-in data;
                    // therefore, the memory of the receiving variable would
                    // remain as the previously stored value where the read-in
                    // data type is the same as the receiving varaible.
                    
                    ++countStepkeyword;
                    //cout << line_index << "\n";
                } else if (strData1=="TIMESTEP") {
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
    } else {
        cout
        << "in AmdatAnalysis::change_atom_order():\n"
        << in << " file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::change_productionAtomTypes(const StructureClass& sysVar,
                                               const int n_trl,
                                               const int n_sys,
                                               const double Temp_d)
{
    /** Modified on 20160127 (SeanJH Hung)
     *  old style was to create a "new" trajectory where new atom_type infos are
     *  put; new style is changing the atom_type in the original trajectory file
     *  and writing new atom_types to a temporary fle to new atom_types, while
     *  keeping all other parts unchanged. The original trajectory then is
     *  replaced by the temporary file while keeping the filename the same.
     *  To recover the original atom_type information, use the function
     *  recover_productionAtomTypes() in the same class.
     *  This changing atom_type process is for the convenience of doing amdat
     *  analyses where if we use original atom_types, the generic scripting
     *  design would be much more difficult than using the current work around.
     **/
    
    string in;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/trajectory/");
    in.append("trajectory_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.custom");
    ifstream readFile(in.c_str());
    
    string out;
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/trajectory/");
    out.append("new_");//NOTE
    out.append("trajectory_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.custom");
    ofstream writeFile(out.c_str());
    
    vector<int> backbone  = sysVar.get_backboneTypes();
    vector<int> sideGroup = sysVar.get_sideGroupTypes();
    vector<int>::iterator itr0,itr1;
    
    int n_columns_of_data=7; //# of columns in custom file
    
    if (readFile.is_open())
    {
        bool   is_lineForAtomNum=false;
        int    intData0=0;
        int    countStepkeyword=0;
        int    n_atoms=0;
        int    line_index=0;
        string lineContent;
        string strData0,strData1,strData2,strData3;
        
        vector<double> vecData; // used for storing line content
        
        while (getline(readFile,lineContent))
        {
            ++line_index;
            istringstream iss(lineContent);
            
            // write out "type x y z ix iy iz" in custom file
            //------------------------------------------------------------------
            if ((countStepkeyword!=0)&&(countStepkeyword<=n_atoms)) {
                
                for (int i=0; i<n_columns_of_data; ++i) {
                    vecData.push_back(0);
                    iss >> vecData[i];
                }
                
                intData0=(int)vecData[0]; /** this is the atom type **/
                
                /** if found matching (backbone) type, change it to the
                 *  assigned type number **/
                itr0=find(backbone.begin(),backbone.end(),intData0);
                if (itr0!=backbone.end()) {
                    vecData[0]=1;
                }
                
                /** if found matching (sideGroup) type, change it to the
                 *  assigned type number **/
                itr1=find(sideGroup.begin(),sideGroup.end(),intData0);
                if (itr1!=sideGroup.end()) {
                    vecData[0]=2;
                }
                
                // write rest of content in the same row
                for (int i=0; i<n_columns_of_data; ++i) {
                    writeFile << vecData[i] << " ";
                } writeFile << "\n";
                
                // clean container content after each use
                vecData.clear();
                
            } else {
                
                writeFile << lineContent << "\n";
            }
            //------------------------------------------------------------------
            
            iss >> strData0;
            
            if (is_lineForAtomNum) {
                n_atoms=stoi(strData0); // use # of atoms defined in file
                //cout << n_atoms << "\n";
                is_lineForAtomNum=false;
            }
            
            if (strData0=="ITEM:") {
                
                iss >> strData1;
                
                if (strData1=="NUMBER") {
                    is_lineForAtomNum=true;
                    //cout << line_index << "\n";
                } else if (strData1=="ATOMS") {
                    ++countStepkeyword;
                    //cout << line_index << "\n";
                } else if (strData1=="TIMESTEP") {
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
    } else {
        cout
        << "in AmdatAnalysis::change_productionAtomTypes():\n"
        << in << " file cannot open.\n";
        exit(EXIT_FAILURE);
    }
    
    if (!is_keep_new_trajectory) {
        string bash="rm "+in+";"+"mv "+out+" "+in;
        system(bash.c_str());
    }
}





void AmdatAnalysis::recover_productionAtomTypes(const StructureClass& sysVar,
                                                const int n_trl,
                                                const int n_sys,
                                                const double Temp_d)
{
    string in;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/trajectory/");
    in.append("trajectory_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.custom");
    ifstream readFile(in.c_str());
    
    string out;
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/trajectory/");
    out.append("recovered_");//NOTE
    out.append("trajectory_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.custom");
    ofstream writeFile(out.c_str());
    
    /** path to the data file
     ** the idea behind this recovery function is that the order of atoms in
     ** the sorted custom file is the same as that in the data file; so we can
     ** use the atom type info in data file to recover the modified ones in
     ** the custom file. **/
    string data;
    data.append(return_SimulationFolderPath(sysVar));
    data.append("/lammps_inputs/start_data/");
    data.append("input");
    data.append("_00"+to_string((long long int)n_trl));
    data.append(".data");
    ifstream datafile(data.c_str());
    
    bool   is_inblock=false;
    int    intvar=0;
    int    n_total=0;
    string lineContent;
    string strvar,strvar1;
    
    vector<int> atomtypeVec;
    
    /** store atom type info in container **/
    if (datafile.is_open())
    {
        string strvar,str_now,str_pre;
        while (getline(datafile,lineContent))
        {
            istringstream iss(lineContent);
            iss >> strvar;
            /** get the total number of atom info from data file **/
            if (strvar=="LAMMPS") {
                while(getline(datafile,lineContent)) {
                    istringstream iss1(lineContent);
                    iss1 >> strvar;
                    iss1 >> strvar1;
                    if (strvar1=="atoms") {
                        n_total=atoi(strvar.c_str());
                        cout << "n_total " << n_total << "\n";
                        break;
                    }
                }
            }
            /** get inside the Atoms section and collect original atom types **/
            str_now=strvar;
            if (strvar=="Atoms") {
                is_inblock=true;
                continue;
            } else if (strvar=="Bonds") {
                break;
            }
            if (is_inblock) {
                if (str_now==str_pre) break;
                iss >> intvar; // molid
                iss >> intvar; // atomType
                atomtypeVec.push_back(intvar);
            } str_pre=str_now;
        } datafile.close();
        
    } else {
        cout
        << "in AmdatAnalysis::recover_productionAtomTypes():\n"
        << data << " data file cannot open.\n";
        exit(EXIT_FAILURE);
    }
    
    if (n_total!=(int)atomtypeVec.size())
    {
        cout << "\n"
        << "in AmdatAnalysis::recover_productionAtomTypes():" << "\n"
        << "read-in number of atoms from data file = " << n_total << "\n"
        << "size of container carrying atom types  = " << atomtypeVec.size() << "\n"
        << "Please check and make sure they are equal." << "\n"
        << "Program aborted.\n\n"; exit(EXIT_FAILURE);
    }
    
    int countStepkeyword=0,n_atoms=0,line_index=0;
    int n_columns_of_data=7; // # of columns in custom file
    string strData0,strData1,strData2,strData3;
    vector<double> vecData;
    bool is_lineForAtomNum=false;
    
    /** recover the atom type in the modified custom file using the stored
     atom types from the data file **/
    if (readFile.is_open())
    {
        lineContent.clear();
        while (getline(readFile,lineContent))
        {
            ++line_index;
            istringstream iss(lineContent);
            if ((countStepkeyword!=0)&&(countStepkeyword<=n_atoms)) {
                for (indexi=0; indexi<n_columns_of_data; ++indexi) {
                    vecData.push_back(0);
                    iss >> vecData.at(indexi);
                }
                vecData.at(0)=(double)atomtypeVec.at(countStepkeyword-1);
                /** write line content to the recovered file **/
                for (indexi=0; indexi<n_columns_of_data; ++indexi) {
                    writeFile << vecData.at(indexi) << " ";
                } writeFile << "\n";
                vecData.clear();
            } else {
                writeFile << lineContent << "\n";
            }
            /** NOTE:
             **  the line below is where the first iss read-in is put, in order
             **  to carefully trigger the proper handling section **/
            iss >> strData0;
            if (is_lineForAtomNum) {
                n_atoms=stoi(strData0); // use # of atoms defined in file
                //cout << n_atoms << "\n";
                is_lineForAtomNum=false;
            }
            if (strData0=="ITEM:") {
                iss >> strData1;
                if (strData1=="NUMBER") {
                    is_lineForAtomNum=true;
                    //cout << line_index << "\n";
                } else if (strData1=="ATOMS") {
                    ++countStepkeyword;
                    //cout << line_index << "\n";
                } else if (strData1=="TIMESTEP") {
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
        } readFile.close(); writeFile.close();
    } else {
        cout
        << "in AmdatAnalysis::recover_productionAtomTypes():\n"
        << "recover_productionAtomTypes: custom trajectory file cannot open.\n";
        exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::change_prdcustomfilesteps(const StructureClass& sysVar,
                                              const int n_trl,
                                              const int n_sys,
                                              const double Temp_d)
{
    string in;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/trajectory/");
    in.append("trajectory_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.custom");
    ifstream readFile(in.c_str());
    
    string out;
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/trajectory/");
    out.append("newStep_");//NOTE
    out.append("trajectory_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.custom");
    ofstream writeFile(out.c_str());
    
    if (readFile.is_open())
    {
        int steps_zero=0;
        int steps_corr=0;
        int countStepkeyword=0;
        int line_index=0;
        int block_index=0;
        int frame_index=0;
        int countaccess=0;
        
        string lineContent;
        string strData0,strData1;
        
        while (getline(readFile,lineContent))
        {
            ++line_index;
            istringstream iss(lineContent);
            
            // NOTE:
            // the very first frame is at timestep=0;
            // in the using exponential timestepping scheme,
            // 0th block has N+1 frames, 1st-final blocks
            // each has N frames.
            if ((countStepkeyword==0)&&(countaccess>1)) {
                
                if (block_index>0) {
                    iss >> strData0;
                    steps_zero=stoi(strData0);
                    steps_corr=steps_zero+block_index*(int)pow(prd_exp_base,prd_blocksize-1);
                    writeFile << steps_corr << "\n";
                } else {
                    writeFile << lineContent << "\n";
                }
                
                ++frame_index;
                if (frame_index==prd_blocksize) {
                    ++block_index;
                    frame_index=0;
                }
            } else {
                writeFile << lineContent << "\n";
            }
            
            ++countStepkeyword;
            
            iss >> strData0;
            
            if (strData0=="ITEM:") {
                
                iss >> strData1;
                
                if (strData1=="TIMESTEP") {
                    ++countaccess;
                    countStepkeyword=0;
                    //cout << line_index << "\n";
                }
            }
        } readFile.close(); writeFile.close();
        
    } else {
        cout
        << "in AmdatAnalysis::change_prdcustomfilesteps():\n"
        << in << " file cannot open.\n"; exit(EXIT_FAILURE);
    }
}





void AmdatAnalysis::delete_zerothFrames(const StructureClass& sysVar,
                                        const int n_trl,
                                        const int n_sys,
                                        const double Temp_d)
{
    string in;
    in.append(return_SimulationFolderPath(sysVar));
    in.append("/production/trajectory/");
    in.append("trajectory_");
    in.append(sysVar.get_usic());
    in.append("_00"+to_string((long long int)n_trl));
    in.append("_"+sysVar.get_nameString(n_sys));
    in.append("_T"+to_string((long long int)Temp_d));
    in.append(".prd.custom");
    ifstream readFile(in.c_str());
    
    string out;
    out.append(return_SimulationFolderPath(sysVar));
    out.append("/production/trajectory/");
    out.append("new_");//NOTE
    out.append("trajectory_");
    out.append(sysVar.get_usic());
    out.append("_00"+to_string((long long int)n_trl));
    out.append("_"+sysVar.get_nameString(n_sys));
    out.append("_T"+to_string((long long int)Temp_d));
    out.append(".prd.custom");
    ofstream writeFile(out.c_str());
    
    if (readFile.is_open())
    {
        bool is_readinframe=true;
        int  line_index=0;
        int  countaccess=0;
        string lineContent,lineContent1;
        string strData0,strData1;
        
        // Read every line one at a time of the file
        while (getline(readFile,lineContent))
        {
            ++line_index;
            istringstream iss(lineContent);
            iss >> strData0;
            if (strData0=="ITEM:") {
                iss >> strData1;
                if (strData1=="TIMESTEP") {
                    ++countaccess;
                    //cout << line_index << "\n";
                    getline(readFile,lineContent1);
                    /** only delete excess 0th frames **/
                    if ((stoi(lineContent1)==0)&&(countaccess>1)) {
                        is_readinframe=false;
                    } else {
                        is_readinframe=true;
                    }
                    if (is_readinframe) {
                        writeFile
                        << lineContent  << "\n"
                        << lineContent1 << "\n";
                    }
                } else {
                    if (is_readinframe) writeFile << lineContent << "\n";
                }
            } else {
                if (is_readinframe) writeFile << lineContent << "\n";
            }
        } readFile.close();
        
    } else {
        cout
        << "in AmdatAnalysis::delete_zerothFrames():\n"
        << in << " cannot open." << "\n";
        exit(EXIT_FAILURE);
    } writeFile.close();
    
    /** mv tmp file name to original **/
    string bash="rm "+in+";"+"mv "+out+" "+in;
    system(bash.c_str());
}





const string AmdatAnalysis::amdat_target(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys,
                                         const double Temp_d,
                                         const string& phase) const
{
    string target;
    target.append(return_AnalysisFolderPath(sysVar));
    target.append("/statistics");
    target.append("/"+phase);
    target.append("/"+phase+"_");
    target.append(sysVar.get_usic());
    target.append("_00"+to_string((long long int)n_trl));
    target.append("_"+sysVar.get_nameString(n_sys));
    target.append("_T"+to_string((long long int)Temp_d));
    target.append("."+analysispart+".dat");
    return target;
}
const string AmdatAnalysis::amdat_target(const StructureClass& sysVar,
                                         const int n_trl,
                                         const int n_sys,
                                         const double Temp_d,
                                         const string& phase,
                                         const int frame) const
{
    string target;
    target.append(return_AnalysisFolderPath(sysVar));
    target.append("/statistics");
    target.append("/"+phase);
    target.append("/"+phase+"_");
    target.append(sysVar.get_usic());
    target.append("_00"+to_string((long long int)n_trl));
    target.append("_"+sysVar.get_nameString(n_sys));
    target.append("_T"+to_string((long long int)Temp_d));
    target.append("_frame"+to_string((long long int)frame)); // NOTE
    target.append("."+analysispart+".dat");
    return target;
}




/** public setters **/
//------------------------------------------------------------------------------
/* bool */
void AmdatAnalysis::set_is_changeAtomOrder(const bool b){is_changeAtomOrder=b;}
void AmdatAnalysis::set_is_changeAtomType(const bool b){is_changeAtomType=b;}
void AmdatAnalysis::set_is_recoverAtomType(const bool b){is_recoverAtomType=b;}
void AmdatAnalysis::set_is_changeFrameSteps(const bool b){is_changeFrameSteps=b;}
void AmdatAnalysis::set_is_deleteZerothFrames(const bool b){is_deleteZerothFrames=b;}
void AmdatAnalysis::set_is_NPT(const bool b){is_NPT=b;}
void AmdatAnalysis::set_is_fixed_q(const bool b){is_fixed_q=b;}
void AmdatAnalysis::set_is_fixed_maxL(const bool b){is_fixed_maxL=b;}
void AmdatAnalysis::set_is_auto_qvalue(const bool b){is_auto_qvalue=b;}
void AmdatAnalysis::set_is_strFac(const bool b){is_strFac=b;}
void AmdatAnalysis::set_is_rdf(const bool b){is_rdf=b;}
void AmdatAnalysis::set_is_msd(const bool b){is_msd=b;}
void AmdatAnalysis::set_is_ngp(const bool b){is_ngp=b;}
void AmdatAnalysis::set_is_isfs(const bool b){is_isfs=b;}
void AmdatAnalysis::set_is_baf(const bool b){is_baf=b;}
void AmdatAnalysis::set_is_composition(const bool b){is_composition=b;}
void AmdatAnalysis::set_is_u2dist(const bool b){is_u2dist=b;}
void AmdatAnalysis::set_is_stiffness_dist(const bool b){is_stiffness_dist=b;}
void AmdatAnalysis::set_is_isf(const bool b){is_isf=b;}
void AmdatAnalysis::set_is_mbodies(const bool b){is_mbodies=b;}
void AmdatAnalysis::set_is_strings(const bool b){is_strings=b;}
void AmdatAnalysis::set_is_peak_frame(const bool b){is_peak_frame=b;}
void AmdatAnalysis::set_is_keep_new_trajectory(const bool b){is_keep_new_trajectory=b;}
void AmdatAnalysis::set_is_streamlined_strings(const bool b){is_streamlined_strings=b;}
void AmdatAnalysis::set_is_monoStruct(const bool b){is_monoStruct=b;}
void AmdatAnalysis::set_is_use_voroNeighbors(const bool b){is_use_voroNeighbors=b;}
void AmdatAnalysis::set_is_binning(const bool b){is_binning=b;}
/* int */
void AmdatAnalysis::set_amdat_numCores(const int i){amdat_numCores=i;}
void AmdatAnalysis::set_amdat_run_cores(const int i){amdat_run_cores=i;}
void AmdatAnalysis::set_amdat_priority(const int i){amdat_priority=i;}
void AmdatAnalysis::set_fullblock(const int i){fullblock=i;}
/* double */
void AmdatAnalysis::set_logtaubeta(const double d)
{
    logtaubeta=d;
    DWF_time=pow(10,logtaubeta);
}
void AmdatAnalysis::set_strings_threshold(const double d){strings_threshold=d;}
/* string */
void AmdatAnalysis::set_amdat_exe(const string& s){amdat_exe=s;}
void AmdatAnalysis::set_symmetry(const std::string& str){symmetry=str;}
void AmdatAnalysis::set_geometry(const std::string& str){geometry=str;}
void AmdatAnalysis::set_relaxation_target(const std::string& str){relaxation_target=str;}
void AmdatAnalysis::set_centertype(const std::string& str){centertype=str;}
void AmdatAnalysis::set_analysispart(const std::string& str){analysispart=str;}
void AmdatAnalysis::set_analysispart_bin(const std::string& str){analysispart_bin=str;}
void AmdatAnalysis::set_segmode(const std::string& str){segmode=str;}
/* STL */
void AmdatAnalysis::set_n_typeSet(const vector<vector<int>>& vvi){n_typeSet=vvi;}




