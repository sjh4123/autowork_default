//
//  Created by Sean Jh H. on 6/24/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef AMDATANALYSIS_H
#define AMDATANALYSIS_H

#include "structureclass.h"
#include "lmpscripts.h"
#include "workscripts.h"
#include <algorithm>

namespace autoWork {
    
    class AmdatAnalysis
    {
        bool is_changeAtomOrder;
        bool is_changeAtomType;
        bool is_changeFrameSteps;
        bool is_deleteZerothFrames;
        bool is_NPT;
        bool is_fixed_q;
        bool is_fixed_maxL;
        
        bool is_strFac;
        bool is_msd;
        bool is_ngp;
        bool is_isfs;
        
        bool is_composition;
        bool is_u2dist;
        bool is_stiffness_dist;
        bool is_isf;
        bool is_mbodies;
        
        bool is_strings;
        bool is_peak_frame;
        
        int amdat_cores;    // # of cores used for AMDAT jobs (def. = 1)
        int amdat_priority; // priority of AMDAT jobs (def. = 0)
        int fullblock;      // default = 0
        int vhs_nbins;
        int u2dist_nbins;
        int stiffness_nbins;
        
        int prd_blocksize;
        double prd_exp_base;
        double n_prd_blocks;
        double DWF_time;
        
        double vhs_rangebinned;
        double strings_threshold;
        double max_u2;
        double max_stiffness;
        
        std::string threshold_keyword; // greater,less,between
        double threshold_percentile;
        double threshold_percentile_lo;
        double threshold_percentile_hi;
        
        std::string amdat_exe;
        std::string symmetry;
        std::string geometry;
        std::string analysispart;
        std::string centertypename;
        
    public:
        
        /* Constructors */
        AmdatAnalysis()=delete;
        AmdatAnalysis(const StructureClass&,
                      const WorkScripts&);
        AmdatAnalysis(const AmdatAnalysis&)=default;
        
        /* Assignment */
        AmdatAnalysis& operator= (const AmdatAnalysis&)=default;
        
        /* Destructor */
        ~AmdatAnalysis()=default;
        
        
        //===============================================
        // make AMDAT inputFile
        //===============================================
        void systemBlock_definition(std::ofstream& ai,
                                    const std::string& path);
        
        void make_amdatInputFile(const StructureClass&,
                                 const WorkScripts&);
        
        void make_amdatInputFile_strings_loop(const StructureClass&,
                                              const WorkScripts&,
                                              const int blockframe);
        
        std::vector<int> make_amdatInputFile_strings_loop(const StructureClass&,
                                                          const WorkScripts&,
                                                          const int n_trl,
                                                          const int n_sys,
                                                          const double T,
                                                          const int peakframe,
                                                          const int blockframe);
        
        //===============================================
        // make AMDAT submission files
        //===============================================
        void make_amdatSubFile(const StructureClass&,
                               const WorkScripts&,
                               const int n_trl,
                               const int n_sys,
                               const double T,
                               const int frame=0);
        
        
        /* find the first waveindex from structure factor file */
        int find_waveindex(const StructureClass&,
                           const int n_trl,
                           const int n_sys,
                           const double T);
        
        
        int find_ngp_frame(const StructureClass&,
                           const int n_trl,
                           const int n_sys,
                           const double T);
        
        
        //===============================================
        // make AMDAT job submission scripts
        //===============================================
        void make_amdatSubScript(const StructureClass&,
                                 const std::vector<std::vector<double>>) const;
        
        
        //===============================================
        // Echo the paths of AMDAT target files
        //===============================================
        const std::string amdat_target(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d,
                                       const std::string& phase) const;
        
        const std::string amdat_target(const StructureClass& sysVar,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d,
                                       const std::string& phase,
                                       const int frame) const;
        
        /* change the order of atom types in the custom trajectory file */
        void change_atom_order(const StructureClass&,
                               const int n_trl,
                               const int n_sys,
                               const double Temp_d);
        
        /* change atom types in the custom trajectory file */
        void change_productionAtomTypes(const StructureClass&,
                                        const int n_trl,
                                        const int n_sys,
                                        const double Temp_d);
        
        /* change frame timesteps in the custom trajectory file */
        void change_prdcustomfilesteps(const StructureClass&,
                                       const int n_trl,
                                       const int n_sys,
                                       const double Temp_d);
        
        /* change frame timesteps in the custom trajectory file */
        void delete_zerothFrames(const StructureClass&,
                                 const int n_trl,
                                 const int n_sys,
                                 const double Temp_d);
        
        /*==( public setters )==*/
        
        /* bool */
        void set_is_changeAtomOrder(const bool);
        void set_is_changeAtomType(const bool);
        void set_is_changeFrameSteps(const bool);
        void set_is_deleteZerothFrames(const bool);
        void set_is_NPT(const bool);
        void set_is_fixed_q(const bool);
        void set_is_fixed_maxL(const bool);
        void set_is_strFac(const bool);
        void set_is_msd(const bool);
        void set_is_ngp(const bool);
        void set_is_isfs(const bool);
        void set_is_composition(const bool);
        void set_is_u2dist(const bool);
        void set_is_stiffness_dist(const bool);
        void set_is_isf(const bool);
        void set_is_mbodies(const bool);
        void set_is_strings(const bool);
        void set_is_peak_frame(const bool);
        
        /* int */
        void set_amdat_cores(const int);
        void set_amdat_priority(const int);
        void set_fullblock(const int);
        void set_vhs_nbins(const int);
        
        /* double */
        void set_vhs_rangebinned(const double);
        void set_strings_threshold(const double);
        
        /* string */
        void set_amdat_exe(const std::string&);
        void set_symmetry(const std::string&);
        void set_geometry(const std::string&);
        void set_analysispart(const std::string&);
        void set_centertypename(const std::string&);
        
        
        
        /*==( public getters )==*/
        
        /* bool */
        bool get_is_changeAtomOrder() const {return is_changeAtomOrder;}
        bool get_is_changeAtomType() const {return is_changeAtomType;}
        bool get_is_changeFrameSteps() const {return is_changeFrameSteps;}
        bool get_is_deleteZerothFrames() const {return is_deleteZerothFrames;}
        bool get_is_NPT() const {return is_NPT;}
        bool get_is_fixed_q() const {return is_fixed_q;}
        bool get_is_fixed_maxL() const {return is_fixed_maxL;}
        bool get_is_strFac() const {return is_strFac;}
        bool get_is_msd() const {return is_msd;}
        bool get_is_ngp() const {return is_ngp;}
        bool get_is_isfs() const {return is_isfs;}
        bool get_is_composition() const {return is_composition;}
        bool get_is_u2dist() const {return is_u2dist;}
        bool get_is_stiffness_dist() const {return is_stiffness_dist;}
        bool get_is_isf() const {return is_isf;}
        bool get_is_mbodies() const {return is_mbodies;}
        bool get_is_strings() const {return is_strings;}
        bool get_is_peak_frame() const {return is_peak_frame;}
        
        /* int */
        int get_amdat_cores() const {return amdat_cores;}
        int get_amdat_priority() const {return amdat_priority;}
        int get_fullblock() const {return fullblock;}
        int get_vhs_nbins() const {return vhs_nbins;}
        
        /* double */
        double get_vhs_rangebinned() const {return vhs_rangebinned;}
        double get_strings_threshold() const {return strings_threshold;}
        double get_DWF_time() const {return DWF_time;}
        
        /* string */
        const std::string get_amdat_exe() const {return amdat_exe;}
        const std::string get_symmetry() const {return symmetry;}
        const std::string get_geometry() const {return geometry;}
        const std::string get_analysispart() const {return analysispart;}
        const std::string get_centertypename() const {return centertypename;}
    };
    
}
#endif
