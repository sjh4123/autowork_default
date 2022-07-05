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

#ifndef MOLTEMPLATELMPDATA_H
#define MOLTEMPLATELMPDATA_H

#include "structureclass.h"

namespace autoWork
{
    class MoltemplateLmpData
    {
        bool is_restrictbynAtoms;
        bool is_optimizeCharge;
        bool is_mass_normalization;
        bool is_correctAlkylDihedral;
        bool is_HCHexisting;
        
        int dop;
        int indexi;
        int indexii;
        int n_trials;
        int n_total_atoms;
        int sequenceNum;
        int sequenceLen;
        int n_light;
        int n_heavy;
        
        int n_HCHbondType;
        int n_totalatomTypes;
        int atomType_CH3;
        int atomType_CH2;
        int atomType_CH;
        int atomType_H;
        
        double FFmodify_CH3_charge;
        double FFmodify_CH2_charge;
        double FFmodify_CH_charge;
        double FFmodify_H_charge;
        double FFmodify_HCH_charge;
        
        double moltemplateBoxSize;
        double offset;
        double rotate;
        double packingL;
        double offset_spacing;
        double packingL_spacing;
        
        std::string path_cwd;
        std::string path_master;
        std::string path_moltemplatesrc;
        std::string path_oplsaaprm;
        std::string path_monomerBank;
        std::string path_sequenceBank;
        std::string path_scripts;
        std::string merA;
        std::string merB;
        std::string copolymerType;
        std::string tacticity;
        
        /** STL containers **/
        //----------------------------------------------------------------------
        std::vector<std::string> merSet;
        std::vector<std::string> chem_mono;
        std::vector<int> merComp;
        std::vector<int> mono_beads;
        std::vector<std::vector<std::string>> sequenceSet;
        std::vector<int> types_all;
        std::vector<int> types_light;
        std::vector<int> types_heavy;
        std::vector<int> HCHatomIDs;
        std::vector<double> FFmodify_alkylDihedral;
        std::vector<std::vector<int>> n_types_all;
        std::vector<std::vector<int>> n_types_light;
        std::vector<std::vector<int>> n_types_heavy;
        std::vector<std::vector<int>> HCHbondType;
        
        /* internal utility functions */
        int  show_positinindex(const int, const std::vector<int>&);
        int  n_monomerAtoms(const std::string&);
        bool check_monomerbank(const std::string&);
        bool check_is_heavyAtom(const int);
        void check_cwd(const StructureClass&);
        void check_path(const std::string&);
        void check_lasttwo(std::vector<int>&);
        void check_lasttwo(std::vector<double>&);
        void check_lasttwo(std::vector<std::vector<int>>&);
        void check_lasttwo(std::vector<std::vector<double>>&);
        void copy_to_cwd(const std::string&);
        void make_oplsaa_subset(const std::vector<std::string>&);
        void make_oplsaalt(const std::vector<std::string>&);
        void make_polylt(const int,const std::vector<std::string>&);
        void evaluate_offset(const std::string&);
        void evaluate_boxLen();
        void make_systemlt();
        void invoke_moltemplate();
        void optimize_charge();
        void mv_files();
        void initialize_n_typesVec();
        void finalize_n_typesVec(const int);
        void add_n_types_all(const int);
        void add_n_types_heavy(const int);
        void add_n_types_light(const int);
        void getridof_full();
        void getridof_ljcutcoullong();
        void mass_normalization();
        /* FF correction for long alkyl backbone */
        void FFmodify_alkyl_dihedral_oplsaa();
        void FFmodify_alkyl_charge_oplsaa();
        void FFmodify_alkyl_findHCHbond_oplsaa();
        void FFmodify_alkyl_findHCHinBonds_systemdata();
        void FFmodify_alkyl_writeHCHinMasses_systemdata();
        void FFmodify_alkyl_writeHCHinAtoms_systemdata();
        void FFmodify_alkyl_writeHCHincharges_systemcharges();
        void FFmodify_alkyl_writeHCHinsettings_systemsettings();
        
    public:
        
        MoltemplateLmpData()=default;
        MoltemplateLmpData(const StructureClass&);
        MoltemplateLmpData(const MoltemplateLmpData&)=default;
        MoltemplateLmpData& operator= (const MoltemplateLmpData&)=default;
        ~MoltemplateLmpData()=default;
        
        /** make sequence set where atacticity is applied **/
        std::vector<std::vector<std::string>> make_atacticSequenceSet();
        std::vector<std::vector<std::string>> make_SequenceSet_adhoc();
        
        /** Make LAMMPS data file by Moltemplate **/
        void make_LmpDataFilebyMoltemplate(StructureClass&);
        
        /** Set light and heavy atom types and numbers **/
        void set_atomTypes(StructureClass&);
        void set_atomNumbers(StructureClass&);
        
        /** utility function for reading source sequence set from file **/
        std::vector<std::vector<std::string>> read_sequenceSetfromfile(const StructureClass&);
        
        /** public setters **/
        //----------------------------------------------------------------------
        /* bool */
        void set_is_restrictbynAtoms(const bool);
        void set_is_optimizeCharge(const bool);
        void set_is_mass_normalization(const bool);
        void set_is_correctAlkylDihedral(const bool);
        /* int */
        void set_dop(const int);
        void set_n_total_atoms(const int);
        void set_sequenceNum(const int);
        void set_sequenceLen(const int);
        /* double */
        void set_moltemplateBoxSize(const double);
        /* string */
        void set_merA(const std::string&);
        void set_merB(const std::string&);
        void set_copolymerType(const std::string&);
        void set_tacticity(const std::string&);
        /* vector<string> */
        void set_merSet(const std::vector<std::string>&);
        void set_chem_mono(const std::vector<std::string>&);
        /* vector<int> */
        void set_merComp(const std::vector<int>&);
        void set_mono_beads(const std::vector<int>&);
        /* vector<vector<string>> */
        void set_sequenceSet(const std::vector<std::vector<std::string>>&);
        
        /** public getters **/
        //----------------------------------------------------------------------
        const std::string get_path_master() const {return path_master;}
        const std::string get_path_sequenceBank() const {return path_sequenceBank;}
        const std::string get_path_monomerBank() const {return path_monomerBank;}
        const std::string get_copolymerType() const {return copolymerType;}
        const std::string get_tacticity() const {return tacticity;}
        const std::vector<std::string>& get_merSet() const {return merSet;}
        const std::vector<int> get_merComp() const {return merComp;}
        bool get_is_restrictbynAtoms() const {return is_restrictbynAtoms;}
        int get_dop() const {return dop;}
        int get_n_total_atoms() const {return n_total_atoms;}
        int get_sequenceNum() const {return sequenceNum;}
        int get_sequenceLen() const {return sequenceLen;}
    };
}
#endif /* MOLTEMPLATELMPDATA_H */




