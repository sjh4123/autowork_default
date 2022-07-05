//
//  moltemplatelmpdata.h
//  cppWork
//
//  Created by SJH on 12/18/15.
//  Copyright © 2015 Sean Jh H. All rights reserved.
//

#ifndef MOLTEMPLATELMPDATA_H
#define MOLTEMPLATELMPDATA_H

#include "structureclass.h"

namespace autoWork {
    
    class MoltemplateLmpData
    {
        int n_trials;
        int indexi,indexii;
        int sequenceNum;
        int sequenceLen;
        int n_light;
        int n_heavy;
        
        double moltemplateBoxSize;
        double offset;
        double rotate;
        
        std::string path_master;
        std::string path_moltemplatesrc;
        std::string path_oplsaaprm;
        std::string path_monomerBank;
        std::string path_cwd;
        
        std::vector<std::vector<std::string>> sequenceSet;
        std::vector<int> types_all;
        std::vector<int> types_light;
        std::vector<int> types_heavy;
        std::vector<std::vector<int>> n_types_all;
        std::vector<std::vector<int>> n_types_light;
        std::vector<std::vector<int>> n_types_heavy;
        
        /* internal utility functions */
        bool check_monomerbank(const std::string&);
        bool check_is_heavyAtom(const int);
        void check_path(const std::string&);
        void check_lasttwo(std::vector<int>&);
        void copy_to_cwd(const std::string&);
        void make_oplsaa_subset(const std::vector<std::string>&);
        void make_oplsaalt(const std::vector<std::string>&);
        void make_polylt(const int,const std::vector<std::string>&);
        void make_systemlt();
        void invoke_moltemplate();
        void optimize_charge();
        void mv_files();
        int  show_positinindex(const int, const std::vector<int>&);
        void initialize_n_typesVec();
        void finalize_n_typesVec(const int);
        void add_n_types_all(const int);
        void add_n_types_heavy(const int);
        void add_n_types_light(const int);
        
    public:
        
        /* Constructors */
        MoltemplateLmpData()=delete;
        MoltemplateLmpData(const StructureClass&);
        MoltemplateLmpData(const MoltemplateLmpData&)=default;
        /* Assignment */
        MoltemplateLmpData& operator= (const MoltemplateLmpData&)=default;
        /* Destructor */
        ~MoltemplateLmpData()=default;
        
        
        /** Make LAMMPS data file by Moltemplate **/
        void make_LmpDataFilebyMoltemplate(const StructureClass&);
        
        /** Set light and heavy atom types and numbers **/
        void set_atomTypes(StructureClass&);
        void set_atomNumbers(StructureClass&);
        
        /** public utility functions **/
        std::vector<std::vector<std::string>> read_sequenceSetfromfile(const StructureClass&);
        
        
        
        /*==( public setters )==*/
        /* int */
        void set_sequenceNum(const int);
        void set_sequenceLen(const int);
        /* double */
        void set_moltemplateBoxSize(const double);
        /* vector<vector<string>> */
        void set_sequenceSet(const std::vector<std::vector<std::string>>&);
    };
}
#endif /* moltemplatelmpdata_h */
