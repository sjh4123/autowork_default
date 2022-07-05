//
//  Created by Sean Jh H. on 6/3/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef STRUCTURECLASS_H
#define STRUCTURECLASS_H

#include "sysinfo.h"

namespace autoWork {
    
    class StructureClass: public SysInfo
    {
        std::string nameString;
        std::string backboneInWords;
        std::string sideGroupInWords;
        
        int typeB;         // type of backbonedouble
        int typeS;         // type of side group
        int	typeBackbone;  // 10^(typeB)
        int typeSideGroup; // 2^(typeS)
        
        int backLen;       // length of backbone
        int sideLen;       // number of beads per side group
        int n_poly;        // total number of polymer beads
        int numBack;       // total number of backbone beads
        int numSide;       // total number of sideGroup beads
        int ringBeads;     // number of beads per ring structure
        int polySize;      // total number of beads per polymer
        int polyBeads;     // total number of beads in the system
        int totalBeads;    //
        int	atomType;      // Nominal # of atom types
        int bondType;      // Nominal # of bond types
        int angleType;     // Nominal # of angle types
        
        // NOTE:
        // real # of atom,bond,angle types are stored in
        // AtomType,BondType,AngleTypes as local variables
        // in corresponding block;
        // and at the end of class, strore the real values into
        // trueAtoms,trueBonds,trueAngles.
        
        int trueAtoms;     // real # of atom types
        int trueBonds;     // real # of bond types
        int trueAngles;    // real # of angle types
        int AtomID;        //
        int BondID;        //
        int AngleID;       //
        double rcutl;      //
        double rcuth;      //
        double rcut;
        double fenebondK;
        double maxfenebondLen;
        double fenebond_eps;
        double fenebond_sig;
        
        std::vector<int> backLenVec;        //
        std::vector<int> sideLenVec;        //
        std::vector<int> n_polyVec;         //
        std::vector<int> waveindex;         //
        std::vector<int> backboneTypes;
        std::vector<int> sideGroupTypes;
        
        std::vector<double> boxL,pkmL;      //
        std::vector<double> mass;           //
        std::vector<double> bondLen;        //
        std::vector<double> bondCoeffs;     //
        std::vector<double> theta;          //
        std::vector<double> angCoeffs;      //
        std::vector<double> lj_eps;         //
        std::vector<double> lj_sig;         //
        std::vector<double> maxLenScale;    //
        std::vector<double> composition;
        
        std::vector<bool> is_backboneType;  //
        std::vector<bool> is_sidegroupType; //
        std::vector<bool> is_shakeBonds;    //
        std::vector<bool> is_shakeAngles;   //
        
        /** atom types from moltemplate **/
        int n_light;
        int n_heavy;
        std::vector<int> types_all;
        std::vector<int> types_light;
        std::vector<int> types_heavy;
        std::vector<std::vector<int>> n_types_all;
        std::vector<std::vector<int>> n_types_light;
        std::vector<std::vector<int>> n_types_heavy;
        
    public:
        
        /* Constructors */
        StructureClass()=delete; // must have 2 arguments
        StructureClass(const int,const int); // args: (typeB, typeS)
        StructureClass(const StructureClass&)=default;
        /* Assignment */
        StructureClass& operator= (const StructureClass&)=default;
        /* Destructor */
        ~StructureClass()=default;
        
        
        // Unique Numbering of (Backbone + SideGroup)
        //============================================
        // Molecular Structure Definition:
        // The backbones are represented as constant
        // integers with a base of 10;
        // The side groups are with a base of 2;
        // This basically ensures that the additions of
        // the backbone and side group constants
        // are unique. (For use of 'Inter' relations)
        // There's a function in 'main.cpp' that checks
        // the repetition of numbers based on the
        // two base numbers given:
        // 'void testRepeat(const int&,const int&)'
        //============================================
        enum typeB_pool {
            
            linear	= 1,
            typeB1	= 10,
            typeB2	= 100,
            typeB3	= 1000,
            typeB4	= 10000,
            LJliqB  = 100000
        };
        enum typeS_pool {
            
            phenyl	= 1,
            alkyl	= 2,
            tButyl	= 4,
            clProp	= 8,
            C3H8O	= 16,
            stdFENE = 32,
            LJliqS  = 64
        };
        
        
        /* Setup parameters relevant to the structure */
        void set_structureParam();
        
        /* Generate single chain XYZ file */
        void make_SingleChainXYZ(const int);
        
        /* Generate Bond file to the structure */
        void make_BondFile(const int);
        
        /* Generate Angle file to the structure */
        void make_AngleFile(const int);
        
        /* Write Force Fields to file */
        void write_ForcefieldToFile();
        
        /* use PACKMOL for making system configuration */
        void invoke_PACKMOL(const int);
        
        /* Make LAMMPS Data File */
        void make_LammpsDataFile(const int);
        
        
        /* Convenience Functions */
        double return_PI() const {return 3.1415926535897932;}
        double return_JtoCal() const {return 0.239005736;}
        
        void echoInformation();
        
        void rotate(const double &theta_d,std::vector<double> &xy);
        
        void write_MassToFile(const int,const std::vector<double>&,
                              std::ofstream&);
        
        void write_BondCoeffsToFile(const std::vector<std::string>&,
                                    const StructureClass&,
                                    std::ofstream&);
        
        void write_AngleCoeffsToFile(const std::vector<std::string>&,
                                     const StructureClass&,
                                     std::ofstream&);
        
        void write_PairCoeffsToFile(const std::vector<std::string>&,
                                    const StructureClass&,
                                    std::ofstream&);
        
        const std::vector<double> calc_JtoCal(const int,
                                              const std::vector<double>&,
                                              std::string s="");
        
        const std::vector<double> calc_JtoCal_half(const int,
                                                   const std::vector<double>&,
                                                   std::string s="");
        
        const std::vector<double> setlj(const int,const std::vector<double>&);
        
        
        /*==( public setters )==*/
        
        /* int */
        void set_typeB(const int);
        void set_typeS(const int);
        
        void set_n_poly(const int);
        void set_n_light(const int);
        void set_n_heavy(const int);
        
        /* vector<int> */
        void set_backLenVec(const std::vector<int>&);
        void set_sideLenVec(const std::vector<int>&);
        void set_n_polyVec(const std::vector<int>&);
        void set_waveindex(const std::vector<int>&);
        void set_backboneTypes(const std::vector<int>&);
        void set_sideGroupTypes(const std::vector<int>&);
        void set_types_all(const std::vector<int>&);
        void set_types_light(const std::vector<int>&);
        void set_types_heavy(const std::vector<int>&);
        
        /* vector<vector<int>> */
        void set_n_types_all(const std::vector<std::vector<int>>&);
        void set_n_types_light(const std::vector<std::vector<int>>&);
        void set_n_types_heavy(const std::vector<std::vector<int>>&);
        
        /* vector<double> */
        void set_boxL(const std::vector<double>&);
        void set_pkmL(const std::vector<double>&);
        void set_maxLenScale(const std::vector<double>&);
        
        /*==( public getters )==*/
        
        /* The Naming String */
        const std::string get_nameString(const int) const;
        const std::string get_backboneInWords() const {return backboneInWords;}
        const std::string get_sideGroupInWords() const {return sideGroupInWords;}
        
        /* int */
        int get_typeB() const {return typeB;}
        int get_typeS() const {return typeS;}
        int get_typeBackbone() const {return (int)pow(10,typeB);}
        int get_typeSideGroup() const {return (int)pow(2,typeS);}
        int get_numBack() const {return numBack;}
        int get_numSide() const {return numSide;}
        int get_ringBeads() const {return ringBeads;}
        int get_polySize() const {return polySize;}
        int get_polyBeads() const {return polyBeads;}
        int get_totalBeads() const {return totalBeads;}
        int get_atomType() const {return atomType;}
        int get_bondType() const {return bondType;}
        int get_angleType() const {return angleType;}
        int get_trueAtoms() const {return trueAtoms;}
        int get_trueBonds() const {return trueBonds;}
        int get_trueAngles() const {return trueAngles;}
        
        int get_n_poly() const {return n_poly;}
        int get_n_light() const {return n_light;}
        int get_n_heavy() const {return n_heavy;}
        
        /* double */
        double get_rcutl() const {return rcutl;}
        double get_rcuth() const {return rcuth;}
        double get_rcut() const {return rcut;}
        
        /* vector<int> */
        const std::vector<int> get_backLenVec() const {return backLenVec;}
        const std::vector<int> get_sideLenVec() const {return sideLenVec;}
        const std::vector<int> get_n_polyVec() const {return n_polyVec;}
        const std::vector<int> get_waveindex() const {return waveindex;}
        const std::vector<int> get_backboneTypes() const {return backboneTypes;}
        const std::vector<int> get_sideGroupTypes() const {return sideGroupTypes;}
        const std::vector<int> get_types_all() const {return types_all;}
        const std::vector<int> get_types_light() const {return types_light;}
        const std::vector<int> get_types_heavy() const {return types_heavy;}
        
        /* vector<vector<int>> */
        const std::vector<std::vector<int>> get_n_types_all() const
        {return n_types_all;}
        const std::vector<std::vector<int>> get_n_types_light() const
        {return n_types_light;}
        const std::vector<std::vector<int>> get_n_types_heavy() const
        {return n_types_heavy;}
        
        /* vector<double> */
        const std::vector<double> get_boxL() const {return boxL;}
        const std::vector<double> get_pkmL() const {return pkmL;}
        const std::vector<double> get_mass() const {return mass;}
        const std::vector<double> get_bondLen() const {return bondLen;}
        const std::vector<double> get_bondCoeffs() const {return bondCoeffs;}
        const std::vector<double> get_theta() const {return theta;}
        const std::vector<double> get_angCoeffs() const {return angCoeffs;}
        const std::vector<double> get_lj_eps() const {return lj_eps;}
        const std::vector<double> get_lj_sig() const {return lj_sig;}
        const std::vector<double> get_maxLenScale() const {return maxLenScale;}
        
        /* vector<bool> */
        const std::vector<bool> get_is_backboneType() const {return is_backboneType;}
        const std::vector<bool> get_is_sidegroupType() const {return is_sidegroupType;}
        const std::vector<bool> get_is_shakeBonds() const {return is_shakeBonds;}
        const std::vector<bool> get_is_shakeAngles() const {return is_shakeAngles;}
    };
}
#endif