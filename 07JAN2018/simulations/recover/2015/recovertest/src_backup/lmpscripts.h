//
//  Created by Sean Jh H. on 6/18/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#ifndef LMPSCRIPTS_H
#define LMPSCRIPTS_H

#include "structureclass.h"

namespace autoWork {
    
    class LmpScripts
    {
        std::string GPU_mode; // default = "force/neigh"
        
        bool is_tampering;
        bool is_respa;
        bool is_resize;
        bool is_fixMomentum;
        bool is_fixPrdMomentum;
        bool is_fixShake;
        bool is_bondswap;     // default = false
        
        int fixpsteps;
        int fixpsteps_prd;
        int delay;            // default = 5
        int shakeIter;        // default = 20
        int print_shake;      // default = 5000
        
        double shakeTol;      // default = 0.0001
        double Tdamp;         // default = 1000.0
        double Pdamp;         // default = 4000.0
        double pressure;
        
    public:
        
        /* Constructors */
        LmpScripts()=delete;
        LmpScripts(const StructureClass&);
        LmpScripts(const LmpScripts&)=default;
        /* Assignment */
        LmpScripts& operator= (const LmpScripts&)=default;
        /* Destructor */
        ~LmpScripts()=default;
        
        void make_GenerationInputFile(const StructureClass&) const;
        
        void make_QuenchInputFile(const StructureClass&) const;
        
        void make_EquilibrationInputFile(const StructureClass&) const;
        
        void make_TamperingInputFile(const StructureClass&) const;
        
        void make_ResizeInputFile(const StructureClass&) const;
        
        void make_ProductionInputFile(const StructureClass&,
                                      const int n_sys) const;
        
        const std::string work_target(const StructureClass&,
                                      const std::string&,
                                      const int n_trl,
                                      const int n_sys,
                                      const double T) const;
        
        
        /*==( public setters )==*/
        
        /* string */
        void set_GPU_mode(const std::string&);
        
        /* bool */
        void set_is_tampering(const bool);
        void set_is_resize(const bool);
        void set_is_fixMomentum(const bool);
        void set_is_fixPrdMomentum(const bool);
        void set_is_fixShake(const bool);
        void set_is_respa(const bool);
        void set_is_bondswap(const bool);
        
        /* int */
        void set_fixpsteps(const int);
        void set_fixpsteps_prd(const int);
        void set_delay(const int);
        void set_shakeIter(const int);
        void set_print_shake(const int);
        
        /* double */
        void set_shakeTol(const double);
        void set_Tdamp(const double);
        void set_Pdamp(const double);
        void set_pressure(const double);
        
        
        /*==( public getters )==*/
        
        /* string */
        const std::string get_GPU_mode() const {return GPU_mode;}
        
        /* bool */
        bool get_is_tampering() const {return is_tampering;}
        bool get_is_resize() const {return is_resize;}
        bool get_is_fixMomentum() const {return is_fixMomentum;}
        bool get_is_fixPrdMomentum() const {return is_fixPrdMomentum;}
        bool get_is_fixShake() const {return is_fixShake;}
        bool get_is_respa() const {return is_respa;}
        bool get_is_bondswap() const {return is_bondswap;}
        
        /* int */
        int get_fixpsteps() const {return fixpsteps;}
        int get_fixpsteps_prd() const {return fixpsteps_prd;}
        int get_delay() const {return delay;}
        
        /* double */
        double get_Tdamp() const {return Tdamp;};
        double get_Pdamp() const {return Pdamp;};
        double get_pressure() const {return pressure;}
    };
    
}
#endif