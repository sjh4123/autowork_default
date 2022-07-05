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

#ifndef LMPSCRIPTS_H
#define LMPSCRIPTS_H

#include "structureclass.h"

namespace autoWork
{
    class LmpScripts
    {
        bool is_tampering;
        bool is_respa;
        bool is_resize;
        bool is_fixMomentum;
        bool is_fixPrdMomentum;
        bool is_fixShake;
        bool is_bondswap;
        
        int fixpsteps;
        int fixpsteps_prd;
        int delay;
        int shakeIter;
        int print_shake;
        
        double shakeTol;
        double Tdamp;
        double Pdamp;
        double pressure;
        
        std::string GPU_mode;
        
    public:
        
        LmpScripts()=delete;
        LmpScripts(const StructureClass&);
        LmpScripts(const LmpScripts&)=default;
        LmpScripts& operator= (const LmpScripts&)=default;
        ~LmpScripts()=default;
        
        /** LAMMPS input scripts generation **/
        void make_GenerationInputFile(const StructureClass&) const;
        void make_QuenchInputFile(const StructureClass&) const;
        void make_EquilibrationInputFile(const StructureClass&) const;
        void make_TamperingInputFile(const StructureClass&) const;
        void make_ResizeInputFile(const StructureClass&) const;
        void make_ProductionInputFile(const StructureClass&,const int) const;
        
        /** simulation target paths **/
        const std::string work_target(const StructureClass&,const std::string&,const int,const int,const double) const;
        
        /** public setters **/
        //----------------------------------------------------------------------
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
        /* string */
        void set_GPU_mode(const std::string&);
        
        /** public getters **/
        //----------------------------------------------------------------------
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
        /* string */
        const std::string get_GPU_mode() const {return GPU_mode;}
    };
}
#endif /* LMPSCRIPTS_H */




