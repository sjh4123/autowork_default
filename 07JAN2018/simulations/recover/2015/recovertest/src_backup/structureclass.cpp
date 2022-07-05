//
//  Created by Sean Jh H. on 7/4/14.
//  Copyright (c) 2014 Sean Jh H. All rights reserved.
//

#include "structureclass.h"
#include "functions.h"

using namespace std;
using namespace autoWork;

StructureClass::StructureClass(const int typeB_d,const int typeS_d):

/* string */
nameString(),
backboneInWords(),
sideGroupInWords(),

/* int */
typeB(typeB_d),
typeS(typeS_d),
typeBackbone((int)pow(10,typeB_d)),
typeSideGroup((int)pow(2,typeS_d)),

backLen(),
sideLen(),
n_poly(),
numBack(),
numSide(),
ringBeads(),
polySize(),
polyBeads(),
totalBeads(),
atomType(),
bondType(),
angleType(),
trueAtoms(),
trueBonds(),
trueAngles(),
AtomID(),
BondID(),
AngleID(),

/* double */
rcutl(),
rcuth(),
rcut(),
fenebondK(),
maxfenebondLen(),
fenebond_eps(),
fenebond_sig(),

/* vector<int> */
backLenVec(),
sideLenVec(),
n_polyVec(),
waveindex(),
backboneTypes({1}),
sideGroupTypes({2}),

/* vector<double> */
boxL(),
pkmL(),
mass(),
bondLen(),
bondCoeffs(),
theta(),
angCoeffs(),
lj_eps(),
lj_sig(),
maxLenScale(),
composition({1.0}),

/* vector<bool> */
is_backboneType(),
is_sidegroupType(),
is_shakeBonds(),
is_shakeAngles()

{
    //===============================================
    // Backbone Initiation
    //===============================================
    switch (typeBackbone)
    {
        case (linear):
        {
            backboneInWords="linear";
            atomType += 1; // 1 backbone atom:  B
            bondType += 1; // 1 backbone bond:  B-B
            angleType+= 1; // 1 backbone angle: B-B-B
        }
            break;
            
        case (typeB1):
        {
            backboneInWords="Methyl Acrylate";
            atomType += 2; // 2 backbone atoms:  B,Bs
            bondType += 2; // 2 backbone bonds:  B-B,B-Bs
            angleType+= 2; // 2 backbone angles: B-B-B,B-B-Bs
        }
            break;
            
        case (typeB2):
        {
            backboneInWords="Methyl MethAcrylate";
            atomType += 2; // 2 backbone atoms:  B,Bs
            bondType += 2; // 2 backbone bonds:  B-B,B-Bs
            angleType+= 2; // 2 backbone angles: B-B-B,B-B-Bs
        }
            break;
            
        case (typeB3):
        {
            backboneInWords="Methyl Acrylamide";
            atomType += 2; // 2 backbone atoms:  B,Bs
            bondType += 2; // 2 backbone bonds:  B-B,B-Bs
            angleType+= 2; // 2 backbone angles: B-B-B,B-B-Bs
        }
            break;
            
        case (typeB4):
        {
            backboneInWords="1-ethyl-4-vinylbenzene";
            atomType += 2; // 2 backbone atoms:  B,Bs
            bondType += 3; // 3 backbone bonds:  B-B,B-Bs,Bs-Bs
            angleType+= 4; // 3 backbone angles: B-B-B,Bs-B-Bs,B-B-Bs,B-Bs-Bs
        }
            break;
            
        case (LJliqB):
        {
            backboneInWords="LJliqB";
            atomType += 2; // 2 backbone atom (binary LJ)
            bondType += 0; // 0 backbone bond
            angleType+= 0; // 0 backbone angle
        }
            break;
            
        default:
        {
            cout
            << "StructureClass: 'typeBackbone' NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    //===============================================
    // End of Backbone Initiation
    
    
    //===============================================
    // SideGroup Initiation
    //===============================================
    switch (typeSideGroup)
    {
        case (phenyl):
        {
            sideGroupInWords="phenyl";
            atomType += 1; // 1 side group atom:  R
            bondType += 1; // 1 side group bond:  R-R
            angleType+= 1; // 1 side group angle: R-R-R
        }
            break;
            
        case (alkyl):
        {
            sideGroupInWords="alkyl";
            atomType += 1; // 1 side group atom:  R
            bondType += 1; // 1 side group bond:  R-R
            angleType+= 1; // 1 side group angle: R-R-R
        }
            break;
            
        case (tButyl):
        {
            sideGroupInWords="tert-butyl";
            atomType += 1; // 1 side group atom:  R
            bondType += 0; // 0 side group bond
            angleType+= 0; // 0 side group angle
        }
            break;
            
        case (clProp):
        {
            sideGroupInWords="1-chloropropane";
            atomType += 1; // 1 side group atom:  R
            bondType += 0; // 0 side group bond
            angleType+= 0; // 0 side group angle
        }
            break;
            
        case (C3H8O):
        {
            sideGroupInWords="methoxyethane";
            atomType += 1; // 1 side group atom:  R
            bondType += 0; // 0 side group bond
            angleType+= 0; // 0 side group angle
        }
            break;
            
        case (stdFENE):
        {
            sideGroupInWords="stdFENE";
            angleType=0;
        }
            break;
            
        case (LJliqS):
        {
            sideGroupInWords="LJliqS";
            angleType=0;
        }
            break;
            
        default:
        {
            cout
            << "StructureClass: 'typeSideGroup' NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    //===============================================
    // End of SideGroup Initiation
    
    
    //===============================================
    // 'Inter' Relation
    //===============================================
    // NOTE:
    // The order of structure building is:
    //
    // For Bonds and Angles:
    // always 'Backbone'-->'Side'-->'Inter'
    //
    // For Nonbonded Pairs:
    // follow N*(N+1)/2 convention(11,12,22)
    //===============================================
    
    
    /* Pair Interactions (LJ/GROMACS) */
    //-----------------------------------
    rcutl = 9.0;
    rcuth = 12.0;
    //-----------------------------------
    
    
    switch (typeBackbone+typeSideGroup)
    {
        case (linear+stdFENE):
        {
            mass = {1.0};
            
            composition = {1.0};
            
            // FENE Bonding
            fenebondK      = 30.0;
            maxfenebondLen = 1.5;
            fenebond_eps   = 1.0;
            fenebond_sig   = 1.0;
            
            // FENE bond_style: K, R0, epsilon, sigma
            bondCoeffs     = {30.0,1.5,1.0,1.0};
            
            for (int i=0; i<bondType; ++i) {
                bondLen.push_back((2.0/3.0)*maxfenebondLen);
            }
            
            // No angle forcefield
            // angleType=0;
            
            // LJ Pair Interaction
            rcut   = 2.5;
            lj_eps = setlj(atomType,{1.0});
            lj_sig = setlj(atomType,{1.0});
            
            is_backboneType  = {1};
            is_sidegroupType = {0};
            
        }
            break;
            
        case (LJliqB+LJliqS):
        {
            
            /*******************************************************************
             Binary LJ Forcefield from:
             Kob, W., Andersen, H. C. Phys. Rev. Lett. 1994, 73 (10), 1376–1379
             ******************************************************************/
            
            mass = {1.0,1.0};
            
            composition = {0.8,0.2};
            
            // LJ Pair Interaction
            rcut   = 2.5;
            lj_eps = setlj(atomType,{1.0,1.5,0.5});  // AA AB BB
            lj_sig = setlj(atomType,{1.0,0.8,0.88}); // AA AB BB
            
            is_backboneType  = {1,1};
            is_sidegroupType = {0,0};
            
        }
            break;
            
        case (linear+phenyl):
        {
            bondType += 1; // 1 inter bond:   B-R
            angleType+= 2; // 2 inter angles: R-B-B, B-R-R
            
            /************************************************
             linear + phenyl = PS
             Ref: Rossi et al., 2011 Soft Matter
             ************************************************/
            mass
            = {54.0,54.0};
            
            // ## Bonded (units: Angstrom, kJ/mol, degree)
            bondLen
            = {2.5,2.7,2.17};
            
            bondCoeffs
            //= calc_JtoCal_half(bondType,{8000.0,8000.0,8000.0});
            //= calc_JtoCal_half(bondType,{5000.0,5000.0,5000.0});
            = calc_JtoCal_half(bondType,{2000.0,2000.0,2000.0});
            
            theta
            = {170.0,60.0,125.0,156.0};
            
            angCoeffs
            //= calc_JtoCal_half(angleType,{25.0,200.0,45.0,200.0});
            = calc_JtoCal_half(angleType,{25.0,0.0,45.0,200.0});
            
            // ## NonBonded (units: Angstrom, kJ/mol, degree)
            lj_eps
            = calc_JtoCal(atomType,{2.625,2.325,2.4},"pair");
            
            lj_sig
            = setlj(atomType,{4.3,4.3,4.1});
            
            //===============================================
            // Specify backbone and sidegroup types
            //===============================================
            is_backboneType  = {1,0};
            is_sidegroupType = {0,1};
            //===============================================
            // Specify SHAKE Bonds and Angles
            //===============================================
            //is_shakeBonds  = {0,1,1};
            //is_shakeAngles = {0,1,0,0};
            is_shakeBonds  = {0,0,0};
            is_shakeAngles = {0,0,0,0};
        }
            break;
            
        case (typeB1+alkyl):
        {
            bondType += 1; // 1 inter bond:  Bs-R
            angleType+= 1; // Unsure; use 1: B-Bs-R
            
            /************************************************
             // typeB1 + alkyl
             // Unsure!!
             ************************************************/
            mass
            = {54.0,72.0,54.0};
            
            composition
            = {1.0};
            
            // ## Bonded (units: Angstrom, kJ/mol, degree)
            bondLen
            = {2.82,2.82,2.82,2.82};
            
            bondCoeffs
            = calc_JtoCal_half(bondType,{8000.0,8000.0,8000.0,8000.0});
            
            theta
            = {131.0,71.0,180.0,180.0};
            
            angCoeffs
            = calc_JtoCal_half(angleType,{25.0,85.0,25.0,25.0});
            
            // ## NonBonded (units: Angstrom, kJ/mol, degree)
            lj_eps
            = calc_JtoCal(atomType,{3.56,2.7,3.56,3.7,2.7,3.56},"pair");
            
            lj_sig
            = setlj(atomType,{4.25,4.73,4.25,4.73,4.73,4.25});
            
            //===============================================
            // Specify backbone and sidegroup types
            //===============================================
            is_backboneType  = {1,1,0};
            is_sidegroupType = {0,0,1};
            //===============================================
            // Specify SHAKE Bonds and Angles
            //===============================================
            is_shakeBonds  = {0,1,0,1};
            is_shakeAngles = {0,0,0,0};
        }
            break;
            
        case (typeB2+alkyl):
        {
            bondType += 1; // 1 inter bond:  Bs-R
            angleType+= 1; // Unsure; use 1: B-Bs-R
            
            /************************************************
             // typeB2 + alkyl = PMMA
             // Ref: Uttarwar et al., 2012 I&ECR
             ************************************************/
            mass
            = {54.0,72.0,54.0};
            
            composition
            = {1.0};
            
            // ## Bonded (units: Angstrom, kJ/mol, degree)
            bondLen
            = {2.82,2.82,2.82,2.82};
            
            bondCoeffs
            //= calc_JtoCal_half(bondType,{8000.0,8000.0,8000.0,8000.0});
            = calc_JtoCal_half(bondType,{2000.0,2000.0,2000.0,2000.0});
            
            theta
            = {131.0,71.0,180.0,180.0};
            
            angCoeffs
            = calc_JtoCal_half(angleType,{25.0,85.0,25.0,25.0});
            
            // ## NonBonded (units: Angstrom, kJ/mol, degree)
            lj_eps
            = calc_JtoCal(atomType,{3.56,2.7,3.56,3.7,2.7,3.56},"pair");
            
            lj_sig
            = setlj(atomType,{4.25,4.73,4.25,4.73,4.73,4.25});
            
            //===============================================
            // Specify backbone and sidegroup types
            //===============================================
            is_backboneType  = {1,1,0};
            is_sidegroupType = {0,0,1};
            //===============================================
            // Specify SHAKE Bonds and Angles
            //===============================================
            //is_shakeBonds  = {0,1,0,1};
            //is_shakeAngles = {0,0,0,0};
            is_shakeBonds  = {0,0,0,0};
            is_shakeAngles = {0,0,0,0};
        }
            break;
            
        case (typeB3+alkyl):
        {
            bondType += 1; // 1 inter bond:  Bs-R
            angleType+= 1; // Unsure; use 1: B-Bs-R
            
            /************************************************
             // typeB3 + alkyl
             // Unsure!!
             ************************************************/
            mass
            = {54.0,72.0,54.0};
            
            composition
            = {1.0};
            
            // ## Bonded (units: Angstrom, kJ/mol, degree)
            bondLen
            = {2.82,2.82,2.82,2.82};
            
            bondCoeffs
            //= calc_JtoCal_half(bondType,{8000.0,8000.0,8000.0,8000.0});
            = calc_JtoCal_half(bondType,{2000.0,2000.0,2000.0,2000.0});
            
            theta
            = {131.0,71.0,180.0,180.0};
            
            angCoeffs
            = calc_JtoCal_half(angleType,{25.0,85.0,25.0,25.0});
            
            // ## NonBonded (units: Angstrom, kJ/mol, degree)
            lj_eps
            = calc_JtoCal(atomType,{3.56,2.7,3.56,3.7,2.7,3.56},"pair");
            
            lj_sig
            = setlj(atomType,{4.25,4.73,4.25,4.73,4.73,4.25});
            
            //===============================================
            // Specify backbone and sidegroup types
            //===============================================
            is_backboneType  = {1,1,0};
            is_sidegroupType = {0,0,1};
            //===============================================
            // Specify SHAKE Bonds and Angles
            //===============================================
            is_shakeBonds  = {0,1,0,1};
            is_shakeAngles = {0,0,0,0};
        }
            break;
            
        case (typeB4+tButyl):
        {
            bondType += 1; // 1 inter bond:  Bs-R
            angleType+= 1; // Unsure; use 1: Bs-R-Bs
            
            /************************************************
             // typeB4 + alkyl
             // Ref: Rossi et al., 2011 Soft Matter, A-mapping
             ************************************************/
            mass
            = {54.0,54.0,54.0};
            
            composition
            = {1.0};
            
            // ## Bonded (units: Angstrom, kJ/mol, degree)
            bondLen
            = {2.5,2.7,2.17,2.5};
            
            bondCoeffs
            //= calc_JtoCal_half(bondType,{8000.0,8000.0,8000.0,8000.0});
            = calc_JtoCal_half(bondType,{2000.0,2000.0,2000.0,2000.0});
            
            theta
            = {170.0,60.0,125.0,156.0,52.0};
            
            angCoeffs
            = calc_JtoCal_half(angleType,{25.0,200.0,45.0,200.0,550.0});
            
            // ## NonBonded (units: Angstrom, kJ/mol, degree)
            lj_eps
            = calc_JtoCal(atomType,{2.625,2.325,2.625,2.4,2.325,2.625},"pair");
            
            lj_sig
            = setlj(atomType,{4.3,4.3,4.3,4.1,4.3,4.3});
            
            //===============================================
            // Specify backbone and sidegroup types
            //===============================================
            is_backboneType  = {1,1,0};
            is_sidegroupType = {0,0,1};
            //===============================================
            // Specify SHAKE Bonds and Angles
            //===============================================
            is_shakeBonds  = {0,1,1,1};
            is_shakeAngles = {0,1,0,0,1};
        }
            break;
            
        default:
        {
            cout
            << "StructureClass: 'Inter' relation NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    //===============================================
    // End of 'Inter'
}
//===============================================
// End of Constructor





void StructureClass::set_structureParam()
{
    
    //----------------------------------------------------------
    // structural variables are vectors with size of "n_system"
    //----------------------------------------------------------
    
    switch (typeBackbone+typeSideGroup) {
            
        case (linear+stdFENE): // Kremer-Grest Bead-Spring
        {
            set_backLenVec({20});
            set_sideLenVec({0});
            set_n_polyVec({80});
            set_boxL({25.0});
            set_pkmL({20.0});
            
            set_maxLenScale({12.43848});
            set_waveindex({26});
        }
            break;
            
        case (LJliqB+LJliqS):
        {
            set_backLenVec({1});
            set_sideLenVec({0});
            set_n_polyVec({1000});
            set_boxL({25.0});
            set_pkmL({20.0});
            
            set_maxLenScale({12.43848});
            set_waveindex({26});
        }
            break;
            
        case (linear+phenyl): // Martini PS
        {
            set_backLenVec({30,10,60,100,150,200});
            set_sideLenVec({3,3,3,3,3,3});
            set_n_polyVec({32,75,18,18,15,15});
            set_boxL({300.0,200.0,400.0,600.0,900.0,1000.0});
            set_pkmL({150.0,100.0,200.0,300.0,450.0,600.0});
            
            set_maxLenScale({55,55,55,55,55,55});
            set_waveindex({23,23,23,23,23,23});
        }
            break;
            
        case (typeB2+alkyl): // PMMA, PMA
        {
            set_backLenVec({30,60});
            set_sideLenVec({0,0});
            set_n_polyVec({32,18});
            set_boxL({400.0,600.0});
            set_pkmL({200.0,300.0});
            
            set_maxLenScale({55,55}); // note: use those from PS as standard
            set_waveindex({23,23});   // note: use those from PS as standard
        }
            break;
            
        default:
            break;
    }
}





///////////////////////////////////////////////////////
//                                                   //
//  ATOMS                                            //
// ================================================= //
//  Generate single chain topology (singleChain.xyz) //
//  Building Order:                                  //
//  Backbone --> Side Chain                          //
//  '0' based counting                               //
//  Independent of BONDS and ANGLES sections         //
//                                                   //
///////////////////////////////////////////////////////

void StructureClass::make_SingleChainXYZ(const int n_sys)
{
    n_poly  = n_polyVec[n_sys];
    sideLen = sideLenVec[n_sys];
    
    switch (typeBackbone)
    {
        case (linear):
        {
            backLen = backLenVec[n_sys];
            numBack = backLen;
        }
            break;
        case (typeB1):
        {
            backLen = backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
        case (typeB2):
        {
            backLen = backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
        case (typeB3):
        {
            backLen = backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
        case (typeB4):
        {
            ringBeads= 3;
            backLen = backLenVec[n_sys];
            numBack = backLen*(1+ringBeads);
        }
            break;
        case (LJliqB):
        {
            backLen = backLenVec[n_sys];
            numBack = backLen;
        }
            break;
        default:
        {
            cout
            << "makeStructure: 'typeBackbone' NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    
    numSide	 = sideLen*backLen; // total number of beads of the side groups
    polySize = numSide+numBack; // total number of beads per molecule
    polyBeads= n_poly*polySize;
    
    /*=====================================================
     allocate vector size according to the size of chain
     ====================================================*/
    vector<double> x(polySize);
    vector<double> y(polySize);
    vector<double> z(polySize);
    
    string singleChainFiles="singleChain.xyz";
    ofstream AtomFile(singleChainFiles.c_str());
    
    if(AtomFile.is_open())
    {
        AtomFile
        << polySize << "\n"
        << "poly_B" << typeB << "_" << backLen
        << "_S" << typeS << "_" << sideLen
        << "\n";
        
        
        
        /* INITIALIZE ATOM TYPE TO 1 */
        int AtomType = 1;
        
        
        
        /*
         used as the index of the 'current' bond type
         i.e. it's in 'backbone', or 'side', or 'inter'. */
        int BondType = 0;
        
        
        
        
        //===========================================
        // Backbone Beads
        //===========================================
        // beads are allocated on the x_direction;
        // backbone beads are symmetric to the origin
        //===========================================
        switch (typeBackbone)
        {
                //===============================================
                // Backbone: ## linear
                //===============================================
            case (linear):
            {
                
                x[0] = -bondLen[BondType]*((double)backLen/2.0 - 0.5);
                for (int i=1; i<backLen; ++i)
                    x[i] = x[i-1] + bondLen[BondType];
            }
                break;
                
                
                
                //===============================================
                // Backbone: ## typeB1
                //
                // (Methyl Acrylate)n; (C4H6O2)n
                //
                //===============================================
            case (typeB1):
            {
                // backbone atom type 1 (linear chain)
                x[0] = -bondLen[BondType]*((double)backLen/2.0 - 0.5);
                for (int i=1; i<backLen; ++i)
                    x[i] = x[i-1] + bondLen[BondType];
                
                AtomType += 1; // there are 2 atom types in the backbone
                BondType += 1; // the current bond type index is 1
                
                // backbone atom type 2
                for (int i=0; i<backLen; ++i)
                {
                    x[i+backLen] = x[i];
                    if ( (i % 2) == 0 )
                        y[i+backLen] = y[i] - bondLen[BondType];
                    else
                        y[i+backLen] = y[i] + bondLen[BondType];
                }
            }
                break;
                
                
                
                //===============================================
                // Backbone: ## typeB2
                //
                // (Methyl MethAcrylate)n; (C5H8O2)n
                //
                //===============================================
            case (typeB2):
            {
                // backbone atom type 1 (linear chain)
                x[0] = -bondLen[BondType]*((double)backLen/2.0 - 0.5);
                for (int i=1; i<backLen; ++i)
                    x[i] = x[i-1] + bondLen[BondType];
                
                AtomType += 1; // there are 2 atom types in the backbone
                BondType += 1; // the current bond type index is 1
                
                // backbone atom type 2
                for (int i=0; i<backLen; ++i)
                {
                    x[i+backLen] = x[i];
                    if ( (i % 2) == 0 )
                        y[i+backLen] = y[i] - bondLen[BondType];
                    else
                        y[i+backLen] = y[i] + bondLen[BondType];
                }
            }
                break;
                
                
                
                //===============================================
                // Backbone: ## typeB3
                //
                // (Methyl Acrylamide)n; (C4H7NO)n
                //
                //===============================================
            case (typeB3):
            {
                // backbone atom type 1 (linear chain)
                x[0] = -bondLen[BondType]*((double)backLen/2.0 - 0.5);
                for (int i=1; i<backLen; ++i)
                    x[i] = x[i-1] + bondLen[BondType];
                
                AtomType += 1; // there are 2 atom types in the backbone
                BondType += 1; // the current bond type index is 1
                
                // backbone atom type 2
                for (int i=0; i<backLen; ++i)
                {
                    x[i+backLen] = x[i];
                    if ( (i % 2) == 0 )
                        y[i+backLen] = y[i] - bondLen[BondType];
                    else
                        y[i+backLen] = y[i] + bondLen[BondType];
                }
            }
                break;
                
                
                
                //===============================================
                // Backbone: ## typeB4
                //
                // (1-ethyl-4-vinylbenzene)n; (C10H12)n
                //
                // This backbone structure is the same with a
                // linear+phenyl chain, if there's no side beads
                //===============================================
            case (typeB4):
            {
                // backbone atom type 1 (linear chain)
                x[0] = -bondLen[BondType]*((double)backLen/2.0 - 0.5);
                for (int i=1; i<backLen; ++i)
                    x[i] = x[i-1] + bondLen[BondType];
                
                AtomType += 1;
                BondType += 1;
                
                // backbone atom type 2 (phenyl ring)
                for (int i=0; i<backLen; ++i)
                {
                    int sum=0;
                    if ( i>0 )
                        for (int ii=0; ii<i; ++ii)
                            sum += ringBeads;
                    
                    int first = backLen + sum; // NOTE: here is zero-based.
                    
                    x[first]=x[i];
                    if( (i % 2) == 0 )
                        y[first]=y[i]-bondLen[BondType+1]; // NOTE: BondType+1
                    else
                        y[first]=y[i]+bondLen[BondType+1];
                    
                    //===============================================
                    // Rotation of the ring bonds
                    //===============================================
                    double vec_xdir=x[first]-x[i]; // x_dir
                    double vec_ydir=y[first]-y[i]; // y_dir
                    
                    /*
                     GET Y_DIR UNIT VECTOR */
                    if ( sqrt(vec_ydir*vec_ydir)!=0.0 ) {
                        vec_ydir=vec_ydir/sqrt(vec_ydir*vec_ydir);
                    }
                    else {
                        cout
                        << "WARNING: typeB4: /0 error!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    vector<double> xy(2);
                    xy[0]=bondLen[BondType]*vec_xdir; // x_dir
                    xy[1]=bondLen[BondType]*vec_ydir; // y_dir
                    
                    for (int ii=1; ii<ringBeads; ++ii)
                    {
                        double phi0=(2.0*return_PI()/(double)ringBeads);
                        if ( ii == 1 )
                        {
                            double phi1 = (return_PI()-phi0)/2.0;
                            rotate(phi1, xy);
                        }
                        else
                            rotate(-phi0, xy);
                        
                        x[first+ii]=x[first+ii-1];
                        y[first+ii]=y[first+ii-1]+xy[1];
                        z[first+ii]=z[first+ii-1]+xy[0];
                    }
                    //===============================================
                }
            }
                break;
                
                //===============================================
                // Backbone: ## LJliqB
                //===============================================
            case (LJliqB):
            {
                
                x[0] = 0;
                y[0] = 0;
                z[0] = 0;
                
                AtomType = composition.size();
                
            }
                break;
                
            default:
            {
                cout
                << "ATOMS: The 'backbone' is NOT defined."
                << "\n";
                system("read -t5");exit(EXIT_FAILURE);
            }
                break;
        }
        
        
        
        //===============================================
        // Write to File
        //===============================================
        // Note:
        // For linear, AtomType should only be 1;
        // For other backbones, AtomType = 1, 2;
        //===============================================
        for (int i=0; i<numBack; ++i)
        {
            if ((typeBackbone==linear)||(typeBackbone==LJliqB))
                // type 1 backbone bead
                AtomFile
                << AtomType << " "
                << x[i] << " " << y[i] << " " << z[i] << "\n";
            else
            {
                if (i<backLen)
                    // type 1 backbone bead
                    AtomFile
                    << AtomType-1 << " "
                    << x[i] << " " << y[i] << " " << z[i] << "\n";
                else
                    // type 2 backbone bead
                    AtomFile
                    << AtomType << " "
                    << x[i] << " " << y[i] << " " << z[i] << "\n";
            }
        }
        //===============================================
        // End of 'Backbone' Beads
        
        
        
        
        
        /* ===============================================
         index_FinalBack: index of the final backbone bead
         
         NOTE:
         here is 0-based counting
         =============================================== */
        int index_FinalBack = numBack - 1;
        
        
        
        
        //===============================================
        // Side Chain Beads
        //===============================================
        // beads are allocated on the x-y plane
        // originally in a 'syndiotactic' configuration,
        // i.e.
        // neighbor groups point in the opposite directions
        //===============================================
        
        if (sideLen>0) // should have at least 1 side bead
        {
            
            
            AtomType += 1;
            
            
            
            if (typeBackbone==typeB4)
                BondType += 2;
            else
                BondType += 1;
            
            
            
            switch (typeSideGroup)
            {
                    //===============================================
                    // Side Chain: ## stdFENE(5)
                    // Can pair with linear(0)
                    //===============================================
                case (stdFENE):
                {
                    // no side group atoms defined
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Side Chain: ## phenyl(0)
                    // Can pair with: linear(0),B(3)
                    //===============================================
                case (phenyl):
                {
                    
                    //===============================================
                    // Restrictions
                    //===============================================
                    if ((typeBackbone==typeB1) ||
                        (typeBackbone==typeB2) ||
                        (typeBackbone==typeB4))
                    {
                        cout
                        << "WARNING: backbone type can't pair with phenyl!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    
                    
                    if (sideLen!=3)
                    {
                        cout
                        << "WARNING: phenyl sideLen != 3; make it = 3 or 0"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    //===============================================
                    
                    //===============================================
                    // Build phenyl side groups
                    //===============================================
                    for (int i=0; i<backLen; ++i)
                    {
                        //===============================================
                        int sideSum = 0;
                        if(i>0)
                            for(int ii=0; ii<i; ++ii)
                                sideSum += sideLen;
                        
                        int index_root;
                        if (typeBackbone==linear)
                            index_root = i;
                        else
                            index_root = backLen + i;
                        
                        int index_FirstSide = index_FinalBack + sideSum + 1;
                        //===============================================
                        
                        // Inter (NOTE: BondType)
                        x[index_FirstSide]=x[index_root];
                        if ( (i%2)==0 )
                            y[index_FirstSide]=y[index_root]-bondLen[BondType+1];
                        else
                            y[index_FirstSide]=y[index_root]+bondLen[BondType+1];
                        
                        /*===============================================
                         Rotate the bond vectors
                         in order to make the side group geometry
                         
                         NOTE:
                         side beads will be formed on the y-z plane
                         whereas the backbone is on the x-y plane;
                         this way strong repulsion between side beads
                         can be prevented and thus better initial config
                         ===============================================*/
                        //double tmp[2];
                        //tmp[0]=x[index_FirstSide]-x[index_root]; // x_dir
                        //tmp[1]=y[index_FirstSide]-y[index_root]; // y_dir
                        
                        double vec_xdir=x[index_FirstSide]-x[index_root];
                        double vec_ydir=y[index_FirstSide]-y[index_root];
                        
                        /*
                         GET Y_DIR UNIT VECTOR */
                        if ( sqrt(vec_ydir*vec_ydir)!=0.0 ) {
                            vec_ydir=vec_ydir/sqrt(vec_ydir*vec_ydir);
                        }
                        else {
                            cout
                            << "WARNING: phenyl sideGroup: /0 error!"
                            << "\n";
                            system("read -t5");exit(EXIT_FAILURE);
                        }
                        vector<double> xy(2);
                        xy[0]=bondLen[BondType]*vec_xdir; // x_dir
                        xy[1]=bondLen[BondType]*vec_ydir; // y_dir
                        
                        /*===============================================
                         loops through rest of the side beads
                         
                         NOTE:
                         rotate in the "y-z" plane
                         ===============================================*/
                        for (int ii=1; ii<sideLen; ++ii)
                        {
                            /*
                             The angle that equally separates side beads */
                            double tmptheta=(2.0*return_PI()/(double)sideLen);
                            if ( ii==1 ) {
                                double theta1=(return_PI()-tmptheta)/2.0;
                                rotate(theta1,xy);
                            }
                            else
                                rotate(-tmptheta,xy);
                            
                            x[index_FirstSide+ii]=x[index_FirstSide+ii-1];
                            y[index_FirstSide+ii]=y[index_FirstSide+ii-1]+xy[1];
                            z[index_FirstSide+ii]=z[index_FirstSide+ii-1]+xy[0];
                        }
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Side Chain: ## alkyl(1)
                    // Can pair with: linear(0),B(1),B(2),B(3)
                    //===============================================
                case (alkyl):
                {
                    //===============================================
                    // Restrictions
                    //===============================================
                    if (typeBackbone==typeB4)
                    {
                        cout
                        << "WARNING: typeB4 can't pair with alkyl!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    //===============================================
                    
                    //===============================================
                    // Build alkyl side groups
                    //===============================================
                    for (int i=0; i<backLen; ++i)
                    {
                        //===============================================
                        int sideSum = 0;
                        if( i > 0 )
                            for(int ii=0; ii<i; ++ii)
                                sideSum += sideLen;
                        
                        int index_root;
                        if (typeBackbone==linear)
                            index_root = i;
                        else
                            index_root = backLen + i;
                        
                        int index_FirstSide = index_FinalBack + sideSum + 1;
                        //===============================================
                        
                        x[index_FirstSide] = x[index_root];
                        if ( (i%2)==0 )
                        {
                            // Inter (NOTE: BondType)
                            y[index_FirstSide] = y[index_root] - bondLen[BondType+1];
                            if (sideLen>1)
                                for (int ii=1; ii<sideLen; ++ii)
                                {
                                    x[index_FirstSide+ii] = x[index_FirstSide+ii-1];
                                    y[index_FirstSide+ii] = y[index_FirstSide+ii-1] - bondLen[BondType];
                                }
                        }
                        else
                        {
                            // Inter (NOTE: BondType)
                            y[index_FirstSide] = y[index_root] + bondLen[BondType+1];
                            if (sideLen>1)
                                for (int ii=1; ii<sideLen; ++ii)
                                {
                                    x[index_FirstSide+ii] = x[index_FirstSide+ii-1];
                                    y[index_FirstSide+ii] = y[index_FirstSide+ii-1] + bondLen[BondType];
                                }
                        }
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Side Chain: ## tert-butyl(2)
                    // Can pair with: linear(0),B(1),B(2),B(4)
                    // sideLen=1
                    //===============================================
                case (tButyl):
                {
                    //===============================================
                    // Restrictions
                    //===============================================
                    if (typeBackbone==typeB3)
                    {
                        cout
                        << "WARNING: typeB3 can't pair with tButyl!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    if ( sideLen!=1 )
                    {
                        cout
                        << "WARNING: tButyl can only have 1 bead!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    //===============================================
                    
                    //===============================================
                    // build tButyl side groups
                    //===============================================
                    for (int i=0; i<backLen; ++i)
                    {
                        //===============================================
                        int sideSum = 0;
                        int ringSum = 0;
                        if( i > 0 )
                            for(int ii=0; ii<i; ++ii)
                            {
                                sideSum += sideLen;
                                if (typeBackbone==typeB4)
                                    ringSum += ringBeads;
                            }
                        
                        int index_root;
                        if (typeBackbone==linear)
                            index_root = i;
                        else if (typeBackbone==typeB4)
                            index_root = backLen + ringSum;
                        else
                            index_root = backLen + i;
                        
                        int index_FirstSide = index_FinalBack + sideSum + 1;
                        //===============================================
                        double tmp = y[index_root] - y[i]; // y_dir vector
                        if ( tmp != 0.0 ) // denominator can't be zero
                            tmp = tmp / sqrt(tmp*tmp); // get unit vectors
                        
                        // phi1: Bs-(Bs)-R
                        // the0: (Bs)-R-Bs
                        double phi0 = (2.0*return_PI()/(double)ringBeads);
                        double phi1 = (return_PI()-phi0)/2.0;
                        double the0 = theta[angleType-2]*0.5*(2.0*return_PI()/180.0);
                        
                        double xy;
                        xy = fabs(tmp)*(bondLen[BondType]*cos(the0)+bondLen[BondType-2]*cos(phi1));
                        
                        x[index_FirstSide] = x[index_root];
                        if ( (i%2)==0 )
                            y[index_FirstSide] = y[index_root] - xy;
                        else
                            y[index_FirstSide] = y[index_root] + xy;
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Side Chain: ## 1-chloropropane(3)
                    // Can pair with: linear(0),B(1)
                    // sideLen=1
                    //===============================================
                case (clProp):
                {
                    //===============================================
                    // Restrictions
                    //===============================================
                    if ((typeBackbone==typeB2) ||
                        (typeBackbone==typeB3) ||
                        (typeBackbone==typeB4))
                    {
                        cout
                        << "WARNING: backbone type can't pair with 1-chloropropane!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    
                    if ( sideLen!=1 )
                    {
                        cout
                        << "WARNING: 1-chloropropane can only have 1 bead!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    //===============================================
                    
                    // build 1-chloropropane side groups
                    //===============================================
                    for (int i=0; i<backLen; ++i)
                    {
                        
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Side Chain: ## methoxyethane(4)
                    // Can pair with: linear(0),B(3)
                    // sideLen=1
                    //===============================================
                case (C3H8O):
                {
                    //===============================================
                    // Restrictions
                    //===============================================
                    if ((typeBackbone==typeB1) ||
                        (typeBackbone==typeB2) ||
                        (typeBackbone==typeB4))
                    {
                        cout
                        << "WARNING: backbone type can't pair with methoxyethane!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    if ( sideLen!=1 )
                    {
                        cout
                        << "WARNING: methoxyethane can only have 1 bead!"
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                    //===============================================
                    
                    //===============================================
                    // build methoxyethane side groups
                    //===============================================
                    for (int i=0; i<backLen; ++i)
                    {
                        
                    }
                }
                    break;
                    
                    
                    
                default:
                {
                    cout
                    << "ATOMS: The 'side chain' is NOT defined."
                    << "\n";
                    system("read -t5");exit(EXIT_FAILURE);
                }
                    break;
            }
            
            
            
            //===============================================
            // Write to File
            //===============================================
            for (int i=numBack; i<polySize; ++i)
            {
                AtomFile << AtomType << " "
                << x[i] << " " << y[i] << " " << z[i] << "\n";
            }
            //===============================================
            // End of 'Side Chain' Beads
        }
        else
        {
            /*
             NOTE:
             
             Here assuming side chain only has 'one' atom type
             in side chain itself */
            
            static int tmp=0; // static variable here to ensure only do once
            if(tmp==0)
            {
                if (sideLen==0)
                {
                    if ((typeBackbone==linear)&&(typeSideGroup!=stdFENE))
                    {
                        /*
                         for linear backbone, index of side chain beads should
                         start at '1' so, it's the 2nd element in the arrays use
                         'begin()+1' to point to the 2nd element and erase it. */
                        mass.erase(mass.begin()+1);
                        is_backboneType.erase(is_backboneType.begin()+1);
                        is_sidegroupType.erase(is_sidegroupType.begin()+1);
                        
                        cout << "**************************\n";
                        cout << "size of mass array= " << mass.size() << "\n";
                        for (size_t i=0; i<mass.size(); ++i) cout << " " << mass[i]; cout << "\n";
                        cout << "size of is_backboneType array= " << is_backboneType.size() << "\n";
                        for (size_t i=0; i<is_backboneType.size(); ++i) cout << " " << is_backboneType[i]; cout << "\n";
                        cout << "**************************\n";
                    }
                    else if ((typeSideGroup!=stdFENE)&&(typeSideGroup!=LJliqS))
                    {
                        /*
                         for typeB1~4, index of side chain beads should start at
                         '2' so, it's the 3rd element in the arrays use
                         'begin()+2' to point to the 3rd element and erase it. */
                        mass.erase(mass.begin()+2);
                        is_backboneType.erase(is_backboneType.begin()+2);
                        is_sidegroupType.erase(is_sidegroupType.begin()+2);
                        
                        cout << "**************************\n";
                        cout << "size of mass array= " << mass.size() << "\n";
                        for (size_t i=0; i<mass.size(); ++i) cout << " " << mass[i]; cout << "\n";
                        cout << "size of is_backboneType array= " << is_backboneType.size() << "\n";
                        for (size_t i=0; i<is_backboneType.size(); ++i) cout << " " << is_backboneType[i]; cout << "\n";
                        cout << "**************************\n";
                    }
                    ++tmp;
                }
                else
                {
                    cout
                    << "Error: Number of Side Beads < 0 !" << "\n";
                    system("read -t5");exit(EXIT_FAILURE);
                }
            }
        }
        
        // 'trueAtoms' records the number of beads in a molecule
        // taking into account that 'sideLen' might = 0;
        
        trueAtoms = AtomType;
        AtomFile.close();
    }
    else
    {
        cout
        << "'singleChain.xyz' cannot open."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}
//===============================================
// End of ATOMS





///////////////////////////////////////////////////////
//                                                   //
//  BONDS                                            //
// ===============================================   //
//  Generate 'Bond.txt'                              //
//  Building Order:                                  //
//  Backbone --> Side Chain --> Inter                //
//  '1' based counting (rahter than '0' based)       //
//  Independent of ATOMS and ANGLES sections         //
//                                                   //
//  Note:                                            //
//  By default, typeB1, typeB2, typeB3               //
//  have the same bond configuration                 //
//                                                   //
///////////////////////////////////////////////////////

void StructureClass::make_BondFile(const int n_sys)
{
    n_poly = n_polyVec[n_sys];
    sideLen = sideLenVec[n_sys];
    
    switch (typeBackbone)
    {
        case (linear):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen;
        }
            break;
            
        case (typeB1):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB2):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB3):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB4):
        {
            ringBeads= 3;
            backLen	 = backLenVec[n_sys];
            numBack  = backLen*(1+ringBeads);
        }
            break;
            
        case (LJliqB):
        {
            // no bonds
            return;
        }
            break;
            
        default:
        {
            cout
            << "BONDS: 'typeBackbone' NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    
    numSide	= sideLen * backLen;
    polySize= numSide + numBack;
    
    ofstream BondFile("Bonds.txt");
    if( BondFile.is_open() )
    {
        BondID = 0;
        
        //===============================================
        // Generic Looping Structure:
        // Iterate through all the chains in the system
        //===============================================
        for (int i=1; i<=n_poly; ++i)
        {
            
            
            
            /* INITIALIZE BOND TYPE TO 1 */
            int BondType = 1;
            
            
            
            //===============================================
            // 1. Backbone Bonds
            //===============================================
            switch (typeBackbone)
            {
                    //===============================================
                    // Backbone: ## linear
                    //===============================================
                case (linear):
                {
                    // Type (B-B)
                    for (int ii=1; ii<backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+1) + (i-1)*polySize
                        << "\n";
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Backbone: ## typeB1
                    // (Methyl Acrylate)n; (C4H6O2)n
                    //===============================================
                case (typeB1):
                {
                    // backbone bond type 1
                    // Type (B-B)
                    for (int ii=1; ii<backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+1) + (i-1)*polySize
                        << "\n";
                    }
                    
                    BondType += 1;
                    
                    // backbone bond type 2
                    // Type (B-Bs)
                    for (int ii=1; ii<=backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+backLen) + (i-1)*polySize
                        << "\n";
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Backbone: ## typeB2
                    // (Methyl MethAcrylate)n; (C5H8O2)n
                    //===============================================
                case (typeB2):
                {
                    // backbone bond type 1
                    // Type (B-B)
                    for (int ii=1; ii<backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+1) + (i-1)*polySize
                        << "\n";
                    }
                    
                    BondType += 1;
                    
                    // backbone bond type 2
                    // Type (B-Bs)
                    for (int ii=1; ii<=backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+backLen) + (i-1)*polySize
                        << "\n";
                    }
                }
                    break;
                    
                    
                    
                    //===============================================
                    // Backbone: ## typeB3
                    // (Methyl Acrylamide)n; (C4H7NO)n
                    //===============================================
                case (typeB3):
                {
                    // backbone bond type 1
                    // Type (B-B)
                    for (int ii=1; ii<backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+1) + (i-1)*polySize
                        << "\n";
                    }
                    
                    BondType += 1;
                    
                    // backbone bond type 2
                    // Type (B-Bs)
                    for (int ii=1; ii<=backLen; ++ii)
                    {
                        BondID += 1;
                        BondFile
                        << BondID << " "
                        << BondType << " "
                        << (ii) + (i-1)*polySize
                        << " "
                        << (ii+backLen) + (i-1)*polySize
                        << "\n";
                    }
                }
                    break;
                    
                    
                    
                default:
                {
                    cout
                    << "BONDS: 'backbone' relation NOT defined."
                    << "\n";
                    system("read -t5");exit(EXIT_FAILURE);
                }
                    break;
            }
            //===============================================
            // End of 'Backbone' Bonds
            
            
            
            /*
             index_FinalBack: index of the final backbone bead
             
             Note:
             NOT 'numBack-1' (here is 1-based counting) */
            int index_FinalBack = numBack;
            
            
            
            //===============================================
            // 2. Side Chain Bonds
            //===============================================
            // side chains should have at least 2 beads
            // to have side chain bonds
            //===============================================
            if (sideLen>1)
            {
                
                
                
                BondType += 1;
                
                
                
                switch (typeSideGroup)
                {
                        //===============================================
                        // Side Chain: ## phenyl
                        //===============================================
                    case (phenyl):
                    {
                        // Type (R-R)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            /*
                             from 2nd side group on
                             add up all previous side beads */
                            int sideSum = 0;
                            if( ii > 1 )
                                for (int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            
                            /*
                             INTRA-RING BONDS */
                            for (int iii=1; iii<=sideLen; ++iii)
                            {
                                /*
                                 THE FINAL BEAD IN THE RING LINKS TO
                                 THE FIRST BEAD OF THE RING */
                                if ( iii==sideLen ) {
                                    
                                    BondID += 1;
                                    BondFile
                                    << BondID << " "
                                    << BondType << " "
                                    << (index_FinalBack+sideSum+iii)+(i-1)*polySize
                                    << " "
                                    << (index_FinalBack+sideSum+1)+(i-1)*polySize
                                    << "\n";
                                }
                                else {
                                    
                                    /*
                                     WHEN SHAKE_BOND IS TURNED OFF FOR RINGS */
                                    if (get_is_shakeBonds()[1]==false) {
                                        
                                        BondID += 1;
                                        BondFile
                                        << BondID << " "
                                        << BondType << " "
                                        << (index_FinalBack+sideSum+iii)+(i-1)*polySize
                                        << " "
                                        << (index_FinalBack+sideSum+iii+1)+(i-1)*polySize
                                        << "\n";
                                    }
                                    /*
                                     WHEN SHAKE_BOND IS TURNED ON FOR RINGS */
                                    else {
                                        
                                        // NOTE:
                                        //===============================================
                                        // The special bonding configuration
                                        // is to fit the LAMMPS fix 'SHAKE' setting
                                        // should take note if ${sideBead} != 2
                                        //===============================================
                                        if ( iii != 2 ) {
                                            BondID += 1;
                                            BondFile
                                            << BondID << " "
                                            << BondType << " "
                                            << (index_FinalBack+sideSum+iii)+(i-1)*polySize
                                            << " "
                                            << (index_FinalBack+sideSum+iii+1)+(i-1)*polySize
                                            << "\n";
                                        }
                                    }
                                }
                            }
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Side Chain: ## alkyl
                        //===============================================
                    case (alkyl):
                    {
                        // Type (R-R)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            int sideSum = 0;
                            if( ii > 1 )
                                for(int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            
                            for (int iii=1; iii<sideLen; ++iii)
                            {
                                BondID += 1;
                                BondFile
                                << BondID << " "
                                << BondType << " "
                                << (index_FinalBack+sideSum+iii)+(i-1)*polySize
                                << " "
                                << (index_FinalBack+sideSum+iii+1)+(i-1)*polySize
                                << "\n";
                            }
                        }
                    }
                        break;
                        
                        
                        
                    default:
                    {
                        cout
                        << "BONDS: 'side chain' relation NOT defined."
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                        break;
                }
            }
            else
            {
                // NOTE:
                // Here assuming side chain only has 'one' bond type
                // in side chain itself
                
                static int tmp=0; // static variable here to ensure only do once
                if(tmp==0)
                {
                    if (sideLen==1) // if sideLen=1, erase side bonds of typeB1~3
                    {
                        if((typeBackbone==linear)&&(typeSideGroup!=stdFENE))
                        {
                            bondLen.erase(bondLen.begin()+1);
                            bondCoeffs.erase(bondCoeffs.begin()+1);
                            is_shakeBonds.erase(is_shakeBonds.begin()+1);
                            
                            cout << "**************************\n";
                            cout << "bondLen size= " << bondLen.size() << "\n";
                            for (size_t i=0; i<bondLen.size(); ++i) cout << " " << bondLen[i]; cout << "\n";
                            cout << "is_shakeBonds size= " << is_shakeBonds.size() << "\n";
                            for (size_t i=0; i<is_shakeBonds.size(); ++i) cout << " " << is_shakeBonds[i]; cout << "\n";
                            cout << "**************************\n";
                        }
                        else if ((typeBackbone!=linear)&&(typeBackbone!=typeB4))
                        {
                            bondLen.erase(bondLen.begin()+2);
                            bondCoeffs.erase(bondCoeffs.begin()+2);
                            is_shakeBonds.erase(is_shakeBonds.begin()+2);
                            
                            cout << "**************************\n";
                            cout << "bondLen size= " << bondLen.size() << "\n";
                            for (size_t i=0; i<bondLen.size(); ++i) cout << " " << bondLen[i]; cout << "\n";
                            cout << "is_shakeBonds size= " << is_shakeBonds.size() << "\n";
                            for (size_t i=0; i<is_shakeBonds.size(); ++i) cout << " " << is_shakeBonds[i]; cout << "\n";
                            cout << "**************************\n";
                        }
                        ++tmp;
                    }
                }
            }
            //===============================================
            // End of 'Side Chain' Bonds
            
            
            //===============================================
            // 3. Inter Bonds
            //===============================================
            // there should have at least 1 side bead
            // to have inter bonds
            //===============================================
            if (sideLen>0)
            {
                
                
                
                BondType += 1;
                
                
                
                switch (typeBackbone+typeSideGroup)
                {
                        //===============================================
                        // Inter: ## linear + phenyl (PS)
                        //===============================================
                    case (linear+phenyl):
                    {
                        // Type (B-R)
                        // each linear backbone bead has one phenyl side group
                        // each phenyl side group has #{sideLen_d} beads
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            int index_root = ii;
                            int sideSum = 0;
                            if( ii > 1 )
                                for(int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            BondID += 1;
                            BondFile << BondID << " " << BondType << " "
                            << (index_root) + (i-1)*polySize << " "
                            << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Inter: ## typeB1 + alkyl
                        //===============================================
                    case (typeB1+alkyl):
                    {
                        // Type (Bs-R)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            int index_root = ii + backLen;
                            int sideSum = 0;
                            if( ii > 1 )
                                for(int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            BondID += 1;
                            BondFile << BondID << " " << BondType << " "
                            << (index_root) + (i-1)*polySize << " "
                            << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Inter: ## typeB2 + alkyl
                        //===============================================
                    case (typeB2+alkyl):
                    {
                        // Type (Bs-R)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            int index_root = ii + backLen;
                            int sideSum = 0;
                            if( ii > 1 )
                                for(int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            BondID += 1;
                            BondFile << BondID << " " << BondType << " "
                            << (index_root) + (i-1)*polySize << " "
                            << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Inter: ## typeB3 + alkyl
                        //===============================================
                    case (typeB3+alkyl):
                    {
                        // Type (Bs-R)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            int index_root = ii + backLen;
                            int sideSum = 0;
                            if( ii > 1 )
                                for(int iii=1; iii<ii; ++iii)
                                    sideSum += sideLen;
                            BondID += 1;
                            BondFile << BondID << " " << BondType << " "
                            << (index_root) + (i-1)*polySize << " "
                            << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                        }
                    }
                        break;
                        
                        
                        
                    default:
                    {
                        cout
                        << "BONDS: 'Inter' relation NOT defined."
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                        break;
                }
                //===============================================
                // End of 'Inter' Bonds
            }
            // 'trueBonds' records the number of bonds
            // taking into account that 'sideLen' might < 2
            // so that side chain bonds may reduce;
            // also, if 'sideLen' > 0,
            // there should be inter bonds
            trueBonds = BondType;
        }
        //===============================================
        // End of looping thorough all the chains
        BondFile.close();
    }
    else
    {
        cout
        << "'Bonds.txt' cannot open."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}
//===============================================
// End of BONDS





///////////////////////////////////////////////////////
//                                                   //
//  ANGLES                                           //
// ===============================================   //
//  Generate 'Angle.txt'                             //
//  Building Order:                                  //
//  Backbone --> Side Chain --> Inter                //
//  '1' based counting (rahter than '0' based)       //
//  Independent of ATOMS and BONDS sections          //
//                                                   //
//  Note:                                            //
//  By default, typeB1, typeB2, typeB3               //
//  have the same angle configuration                //
//                                                   //
///////////////////////////////////////////////////////

void StructureClass::make_AngleFile(const int n_sys)
{
    n_poly = n_polyVec[n_sys];
    sideLen = sideLenVec[n_sys];
    
    switch (typeBackbone)
    {
        case (linear):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen;
        }
            break;
            
        case (typeB1):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB2):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB3):
        {
            backLen	= backLenVec[n_sys];
            numBack = backLen*2;
        }
            break;
            
        case (typeB4):
        {
            ringBeads= 3;
            backLen	 = backLenVec[n_sys];
            numBack  = backLen*(1+ringBeads);
        }
            break;
            
        case (LJliqB):
        {
            // no angles
            return;
        }
            break;
            
        default:
        {
            cout
            << "ANGLES: 'typeBackbone' NOT defined!"
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
            break;
    }
    
    numSide	= sideLen * backLen;
    polySize= numSide + numBack;
    
    if (angleType>0) {
        
        ofstream AngleFile("Angles.txt");
        if ( AngleFile.is_open() )
        {
            AngleID = 0;
            
            //===============================================
            // Generic Looping Structure:
            // Iterate through all the chains in the system
            //===============================================
            for (int i=1; i<=n_poly; ++i)
            {
                
                
                
                /* INITIALIZE ATOM TYPE TO 1 */
                int AngleType = 1;
                
                
                
                
                //===============================================
                // 1. Backbone Angles
                //===============================================
                switch (typeBackbone)
                {
                        //===============================================
                        // Backbone: ## linear
                        //===============================================
                    case (linear):
                    {
                        // Type (B-B-B)
                        for (int ii=1; ii<backLen-1; ++ii)
                        {
                            AngleID += 1;
                            AngleFile << AngleID << " " << AngleType << " "
                            << (ii) + (i-1)*polySize << " "
                            << (ii+1) + (i-1)*polySize << " "
                            << (ii+2) + (i-1)*polySize << "\n";
                        }
                        
                        if (typeSideGroup==stdFENE) {
                            
                            // no angle in standard FENE structure by default
                            
                            --AngleType;
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Backbone: ## typeB1
                        //===============================================
                    case (typeB1):
                    {
                        // Type (B-B-B)
                        for (int ii=1; ii<backLen-1; ++ii)
                        {
                            AngleID += 1;
                            AngleFile << AngleID << " " << AngleType << " "
                            << (ii) + (i-1)*polySize << " "
                            << (ii+1) + (i-1)*polySize << " "
                            << (ii+2) + (i-1)*polySize << "\n";
                        }
                        
                        AngleType += 1;
                        
                        // Type (B-B-Bs)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            // frist bead: only backward counting
                            if ( ii == 1 )
                            {
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+1) + (i-1)*polySize << "\n";
                            }
                            else
                            {
                                // final bead: only forward counting
                                if ( ii == backLen )
                                {
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                }
                                else
                                {
                                    // backward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                    // forward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii+1) + (i-1)*polySize << "\n";
                                }
                            }
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Backbone: ## typeB2
                        //===============================================
                    case (typeB2):
                    {
                        // Type (B-B-B)
                        for (int ii=1; ii<backLen-1; ++ii)
                        {
                            AngleID += 1;
                            AngleFile << AngleID << " " << AngleType << " "
                            << (ii) + (i-1)*polySize << " "
                            << (ii+1) + (i-1)*polySize << " "
                            << (ii+2) + (i-1)*polySize << "\n";
                        }
                        
                        AngleType += 1;
                        
                        // Type (B-B-Bs)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            // frist bead: only backward counting
                            if ( ii == 1 )
                            {
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+1) + (i-1)*polySize << "\n";
                            }
                            else
                            {
                                // final bead: only forward counting
                                if ( ii == backLen )
                                {
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                }
                                else
                                {
                                    // backward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                    // forward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii+1) + (i-1)*polySize << "\n";
                                }
                            }
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Backbone: ## typeB3
                        //===============================================
                    case (typeB3):
                    {
                        // Type (B-B-B)
                        for (int ii=1; ii<backLen-1; ++ii)
                        {
                            AngleID += 1;
                            AngleFile << AngleID << " " << AngleType << " "
                            << (ii) + (i-1)*polySize << " "
                            << (ii+1) + (i-1)*polySize << " "
                            << (ii+2) + (i-1)*polySize << "\n";
                        }
                        
                        AngleType += 1;
                        
                        // Type (B-B-Bs)
                        for (int ii=1; ii<=backLen; ++ii)
                        {
                            // frist bead: only backward counting
                            if ( ii == 1 )
                            {
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+1) + (i-1)*polySize << "\n";
                            }
                            else
                            {
                                // final bead: only forward counting
                                if ( ii == backLen )
                                {
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                }
                                else
                                {
                                    // backward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii-1) + (i-1)*polySize << "\n";
                                    // forward counting
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (ii+backLen) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii+1) + (i-1)*polySize << "\n";
                                }
                            }
                        }
                    }
                        break;
                        
                        
                        
                        //===============================================
                        // Backbone: ## typeB4
                        //===============================================
                    case (typeB4):
                    {
                        //
                    }
                        break;
                        
                    default:
                    {
                        cout
                        << "ANGLES: 'backbone' relation NOT defined."
                        << "\n";
                        system("read -t5");exit(EXIT_FAILURE);
                    }
                        break;
                }
                //===============================================
                // End of 'Backbone' Angles
                
                
                
                /* index_FinalBack: index of the final backbone bead */
                // Note: NOT 'numBack-1' (here is 1-based counting)
                int index_FinalBack = numBack;
                
                
                
                //===============================================
                // 2. Side Chain Angles
                //===============================================
                // side chains should have at least 3 beads
                // to have side chain angles
                //===============================================
                if (sideLen>2)
                {
                    
                    
                    
                    AngleType += 1;
                    
                    
                    
                    switch (typeSideGroup)
                    {
                            //===============================================
                            // Side Chain: ## phenyl
                            //===============================================
                        case (phenyl):
                        {
                            // Type (R-R-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                // from 2nd side group on
                                // add up all previous side beads
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                // (left)-(center)-(right)
                                << (index_FinalBack+sideSum+2) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+sideLen) + (i-1)*polySize << "\n";
                            }
                        }
                            break;
                            
                            
                            
                            //===============================================
                            // Side Chain: ## alkyl
                            //===============================================
                        case (alkyl):
                        {
                            // Type (R-R-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                
                                for (int iii=1; iii<sideLen-1; ++iii)
                                {
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (index_FinalBack+sideSum+iii) + (i-1)*polySize << " "
                                    << (index_FinalBack+sideSum+iii+1) + (i-1)*polySize << " "
                                    << (index_FinalBack+sideSum+iii+2) + (i-1)*polySize << "\n";
                                }
                            }
                        }
                            break;
                            
                            
                            
                        default:
                        {
                            cout
                            << "ANGLES: 'side chain' relation NOT defined."
                            << "\n";
                            system("read -t5");exit(EXIT_FAILURE);
                        }
                            break;
                    }
                }
                else
                {
                    // NOTE:
                    // Here assuming side chain only has 'one' angle type
                    // in side chain itself
                    
                    static int tmp=0; // static variable here to ensure only do once
                    if(tmp==0)
                    {
                        // if sideLen=1, erase side bonds of typeB1~3
                        if ((typeBackbone==linear)&&(typeSideGroup!=stdFENE))
                        {
                            theta.erase(theta.begin()+1);
                            angCoeffs.erase(angCoeffs.begin()+1);
                            is_shakeAngles.erase(is_shakeAngles.begin()+1);
                            
                            cout << "**************************\n";
                            cout << "theta size= " << theta.size() << "\n";
                            for (size_t i=0; i<theta.size(); ++i) cout << " " << theta[i]; cout << "\n";
                            cout << "is_shakeAngles size= " << is_shakeAngles.size() << "\n";
                            for (size_t i=0; i<is_shakeAngles.size(); ++i) cout << " " << is_shakeAngles[i]; cout << "\n";
                            cout << "**************************\n";
                        }
                        else if ( (typeBackbone!=linear) && (typeBackbone!=typeB4) )
                        {
                            theta.erase(theta.begin()+2);
                            angCoeffs.erase(angCoeffs.begin()+2);
                            is_shakeAngles.erase(is_shakeAngles.begin()+2);
                            
                            cout << "**************************\n";
                            cout << "theta size= " << theta.size() << "\n";
                            for (size_t i=0; i<theta.size(); ++i) cout << " " << theta[i]; cout << "\n";
                            cout << "is_shakeAngles size= " << is_shakeAngles.size() << "\n";
                            for (size_t i=0; i<is_shakeAngles.size(); ++i) cout << " " << is_shakeAngles[i]; cout << "\n";
                            cout << "**************************\n";
                        }
                        ++tmp;
                    }
                }
                //===============================================
                // End of 'Side Chain' Angles
                
                
                //===============================================
                // 3. Inter Angles
                //===============================================
                // there should have at least 1 side bead
                // to have inter angles
                //===============================================
                if (sideLen>0)
                {
                    
                    
                    
                    AngleType += 1;
                    
                    
                    
                    
                    switch (typeBackbone+typeSideGroup)
                    {
                            //===============================================
                            // Inter: ## linear + phenyl = PS
                            //===============================================
                        case (linear+phenyl):
                        {
                            // Type (R-B-B)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                //-------------------------------------
                                // NOTE:
                                // At chain ends, i.e.
                                // ii = 1 and backLen
                                // only count one side of R-B-B;
                                // otherwise, count two sides of R-B-B.
                                //-------------------------------------
                                // frist bead: only backward counting
                                if ( ii == 1 )
                                {
                                    AngleID += 1;
                                    AngleFile << AngleID << " " << AngleType << " "
                                    << (index_FinalBack+1) + (i-1)*polySize << " "
                                    << (ii) + (i-1)*polySize << " "
                                    << (ii+1) + (i-1)*polySize << "\n";
                                }
                                else
                                {
                                    int sideSum = 0;
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                    // final bead: only forward counting
                                    if ( ii == backLen )
                                    {
                                        AngleID += 1;
                                        AngleFile << AngleID << " " << AngleType << " "
                                        << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                        << (ii) + (i-1)*polySize << " "
                                        << (ii-1) + (i-1)*polySize << "\n";
                                    }
                                    else
                                    {
                                        // forward counting
                                        AngleID += 1;
                                        AngleFile << AngleID << " " << AngleType << " "
                                        << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                        << (ii) + (i-1)*polySize << " "
                                        << (ii-1) + (i-1)*polySize << "\n";
                                        // backward counting
                                        AngleID += 1;
                                        AngleFile << AngleID << " " << AngleType << " "
                                        << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                        << (ii) + (i-1)*polySize << " "
                                        << (ii+1) + (i-1)*polySize << "\n";
                                    }
                                }
                            }
                            
                            AngleType += 1;
                            
                            // Type (B-R-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                //-------------------------------------
                                // NOTE:
                                // Every angle should count two sides
                                // of B-R-R.
                                //-------------------------------------
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+2) + (i-1)*polySize << "\n";
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+sideLen) + (i-1)*polySize << "\n";
                            }
                            
                        }
                            break;
                            
                            
                            
                            //===============================================
                            // Inter: ## typeB1 + alkyl
                            //===============================================
                        case (typeB1+alkyl):
                        {
                            // Type (B-Bs-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                            }
                        }
                            break;
                            
                            //===============================================
                            // Inter: ## typeB2 + alkyl
                            //===============================================
                        case (typeB2+alkyl):
                        {
                            // Type 4 (B-Bs-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                            }
                        }
                            break;
                            
                            
                            
                            //===============================================
                            // Inter: ## typeB3 + alkyl
                            //===============================================
                        case (typeB3+alkyl):
                        {
                            // Type 4 (B-Bs-R)
                            for (int ii=1; ii<=backLen; ++ii)
                            {
                                int sideSum = 0;
                                if( ii > 1 )
                                    for (int iii=1; iii<ii; ++iii)
                                        sideSum += sideLen;
                                AngleID += 1;
                                AngleFile << AngleID << " " << AngleType << " "
                                << (ii) + (i-1)*polySize << " "
                                << (ii+backLen) + (i-1)*polySize << " "
                                << (index_FinalBack+sideSum+1) + (i-1)*polySize << "\n";
                            }
                        }
                            break;
                            
                            
                            
                        default:
                        {
                            cout
                            << "ANGLES: 'Inter' relation NOT defined."
                            << "\n";
                            system("read -t5");exit(EXIT_FAILURE);
                        }
                            break;
                    }
                    //===============================================
                    // End of 'Inter' Angles
                }
                // 'trueAngles' records the number of angles
                // taking into account that 'sideLen' might < 3
                // so that side chain angles may reduce;
                // also, if 'sideLen' > 0,
                // there should be inter angles.
                trueAngles = AngleType;
            }
            //===============================================
            // End of looping thorough all the chains
            AngleFile.close();
        }
        else
        {
            cout
            << "'Angles.txt' cannot open."
            << "\n";
            system("read -t5");exit(EXIT_FAILURE);
        }
    }
}
//===============================================
// End of ANGLES





//===============================================
// Function: Write FF Parameters to file
//===============================================
void StructureClass::write_ForcefieldToFile()
{
    ofstream ffParam("sys_Params.txt");
    ofstream lmphead("lmp_header.txt");
    ofstream tmpMass("tmp_mass.txt");
    ofstream bondPar("bonds_Params.txt");
    ofstream anglePa("angles_Params.txt");
    ofstream pairPar("pairs_Params.txt");
    
    if (ffParam.is_open() &&
        lmphead.is_open() &&
        tmpMass.is_open() &&
        bondPar.is_open() &&
        anglePa.is_open() &&
        pairPar.is_open())
    {
        //===============================================
        // LAMMPS Data File Header
        //===============================================
        if (typeSideGroup==stdFENE) {
            lmphead << get_polyBeads()  << " atoms"       << "\n";
            lmphead << BondID           << " bonds"       << "\n";
            lmphead << get_trueAtoms()  << " atom types"  << "\n";
            lmphead << get_trueBonds()  << " bond types"  << "\n";
            // Keep a copy in the archive
            ffParam << n_poly*polySize  << " atoms"       << "\n";
            ffParam << BondID           << " bonds"       << "\n";
            ffParam << get_trueAtoms()  << " atom types"  << "\n";
            ffParam << get_trueBonds()  << " bond types"  << "\n";
        }
        else if (typeSideGroup==LJliqB) {
            lmphead << get_polyBeads()  << " atoms"       << "\n";
            lmphead << get_trueAtoms()  << " atom types"  << "\n";
            // Keep a copy in the archive
            ffParam << n_poly*polySize  << " atoms"       << "\n";
            ffParam << get_trueAtoms()  << " atom types"  << "\n";
        }
        else {
            lmphead << get_polyBeads()  << " atoms"       << "\n";
            lmphead << BondID           << " bonds"       << "\n";
            lmphead << AngleID          << " angles"      << "\n\n";
            lmphead << get_trueAtoms()  << " atom types"  << "\n";
            lmphead << get_trueBonds()  << " bond types"  << "\n";
            lmphead << get_trueAngles() << " angle types" << "\n";
            // Keep a copy in the archive
            ffParam << n_poly*polySize  << " atoms"       << "\n";
            ffParam << BondID           << " bonds"       << "\n";
            ffParam << AngleID          << " angles"      << "\n\n";
            ffParam << get_trueAtoms()  << " atom types"  << "\n";
            ffParam << get_trueBonds()  << " bond types"  << "\n";
            ffParam << get_trueAngles() << " angle types" << "\n";
        }
        
        
        //===============================================
        // Specify the FF styles
        //===============================================
        
        tmpMass << "Masses" << "\n\n";
        write_MassToFile(trueAtoms,mass,tmpMass);
        tmpMass << "\n";
        
        ffParam << "\n" << "Masses" << "\n\n";
        write_MassToFile(trueAtoms,mass,ffParam);
        
        vector<string> strVar;
        ostringstream oss;
        
        if (typeSideGroup==stdFENE) {
            
            const string bond_style ="bond_style fene";
            const string pair_style ="pair_style lj/cut";
            
            // Bond Coeffs
            ffParam << "\n" << "Bond Coeffs" << "\n\n";
            strVar.push_back(bond_style);
            strVar.push_back("bond_coeff");
            write_BondCoeffsToFile(strVar,*this,ffParam);
            write_BondCoeffsToFile(strVar,*this,bondPar);
            
            // Pair LJ Coeffs
            ffParam << "\n" << "PairIJ Coeffs" << "\n\n";
            strVar.clear();
            strVar.push_back(pair_style);
            strVar.push_back("pair_coeff");
            write_PairCoeffsToFile(strVar,*this,ffParam);
            write_PairCoeffsToFile(strVar,*this,pairPar);
        }
        else if (typeSideGroup==LJliqS) {
            
            const string pair_style ="pair_style lj/cut";
            
            // Pair LJ Coeffs
            ffParam << "\n" << "PairIJ Coeffs" << "\n\n";
            strVar.clear();
            strVar.push_back(pair_style);
            strVar.push_back("pair_coeff");
            write_PairCoeffsToFile(strVar,*this,ffParam);
            write_PairCoeffsToFile(strVar,*this,pairPar);
        }
        else {
            
            const string bond_style ="bond_style harmonic";
            const string anlge_style="angle_style harmonic";
            const string pair_style ="pair_style lj/gromacs";
            
            // Bond Coeffs
            ffParam << "\n" << "Bond Coeffs" << "\n\n";
            strVar.push_back(bond_style);
            strVar.push_back("bond_coeff");
            write_BondCoeffsToFile(strVar,*this,ffParam);
            write_BondCoeffsToFile(strVar,*this,bondPar);
            
            // Angle Coeffs
            ffParam << "\n" << "Angle Coeffs" << "\n\n";
            strVar.clear();
            strVar.push_back(anlge_style);
            strVar.push_back("angle_coeff");
            write_AngleCoeffsToFile(strVar,*this,ffParam);
            write_AngleCoeffsToFile(strVar,*this,anglePa);
            
            // Pair LJ Coeffs
            ffParam << "\n" << "PairIJ Coeffs" << "\n\n";
            strVar.clear();
            strVar.push_back(pair_style);
            strVar.push_back("pair_coeff");
            write_PairCoeffsToFile(strVar,*this,ffParam);
            write_PairCoeffsToFile(strVar,*this,pairPar);
        }
        
        
        ffParam.close();lmphead.close();tmpMass.close();
        bondPar.close();anglePa.close();pairPar.close();
    }
    else
    {
        cout
        << "'sys_Params.txt' cannot open."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}





void StructureClass::invoke_PACKMOL(const int n_sys)
{
    int    n_poly_d = n_polyVec[n_sys];
    double pkmL_d   = pkmL[n_sys];
    double hpkmL    = pkmL_d*0.5;
    double TOL      = lj_sig[0];  // TOL uses lj sigma distance
    
    vector<int> n_comp;
    
    static time_t rndseed=time(NULL);
    srand((unsigned int)rndseed);
    rndseed+=3333;
    
    int pseed=0; // in (1e5,1e6)
    do {
        pseed=abs((333333*rand())%1000000);
    } while (pseed<1e+5);
    
    ofstream pkmInput("pkm.inp");
    if ( pkmInput.is_open() ) {
        
        pkmInput
        << "tolerance "
        << fixed << setprecision(1) << TOL                              << "\n"
        << "output pkm_output.xyz"                                      << "\n"
        << "filetype xyz"                                               << "\n"
        << "seed " << pseed                                             << "\n";
        
        if (composition.size()==1) {
            
            pkmInput
            << "structure singleChain.xyz"                              << "\n"
            << "  number " << n_poly_d                                  << "\n"
            << "  inside box "
            << -hpkmL << " "
            << -hpkmL << " "
            << -hpkmL << " "
            <<  hpkmL << " "
            <<  hpkmL << " "
            <<  hpkmL                                                   << "\n"
            << "end structure"                                          << "\n";
        }
        else {
            
            for (int i=0; i<composition.size(); i++) {
                
                n_comp.push_back((int)(n_poly_d*composition[i]));
                
                pkmInput
                << "structure singleChain.xyz"                          << "\n"
                << "  number " << n_comp[i]                             << "\n"
                << "  inside box "
                << -hpkmL << " "
                << -hpkmL << " "
                << -hpkmL << " "
                <<  hpkmL << " "
                <<  hpkmL << " "
                <<  hpkmL                                               << "\n"
                << "end structure"                                      << "\n";
            }
        }
        
        pkmInput.close();
    }
    else cout << "pkm_inp: 'pkm.inp' cannot open." << "\n";
    
    
    string str="./packmol < pkm.inp";
    int ret=system(str.c_str());
    if (ret!=0) {
        cout << "Error: ./packmol << pkm.inp" << "\n";
        exit(EXIT_FAILURE);
    }
    
    
    if (composition.size()>1) {
        
        int accumu=0;
        
        string input="pkm_output.xyz";
        string output="pkm_tmp.txt";
        
        for (int i=0; i<composition.size(); ++i) {
            
            ifstream input_file(input.c_str());
            ofstream output_file(output.c_str());
            
            string lineContent;
            string tmpStr;
            vector<string> rowdata;
            int line_index=0;
            
            // 2 trash lines
            getline(input_file,lineContent);
            output_file << lineContent << "\n";
            getline(input_file,lineContent);
            output_file << lineContent << "\n";
            
            while (getline(input_file,lineContent)) {
                
                ++line_index;
                
                if ((line_index>accumu) &&
                    (line_index<=(n_comp[i]+accumu))) {
                    
                    istringstream iss(lineContent);
                    
                    rowdata.clear();
                    while(iss>>tmpStr) rowdata.push_back(tmpStr);
                    
                    output_file << to_string((long long int)(i+1)) << "  ";
                    
                    for (int ii=1; ii<(int)rowdata.size(); ++ii) {
                        output_file << rowdata[ii] << "  ";
                    }
                    output_file << "\n";
                    
                }
                else {
                    output_file << lineContent << "\n";
                }
            }
            
            accumu += n_comp[i];
            
            
            string rm="rm "+input;
            string mv="mv "+output+" "+input;
            system(rm.c_str());
            system(mv.c_str());
        }
    }
}





void StructureClass::make_LammpsDataFile(const int n_sys)
{
    const double boxL_d   = get_boxL()[n_sys];
    const double n_poly_d = get_n_polyVec()[n_sys];
    const double hboxL    = 0.5*boxL_d;
    
    ofstream writeData("start.data");
    if( writeData.is_open() ) {
        
        writeData
        << "# LAMMPS Data File"
        << "\n\n";
        
        
        /* The header section of the data file */
        //----------------------------------------------------------------------
        ifstream readFF("lmp_header.txt");
        if ( readFF.is_open() ) {
            string lineContent;
            while ( getline(readFF,lineContent) )
                writeData << lineContent << "\n";
            readFF.close();
        }
        else cout << "lammpsData: 'lmp_header.txt' cannot open." << "\n";
        
        writeData
        << "\n"
        << fixed <<  setprecision(1)
        << -hboxL << " " << hboxL << " xlo xhi" << "\n"
        << -hboxL << " " << hboxL << " ylo yhi" << "\n"
        << -hboxL << " " << hboxL << " zlo zhi" << "\n\n";
        system("rm ./lmp_header.txt");
        
        writeData.unsetf(ios_base::floatfield); // reset float precision
        
        /* Mass from tmp_mass.txt temporary file */
        ifstream tmpMass("tmp_mass.txt");
        if ( tmpMass.is_open() ) {
            string lineContent;
            while ( getline(tmpMass,lineContent) )
                writeData << lineContent << "\n";
            tmpMass.close();
        }
        else cout << "lammpsData: 'tmp_mass.txt' cannot open." << "\n";
        system("rm ./tmp_mass.txt");
        //----------------------------------------------------------------------
        
        
        
        /* Atoms Section (from PACKMOL output) */
        //----------------------------------------------------------------------
        writeData
        << "Atoms" << "\n\n";
        
        ifstream readAtom("pkm_output.xyz");
        if ( readAtom.is_open() ) {
            string lineContent;
            int intData;
            double doubleData;
            
            /* First line of the packmol output file */
            getline(readAtom,lineContent);
            istringstream lineStream(lineContent);
            lineStream >> intData;
            
            int polyBeads_d = intData;
            int polySize_d	= polyBeads_d / n_poly_d;
            
            /* Second line */
            getline(readAtom,lineContent);
            
            /*
             Body of the packmol output;
             format: atomType,x,y,z */
            int atomID = 1, molID = 1;
            while ( getline(readAtom,lineContent) )	{
                
                writeData
                << atomID << " " << molID  << " ";
                
                istringstream lineStream(lineContent);
                lineStream >> intData;
                writeData << intData << " ";
                
                while (lineStream>>doubleData) {
                    writeData
                    << fixed << setprecision(6)
                    << doubleData << " "; // x,y,z
                }
                writeData << "\n";
                
                writeData.unsetf(ios_base::floatfield); // reset float precision
                
                if ( (atomID%polySize_d)==0 ){ molID += 1; }
                atomID += 1;
            }
            readAtom.close();
        }
        else cout << "lammpsData: 'pkm_output.xyz' cannot open." << "\n";
        //----------------------------------------------------------------------
        
        
        
        /* Bonds Section */
        //----------------------------------------------------------------------
        if (typeBackbone!=LJliqB) {
            
            writeData
            << "\n"
            << "Bonds" << "\n\n";
            
            ifstream readBond("Bonds.txt");
            if ( readBond.is_open() )	{
                string lineContent;
                while ( getline(readBond,lineContent) )
                    writeData << lineContent << "\n";
                readBond.close();
            }
            else cout << "lammpsData: 'Bond.txt' cannot open." << "\n";
        }
        //----------------------------------------------------------------------
        
        
        
        /* Angles Section */
        //----------------------------------------------------------------------
        if ((typeSideGroup!=stdFENE)&&(typeSideGroup!=LJliqS)) {
            
            writeData
            << "\n"
            << "Angles" << "\n\n";
            
            ifstream readAngle("Angles.txt");
            if ( readAngle.is_open() ) {
                string lineContent;
                while ( getline(readAngle,lineContent) )
                    writeData << lineContent << "\n";
                readAngle.close();
            }
            else cout << "lammpsData: 'Angle.txt' cannot open." << "\n";
        }
        //----------------------------------------------------------------------
        
        
        
        writeData.close();
    }
    else {
        cout << "lammpsData: 'start.data' cannot open." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}





void StructureClass::echoInformation()
{
    cout << "******************************\n";
    cout << "USIC          = " << get_usic() << "\n";
    cout << "typeBackbone  = # " << get_typeB() << "\n";
    cout << "typeSideGroup = # " << get_typeS() << "\n";
    cout << "\n";
    cout << "atomType      = " << atomType << "\n";
    cout << "trueAtoms     = " << trueAtoms << "\n";  // note: 'trueAtoms'
    cout << "bondType      = " << bondType << "\n";
    cout << "trueBonds     = " << trueBonds << "\n";  // note: 'trueBonds'
    cout << "angleType     = " << angleType << "\n";
    cout << "trueAngles    = " << trueAngles << "\n"; // note: 'trueAngles'
    cout << "\n";
    
    int tmp=0;
    // check backbone atom types
    for (int i=0; i<trueAtoms; ++i) if(is_backboneType[i]) ++tmp;
    if (tmp>0)
    {
        cout << "Backbone ID   = #";
        for (int i=0; i<trueAtoms; ++i) // NOTE: use trueBonds
            if(is_backboneType[i]) cout << " " << i+1;
        cout << "\n";
    }
    else
    {
        cout << "Error: No backbone atoms!" << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    
    tmp=0;
    // check side group atom types
    for (int i=0; i<trueAtoms; ++i) if(is_sidegroupType[i]) ++tmp;
    if (tmp>0)
    {
        cout << "SideGroup ID  = #";
        for (int i=0; i<trueAtoms; ++i) // NOTE: use trueBonds
            if(is_sidegroupType[i]) cout << " " << i+1;
        cout << "\n";
    }
    else
    {
        cout << "No SideGroups" << "\n";
        // exception (<0) alreay handled in makeStructure();
    }
    cout << "\n";
    
    tmp=0;
    // check if there's SHAKE bonds
    for (int i=0; i<trueBonds; ++i) if(is_shakeBonds[i]) ++tmp;
    if (tmp>0)
    {
        cout << "SHAKE Bonds   = #";
        for (int i=0; i<trueBonds; ++i) // NOTE: use trueBonds
            if(is_shakeBonds[i]) cout << " " << i+1;
        cout << "\n";
    }
    else cout << "No SHAKE Bonds" << "\n";
    
    tmp=0;
    // check if there's SHAKE angles
    for (int i=0; i<trueAngles; ++i) if(is_shakeAngles[i]) ++tmp;
    if (tmp>0)
    {
        cout << "SHAKE Angles  = #";
        for (int i=0; i<trueAngles; ++i) // NOTE: use trueAngles
            if(is_shakeAngles[i]) cout << " " << i+1;
        cout << "\n";
    }
    else cout << "No SHAKE Angles" << "\n";
    cout << "******************************\n";
    cout << "\n";
}





void StructureClass::rotate(const double &theta_d,
                            vector<double> &xy)
{
    double rm[2][2]={0.0};
    double  a[2][1]={0.0};
    double  b[2][1]={0.0};
    
    rm[0][0] = cos(theta_d);
    rm[1][0] = sin(theta_d);
    rm[0][1] = -sin(theta_d);
    rm[1][1] = cos(theta_d);
    
    b[0][0] = xy[0];
    b[1][0] = xy[1];
    
    // Rotation of a space vector
    // xy = rm .x. xy
    for (int i=0; i<2; ++i)
        for (int j=0; j<1; ++j)
            for (int k=0; k<2; ++k)
                a[i][j] += rm[i][k]*b[k][j];
    
    xy[0] = a[0][0];
    xy[1] = a[1][0];
}





void StructureClass::write_MassToFile(const int trueAtoms,
                                      const vector<double>& vD,
                                      ofstream& writeFile)
{
    for (int i=0; i<trueAtoms; ++i)
        writeFile << i+1 << " " << fixed << setprecision(2)	<< vD[i] << "\n";
}





void StructureClass::write_BondCoeffsToFile(const vector<string>& strVar,
                                            const StructureClass& sysVar,
                                            ofstream &writeFile)
{
    const int trueBonds_d = sysVar.get_trueBonds();
    
    const string bond_style = strVar[0]; // bond_style
    const string bond_coeff = strVar[1]; // bond_coeff
    
    const vector<double>& bondCoeffs_d = sysVar.get_bondCoeffs();
    const vector<double>& bondLen_d = sysVar.get_bondLen();
    
    if (typeSideGroup==stdFENE) {
        
        writeFile << bond_style << "\n";
        
        for (int i=0; i<trueBonds_d; ++i) {
            
            writeFile
            << bond_coeff
            << " " << i+1
            << " "
            << fixed << setprecision(2);
            
            for (size_t ii=0; ii<bondCoeffs_d.size(); ++ii) {
                writeFile << bondCoeffs_d[ii] << " ";
            }
            writeFile << "\n";
        }
    }
    else {
        
        writeFile << bond_style << "\n";
        
        for (int i=0; i<trueBonds_d; ++i) {
            
            writeFile
            << bond_coeff << " " << i+1 << " "
            << fixed << setprecision(2)
            << bondCoeffs_d[i] << " "
            << bondLen_d[i]
            << "\n";
        }
    }
}





void StructureClass::write_AngleCoeffsToFile(const vector<string>& strVar,
                                             const StructureClass& sysVar,
                                             ofstream &writeFile)
{
    const int trueAngles_d = sysVar.get_trueAngles();
    
    const string anlge_style = strVar[0]; // anlge_style
    const string angle_coeff = strVar[1]; // angle_coeff
    
    const vector<double>& angCoeffs_d = sysVar.get_angCoeffs();
    const vector<double>& theta_d     = sysVar.get_theta();
    
    if ( angCoeffs_d.size() == theta_d.size() )
    {
        writeFile << anlge_style << "\n";
        
        for (int i=0; i<trueAngles_d; ++i)
            
            writeFile
            << angle_coeff
            << " " << i+1
            << " "
            << fixed << setprecision(2)
            << angCoeffs_d[i] << " "
            << theta_d[i]
            << "\n";
    }
    else
    {
        cout
        << "write_AngleCoeffsToFile: Error! array sizes are not consistent."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}





void StructureClass::write_PairCoeffsToFile(const vector<string>& strVar,
                                            const StructureClass& sysVar,
                                            ofstream &writeFile)
{
    const int trueAtoms_d     = sysVar.get_trueAtoms();
    
    const string pair_style = strVar[0]; // pair_style
    const string pair_coeff = strVar[1]; // pair_coeff
    
    const vector<double> epsilon_d = sysVar.get_lj_eps();
    const vector<double> sigma_d   = sysVar.get_lj_sig();
    
    
    if ( epsilon_d.size() == sigma_d.size() )
    {
        
        int tmp = 0;
        
        if (typeSideGroup==stdFENE) {
            
            // pair_style lj/cut
            
            const double rcut_d = sysVar.get_rcut();
            
            writeFile
            << pair_style
            << fixed << setprecision(1)
            << " " << rcut_d
            << "\n";
            
            for (int i=1; i<=trueAtoms_d; ++i) {
                for (int j=i; j<=trueAtoms_d; ++j)
                {
                    writeFile
                    << pair_coeff
                    << " "
                    << i << " "
                    << j << " "
                    
                    << setprecision(6)
                    << epsilon_d[tmp] << " "
                    << sigma_d[tmp]   << " ";
                    
                    ++tmp;
                }
            }
        }
        else if (typeSideGroup==LJliqS) {
            
            // pair_style lj/cut
            
            const double rcut_d = sysVar.get_rcut();
            
            writeFile
            << pair_style
            << fixed << setprecision(1)
            << " " << rcut_d
            << "\n";
            
            for (int i=1; i<=trueAtoms_d; ++i) {
                for (int j=i; j<=trueAtoms_d; ++j)
                {
                    writeFile
                    << pair_coeff
                    << " "
                    << i << " "
                    << j << " "
                    
                    << setprecision(6)
                    << epsilon_d[tmp] << " "
                    << sigma_d[tmp]   << " ";
                    
                    ++tmp;
                    
                    writeFile << "\n";
                }
            }
        }
        else {
            
            // pair_style lj/gromacs
            
            const double rcutl_d = sysVar.get_rcutl();
            const double rcuth_d = sysVar.get_rcuth();
            
            writeFile
            << pair_style
            << fixed << setprecision(1)
            << " " << rcutl_d
            << " " << rcuth_d
            << "\n";
            
            for (int i=1; i<=trueAtoms; ++i) {
                for (int j=i; j<=trueAtoms; ++j)
                {
                    writeFile
                    << pair_coeff << " "
                    << i << " "
                    << j << " "
                    
                    << setprecision(6)
                    << epsilon_d[tmp] << " "
                    << sigma_d[tmp]   << " "
                    
                    << setprecision(2)
                    << rcutl_d << " "
                    << rcuth_d << "\n";
                    
                    ++tmp;
                }
            }
        }
    }
    else
    {
        cout
        << "write_PairCoeffsToFile: Error! array sizes are not consistent."
        << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
}





const vector<double> StructureClass::calc_JtoCal(const int inputInt,
                                                 const vector<double>& inputVector,
                                                 string inputString)
{
    size_t vecsize = size_t(inputInt);
    if ( inputString == "pair" )
        vecsize = ((inputInt+1)*inputInt)/2;
    
    vector<double> returnVector;
    if ( vecsize != inputVector.size() )
    {
        cout << "calc_JtoCal: Check Array Size and Element." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
        returnVector.push_back(0);
    }
    else
        for (size_t i=0; i<vecsize; ++i)
        {
            returnVector.push_back(inputVector[i]);
            returnVector[i] *= return_JtoCal();
        }
    
    return returnVector;
}





const vector<double> StructureClass::calc_JtoCal_half(const int inputInt,
                                                      const vector<double>& inputVector,
                                                      string inputString)
{
    size_t vecsize = size_t(inputInt);
    if ( inputString == "pair" )
        vecsize = ((inputInt+1)*inputInt)/2;
    
    vector<double> returnVector;
    if ( vecsize != inputVector.size() )
    {
        cout << "calc_JtoCal_half: Check Array Size and Element." << "\n";
        // assign null values
        for ( size_t i=0; i<vecsize; ++i)
            returnVector.push_back(0);
    }
    else
        for (size_t i=0; i<vecsize; ++i)
        {
            returnVector.push_back(inputVector[i]);
            returnVector[i] *= return_JtoCal()*0.5;
        }
    
    return returnVector;
}





const vector<double> StructureClass::setlj(const int inputInt,const vector<double>& vD)
{
    vector<double> returnVector;
    int tmp = ((inputInt+1)*inputInt)/2;
    if (size_t(tmp)!=vD.size()){
        cout << "setlj: Check Array Size and Elements." << "\n";
        system("read -t5");exit(EXIT_FAILURE);
    }
    else {
        returnVector=vD;
    }
    return returnVector;
}





/*==( public setters )==*/

/* int */
void StructureClass::set_typeB(const int i){typeB=i;}
void StructureClass::set_typeS(const int i){typeS=i;}
void StructureClass::set_n_poly(const int i){n_poly=i;}
void StructureClass::set_n_light(const int i){n_light=i;}
void StructureClass::set_n_heavy(const int i){n_heavy=i;}
/* vector<int> */
void StructureClass::set_backLenVec(const vector<int>& vi){backLenVec=vi;}
void StructureClass::set_sideLenVec(const vector<int>& vi){sideLenVec=vi;}
void StructureClass::set_n_polyVec(const vector<int>& vi){n_polyVec=vi;}
void StructureClass::set_waveindex(const vector<int>& vi){waveindex=vi;}
void StructureClass::set_backboneTypes(const std::vector<int>& vi){backboneTypes=vi;}
void StructureClass::set_sideGroupTypes(const std::vector<int>& vi){sideGroupTypes=vi;}
void StructureClass::set_types_all(const std::vector<int>& vi){types_all=vi;}
void StructureClass::set_types_light(const std::vector<int>& vi){types_light=vi;}
void StructureClass::set_types_heavy(const std::vector<int>& vi){types_heavy=vi;}
/* vector<vector<int>> */
void StructureClass::set_n_types_all(const std::vector<std::vector<int>>& vvi)
{n_types_all=vvi;}
void StructureClass::set_n_types_light(const std::vector<std::vector<int>>& vvi)
{n_types_light=vvi;}
void StructureClass::set_n_types_heavy(const std::vector<std::vector<int>>& vvi)
{n_types_heavy=vvi;}
/* vector<double> */
void StructureClass::set_boxL(const vector<double>& vd){boxL=vd;}
void StructureClass::set_pkmL(const vector<double>& vd){pkmL=vd;}
void StructureClass::set_maxLenScale(const vector<double>& vd){maxLenScale=vd;}

/*==( public getters )==*/

/* The Naming String */
const string StructureClass::get_nameString(const int n_sys) const
{
    string namestring=
    "B"+to_string((long long int)get_typeB())+"_"+
    to_string((long long int)get_backLenVec()[n_sys])+"_"+
    "S"+to_string((long long int)get_typeS())+"_"+
    to_string((long long int)get_sideLenVec()[n_sys]);
    
    return namestring;
}
