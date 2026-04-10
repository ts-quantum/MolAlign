// created with MolAlign 1.1 (C)2026 Dr. Tobias Schulz
// IRC Trajectory as Array of Molecules
// #include "molecule.inc" into povray
// ---- use "object{molecule[i]}" in code
//declare molecule object array
#declare molecule = array[98];
//
// ---- Atom and Bond Section
//transparency
#declare trans_bd = 0;
#declare trans_atom = 0;
//atom radius
#declare atom_rad_h = 0.24;
#declare atom_rad_2 = 0.35;
#declare atom_rad_3 = 0.42;
#declare atom_rad_def = 0.5;
#declare bond_rad = 0.08;
// predefined finishes:
#declare Fin_Glassy   = finish { phong 0.9 specular 0.8 reflection 0.1 roughness 0.001 }
#declare Fin_Metallic = finish { phong 0.5 metallic 0.7 brilliance 2.0 diffuse 0.3 }
#declare Fin_Matte    = finish { phong 0.0 ambient 0.1 diffuse 0.8 }
// Define Bond Finishes
#declare Fin_Bd_Std = finish { phong 0.2 ambient 0.2 }
// Select bond Finish
#declare BdFinish = Fin_Bd_Std;
// Definde Atom finishes
#declare Fin_Atom_Std = finish { phong 0.6 specular 0.4 ambient 0.2 }
// select Atom Finish
#declare AtomFinish = Fin_Atom_Std;
///////////////////////////////////
        //Begin of Section 1  
#declare AtomsGroup_1 = union {
  sphere { <1.0516, 1.3221, -0.1095>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4355, 1.1499, -0.0525>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1222, -0.0115, 0.0105>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5062, -1.3512, 0.0412>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6892, -1.5770, 0.1199>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5405, 0.5136, -0.6688>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3032, -0.5302>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4720, 1.2743, 0.9122>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2357, -2.2010, -0.0112>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2162, 0.0161, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0271, 2.0736, -0.0516>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_1 = union {
  cylinder { <1.0516, 1.3221, -0.1095>, <-0.4355, 1.1499, -0.0525>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0516, 1.3221, -0.1095>, <1.5405, 0.5136, -0.6688>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0516, 1.3221, -0.1095>, <1.3193, 2.3032, -0.5302>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0516, 1.3221, -0.1095>, <1.4720, 1.2743, 0.9122>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1499, -0.0525>, <-1.1222, -0.0115, 0.0105>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1499, -0.0525>, <-1.0271, 2.0736, -0.0516>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1222, -0.0115, 0.0105>, <-0.5062, -1.3512, 0.0412>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1222, -0.0115, 0.0105>, <-2.2162, 0.0161, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5062, -1.3512, 0.0412>, <0.6892, -1.5770, 0.1199>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5062, -1.3512, 0.0412>, <-1.2357, -2.2010, -0.0112>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[1] = union {
    object {AtomsGroup_1}
    object {BondsGroup_1}
}
// End of Section 1
        //Begin of Section 2  
#declare AtomsGroup_2 = union {
  sphere { <1.0517, 1.3222, -0.1098>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4355, 1.1495, -0.0526>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1222, -0.0114, 0.0105>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5064, -1.3513, 0.0412>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6895, -1.5764, 0.1202>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5407, 0.5134, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3035, -0.5301>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4719, 1.2741, 0.9120>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2357, -2.2014, -0.0113>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2165, 0.0163, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0268, 2.0737, -0.0514>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_2 = union {
  cylinder { <1.0517, 1.3222, -0.1098>, <-0.4355, 1.1495, -0.0526>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0517, 1.3222, -0.1098>, <1.5407, 0.5134, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0517, 1.3222, -0.1098>, <1.3193, 2.3035, -0.5301>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0517, 1.3222, -0.1098>, <1.4719, 1.2741, 0.9120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1495, -0.0526>, <-1.1222, -0.0114, 0.0105>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1495, -0.0526>, <-1.0268, 2.0737, -0.0514>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1222, -0.0114, 0.0105>, <-0.5064, -1.3513, 0.0412>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1222, -0.0114, 0.0105>, <-2.2165, 0.0163, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5064, -1.3513, 0.0412>, <0.6895, -1.5764, 0.1202>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5064, -1.3513, 0.0412>, <-1.2357, -2.2014, -0.0113>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[2] = union {
    object {AtomsGroup_2}
    object {BondsGroup_2}
}
// End of Section 2
        //Begin of Section 3  
#declare AtomsGroup_3 = union {
  sphere { <1.0518, 1.3223, -0.1100>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4355, 1.1494, -0.0527>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1223, -0.0115, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5066, -1.3515, 0.0413>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6896, -1.5759, 0.1204>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5408, 0.5133, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3038, -0.5299>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4717, 1.2739, 0.9118>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2358, -2.2018, -0.0114>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2167, 0.0164, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0265, 2.0737, -0.0512>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_3 = union {
  cylinder { <1.0518, 1.3223, -0.1100>, <-0.4355, 1.1494, -0.0527>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0518, 1.3223, -0.1100>, <1.5408, 0.5133, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0518, 1.3223, -0.1100>, <1.3193, 2.3038, -0.5299>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0518, 1.3223, -0.1100>, <1.4717, 1.2739, 0.9118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1494, -0.0527>, <-1.1223, -0.0115, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1494, -0.0527>, <-1.0265, 2.0737, -0.0512>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1223, -0.0115, 0.0104>, <-0.5066, -1.3515, 0.0413>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1223, -0.0115, 0.0104>, <-2.2167, 0.0164, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5066, -1.3515, 0.0413>, <0.6896, -1.5759, 0.1204>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5066, -1.3515, 0.0413>, <-1.2358, -2.2018, -0.0114>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[3] = union {
    object {AtomsGroup_3}
    object {BondsGroup_3}
}
// End of Section 3
        //Begin of Section 4  
#declare AtomsGroup_4 = union {
  sphere { <1.0519, 1.3223, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4355, 1.1493, -0.0528>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1224, -0.0116, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5067, -1.3516, 0.0413>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6898, -1.5754, 0.1206>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5409, 0.5131, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3041, -0.5299>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4716, 1.2737, 0.9116>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2359, -2.2021, -0.0115>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2170, 0.0165, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0262, 2.0738, -0.0510>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_4 = union {
  cylinder { <1.0519, 1.3223, -0.1101>, <-0.4355, 1.1493, -0.0528>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0519, 1.3223, -0.1101>, <1.5409, 0.5131, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0519, 1.3223, -0.1101>, <1.3193, 2.3041, -0.5299>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0519, 1.3223, -0.1101>, <1.4716, 1.2737, 0.9116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1493, -0.0528>, <-1.1224, -0.0116, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4355, 1.1493, -0.0528>, <-1.0262, 2.0738, -0.0510>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1224, -0.0116, 0.0104>, <-0.5067, -1.3516, 0.0413>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1224, -0.0116, 0.0104>, <-2.2170, 0.0165, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5067, -1.3516, 0.0413>, <0.6898, -1.5754, 0.1206>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5067, -1.3516, 0.0413>, <-1.2359, -2.2021, -0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[4] = union {
    object {AtomsGroup_4}
    object {BondsGroup_4}
}
// End of Section 4
        //Begin of Section 5  
#declare AtomsGroup_5 = union {
  sphere { <1.0520, 1.3224, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4354, 1.1494, -0.0529>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1225, -0.0118, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5066, -1.3517, 0.0414>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6898, -1.5749, 0.1207>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5410, 0.5130, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3043, -0.5298>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4715, 1.2735, 0.9114>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2360, -2.2024, -0.0115>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2172, 0.0166, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0260, 2.0738, -0.0508>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_5 = union {
  cylinder { <1.0520, 1.3224, -0.1101>, <-0.4354, 1.1494, -0.0529>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.5410, 0.5130, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.3193, 2.3043, -0.5298>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.4715, 1.2735, 0.9114>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4354, 1.1494, -0.0529>, <-1.1225, -0.0118, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4354, 1.1494, -0.0529>, <-1.0260, 2.0738, -0.0508>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1225, -0.0118, 0.0104>, <-0.5066, -1.3517, 0.0414>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1225, -0.0118, 0.0104>, <-2.2172, 0.0166, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5066, -1.3517, 0.0414>, <0.6898, -1.5749, 0.1207>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5066, -1.3517, 0.0414>, <-1.2360, -2.2024, -0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[5] = union {
    object {AtomsGroup_5}
    object {BondsGroup_5}
}
// End of Section 5
        //Begin of Section 6  
#declare AtomsGroup_6 = union {
  sphere { <1.0520, 1.3224, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4354, 1.1496, -0.0530>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1226, -0.0121, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5065, -1.3517, 0.0414>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6897, -1.5746, 0.1208>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5411, 0.5128, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3046, -0.5297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4713, 1.2734, 0.9112>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2361, -2.2027, -0.0116>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2173, 0.0167, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0258, 2.0737, -0.0506>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_6 = union {
  cylinder { <1.0520, 1.3224, -0.1101>, <-0.4354, 1.1496, -0.0530>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.5411, 0.5128, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.3193, 2.3046, -0.5297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0520, 1.3224, -0.1101>, <1.4713, 1.2734, 0.9112>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4354, 1.1496, -0.0530>, <-1.1226, -0.0121, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4354, 1.1496, -0.0530>, <-1.0258, 2.0737, -0.0506>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1226, -0.0121, 0.0104>, <-0.5065, -1.3517, 0.0414>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1226, -0.0121, 0.0104>, <-2.2173, 0.0167, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5065, -1.3517, 0.0414>, <0.6897, -1.5746, 0.1208>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5065, -1.3517, 0.0414>, <-1.2361, -2.2027, -0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[6] = union {
    object {AtomsGroup_6}
    object {BondsGroup_6}
}
// End of Section 6
        //Begin of Section 7  
#declare AtomsGroup_7 = union {
  sphere { <1.0521, 1.3224, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4353, 1.1498, -0.0531>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1227, -0.0123, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5063, -1.3518, 0.0415>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6896, -1.5743, 0.1209>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5412, 0.5127, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3047, -0.5297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4712, 1.2733, 0.9110>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2361, -2.2029, -0.0116>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2174, 0.0167, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0256, 2.0737, -0.0505>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_7 = union {
  cylinder { <1.0521, 1.3224, -0.1101>, <-0.4353, 1.1498, -0.0531>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.5412, 0.5127, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.3193, 2.3047, -0.5297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.4712, 1.2733, 0.9110>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4353, 1.1498, -0.0531>, <-1.1227, -0.0123, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4353, 1.1498, -0.0531>, <-1.0256, 2.0737, -0.0505>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1227, -0.0123, 0.0104>, <-0.5063, -1.3518, 0.0415>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1227, -0.0123, 0.0104>, <-2.2174, 0.0167, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5063, -1.3518, 0.0415>, <0.6896, -1.5743, 0.1209>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5063, -1.3518, 0.0415>, <-1.2361, -2.2029, -0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[7] = union {
    object {AtomsGroup_7}
    object {BondsGroup_7}
}
// End of Section 7
        //Begin of Section 8  
#declare AtomsGroup_8 = union {
  sphere { <1.0521, 1.3224, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4352, 1.1499, -0.0532>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1229, -0.0125, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5061, -1.3518, 0.0415>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6894, -1.5741, 0.1210>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5413, 0.5127, -0.6687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3048, -0.5297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4712, 1.2733, 0.9109>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2362, -2.2030, -0.0116>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2175, 0.0168, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0255, 2.0737, -0.0505>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_8 = union {
  cylinder { <1.0521, 1.3224, -0.1101>, <-0.4352, 1.1499, -0.0532>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.5413, 0.5127, -0.6687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.3193, 2.3048, -0.5297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.4712, 1.2733, 0.9109>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4352, 1.1499, -0.0532>, <-1.1229, -0.0125, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4352, 1.1499, -0.0532>, <-1.0255, 2.0737, -0.0505>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, -0.0125, 0.0104>, <-0.5061, -1.3518, 0.0415>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, -0.0125, 0.0104>, <-2.2175, 0.0168, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5061, -1.3518, 0.0415>, <0.6894, -1.5741, 0.1210>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5061, -1.3518, 0.0415>, <-1.2362, -2.2030, -0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[8] = union {
    object {AtomsGroup_8}
    object {BondsGroup_8}
}
// End of Section 8
        //Begin of Section 9  
#declare AtomsGroup_9 = union {
  sphere { <1.0521, 1.3224, -0.1101>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4352, 1.1500, -0.0532>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1229, -0.0126, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5060, -1.3519, 0.0415>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6894, -1.5740, 0.1210>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5413, 0.5126, -0.6688>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3193, 2.3049, -0.5297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4712, 1.2732, 0.9109>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2362, -2.2031, -0.0116>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2176, 0.0168, 0.0369>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0255, 2.0736, -0.0504>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_9 = union {
  cylinder { <1.0521, 1.3224, -0.1101>, <-0.4352, 1.1500, -0.0532>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.5413, 0.5126, -0.6688>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.3193, 2.3049, -0.5297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0521, 1.3224, -0.1101>, <1.4712, 1.2732, 0.9109>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4352, 1.1500, -0.0532>, <-1.1229, -0.0126, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4352, 1.1500, -0.0532>, <-1.0255, 2.0736, -0.0504>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, -0.0126, 0.0104>, <-0.5060, -1.3519, 0.0415>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, -0.0126, 0.0104>, <-2.2176, 0.0168, 0.0369>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5060, -1.3519, 0.0415>, <0.6894, -1.5740, 0.1210>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5060, -1.3519, 0.0415>, <-1.2362, -2.2031, -0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[9] = union {
    object {AtomsGroup_9}
    object {BondsGroup_9}
}
// End of Section 9
        //Begin of Section 10  
#declare AtomsGroup_10 = union {
  sphere { <1.0525, 1.3223, -0.1126>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4347, 1.1487, -0.0531>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1235, -0.0120, 0.0104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5063, -1.3523, 0.0416>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6903, -1.5733, 0.1213>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5409, 0.5131, -0.6681>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3192, 2.3050, -0.5293>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4714, 1.2729, 0.9118>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2367, -2.2037, -0.0118>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2177, 0.0169, 0.0368>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0256, 2.0744, -0.0502>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_10 = union {
  cylinder { <1.0525, 1.3223, -0.1126>, <-0.4347, 1.1487, -0.0531>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0525, 1.3223, -0.1126>, <1.5409, 0.5131, -0.6681>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0525, 1.3223, -0.1126>, <1.3192, 2.3050, -0.5293>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0525, 1.3223, -0.1126>, <1.4714, 1.2729, 0.9118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4347, 1.1487, -0.0531>, <-1.1235, -0.0120, 0.0104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4347, 1.1487, -0.0531>, <-1.0256, 2.0744, -0.0502>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1235, -0.0120, 0.0104>, <-0.5063, -1.3523, 0.0416>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1235, -0.0120, 0.0104>, <-2.2177, 0.0169, 0.0368>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5063, -1.3523, 0.0416>, <0.6903, -1.5733, 0.1213>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5063, -1.3523, 0.0416>, <-1.2367, -2.2037, -0.0118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[10] = union {
    object {AtomsGroup_10}
    object {BondsGroup_10}
}
// End of Section 10
        //Begin of Section 11  
#declare AtomsGroup_11 = union {
  sphere { <1.0528, 1.3223, -0.1141>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4345, 1.1477, -0.0531>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1237, -0.0114, 0.0103>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5067, -1.3525, 0.0417>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6910, -1.5730, 0.1214>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5407, 0.5134, -0.6676>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3191, 2.3050, -0.5290>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4715, 1.2727, 0.9124>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2369, -2.2040, -0.0118>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2177, 0.0169, 0.0368>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0257, 2.0748, -0.0501>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_11 = union {
  cylinder { <1.0528, 1.3223, -0.1141>, <-0.4345, 1.1477, -0.0531>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0528, 1.3223, -0.1141>, <1.5407, 0.5134, -0.6676>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0528, 1.3223, -0.1141>, <1.3191, 2.3050, -0.5290>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0528, 1.3223, -0.1141>, <1.4715, 1.2727, 0.9124>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4345, 1.1477, -0.0531>, <-1.1237, -0.0114, 0.0103>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4345, 1.1477, -0.0531>, <-1.0257, 2.0748, -0.0501>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1237, -0.0114, 0.0103>, <-0.5067, -1.3525, 0.0417>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1237, -0.0114, 0.0103>, <-2.2177, 0.0169, 0.0368>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5067, -1.3525, 0.0417>, <0.6910, -1.5730, 0.1214>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5067, -1.3525, 0.0417>, <-1.2369, -2.2040, -0.0118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[11] = union {
    object {AtomsGroup_11}
    object {BondsGroup_11}
}
// End of Section 11
        //Begin of Section 12  
#declare AtomsGroup_12 = union {
  sphere { <1.0504, 1.3183, -0.1094>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4389, 1.1540, -0.0565>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1244, -0.0090, 0.0197>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5011, -1.3464, 0.0419>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.6966, -1.5620, 0.1159>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5335, 0.5025, -0.6648>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3278, 2.3002, -0.5284>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4607, 1.2612, 0.9143>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2267, -2.2004, -0.0224>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2199, 0.0158, 0.0551>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0285, 2.0778, -0.0686>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_12 = union {
  cylinder { <1.0504, 1.3183, -0.1094>, <-0.4389, 1.1540, -0.0565>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0504, 1.3183, -0.1094>, <1.5335, 0.5025, -0.6648>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0504, 1.3183, -0.1094>, <1.3278, 2.3002, -0.5284>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0504, 1.3183, -0.1094>, <1.4607, 1.2612, 0.9143>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4389, 1.1540, -0.0565>, <-1.1244, -0.0090, 0.0197>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4389, 1.1540, -0.0565>, <-1.0285, 2.0778, -0.0686>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1244, -0.0090, 0.0197>, <-0.5011, -1.3464, 0.0419>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1244, -0.0090, 0.0197>, <-2.2199, 0.0158, 0.0551>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5011, -1.3464, 0.0419>, <0.6966, -1.5620, 0.1159>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5011, -1.3464, 0.0419>, <-1.2267, -2.2004, -0.0224>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[12] = union {
    object {AtomsGroup_12}
    object {BondsGroup_12}
}
// End of Section 12
        //Begin of Section 13  
#declare AtomsGroup_13 = union {
  sphere { <1.0463, 1.3166, -0.1066>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4428, 1.1571, -0.0594>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1250, -0.0051, 0.0289>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4965, -1.3384, 0.0416>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7016, -1.5545, 0.1101>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5261, 0.4932, -0.6621>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3359, 2.2915, -0.5267>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4518, 1.2510, 0.9177>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2156, -2.1952, -0.0320>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2197, 0.0141, 0.0722>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0327, 2.0814, -0.0870>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_13 = union {
  cylinder { <1.0463, 1.3166, -0.1066>, <-0.4428, 1.1571, -0.0594>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0463, 1.3166, -0.1066>, <1.5261, 0.4932, -0.6621>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0463, 1.3166, -0.1066>, <1.3359, 2.2915, -0.5267>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0463, 1.3166, -0.1066>, <1.4518, 1.2510, 0.9177>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4428, 1.1571, -0.0594>, <-1.1250, -0.0051, 0.0289>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4428, 1.1571, -0.0594>, <-1.0327, 2.0814, -0.0870>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1250, -0.0051, 0.0289>, <-0.4965, -1.3384, 0.0416>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1250, -0.0051, 0.0289>, <-2.2197, 0.0141, 0.0722>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4965, -1.3384, 0.0416>, <0.7016, -1.5545, 0.1101>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4965, -1.3384, 0.0416>, <-1.2156, -2.1952, -0.0320>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[13] = union {
    object {AtomsGroup_13}
    object {BondsGroup_13}
}
// End of Section 13
        //Begin of Section 14  
#declare AtomsGroup_14 = union {
  sphere { <1.0429, 1.3109, -0.1083>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4459, 1.1614, -0.0618>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1265, -0.0015, 0.0379>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4921, -1.3314, 0.0411>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7067, -1.5489, 0.1044>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5170, 0.4885, -0.6573>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3435, 2.2826, -0.5247>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4448, 1.2424, 0.9228>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.2055, -2.1894, -0.0407>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2182, 0.0124, 0.0877>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0374, 2.0848, -0.1044>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_14 = union {
  cylinder { <1.0429, 1.3109, -0.1083>, <-0.4459, 1.1614, -0.0618>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0429, 1.3109, -0.1083>, <1.5170, 0.4885, -0.6573>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0429, 1.3109, -0.1083>, <1.3435, 2.2826, -0.5247>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0429, 1.3109, -0.1083>, <1.4448, 1.2424, 0.9228>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4459, 1.1614, -0.0618>, <-1.1265, -0.0015, 0.0379>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4459, 1.1614, -0.0618>, <-1.0374, 2.0848, -0.1044>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1265, -0.0015, 0.0379>, <-0.4921, -1.3314, 0.0411>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1265, -0.0015, 0.0379>, <-2.2182, 0.0124, 0.0877>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4921, -1.3314, 0.0411>, <0.7067, -1.5489, 0.1044>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4921, -1.3314, 0.0411>, <-1.2055, -2.1894, -0.0407>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[14] = union {
    object {AtomsGroup_14}
    object {BondsGroup_14}
}
// End of Section 14
        //Begin of Section 15  
#declare AtomsGroup_15 = union {
  sphere { <1.0401, 1.3062, -0.1062>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4498, 1.1655, -0.0643>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1253, 0.0024, 0.0460>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4865, -1.3264, 0.0409>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7113, -1.5436, 0.0990>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5089, 0.4831, -0.6536>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3505, 2.2747, -0.5228>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4381, 1.2352, 0.9249>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1980, -2.1839, -0.0485>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2182, 0.0108, 0.1014>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0419, 2.0877, -0.1201>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_15 = union {
  cylinder { <1.0401, 1.3062, -0.1062>, <-0.4498, 1.1655, -0.0643>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0401, 1.3062, -0.1062>, <1.5089, 0.4831, -0.6536>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0401, 1.3062, -0.1062>, <1.3505, 2.2747, -0.5228>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0401, 1.3062, -0.1062>, <1.4381, 1.2352, 0.9249>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4498, 1.1655, -0.0643>, <-1.1253, 0.0024, 0.0460>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4498, 1.1655, -0.0643>, <-1.0419, 2.0877, -0.1201>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1253, 0.0024, 0.0460>, <-0.4865, -1.3264, 0.0409>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1253, 0.0024, 0.0460>, <-2.2182, 0.0108, 0.1014>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4865, -1.3264, 0.0409>, <0.7113, -1.5436, 0.0990>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4865, -1.3264, 0.0409>, <-1.1980, -2.1839, -0.0485>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[15] = union {
    object {AtomsGroup_15}
    object {BondsGroup_15}
}
// End of Section 15
        //Begin of Section 16  
#declare AtomsGroup_16 = union {
  sphere { <1.0370, 1.3051, -0.1016>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4530, 1.1700, -0.0668>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1244, 0.0048, 0.0534>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4820, -1.3222, 0.0410>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7176, -1.5389, 0.0940>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5024, 0.4757, -0.6516>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3569, 2.2678, -0.5210>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4319, 1.2290, 0.9255>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1931, -2.1794, -0.0555>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2184, 0.0092, 0.1133>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0458, 2.0905, -0.1341>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_16 = union {
  cylinder { <1.0370, 1.3051, -0.1016>, <-0.4530, 1.1700, -0.0668>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0370, 1.3051, -0.1016>, <1.5024, 0.4757, -0.6516>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0370, 1.3051, -0.1016>, <1.3569, 2.2678, -0.5210>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0370, 1.3051, -0.1016>, <1.4319, 1.2290, 0.9255>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4530, 1.1700, -0.0668>, <-1.1244, 0.0048, 0.0534>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4530, 1.1700, -0.0668>, <-1.0458, 2.0905, -0.1341>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1244, 0.0048, 0.0534>, <-0.4820, -1.3222, 0.0410>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1244, 0.0048, 0.0534>, <-2.2184, 0.0092, 0.1133>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4820, -1.3222, 0.0410>, <0.7176, -1.5389, 0.0940>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4820, -1.3222, 0.0410>, <-1.1931, -2.1794, -0.0555>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[16] = union {
    object {AtomsGroup_16}
    object {BondsGroup_16}
}
// End of Section 16
        //Begin of Section 17  
#declare AtomsGroup_17 = union {
  sphere { <1.0343, 1.3047, -0.0995>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4554, 1.1735, -0.0686>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1247, 0.0064, 0.0600>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4767, -1.3191, 0.0416>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7232, -1.5345, 0.0895>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4962, 0.4686, -0.6495>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3627, 2.2629, -0.5192>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4270, 1.2236, 0.9271>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1902, -2.1756, -0.0617>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2179, 0.0078, 0.1232>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0494, 2.0932, -0.1463>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_17 = union {
  cylinder { <1.0343, 1.3047, -0.0995>, <-0.4554, 1.1735, -0.0686>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0343, 1.3047, -0.0995>, <1.4962, 0.4686, -0.6495>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0343, 1.3047, -0.0995>, <1.3627, 2.2629, -0.5192>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0343, 1.3047, -0.0995>, <1.4270, 1.2236, 0.9271>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4554, 1.1735, -0.0686>, <-1.1247, 0.0064, 0.0600>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4554, 1.1735, -0.0686>, <-1.0494, 2.0932, -0.1463>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0064, 0.0600>, <-0.4767, -1.3191, 0.0416>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0064, 0.0600>, <-2.2179, 0.0078, 0.1232>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.3191, 0.0416>, <0.7232, -1.5345, 0.0895>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.3191, 0.0416>, <-1.1902, -2.1756, -0.0617>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[17] = union {
    object {AtomsGroup_17}
    object {BondsGroup_17}
}
// End of Section 17
        //Begin of Section 18  
#declare AtomsGroup_18 = union {
  sphere { <1.0327, 1.3031, -0.1006>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4580, 1.1757, -0.0696>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1247, 0.0080, 0.0657>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4722, -1.3167, 0.0424>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7279, -1.5303, 0.0854>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4897, 0.4632, -0.6465>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3680, 2.2600, -0.5177>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4233, 1.2188, 0.9299>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1882, -2.1722, -0.0671>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2172, 0.0062, 0.1315>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0525, 2.0957, -0.1569>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_18 = union {
  cylinder { <1.0327, 1.3031, -0.1006>, <-0.4580, 1.1757, -0.0696>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0327, 1.3031, -0.1006>, <1.4897, 0.4632, -0.6465>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0327, 1.3031, -0.1006>, <1.3680, 2.2600, -0.5177>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0327, 1.3031, -0.1006>, <1.4233, 1.2188, 0.9299>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4580, 1.1757, -0.0696>, <-1.1247, 0.0080, 0.0657>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4580, 1.1757, -0.0696>, <-1.0525, 2.0957, -0.1569>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0080, 0.0657>, <-0.4722, -1.3167, 0.0424>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0080, 0.0657>, <-2.2172, 0.0062, 0.1315>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4722, -1.3167, 0.0424>, <0.7279, -1.5303, 0.0854>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4722, -1.3167, 0.0424>, <-1.1882, -2.1722, -0.0671>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[18] = union {
    object {AtomsGroup_18}
    object {BondsGroup_18}
}
// End of Section 18
        //Begin of Section 19  
#declare AtomsGroup_19 = union {
  sphere { <1.0320, 1.3007, -0.1027>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4605, 1.1772, -0.0702>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1240, 0.0092, 0.0707>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4693, -1.3148, 0.0431>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7320, -1.5266, 0.0820>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4833, 0.4590, -0.6431>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3729, 2.2589, -0.5165>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4201, 1.2148, 0.9325>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1864, -2.1690, -0.0716>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2164, 0.0046, 0.1383>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0549, 2.0976, -0.1658>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_19 = union {
  cylinder { <1.0320, 1.3007, -0.1027>, <-0.4605, 1.1772, -0.0702>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0320, 1.3007, -0.1027>, <1.4833, 0.4590, -0.6431>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0320, 1.3007, -0.1027>, <1.3729, 2.2589, -0.5165>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0320, 1.3007, -0.1027>, <1.4201, 1.2148, 0.9325>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4605, 1.1772, -0.0702>, <-1.1240, 0.0092, 0.0707>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4605, 1.1772, -0.0702>, <-1.0549, 2.0976, -0.1658>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1240, 0.0092, 0.0707>, <-0.4693, -1.3148, 0.0431>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1240, 0.0092, 0.0707>, <-2.2164, 0.0046, 0.1383>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4693, -1.3148, 0.0431>, <0.7320, -1.5266, 0.0820>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4693, -1.3148, 0.0431>, <-1.1864, -2.1690, -0.0716>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[19] = union {
    object {AtomsGroup_19}
    object {BondsGroup_19}
}
// End of Section 19
        //Begin of Section 20  
#declare AtomsGroup_20 = union {
  sphere { <1.0316, 1.2991, -0.1038>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4623, 1.1787, -0.0708>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1231, 0.0098, 0.0750>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4673, -1.3130, 0.0436>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7348, -1.5234, 0.0791>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4777, 0.4550, -0.6401>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3771, 2.2585, -0.5155>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4171, 1.2113, 0.9339>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1845, -2.1661, -0.0753>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2157, 0.0030, 0.1439>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0565, 2.0986, -0.1733>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_20 = union {
  cylinder { <1.0316, 1.2991, -0.1038>, <-0.4623, 1.1787, -0.0708>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0316, 1.2991, -0.1038>, <1.4777, 0.4550, -0.6401>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0316, 1.2991, -0.1038>, <1.3771, 2.2585, -0.5155>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0316, 1.2991, -0.1038>, <1.4171, 1.2113, 0.9339>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4623, 1.1787, -0.0708>, <-1.1231, 0.0098, 0.0750>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4623, 1.1787, -0.0708>, <-1.0565, 2.0986, -0.1733>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1231, 0.0098, 0.0750>, <-0.4673, -1.3130, 0.0436>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1231, 0.0098, 0.0750>, <-2.2157, 0.0030, 0.1439>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4673, -1.3130, 0.0436>, <0.7348, -1.5234, 0.0791>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4673, -1.3130, 0.0436>, <-1.1845, -2.1661, -0.0753>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[20] = union {
    object {AtomsGroup_20}
    object {BondsGroup_20}
}
// End of Section 20
        //Begin of Section 21  
#declare AtomsGroup_21 = union {
  sphere { <1.0310, 1.2990, -0.1034>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4631, 1.1803, -0.0715>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1224, 0.0098, 0.0789>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4658, -1.3112, 0.0439>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7363, -1.5207, 0.0766>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4730, 0.4507, -0.6378>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3807, 2.2583, -0.5146>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4141, 1.2084, 0.9339>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1825, -2.1634, -0.0784>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2149, 0.0013, 0.1485>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0576, 2.0989, -0.1794>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_21 = union {
  cylinder { <1.0310, 1.2990, -0.1034>, <-0.4631, 1.1803, -0.0715>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0310, 1.2990, -0.1034>, <1.4730, 0.4507, -0.6378>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0310, 1.2990, -0.1034>, <1.3807, 2.2583, -0.5146>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0310, 1.2990, -0.1034>, <1.4141, 1.2084, 0.9339>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4631, 1.1803, -0.0715>, <-1.1224, 0.0098, 0.0789>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4631, 1.1803, -0.0715>, <-1.0576, 2.0989, -0.1794>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1224, 0.0098, 0.0789>, <-0.4658, -1.3112, 0.0439>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1224, 0.0098, 0.0789>, <-2.2149, 0.0013, 0.1485>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4658, -1.3112, 0.0439>, <0.7363, -1.5207, 0.0766>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4658, -1.3112, 0.0439>, <-1.1825, -2.1634, -0.0784>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[21] = union {
    object {AtomsGroup_21}
    object {BondsGroup_21}
}
// End of Section 21
        //Begin of Section 22  
#declare AtomsGroup_22 = union {
  sphere { <1.0303, 1.3002, -0.1020>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4631, 1.1815, -0.0722>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1219, 0.0101, 0.0822>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4646, -1.3094, 0.0441>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7370, -1.5185, 0.0745>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4691, 0.4465, -0.6362>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3836, 2.2580, -0.5137>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4113, 1.2059, 0.9329>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1806, -2.1611, -0.0809>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2140, -0.0004, 0.1523>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0582, 2.0987, -0.1845>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_22 = union {
  cylinder { <1.0303, 1.3002, -0.1020>, <-0.4631, 1.1815, -0.0722>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0303, 1.3002, -0.1020>, <1.4691, 0.4465, -0.6362>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0303, 1.3002, -0.1020>, <1.3836, 2.2580, -0.5137>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0303, 1.3002, -0.1020>, <1.4113, 1.2059, 0.9329>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4631, 1.1815, -0.0722>, <-1.1219, 0.0101, 0.0822>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4631, 1.1815, -0.0722>, <-1.0582, 2.0987, -0.1845>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1219, 0.0101, 0.0822>, <-0.4646, -1.3094, 0.0441>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1219, 0.0101, 0.0822>, <-2.2140, -0.0004, 0.1523>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4646, -1.3094, 0.0441>, <0.7370, -1.5185, 0.0745>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4646, -1.3094, 0.0441>, <-1.1806, -2.1611, -0.0809>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[22] = union {
    object {AtomsGroup_22}
    object {BondsGroup_22}
}
// End of Section 22
        //Begin of Section 23  
#declare AtomsGroup_23 = union {
  sphere { <1.0296, 1.3021, -0.1005>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4627, 1.1821, -0.0728>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1213, 0.0107, 0.0850>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4640, -1.3076, 0.0441>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7373, -1.5168, 0.0729>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4658, 0.4427, -0.6348>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3859, 2.2575, -0.5128>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4088, 1.2039, 0.9316>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1787, -2.1592, -0.0828>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2132, -0.0020, 0.1553>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0586, 2.0981, -0.1886>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_23 = union {
  cylinder { <1.0296, 1.3021, -0.1005>, <-0.4627, 1.1821, -0.0728>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0296, 1.3021, -0.1005>, <1.4658, 0.4427, -0.6348>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0296, 1.3021, -0.1005>, <1.3859, 2.2575, -0.5128>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0296, 1.3021, -0.1005>, <1.4088, 1.2039, 0.9316>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4627, 1.1821, -0.0728>, <-1.1213, 0.0107, 0.0850>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4627, 1.1821, -0.0728>, <-1.0586, 2.0981, -0.1886>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1213, 0.0107, 0.0850>, <-0.4640, -1.3076, 0.0441>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1213, 0.0107, 0.0850>, <-2.2132, -0.0020, 0.1553>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.3076, 0.0441>, <0.7373, -1.5168, 0.0729>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.3076, 0.0441>, <-1.1787, -2.1592, -0.0828>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[23] = union {
    object {AtomsGroup_23}
    object {BondsGroup_23}
}
// End of Section 23
        //Begin of Section 24  
#declare AtomsGroup_24 = union {
  sphere { <1.0289, 1.3038, -0.0995>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4621, 1.1824, -0.0733>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1208, 0.0119, 0.0873>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4639, -1.3058, 0.0440>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7374, -1.5156, 0.0715>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4629, 0.4398, -0.6335>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3878, 2.2570, -0.5119>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4066, 1.2022, 0.9303>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1769, -2.1578, -0.0843>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2125, -0.0036, 0.1578>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0587, 2.0972, -0.1918>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_24 = union {
  cylinder { <1.0289, 1.3038, -0.0995>, <-0.4621, 1.1824, -0.0733>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0289, 1.3038, -0.0995>, <1.4629, 0.4398, -0.6335>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0289, 1.3038, -0.0995>, <1.3878, 2.2570, -0.5119>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0289, 1.3038, -0.0995>, <1.4066, 1.2022, 0.9303>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4621, 1.1824, -0.0733>, <-1.1208, 0.0119, 0.0873>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4621, 1.1824, -0.0733>, <-1.0587, 2.0972, -0.1918>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1208, 0.0119, 0.0873>, <-0.4639, -1.3058, 0.0440>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1208, 0.0119, 0.0873>, <-2.2125, -0.0036, 0.1578>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4639, -1.3058, 0.0440>, <0.7374, -1.5156, 0.0715>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4639, -1.3058, 0.0440>, <-1.1769, -2.1578, -0.0843>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[24] = union {
    object {AtomsGroup_24}
    object {BondsGroup_24}
}
// End of Section 24
        //Begin of Section 25  
#declare AtomsGroup_25 = union {
  sphere { <1.0283, 1.3049, -0.0993>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4615, 1.1825, -0.0736>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1202, 0.0132, 0.0892>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4640, -1.3042, 0.0438>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7375, -1.5148, 0.0705>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4603, 0.4378, -0.6320>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3892, 2.2566, -0.5111>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4049, 1.2008, 0.9295>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1753, -2.1566, -0.0856>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2118, -0.0050, 0.1597>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0586, 2.0962, -0.1945>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_25 = union {
  cylinder { <1.0283, 1.3049, -0.0993>, <-0.4615, 1.1825, -0.0736>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0283, 1.3049, -0.0993>, <1.4603, 0.4378, -0.6320>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0283, 1.3049, -0.0993>, <1.3892, 2.2566, -0.5111>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0283, 1.3049, -0.0993>, <1.4049, 1.2008, 0.9295>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4615, 1.1825, -0.0736>, <-1.1202, 0.0132, 0.0892>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4615, 1.1825, -0.0736>, <-1.0586, 2.0962, -0.1945>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1202, 0.0132, 0.0892>, <-0.4640, -1.3042, 0.0438>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1202, 0.0132, 0.0892>, <-2.2118, -0.0050, 0.1597>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.3042, 0.0438>, <0.7375, -1.5148, 0.0705>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.3042, 0.0438>, <-1.1753, -2.1566, -0.0856>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[25] = union {
    object {AtomsGroup_25}
    object {BondsGroup_25}
}
// End of Section 25
        //Begin of Section 26  
#declare AtomsGroup_26 = union {
  sphere { <1.0279, 1.3053, -0.0998>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4610, 1.1826, -0.0738>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1197, 0.0146, 0.0908>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4643, -1.3026, 0.0437>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7375, -1.5141, 0.0696>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4579, 0.4367, -0.6305>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3905, 2.2563, -0.5104>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4036, 1.1996, 0.9290>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1740, -2.1558, -0.0866>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2112, -0.0063, 0.1613>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0584, 2.0951, -0.1966>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_26 = union {
  cylinder { <1.0279, 1.3053, -0.0998>, <-0.4610, 1.1826, -0.0738>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0279, 1.3053, -0.0998>, <1.4579, 0.4367, -0.6305>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0279, 1.3053, -0.0998>, <1.3905, 2.2563, -0.5104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0279, 1.3053, -0.0998>, <1.4036, 1.1996, 0.9290>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4610, 1.1826, -0.0738>, <-1.1197, 0.0146, 0.0908>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4610, 1.1826, -0.0738>, <-1.0584, 2.0951, -0.1966>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1197, 0.0146, 0.0908>, <-0.4643, -1.3026, 0.0437>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1197, 0.0146, 0.0908>, <-2.2112, -0.0063, 0.1613>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4643, -1.3026, 0.0437>, <0.7375, -1.5141, 0.0696>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4643, -1.3026, 0.0437>, <-1.1740, -2.1558, -0.0866>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[26] = union {
    object {AtomsGroup_26}
    object {BondsGroup_26}
}
// End of Section 26
        //Begin of Section 27  
#declare AtomsGroup_27 = union {
  sphere { <1.0276, 1.3050, -0.1009>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4605, 1.1829, -0.0740>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1194, 0.0158, 0.0920>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4644, -1.3012, 0.0436>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7374, -1.5137, 0.0689>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4557, 0.4362, -0.6290>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3914, 2.2562, -0.5098>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4027, 1.1987, 0.9288>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1729, -2.1552, -0.0874>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2106, -0.0073, 0.1625>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0582, 2.0940, -0.1983>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_27 = union {
  cylinder { <1.0276, 1.3050, -0.1009>, <-0.4605, 1.1829, -0.0740>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0276, 1.3050, -0.1009>, <1.4557, 0.4362, -0.6290>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0276, 1.3050, -0.1009>, <1.3914, 2.2562, -0.5098>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0276, 1.3050, -0.1009>, <1.4027, 1.1987, 0.9288>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4605, 1.1829, -0.0740>, <-1.1194, 0.0158, 0.0920>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4605, 1.1829, -0.0740>, <-1.0582, 2.0940, -0.1983>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1194, 0.0158, 0.0920>, <-0.4644, -1.3012, 0.0436>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1194, 0.0158, 0.0920>, <-2.2106, -0.0073, 0.1625>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4644, -1.3012, 0.0436>, <0.7374, -1.5137, 0.0689>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4644, -1.3012, 0.0436>, <-1.1729, -2.1552, -0.0874>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[27] = union {
    object {AtomsGroup_27}
    object {BondsGroup_27}
}
// End of Section 27
        //Begin of Section 28  
#declare AtomsGroup_28 = union {
  sphere { <1.0274, 1.3043, -0.1025>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4600, 1.1833, -0.0741>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1192, 0.0169, 0.0932>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4645, -1.2998, 0.0435>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7372, -1.5133, 0.0683>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4535, 0.4363, -0.6273>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3923, 2.2562, -0.5093>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4019, 1.1978, 0.9290>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1719, -2.1548, -0.0881>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2100, -0.0083, 0.1636>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0579, 2.0928, -0.1997>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_28 = union {
  cylinder { <1.0274, 1.3043, -0.1025>, <-0.4600, 1.1833, -0.0741>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0274, 1.3043, -0.1025>, <1.4535, 0.4363, -0.6273>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0274, 1.3043, -0.1025>, <1.3923, 2.2562, -0.5093>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0274, 1.3043, -0.1025>, <1.4019, 1.1978, 0.9290>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4600, 1.1833, -0.0741>, <-1.1192, 0.0169, 0.0932>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4600, 1.1833, -0.0741>, <-1.0579, 2.0928, -0.1997>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0169, 0.0932>, <-0.4645, -1.2998, 0.0435>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0169, 0.0932>, <-2.2100, -0.0083, 0.1636>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4645, -1.2998, 0.0435>, <0.7372, -1.5133, 0.0683>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4645, -1.2998, 0.0435>, <-1.1719, -2.1548, -0.0881>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[28] = union {
    object {AtomsGroup_28}
    object {BondsGroup_28}
}
// End of Section 28
        //Begin of Section 29  
#declare AtomsGroup_29 = union {
  sphere { <1.0273, 1.3031, -0.1044>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4597, 1.1839, -0.0742>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1191, 0.0178, 0.0942>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4644, -1.2986, 0.0434>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7369, -1.5131, 0.0678>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4515, 0.4367, -0.6256>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3931, 2.2563, -0.5088>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4013, 1.1971, 0.9293>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1710, -2.1545, -0.0886>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2095, -0.0092, 0.1644>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0576, 2.0916, -0.2009>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_29 = union {
  cylinder { <1.0273, 1.3031, -0.1044>, <-0.4597, 1.1839, -0.0742>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3031, -0.1044>, <1.4515, 0.4367, -0.6256>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3031, -0.1044>, <1.3931, 2.2563, -0.5088>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3031, -0.1044>, <1.4013, 1.1971, 0.9293>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4597, 1.1839, -0.0742>, <-1.1191, 0.0178, 0.0942>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4597, 1.1839, -0.0742>, <-1.0576, 2.0916, -0.2009>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0178, 0.0942>, <-0.4644, -1.2986, 0.0434>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0178, 0.0942>, <-2.2095, -0.0092, 0.1644>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4644, -1.2986, 0.0434>, <0.7369, -1.5131, 0.0678>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4644, -1.2986, 0.0434>, <-1.1710, -2.1545, -0.0886>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[29] = union {
    object {AtomsGroup_29}
    object {BondsGroup_29}
}
// End of Section 29
        //Begin of Section 30  
#declare AtomsGroup_30 = union {
  sphere { <1.0273, 1.3018, -0.1062>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4593, 1.1846, -0.0743>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1191, 0.0184, 0.0949>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4643, -1.2975, 0.0433>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7365, -1.5129, 0.0674>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4498, 0.4374, -0.6241>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3936, 2.2565, -0.5085>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4009, 1.1966, 0.9297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1704, -2.1543, -0.0891>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2091, -0.0098, 0.1651>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0573, 2.0906, -0.2018>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_30 = union {
  cylinder { <1.0273, 1.3018, -0.1062>, <-0.4593, 1.1846, -0.0743>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3018, -0.1062>, <1.4498, 0.4374, -0.6241>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3018, -0.1062>, <1.3936, 2.2565, -0.5085>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3018, -0.1062>, <1.4009, 1.1966, 0.9297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4593, 1.1846, -0.0743>, <-1.1191, 0.0184, 0.0949>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4593, 1.1846, -0.0743>, <-1.0573, 2.0906, -0.2018>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0184, 0.0949>, <-0.4643, -1.2975, 0.0433>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0184, 0.0949>, <-2.2091, -0.0098, 0.1651>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4643, -1.2975, 0.0433>, <0.7365, -1.5129, 0.0674>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4643, -1.2975, 0.0433>, <-1.1704, -2.1543, -0.0891>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[30] = union {
    object {AtomsGroup_30}
    object {BondsGroup_30}
}
// End of Section 30
        //Begin of Section 31  
#declare AtomsGroup_31 = union {
  sphere { <1.0273, 1.3005, -0.1078>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4591, 1.1852, -0.0743>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1191, 0.0189, 0.0955>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4640, -1.2967, 0.0432>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7362, -1.5127, 0.0671>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4485, 0.4381, -0.6229>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3941, 2.2567, -0.5082>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4007, 1.1962, 0.9302>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1699, -2.1541, -0.0894>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2087, -0.0103, 0.1656>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0571, 2.0898, -0.2025>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_31 = union {
  cylinder { <1.0273, 1.3005, -0.1078>, <-0.4591, 1.1852, -0.0743>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3005, -0.1078>, <1.4485, 0.4381, -0.6229>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3005, -0.1078>, <1.3941, 2.2567, -0.5082>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.3005, -0.1078>, <1.4007, 1.1962, 0.9302>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4591, 1.1852, -0.0743>, <-1.1191, 0.0189, 0.0955>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4591, 1.1852, -0.0743>, <-1.0571, 2.0898, -0.2025>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0189, 0.0955>, <-0.4640, -1.2967, 0.0432>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1191, 0.0189, 0.0955>, <-2.2087, -0.0103, 0.1656>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.2967, 0.0432>, <0.7362, -1.5127, 0.0671>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4640, -1.2967, 0.0432>, <-1.1699, -2.1541, -0.0894>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[31] = union {
    object {AtomsGroup_31}
    object {BondsGroup_31}
}
// End of Section 31
        //Begin of Section 32  
#declare AtomsGroup_32 = union {
  sphere { <1.0273, 1.2996, -0.1089>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4589, 1.1857, -0.0743>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1192, 0.0191, 0.0959>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4638, -1.2962, 0.0432>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7359, -1.5126, 0.0669>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4475, 0.4386, -0.6220>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3944, 2.2568, -0.5080>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4006, 1.1959, 0.9305>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1696, -2.1541, -0.0896>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2084, -0.0107, 0.1659>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0569, 2.0892, -0.2029>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_32 = union {
  cylinder { <1.0273, 1.2996, -0.1089>, <-0.4589, 1.1857, -0.0743>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2996, -0.1089>, <1.4475, 0.4386, -0.6220>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2996, -0.1089>, <1.3944, 2.2568, -0.5080>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2996, -0.1089>, <1.4006, 1.1959, 0.9305>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4589, 1.1857, -0.0743>, <-1.1192, 0.0191, 0.0959>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4589, 1.1857, -0.0743>, <-1.0569, 2.0892, -0.2029>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0191, 0.0959>, <-0.4638, -1.2962, 0.0432>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0191, 0.0959>, <-2.2084, -0.0107, 0.1659>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4638, -1.2962, 0.0432>, <0.7359, -1.5126, 0.0669>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4638, -1.2962, 0.0432>, <-1.1696, -2.1541, -0.0896>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[32] = union {
    object {AtomsGroup_32}
    object {BondsGroup_32}
}
// End of Section 32
        //Begin of Section 33  
#declare AtomsGroup_33 = union {
  sphere { <1.0273, 1.2991, -0.1095>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4588, 1.1860, -0.0743>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1192, 0.0192, 0.0961>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4637, -1.2959, 0.0432>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7357, -1.5126, 0.0668>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4470, 0.4389, -0.6216>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3945, 2.2569, -0.5079>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4005, 1.1957, 0.9307>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1694, -2.1540, -0.0898>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2083, -0.0108, 0.1660>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0568, 2.0889, -0.2031>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_33 = union {
  cylinder { <1.0273, 1.2991, -0.1095>, <-0.4588, 1.1860, -0.0743>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2991, -0.1095>, <1.4470, 0.4389, -0.6216>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2991, -0.1095>, <1.3945, 2.2569, -0.5079>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0273, 1.2991, -0.1095>, <1.4005, 1.1957, 0.9307>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4588, 1.1860, -0.0743>, <-1.1192, 0.0192, 0.0961>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4588, 1.1860, -0.0743>, <-1.0568, 2.0889, -0.2031>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0192, 0.0961>, <-0.4637, -1.2959, 0.0432>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1192, 0.0192, 0.0961>, <-2.2083, -0.0108, 0.1660>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4637, -1.2959, 0.0432>, <0.7357, -1.5126, 0.0668>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4637, -1.2959, 0.0432>, <-1.1694, -2.1540, -0.0898>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[33] = union {
    object {AtomsGroup_33}
    object {BondsGroup_33}
}
// End of Section 33
        //Begin of Section 34  
#declare AtomsGroup_34 = union {
  sphere { <1.0260, 1.2979, -0.0985>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4618, 1.1884, -0.0727>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1170, 0.0239, 0.1028>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4594, -1.2951, 0.0438>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7399, -1.5141, 0.0622>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4427, 0.4446, -0.6295>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4009, 2.2400, -0.5119>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3946, 1.1952, 0.9368>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1646, -2.1512, -0.0938>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2119, -0.0127, 0.1717>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0606, 2.0943, -0.2141>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_34 = union {
  cylinder { <1.0260, 1.2979, -0.0985>, <-0.4618, 1.1884, -0.0727>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0260, 1.2979, -0.0985>, <1.4427, 0.4446, -0.6295>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0260, 1.2979, -0.0985>, <1.4009, 2.2400, -0.5119>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0260, 1.2979, -0.0985>, <1.3946, 1.1952, 0.9368>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4618, 1.1884, -0.0727>, <-1.1170, 0.0239, 0.1028>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4618, 1.1884, -0.0727>, <-1.0606, 2.0943, -0.2141>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1170, 0.0239, 0.1028>, <-0.4594, -1.2951, 0.0438>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1170, 0.0239, 0.1028>, <-2.2119, -0.0127, 0.1717>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4594, -1.2951, 0.0438>, <0.7399, -1.5141, 0.0622>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4594, -1.2951, 0.0438>, <-1.1646, -2.1512, -0.0938>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[34] = union {
    object {AtomsGroup_34}
    object {BondsGroup_34}
}
// End of Section 34
        //Begin of Section 35  
#declare AtomsGroup_35 = union {
  sphere { <1.0240, 1.2987, -0.0860>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4650, 1.1916, -0.0713>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1171, 0.0252, 0.1083>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4554, -1.2953, 0.0446>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7448, -1.5158, 0.0582>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4412, 0.4483, -0.6394>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4064, 2.2246, -0.5164>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3900, 1.1953, 0.9430>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1612, -2.1487, -0.0972>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2146, -0.0133, 0.1761>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0643, 2.1008, -0.2233>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_35 = union {
  cylinder { <1.0240, 1.2987, -0.0860>, <-0.4650, 1.1916, -0.0713>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0240, 1.2987, -0.0860>, <1.4412, 0.4483, -0.6394>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0240, 1.2987, -0.0860>, <1.4064, 2.2246, -0.5164>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0240, 1.2987, -0.0860>, <1.3900, 1.1953, 0.9430>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4650, 1.1916, -0.0713>, <-1.1171, 0.0252, 0.1083>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4650, 1.1916, -0.0713>, <-1.0643, 2.1008, -0.2233>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1171, 0.0252, 0.1083>, <-0.4554, -1.2953, 0.0446>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1171, 0.0252, 0.1083>, <-2.2146, -0.0133, 0.1761>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4554, -1.2953, 0.0446>, <0.7448, -1.5158, 0.0582>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4554, -1.2953, 0.0446>, <-1.1612, -2.1487, -0.0972>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[35] = union {
    object {AtomsGroup_35}
    object {BondsGroup_35}
}
// End of Section 35
        //Begin of Section 36  
#declare AtomsGroup_36 = union {
  sphere { <1.0213, 1.2985, -0.0750>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4691, 1.1964, -0.0698>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1195, 0.0235, 0.1132>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4509, -1.2958, 0.0457>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7498, -1.5176, 0.0546>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4412, 0.4518, -0.6499>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4119, 2.2109, -0.5216>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3868, 1.1957, 0.9512>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1588, -2.1465, -0.1002>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2160, -0.0130, 0.1798>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0678, 2.1076, -0.2317>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_36 = union {
  cylinder { <1.0213, 1.2985, -0.0750>, <-0.4691, 1.1964, -0.0698>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0213, 1.2985, -0.0750>, <1.4412, 0.4518, -0.6499>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0213, 1.2985, -0.0750>, <1.4119, 2.2109, -0.5216>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0213, 1.2985, -0.0750>, <1.3868, 1.1957, 0.9512>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4691, 1.1964, -0.0698>, <-1.1195, 0.0235, 0.1132>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4691, 1.1964, -0.0698>, <-1.0678, 2.1076, -0.2317>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1195, 0.0235, 0.1132>, <-0.4509, -1.2958, 0.0457>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1195, 0.0235, 0.1132>, <-2.2160, -0.0130, 0.1798>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4509, -1.2958, 0.0457>, <0.7498, -1.5176, 0.0546>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4509, -1.2958, 0.0457>, <-1.1588, -2.1465, -0.1002>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[36] = union {
    object {AtomsGroup_36}
    object {BondsGroup_36}
}
// End of Section 36
        //Begin of Section 37  
#declare AtomsGroup_37 = union {
  sphere { <1.0184, 1.2957, -0.0680>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4747, 1.2012, -0.0680>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1226, 0.0208, 0.1173>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4458, -1.2966, 0.0471>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7543, -1.5192, 0.0513>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4414, 0.4562, -0.6593>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4173, 2.1996, -0.5270>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3852, 1.1961, 0.9620>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1573, -2.1447, -0.1029>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2164, -0.0119, 0.1828>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0709, 2.1143, -0.2391>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_37 = union {
  cylinder { <1.0184, 1.2957, -0.0680>, <-0.4747, 1.2012, -0.0680>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0184, 1.2957, -0.0680>, <1.4414, 0.4562, -0.6593>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0184, 1.2957, -0.0680>, <1.4173, 2.1996, -0.5270>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0184, 1.2957, -0.0680>, <1.3852, 1.1961, 0.9620>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4747, 1.2012, -0.0680>, <-1.1226, 0.0208, 0.1173>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4747, 1.2012, -0.0680>, <-1.0709, 2.1143, -0.2391>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1226, 0.0208, 0.1173>, <-0.4458, -1.2966, 0.0471>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1226, 0.0208, 0.1173>, <-2.2164, -0.0119, 0.1828>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4458, -1.2966, 0.0471>, <0.7543, -1.5192, 0.0513>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4458, -1.2966, 0.0471>, <-1.1573, -2.1447, -0.1029>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[37] = union {
    object {AtomsGroup_37}
    object {BondsGroup_37}
}
// End of Section 37
        //Begin of Section 38  
#declare AtomsGroup_38 = union {
  sphere { <1.0159, 1.2904, -0.0657>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4822, 1.2042, -0.0656>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1247, 0.0193, 0.1204>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4410, -1.2974, 0.0486>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7583, -1.5206, 0.0486>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4414, 0.4617, -0.6665>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4225, 2.1911, -0.5322>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3850, 1.1965, 0.9743>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1565, -2.1434, -0.1052>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2164, -0.0106, 0.1851>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0735, 2.1205, -0.2454>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_38 = union {
  cylinder { <1.0159, 1.2904, -0.0657>, <-0.4822, 1.2042, -0.0656>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0159, 1.2904, -0.0657>, <1.4414, 0.4617, -0.6665>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0159, 1.2904, -0.0657>, <1.4225, 2.1911, -0.5322>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0159, 1.2904, -0.0657>, <1.3850, 1.1965, 0.9743>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4822, 1.2042, -0.0656>, <-1.1247, 0.0193, 0.1204>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4822, 1.2042, -0.0656>, <-1.0735, 2.1205, -0.2454>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0193, 0.1204>, <-0.4410, -1.2974, 0.0486>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1247, 0.0193, 0.1204>, <-2.2164, -0.0106, 0.1851>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4410, -1.2974, 0.0486>, <0.7583, -1.5206, 0.0486>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4410, -1.2974, 0.0486>, <-1.1565, -2.1434, -0.1052>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[38] = union {
    object {AtomsGroup_38}
    object {BondsGroup_38}
}
// End of Section 38
        //Begin of Section 39  
#declare AtomsGroup_39 = union {
  sphere { <1.0142, 1.2836, -0.0672>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4904, 1.2045, -0.0629>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1249, 0.0198, 0.1224>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4373, -1.2983, 0.0500>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7617, -1.5217, 0.0464>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4409, 0.4676, -0.6713>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4270, 2.1854, -0.5366>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3858, 1.1968, 0.9863>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1562, -2.1426, -0.1071>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2165, -0.0093, 0.1868>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0754, 2.1257, -0.2505>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_39 = union {
  cylinder { <1.0142, 1.2836, -0.0672>, <-0.4904, 1.2045, -0.0629>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0142, 1.2836, -0.0672>, <1.4409, 0.4676, -0.6713>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0142, 1.2836, -0.0672>, <1.4270, 2.1854, -0.5366>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0142, 1.2836, -0.0672>, <1.3858, 1.1968, 0.9863>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4904, 1.2045, -0.0629>, <-1.1249, 0.0198, 0.1224>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4904, 1.2045, -0.0629>, <-1.0754, 2.1257, -0.2505>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1249, 0.0198, 0.1224>, <-0.4373, -1.2983, 0.0500>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1249, 0.0198, 0.1224>, <-2.2165, -0.0093, 0.1868>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4373, -1.2983, 0.0500>, <0.7617, -1.5217, 0.0464>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4373, -1.2983, 0.0500>, <-1.1562, -2.1426, -0.1071>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[39] = union {
    object {AtomsGroup_39}
    object {BondsGroup_39}
}
// End of Section 39
        //Begin of Section 40  
#declare AtomsGroup_40 = union {
  sphere { <1.0134, 1.2776, -0.0700>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4975, 1.2032, -0.0606>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1239, 0.0216, 0.1234>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4351, -1.2989, 0.0510>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7643, -1.5224, 0.0449>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4402, 0.4725, -0.6738>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4302, 2.1821, -0.5397>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3867, 1.1970, 0.9956>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1560, -2.1421, -0.1084>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2167, -0.0084, 0.1879>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0768, 2.1294, -0.2539>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_40 = union {
  cylinder { <1.0134, 1.2776, -0.0700>, <-0.4975, 1.2032, -0.0606>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0134, 1.2776, -0.0700>, <1.4402, 0.4725, -0.6738>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0134, 1.2776, -0.0700>, <1.4302, 2.1821, -0.5397>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0134, 1.2776, -0.0700>, <1.3867, 1.1970, 0.9956>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4975, 1.2032, -0.0606>, <-1.1239, 0.0216, 0.1234>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4975, 1.2032, -0.0606>, <-1.0768, 2.1294, -0.2539>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1239, 0.0216, 0.1234>, <-0.4351, -1.2989, 0.0510>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1239, 0.0216, 0.1234>, <-2.2167, -0.0084, 0.1879>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4351, -1.2989, 0.0510>, <0.7643, -1.5224, 0.0449>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4351, -1.2989, 0.0510>, <-1.1560, -2.1421, -0.1084>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[40] = union {
    object {AtomsGroup_40}
    object {BondsGroup_40}
}
// End of Section 40
        //Begin of Section 41  
#declare AtomsGroup_41 = union {
  sphere { <1.0132, 1.2741, -0.0720>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.5016, 1.2019, -0.0593>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1229, 0.0231, 0.1238>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4341, -1.2993, 0.0516>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7658, -1.5227, 0.0442>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4397, 0.4753, -0.6748>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4319, 2.1805, -0.5413>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3872, 1.1971, 1.0006>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1560, -2.1419, -0.1090>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2169, -0.0079, 0.1885>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0775, 2.1314, -0.2557>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_41 = union {
  cylinder { <1.0132, 1.2741, -0.0720>, <-0.5016, 1.2019, -0.0593>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0132, 1.2741, -0.0720>, <1.4397, 0.4753, -0.6748>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0132, 1.2741, -0.0720>, <1.4319, 2.1805, -0.5413>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <1.0132, 1.2741, -0.0720>, <1.3872, 1.1971, 1.0006>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5016, 1.2019, -0.0593>, <-1.1229, 0.0231, 0.1238>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.5016, 1.2019, -0.0593>, <-1.0775, 2.1314, -0.2557>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, 0.0231, 0.1238>, <-0.4341, -1.2993, 0.0516>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1229, 0.0231, 0.1238>, <-2.2169, -0.0079, 0.1885>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4341, -1.2993, 0.0516>, <0.7658, -1.5227, 0.0442>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4341, -1.2993, 0.0516>, <-1.1560, -2.1419, -0.1090>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[41] = union {
    object {AtomsGroup_41}
    object {BondsGroup_41}
}
// End of Section 41
        //Begin of Section 42  
#declare AtomsGroup_42 = union {
  sphere { <0.9829, 1.2824, -0.0362>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4999, 1.2178, -0.0554>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1299, 0.0352, 0.1493>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4173, -1.2686, 0.0509>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7813, -1.4831, 0.0068>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3935, 0.4177, -0.6388>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4305, 2.1566, -0.5382>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4007, 1.1634, 0.9792>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1122, -2.1292, -0.1174>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2278, -0.0039, 0.1998>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0730, 2.1233, -0.3035>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_42 = union {
  cylinder { <0.9829, 1.2824, -0.0362>, <-0.4999, 1.2178, -0.0554>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9829, 1.2824, -0.0362>, <1.3935, 0.4177, -0.6388>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9829, 1.2824, -0.0362>, <1.4305, 2.1566, -0.5382>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9829, 1.2824, -0.0362>, <1.4007, 1.1634, 0.9792>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4999, 1.2178, -0.0554>, <-1.1299, 0.0352, 0.1493>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4999, 1.2178, -0.0554>, <-1.0730, 2.1233, -0.3035>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1299, 0.0352, 0.1493>, <-0.4173, -1.2686, 0.0509>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1299, 0.0352, 0.1493>, <-2.2278, -0.0039, 0.1998>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4173, -1.2686, 0.0509>, <0.7813, -1.4831, 0.0068>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4173, -1.2686, 0.0509>, <-1.1122, -2.1292, -0.1174>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[42] = union {
    object {AtomsGroup_42}
    object {BondsGroup_42}
}
// End of Section 42
        //Begin of Section 43  
#declare AtomsGroup_43 = union {
  sphere { <0.9803, 1.2836, -0.0330>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4965, 1.2254, -0.0510>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1371, 0.0481, 0.1703>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4139, -1.2411, 0.0476>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.7976, -1.4412, -0.0252>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3399, 0.3710, -0.5845>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4206, 2.1413, -0.5252>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4132, 1.1309, 0.9550>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0723, -2.1167, -0.1224>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2339, -0.0026, 0.2079>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0693, 2.1127, -0.3429>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_43 = union {
  cylinder { <0.9803, 1.2836, -0.0330>, <-0.4965, 1.2254, -0.0510>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9803, 1.2836, -0.0330>, <1.3399, 0.3710, -0.5845>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9803, 1.2836, -0.0330>, <1.4206, 2.1413, -0.5252>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9803, 1.2836, -0.0330>, <1.4132, 1.1309, 0.9550>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4965, 1.2254, -0.0510>, <-1.1371, 0.0481, 0.1703>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4965, 1.2254, -0.0510>, <-1.0693, 2.1127, -0.3429>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1371, 0.0481, 0.1703>, <-0.4139, -1.2411, 0.0476>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1371, 0.0481, 0.1703>, <-2.2339, -0.0026, 0.2079>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4139, -1.2411, 0.0476>, <0.7976, -1.4412, -0.0252>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4139, -1.2411, 0.0476>, <-1.0723, -2.1167, -0.1224>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[43] = union {
    object {AtomsGroup_43}
    object {BondsGroup_43}
}
// End of Section 43
        //Begin of Section 44  
#declare AtomsGroup_44 = union {
  sphere { <0.9895, 1.2679, -0.0555>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4958, 1.2248, -0.0447>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1410, 0.0597, 0.1859>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4052, -1.2185, 0.0429>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8019, -1.3957, -0.0509>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2857, 0.3360, -0.5228>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4112, 2.1418, -0.5117>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4270, 1.1010, 0.9401>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0411, -2.1071, -0.1252>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2362, -0.0023, 0.2133>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0674, 2.1039, -0.3748>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_44 = union {
  cylinder { <0.9895, 1.2679, -0.0555>, <-0.4958, 1.2248, -0.0447>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9895, 1.2679, -0.0555>, <1.2857, 0.3360, -0.5228>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9895, 1.2679, -0.0555>, <1.4112, 2.1418, -0.5117>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9895, 1.2679, -0.0555>, <1.4270, 1.1010, 0.9401>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4958, 1.2248, -0.0447>, <-1.1410, 0.0597, 0.1859>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4958, 1.2248, -0.0447>, <-1.0674, 2.1039, -0.3748>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0597, 0.1859>, <-0.4052, -1.2185, 0.0429>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0597, 0.1859>, <-2.2362, -0.0023, 0.2133>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4052, -1.2185, 0.0429>, <0.8019, -1.3957, -0.0509>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4052, -1.2185, 0.0429>, <-1.0411, -2.1071, -0.1252>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[44] = union {
    object {AtomsGroup_44}
    object {BondsGroup_44}
}
// End of Section 44
        //Begin of Section 45  
#declare AtomsGroup_45 = union {
  sphere { <0.9918, 1.2508, -0.0719>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4930, 1.2263, -0.0408>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1436, 0.0622, 0.1970>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4010, -1.1997, 0.0390>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8064, -1.3514, -0.0711>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2429, 0.2964, -0.4733>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4074, 2.1570, -0.5045>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4389, 1.0766, 0.9297>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0199, -2.1007, -0.1264>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2357, -0.0022, 0.2169>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0657, 2.0962, -0.3979>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_45 = union {
  cylinder { <0.9918, 1.2508, -0.0719>, <-0.4930, 1.2263, -0.0408>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9918, 1.2508, -0.0719>, <1.2429, 0.2964, -0.4733>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9918, 1.2508, -0.0719>, <1.4074, 2.1570, -0.5045>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9918, 1.2508, -0.0719>, <1.4389, 1.0766, 0.9297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4930, 1.2263, -0.0408>, <-1.1436, 0.0622, 0.1970>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4930, 1.2263, -0.0408>, <-1.0657, 2.0962, -0.3979>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1436, 0.0622, 0.1970>, <-0.4010, -1.1997, 0.0390>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1436, 0.0622, 0.1970>, <-2.2357, -0.0022, 0.2169>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4010, -1.1997, 0.0390>, <0.8064, -1.3514, -0.0711>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4010, -1.1997, 0.0390>, <-1.0199, -2.1007, -0.1264>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[45] = union {
    object {AtomsGroup_45}
    object {BondsGroup_45}
}
// End of Section 45
        //Begin of Section 46  
#declare AtomsGroup_46 = union {
  sphere { <0.9866, 1.2499, -0.0743>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4891, 1.2294, -0.0398>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1449, 0.0599, 0.2054>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4036, -1.1838, 0.0360>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8155, -1.3082, -0.0887>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2110, 0.2418, -0.4376>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4065, 2.1759, -0.5009>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4485, 1.0560, 0.9190>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0046, -2.0962, -0.1266>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2343, -0.0021, 0.2198>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0631, 2.0887, -0.4156>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_46 = union {
  cylinder { <0.9866, 1.2499, -0.0743>, <-0.4891, 1.2294, -0.0398>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9866, 1.2499, -0.0743>, <1.2110, 0.2418, -0.4376>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9866, 1.2499, -0.0743>, <1.4065, 2.1759, -0.5009>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9866, 1.2499, -0.0743>, <1.4485, 1.0560, 0.9190>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4891, 1.2294, -0.0398>, <-1.1449, 0.0599, 0.2054>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4891, 1.2294, -0.0398>, <-1.0631, 2.0887, -0.4156>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1449, 0.0599, 0.2054>, <-0.4036, -1.1838, 0.0360>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1449, 0.0599, 0.2054>, <-2.2343, -0.0021, 0.2198>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4036, -1.1838, 0.0360>, <0.8155, -1.3082, -0.0887>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4036, -1.1838, 0.0360>, <-1.0046, -2.0962, -0.1266>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[46] = union {
    object {AtomsGroup_46}
    object {BondsGroup_46}
}
// End of Section 46
        //Begin of Section 47  
#declare AtomsGroup_47 = union {
  sphere { <0.9784, 1.2659, -0.0658>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4860, 1.2311, -0.0407>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1444, 0.0572, 0.2109>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4064, -1.1729, 0.0338>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8234, -1.2703, -0.1027>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1901, 0.1824, -0.4158>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4067, 2.1909, -0.4994>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4547, 1.0403, 0.9082>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9949, -2.0937, -0.1261>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2333, -0.0021, 0.2221>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0599, 2.0824, -0.4277>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_47 = union {
  cylinder { <0.9784, 1.2659, -0.0658>, <-0.4860, 1.2311, -0.0407>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9784, 1.2659, -0.0658>, <1.1901, 0.1824, -0.4158>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9784, 1.2659, -0.0658>, <1.4067, 2.1909, -0.4994>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9784, 1.2659, -0.0658>, <1.4547, 1.0403, 0.9082>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4860, 1.2311, -0.0407>, <-1.1444, 0.0572, 0.2109>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4860, 1.2311, -0.0407>, <-1.0599, 2.0824, -0.4277>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1444, 0.0572, 0.2109>, <-0.4064, -1.1729, 0.0338>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1444, 0.0572, 0.2109>, <-2.2333, -0.0021, 0.2221>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4064, -1.1729, 0.0338>, <0.8234, -1.2703, -0.1027>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4064, -1.1729, 0.0338>, <-0.9949, -2.0937, -0.1261>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8234, -1.2703, -0.1027>, <1.1901, 0.1824, -0.4158>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[47] = union {
    object {AtomsGroup_47}
    object {BondsGroup_47}
}
// End of Section 47
        //Begin of Section 48  
#declare AtomsGroup_48 = union {
  sphere { <0.9707, 1.2866, -0.0548>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4826, 1.2314, -0.0427>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1435, 0.0547, 0.2142>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4062, -1.1663, 0.0323>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8271, -1.2408, -0.1125>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1773, 0.1334, -0.4031>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4070, 2.1997, -0.4986>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4583, 1.0296, 0.8990>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9898, -2.0927, -0.1255>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2330, -0.0021, 0.2238>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0570, 2.0779, -0.4353>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_48 = union {
  cylinder { <0.9707, 1.2866, -0.0548>, <-0.4826, 1.2314, -0.0427>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9707, 1.2866, -0.0548>, <1.1773, 0.1334, -0.4031>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9707, 1.2866, -0.0548>, <1.4070, 2.1997, -0.4986>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9707, 1.2866, -0.0548>, <1.4583, 1.0296, 0.8990>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4826, 1.2314, -0.0427>, <-1.1435, 0.0547, 0.2142>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4826, 1.2314, -0.0427>, <-1.0570, 2.0779, -0.4353>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1435, 0.0547, 0.2142>, <-0.4062, -1.1663, 0.0323>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1435, 0.0547, 0.2142>, <-2.2330, -0.0021, 0.2238>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4062, -1.1663, 0.0323>, <0.8271, -1.2408, -0.1125>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4062, -1.1663, 0.0323>, <-0.9898, -2.0927, -0.1255>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8271, -1.2408, -0.1125>, <1.1773, 0.1334, -0.4031>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[48] = union {
    object {AtomsGroup_48}
    object {BondsGroup_48}
}
// End of Section 48
        //Begin of Section 49  
#declare AtomsGroup_49 = union {
  sphere { <0.9636, 1.3039, -0.0464>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4770, 1.2319, -0.0449>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1430, 0.0507, 0.2162>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4081, -1.1617, 0.0319>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8309, -1.2207, -0.1193>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1699, 0.1001, -0.3957>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4071, 2.2039, -0.4980>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4604, 1.0228, 0.8928>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9874, -2.0926, -0.1249>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2329, -0.0021, 0.2249>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0551, 2.0752, -0.4397>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_49 = union {
  cylinder { <0.9636, 1.3039, -0.0464>, <-0.4770, 1.2319, -0.0449>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9636, 1.3039, -0.0464>, <1.1699, 0.1001, -0.3957>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9636, 1.3039, -0.0464>, <1.4071, 2.2039, -0.4980>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9636, 1.3039, -0.0464>, <1.4604, 1.0228, 0.8928>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4770, 1.2319, -0.0449>, <-1.1430, 0.0507, 0.2162>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4770, 1.2319, -0.0449>, <-1.0551, 2.0752, -0.4397>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1430, 0.0507, 0.2162>, <-0.4081, -1.1617, 0.0319>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1430, 0.0507, 0.2162>, <-2.2329, -0.0021, 0.2249>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4081, -1.1617, 0.0319>, <0.8309, -1.2207, -0.1193>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4081, -1.1617, 0.0319>, <-0.9874, -2.0926, -0.1249>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8309, -1.2207, -0.1193>, <1.1699, 0.1001, -0.3957>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[49] = union {
    object {AtomsGroup_49}
    object {BondsGroup_49}
}
// End of Section 49
        //Begin of Section 50  
#declare AtomsGroup_50 = union {
  sphere { <0.9573, 1.3151, -0.0418>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4702, 1.2332, -0.0468>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1429, 0.0455, 0.2175>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4142, -1.1579, 0.0327>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8369, -1.2089, -0.1239>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1660, 0.0812, -0.3917>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4070, 2.2054, -0.4976>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4618, 1.0188, 0.8896>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9863, -2.0927, -0.1246>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2329, -0.0022, 0.2255>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0543, 2.0739, -0.4422>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_50 = union {
  cylinder { <0.9573, 1.3151, -0.0418>, <-0.4702, 1.2332, -0.0468>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9573, 1.3151, -0.0418>, <1.1660, 0.0812, -0.3917>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9573, 1.3151, -0.0418>, <1.4070, 2.2054, -0.4976>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9573, 1.3151, -0.0418>, <1.4618, 1.0188, 0.8896>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4702, 1.2332, -0.0468>, <-1.1429, 0.0455, 0.2175>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4702, 1.2332, -0.0468>, <-1.0543, 2.0739, -0.4422>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1429, 0.0455, 0.2175>, <-0.4142, -1.1579, 0.0327>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1429, 0.0455, 0.2175>, <-2.2329, -0.0022, 0.2255>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4142, -1.1579, 0.0327>, <0.8369, -1.2089, -0.1239>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4142, -1.1579, 0.0327>, <-0.9863, -2.0927, -0.1246>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8369, -1.2089, -0.1239>, <1.1660, 0.0812, -0.3917>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[50] = union {
    object {AtomsGroup_50}
    object {BondsGroup_50}
}
// End of Section 50
        //Begin of Section 51  
#declare AtomsGroup_51 = union {
  sphere { <0.9534, 1.3203, -0.0400>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4654, 1.2343, -0.0479>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1429, 0.0416, 0.2182>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4198, -1.1556, 0.0336>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8418, -1.2038, -0.1262>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1645, 0.0732, -0.3900>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4069, 2.2058, -0.4973>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4625, 1.0171, 0.8884>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9858, -2.0928, -0.1244>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2328, -0.0023, 0.2257>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0540, 2.0736, -0.4433>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_51 = union {
  cylinder { <0.9534, 1.3203, -0.0400>, <-0.4654, 1.2343, -0.0479>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9534, 1.3203, -0.0400>, <1.1645, 0.0732, -0.3900>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9534, 1.3203, -0.0400>, <1.4069, 2.2058, -0.4973>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9534, 1.3203, -0.0400>, <1.4625, 1.0171, 0.8884>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4654, 1.2343, -0.0479>, <-1.1429, 0.0416, 0.2182>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4654, 1.2343, -0.0479>, <-1.0540, 2.0736, -0.4433>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1429, 0.0416, 0.2182>, <-0.4198, -1.1556, 0.0336>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1429, 0.0416, 0.2182>, <-2.2328, -0.0023, 0.2257>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4198, -1.1556, 0.0336>, <0.8418, -1.2038, -0.1262>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4198, -1.1556, 0.0336>, <-0.9858, -2.0928, -0.1244>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8418, -1.2038, -0.1262>, <1.1645, 0.0732, -0.3900>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[51] = union {
    object {AtomsGroup_51}
    object {BondsGroup_51}
}
// End of Section 51
        //Begin of Section 52  
#declare AtomsGroup_52 = union {
  sphere { <0.9362, 1.3609, -0.0302>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4588, 1.2530, -0.0470>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1410, 0.0412, 0.2182>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4269, -1.1350, 0.0300>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8615, -1.1715, -0.1383>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1431, -0.0784, -0.3609>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3928, 2.2431, -0.4832>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4726, 1.0032, 0.8621>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9617, -2.0829, -0.1279>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2305, -0.0011, 0.2305>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0590, 2.0788, -0.4566>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_52 = union {
  cylinder { <0.9362, 1.3609, -0.0302>, <-0.4588, 1.2530, -0.0470>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9362, 1.3609, -0.0302>, <1.1431, -0.0784, -0.3609>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9362, 1.3609, -0.0302>, <1.3928, 2.2431, -0.4832>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9362, 1.3609, -0.0302>, <1.4726, 1.0032, 0.8621>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4588, 1.2530, -0.0470>, <-1.1410, 0.0412, 0.2182>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4588, 1.2530, -0.0470>, <-1.0590, 2.0788, -0.4566>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0412, 0.2182>, <-0.4269, -1.1350, 0.0300>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0412, 0.2182>, <-2.2305, -0.0011, 0.2305>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4269, -1.1350, 0.0300>, <0.8615, -1.1715, -0.1383>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4269, -1.1350, 0.0300>, <-0.9617, -2.0829, -0.1279>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8615, -1.1715, -0.1383>, <1.1431, -0.0784, -0.3609>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[52] = union {
    object {AtomsGroup_52}
    object {BondsGroup_52}
}
// End of Section 52
        //Begin of Section 53  
#declare AtomsGroup_53 = union {
  sphere { <0.9301, 1.3741, -0.0240>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4527, 1.2551, -0.0490>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1410, 0.0365, 0.2185>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4338, -1.1326, 0.0310>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8670, -1.1645, -0.1425>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1418, -0.0949, -0.3583>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3927, 2.2449, -0.4828>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4737, 0.9997, 0.8582>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9609, -2.0834, -0.1269>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2305, -0.0012, 0.2312>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0581, 2.0777, -0.4586>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_53 = union {
  cylinder { <0.9301, 1.3741, -0.0240>, <-0.4527, 1.2551, -0.0490>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9301, 1.3741, -0.0240>, <1.1418, -0.0949, -0.3583>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9301, 1.3741, -0.0240>, <1.3927, 2.2449, -0.4828>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9301, 1.3741, -0.0240>, <1.4737, 0.9997, 0.8582>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4527, 1.2551, -0.0490>, <-1.1410, 0.0365, 0.2185>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4527, 1.2551, -0.0490>, <-1.0581, 2.0777, -0.4586>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0365, 0.2185>, <-0.4338, -1.1326, 0.0310>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0365, 0.2185>, <-2.2305, -0.0012, 0.2312>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4338, -1.1326, 0.0310>, <0.8670, -1.1645, -0.1425>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4338, -1.1326, 0.0310>, <-0.9609, -2.0834, -0.1269>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8670, -1.1645, -0.1425>, <1.1418, -0.0949, -0.3583>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[53] = union {
    object {AtomsGroup_53}
    object {BondsGroup_53}
}
// End of Section 53
        //Begin of Section 54  
#declare AtomsGroup_54 = union {
  sphere { <0.9235, 1.3943, -0.0123>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4486, 1.2572, -0.0516>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1412, 0.0338, 0.2189>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4373, -1.1308, 0.0309>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8707, -1.1532, -0.1488>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1392, -0.1238, -0.3539>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3931, 2.2498, -0.4825>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4752, 0.9939, 0.8509>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9594, -2.0843, -0.1250>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2305, -0.0012, 0.2325>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0565, 2.0756, -0.4625>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_54 = union {
  cylinder { <0.9235, 1.3943, -0.0123>, <-0.4486, 1.2572, -0.0516>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9235, 1.3943, -0.0123>, <1.1392, -0.1238, -0.3539>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9235, 1.3943, -0.0123>, <1.3931, 2.2498, -0.4825>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9235, 1.3943, -0.0123>, <1.4752, 0.9939, 0.8509>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4486, 1.2572, -0.0516>, <-1.1412, 0.0338, 0.2189>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4486, 1.2572, -0.0516>, <-1.0565, 2.0756, -0.4625>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1412, 0.0338, 0.2189>, <-0.4373, -1.1308, 0.0309>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1412, 0.0338, 0.2189>, <-2.2305, -0.0012, 0.2325>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4373, -1.1308, 0.0309>, <0.8707, -1.1532, -0.1488>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4373, -1.1308, 0.0309>, <-0.9594, -2.0843, -0.1250>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8707, -1.1532, -0.1488>, <1.1392, -0.1238, -0.3539>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[54] = union {
    object {AtomsGroup_54}
    object {BondsGroup_54}
}
// End of Section 54
        //Begin of Section 55  
#declare AtomsGroup_55 = union {
  sphere { <0.9159, 1.4215, 0.0060>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4472, 1.2601, -0.0541>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1410, 0.0329, 0.2191>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4383, -1.1286, 0.0297>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8730, -1.1409, -0.1572>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1354, -0.1646, -0.3477>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3942, 2.2598, -0.4824>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4777, 0.9853, 0.8399>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9565, -2.0857, -0.1219>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2305, -0.0010, 0.2346>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0543, 2.0725, -0.4691>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_55 = union {
  cylinder { <0.9159, 1.4215, 0.0060>, <-0.4472, 1.2601, -0.0541>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9159, 1.4215, 0.0060>, <1.1354, -0.1646, -0.3477>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9159, 1.4215, 0.0060>, <1.3942, 2.2598, -0.4824>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9159, 1.4215, 0.0060>, <1.4777, 0.9853, 0.8399>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4472, 1.2601, -0.0541>, <-1.1410, 0.0329, 0.2191>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4472, 1.2601, -0.0541>, <-1.0543, 2.0725, -0.4691>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0329, 0.2191>, <-0.4383, -1.1286, 0.0297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1410, 0.0329, 0.2191>, <-2.2305, -0.0010, 0.2346>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4383, -1.1286, 0.0297>, <0.8730, -1.1409, -0.1572>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4383, -1.1286, 0.0297>, <-0.9565, -2.0857, -0.1219>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8730, -1.1409, -0.1572>, <1.1354, -0.1646, -0.3477>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[55] = union {
    object {AtomsGroup_55}
    object {BondsGroup_55}
}
// End of Section 55
        //Begin of Section 56  
#declare AtomsGroup_56 = union {
  sphere { <0.9049, 1.4557, 0.0303>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4461, 1.2655, -0.0558>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1400, 0.0316, 0.2185>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4413, -1.1243, 0.0278>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8741, -1.1390, -0.1670>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1328, -0.2092, -0.3412>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3958, 2.2761, -0.4820>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4822, 0.9735, 0.8254>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9511, -2.0870, -0.1171>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2306, -0.0007, 0.2376>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0524, 2.0693, -0.4797>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_56 = union {
  cylinder { <0.9049, 1.4557, 0.0303>, <-0.4461, 1.2655, -0.0558>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9049, 1.4557, 0.0303>, <1.3958, 2.2761, -0.4820>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.9049, 1.4557, 0.0303>, <1.4822, 0.9735, 0.8254>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4461, 1.2655, -0.0558>, <-1.1400, 0.0316, 0.2185>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4461, 1.2655, -0.0558>, <-1.0524, 2.0693, -0.4797>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1400, 0.0316, 0.2185>, <-0.4413, -1.1243, 0.0278>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1400, 0.0316, 0.2185>, <-2.2306, -0.0007, 0.2376>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4413, -1.1243, 0.0278>, <0.8741, -1.1390, -0.1670>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4413, -1.1243, 0.0278>, <-0.9511, -2.0870, -0.1171>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8741, -1.1390, -0.1670>, <1.1328, -0.2092, -0.3412>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[56] = union {
    object {AtomsGroup_56}
    object {BondsGroup_56}
}
// End of Section 56
        //Begin of Section 57  
#declare AtomsGroup_57 = union {
  sphere { <0.8916, 1.4945, 0.0583>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4444, 1.2739, -0.0559>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1388, 0.0307, 0.2168>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4475, -1.1201, 0.0254>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8716, -1.1606, -0.1758>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1358, -0.2420, -0.3372>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3968, 2.2972, -0.4797>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4886, 0.9595, 0.8080>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9434, -2.0879, -0.1102>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2301, -0.0003, 0.2418>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0519, 2.0667, -0.4948>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_57 = union {
  cylinder { <0.8916, 1.4945, 0.0583>, <-0.4444, 1.2739, -0.0559>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8916, 1.4945, 0.0583>, <1.3968, 2.2972, -0.4797>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8916, 1.4945, 0.0583>, <1.4886, 0.9595, 0.8080>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4444, 1.2739, -0.0559>, <-1.1388, 0.0307, 0.2168>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4444, 1.2739, -0.0559>, <-1.0519, 2.0667, -0.4948>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1388, 0.0307, 0.2168>, <-0.4475, -1.1201, 0.0254>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1388, 0.0307, 0.2168>, <-2.2301, -0.0003, 0.2418>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4475, -1.1201, 0.0254>, <0.8716, -1.1606, -0.1758>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4475, -1.1201, 0.0254>, <-0.9434, -2.0879, -0.1102>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8716, -1.1606, -0.1758>, <1.1358, -0.2420, -0.3372>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[57] = union {
    object {AtomsGroup_57}
    object {BondsGroup_57}
}
// End of Section 57
        //Begin of Section 58  
#declare AtomsGroup_58 = union {
  sphere { <0.8814, 1.5336, 0.0879>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4447, 1.2846, -0.0541>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1377, 0.0338, 0.2139>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4548, -1.1211, 0.0222>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8645, -1.2020, -0.1835>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1457, -0.2611, -0.3358>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3956, 2.3215, -0.4736>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.4957, 0.9457, 0.7877>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9357, -2.0889, -0.1008>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2287, 0.0001, 0.2468>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0529, 2.0651, -0.5140>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_58 = union {
  cylinder { <0.8814, 1.5336, 0.0879>, <-0.4447, 1.2846, -0.0541>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8814, 1.5336, 0.0879>, <1.3956, 2.3215, -0.4736>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8814, 1.5336, 0.0879>, <1.4957, 0.9457, 0.7877>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4447, 1.2846, -0.0541>, <-1.1377, 0.0338, 0.2139>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4447, 1.2846, -0.0541>, <-1.0529, 2.0651, -0.5140>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1377, 0.0338, 0.2139>, <-0.4548, -1.1211, 0.0222>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1377, 0.0338, 0.2139>, <-2.2287, 0.0001, 0.2468>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4548, -1.1211, 0.0222>, <0.8645, -1.2020, -0.1835>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4548, -1.1211, 0.0222>, <-0.9357, -2.0889, -0.1008>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8645, -1.2020, -0.1835>, <1.1457, -0.2611, -0.3358>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[58] = union {
    object {AtomsGroup_58}
    object {BondsGroup_58}
}
// End of Section 58
        //Begin of Section 59  
#declare AtomsGroup_59 = union {
  sphere { <0.8776, 1.5706, 0.1194>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4488, 1.2993, -0.0495>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1342, 0.0426, 0.2085>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4611, -1.1297, 0.0181>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8568, -1.2437, -0.1941>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1609, -0.2860, -0.3347>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3896, 2.3512, -0.4587>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5027, 0.9343, 0.7617>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9331, -2.0936, -0.0872>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2259, 0.0013, 0.2537>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0563, 2.0651, -0.5404>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_59 = union {
  cylinder { <0.8776, 1.5706, 0.1194>, <-0.4488, 1.2993, -0.0495>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8776, 1.5706, 0.1194>, <1.3896, 2.3512, -0.4587>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8776, 1.5706, 0.1194>, <1.5027, 0.9343, 0.7617>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4488, 1.2993, -0.0495>, <-1.1342, 0.0426, 0.2085>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4488, 1.2993, -0.0495>, <-1.0563, 2.0651, -0.5404>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1342, 0.0426, 0.2085>, <-0.4611, -1.1297, 0.0181>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1342, 0.0426, 0.2085>, <-2.2259, 0.0013, 0.2537>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4611, -1.1297, 0.0181>, <0.8568, -1.2437, -0.1941>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4611, -1.1297, 0.0181>, <-0.9331, -2.0936, -0.0872>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8568, -1.2437, -0.1941>, <1.1609, -0.2860, -0.3347>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[59] = union {
    object {AtomsGroup_59}
    object {BondsGroup_59}
}
// End of Section 59
        //Begin of Section 60  
#declare AtomsGroup_60 = union {
  sphere { <0.8773, 1.5694, 0.1186>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4480, 1.3001, -0.0493>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1334, 0.0417, 0.2080>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4609, -1.1289, 0.0182>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8576, -1.2404, -0.1949>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1609, -0.2894, -0.3344>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3888, 2.3520, -0.4571>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5028, 0.9350, 0.7614>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9348, -2.0949, -0.0866>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2256, 0.0015, 0.2539>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0565, 2.0653, -0.5411>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_60 = union {
  cylinder { <0.8773, 1.5694, 0.1186>, <-0.4480, 1.3001, -0.0493>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8773, 1.5694, 0.1186>, <1.3888, 2.3520, -0.4571>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8773, 1.5694, 0.1186>, <1.5028, 0.9350, 0.7614>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4480, 1.3001, -0.0493>, <-1.1334, 0.0417, 0.2080>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4480, 1.3001, -0.0493>, <-1.0565, 2.0653, -0.5411>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1334, 0.0417, 0.2080>, <-0.4609, -1.1289, 0.0182>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1334, 0.0417, 0.2080>, <-2.2256, 0.0015, 0.2539>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4609, -1.1289, 0.0182>, <0.8576, -1.2404, -0.1949>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4609, -1.1289, 0.0182>, <-0.9348, -2.0949, -0.0866>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8576, -1.2404, -0.1949>, <1.1609, -0.2894, -0.3344>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[60] = union {
    object {AtomsGroup_60}
    object {BondsGroup_60}
}
// End of Section 60
        //Begin of Section 61  
#declare AtomsGroup_61 = union {
  sphere { <0.8768, 1.5670, 0.1169>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4468, 1.3016, -0.0489>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1319, 0.0403, 0.2071>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4604, -1.1277, 0.0185>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8589, -1.2345, -0.1962>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1612, -0.2955, -0.3341>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3873, 2.3534, -0.4539>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5029, 0.9363, 0.7608>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9379, -2.0972, -0.0853>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2249, 0.0019, 0.2543>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0569, 2.0656, -0.5424>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_61 = union {
  cylinder { <0.8768, 1.5670, 0.1169>, <-0.4468, 1.3016, -0.0489>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8768, 1.5670, 0.1169>, <1.3873, 2.3534, -0.4539>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8768, 1.5670, 0.1169>, <1.5029, 0.9363, 0.7608>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4468, 1.3016, -0.0489>, <-1.1319, 0.0403, 0.2071>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4468, 1.3016, -0.0489>, <-1.0569, 2.0656, -0.5424>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1319, 0.0403, 0.2071>, <-0.4604, -1.1277, 0.0185>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1319, 0.0403, 0.2071>, <-2.2249, 0.0019, 0.2543>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4604, -1.1277, 0.0185>, <0.8589, -1.2345, -0.1962>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4604, -1.1277, 0.0185>, <-0.9379, -2.0972, -0.0853>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8589, -1.2345, -0.1962>, <1.1612, -0.2955, -0.3341>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[61] = union {
    object {AtomsGroup_61}
    object {BondsGroup_61}
}
// End of Section 61
        //Begin of Section 62  
#declare AtomsGroup_62 = union {
  sphere { <0.8764, 1.5639, 0.1144>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4455, 1.3036, -0.0485>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1300, 0.0389, 0.2057>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4599, -1.1268, 0.0187>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8603, -1.2276, -0.1980>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1621, -0.3031, -0.3338>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3849, 2.3557, -0.4490>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5032, 0.9381, 0.7599>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9423, -2.1000, -0.0833>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2237, 0.0026, 0.2550>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0573, 2.0661, -0.5443>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_62 = union {
  cylinder { <0.8764, 1.5639, 0.1144>, <-0.4455, 1.3036, -0.0485>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8764, 1.5639, 0.1144>, <1.3849, 2.3557, -0.4490>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8764, 1.5639, 0.1144>, <1.5032, 0.9381, 0.7599>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4455, 1.3036, -0.0485>, <-1.1300, 0.0389, 0.2057>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4455, 1.3036, -0.0485>, <-1.0573, 2.0661, -0.5443>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1300, 0.0389, 0.2057>, <-0.4599, -1.1268, 0.0187>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1300, 0.0389, 0.2057>, <-2.2237, 0.0026, 0.2550>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4599, -1.1268, 0.0187>, <0.8603, -1.2276, -0.1980>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4599, -1.1268, 0.0187>, <-0.9423, -2.1000, -0.0833>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8603, -1.2276, -0.1980>, <1.1621, -0.3031, -0.3338>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[62] = union {
    object {AtomsGroup_62}
    object {BondsGroup_62}
}
// End of Section 62
        //Begin of Section 63  
#declare AtomsGroup_63 = union {
  sphere { <0.8764, 1.5606, 0.1111>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4447, 1.3056, -0.0483>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1282, 0.0382, 0.2040>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4593, -1.1270, 0.0186>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8611, -1.2225, -0.2000>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1643, -0.3100, -0.3339>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3818, 2.3590, -0.4423>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5038, 0.9403, 0.7587>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9473, -2.1027, -0.0804>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2218, 0.0034, 0.2559>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0577, 2.0667, -0.5468>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_63 = union {
  cylinder { <0.8764, 1.5606, 0.1111>, <-0.4447, 1.3056, -0.0483>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8764, 1.5606, 0.1111>, <1.3818, 2.3590, -0.4423>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8764, 1.5606, 0.1111>, <1.5038, 0.9403, 0.7587>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4447, 1.3056, -0.0483>, <-1.1282, 0.0382, 0.2040>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4447, 1.3056, -0.0483>, <-1.0577, 2.0667, -0.5468>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1282, 0.0382, 0.2040>, <-0.4593, -1.1270, 0.0186>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1282, 0.0382, 0.2040>, <-2.2218, 0.0034, 0.2559>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4593, -1.1270, 0.0186>, <0.8611, -1.2225, -0.2000>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4593, -1.1270, 0.0186>, <-0.9473, -2.1027, -0.0804>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8611, -1.2225, -0.2000>, <1.1643, -0.3100, -0.3339>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[63] = union {
    object {AtomsGroup_63}
    object {BondsGroup_63}
}
// End of Section 63
        //Begin of Section 64  
#declare AtomsGroup_64 = union {
  sphere { <0.8769, 1.5576, 0.1074>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4448, 1.3074, -0.0484>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1267, 0.0383, 0.2021>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4590, -1.1285, 0.0183>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8609, -1.2211, -0.2019>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1681, -0.3145, -0.3345>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3781, 2.3631, -0.4345>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5045, 0.9426, 0.7574>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9524, -2.1048, -0.0767>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2194, 0.0043, 0.2571>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0579, 2.0672, -0.5496>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_64 = union {
  cylinder { <0.8769, 1.5576, 0.1074>, <-0.4448, 1.3074, -0.0484>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8769, 1.5576, 0.1074>, <1.3781, 2.3631, -0.4345>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8769, 1.5576, 0.1074>, <1.5045, 0.9426, 0.7574>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4448, 1.3074, -0.0484>, <-1.1267, 0.0383, 0.2021>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4448, 1.3074, -0.0484>, <-1.0579, 2.0672, -0.5496>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1267, 0.0383, 0.2021>, <-0.4590, -1.1285, 0.0183>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1267, 0.0383, 0.2021>, <-2.2194, 0.0043, 0.2571>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4590, -1.1285, 0.0183>, <0.8609, -1.2211, -0.2019>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4590, -1.1285, 0.0183>, <-0.9524, -2.1048, -0.0767>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8609, -1.2211, -0.2019>, <1.1681, -0.3145, -0.3345>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[64] = union {
    object {AtomsGroup_64}
    object {BondsGroup_64}
}
// End of Section 64
        //Begin of Section 65  
#declare AtomsGroup_65 = union {
  sphere { <0.8776, 1.5550, 0.1036>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4457, 1.3085, -0.0489>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1254, 0.0390, 0.2001>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4591, -1.1309, 0.0178>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8597, -1.2241, -0.2036>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1733, -0.3159, -0.3357>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3742, 2.3682, -0.4258>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5056, 0.9447, 0.7560>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9572, -2.1060, -0.0725>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2167, 0.0053, 0.2585>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0581, 2.0677, -0.5526>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_65 = union {
  cylinder { <0.8776, 1.5550, 0.1036>, <-0.4457, 1.3085, -0.0489>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8776, 1.5550, 0.1036>, <1.3742, 2.3682, -0.4258>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8776, 1.5550, 0.1036>, <1.5056, 0.9447, 0.7560>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4457, 1.3085, -0.0489>, <-1.1254, 0.0390, 0.2001>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4457, 1.3085, -0.0489>, <-1.0581, 2.0677, -0.5526>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1254, 0.0390, 0.2001>, <-0.4591, -1.1309, 0.0178>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1254, 0.0390, 0.2001>, <-2.2167, 0.0053, 0.2585>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4591, -1.1309, 0.0178>, <0.8597, -1.2241, -0.2036>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4591, -1.1309, 0.0178>, <-0.9572, -2.1060, -0.0725>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8597, -1.2241, -0.2036>, <1.1733, -0.3159, -0.3357>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[65] = union {
    object {AtomsGroup_65}
    object {BondsGroup_65}
}
// End of Section 65
        //Begin of Section 66  
#declare AtomsGroup_66 = union {
  sphere { <0.8783, 1.5530, 0.0998>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4470, 1.3090, -0.0500>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1242, 0.0400, 0.1980>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4599, -1.1336, 0.0170>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8579, -1.2307, -0.2051>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1794, -0.3152, -0.3374>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3703, 2.3743, -0.4167>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5070, 0.9464, 0.7546>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9616, -2.1064, -0.0678>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2139, 0.0063, 0.2600>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0580, 2.0681, -0.5557>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_66 = union {
  cylinder { <0.8783, 1.5530, 0.0998>, <-0.4470, 1.3090, -0.0500>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8783, 1.5530, 0.0998>, <1.3703, 2.3743, -0.4167>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8783, 1.5530, 0.0998>, <1.5070, 0.9464, 0.7546>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4470, 1.3090, -0.0500>, <-1.1242, 0.0400, 0.1980>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4470, 1.3090, -0.0500>, <-1.0580, 2.0681, -0.5557>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1242, 0.0400, 0.1980>, <-0.4599, -1.1336, 0.0170>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1242, 0.0400, 0.1980>, <-2.2139, 0.0063, 0.2600>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4599, -1.1336, 0.0170>, <0.8579, -1.2307, -0.2051>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4599, -1.1336, 0.0170>, <-0.9616, -2.1064, -0.0678>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8579, -1.2307, -0.2051>, <1.1794, -0.3152, -0.3374>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[66] = union {
    object {AtomsGroup_66}
    object {BondsGroup_66}
}
// End of Section 66
        //Begin of Section 67  
#declare AtomsGroup_67 = union {
  sphere { <0.8783, 1.5514, 0.0960>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4485, 1.3087, -0.0521>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1228, 0.0409, 0.1955>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4616, -1.1361, 0.0161>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8560, -1.2406, -0.2070>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1867, -0.3140, -0.3396>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3660, 2.3827, -0.4062>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5090, 0.9479, 0.7531>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9661, -2.1058, -0.0619>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2109, 0.0076, 0.2619>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0577, 2.0686, -0.5592>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_67 = union {
  cylinder { <0.8783, 1.5514, 0.0960>, <-0.4485, 1.3087, -0.0521>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8783, 1.5514, 0.0960>, <1.3660, 2.3827, -0.4062>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8783, 1.5514, 0.0960>, <1.5090, 0.9479, 0.7531>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4485, 1.3087, -0.0521>, <-1.1228, 0.0409, 0.1955>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4485, 1.3087, -0.0521>, <-1.0577, 2.0686, -0.5592>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1228, 0.0409, 0.1955>, <-0.4616, -1.1361, 0.0161>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1228, 0.0409, 0.1955>, <-2.2109, 0.0076, 0.2619>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4616, -1.1361, 0.0161>, <0.8560, -1.2406, -0.2070>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4616, -1.1361, 0.0161>, <-0.9661, -2.1058, -0.0619>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8560, -1.2406, -0.2070>, <1.1867, -0.3140, -0.3396>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[67] = union {
    object {AtomsGroup_67}
    object {BondsGroup_67}
}
// End of Section 67
        //Begin of Section 68  
#declare AtomsGroup_68 = union {
  sphere { <0.8774, 1.5502, 0.0924>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4495, 1.3076, -0.0557>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1209, 0.0408, 0.1925>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4641, -1.1378, 0.0151>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8550, -1.2504, -0.2097>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1937, -0.3159, -0.3419>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3614, 2.3942, -0.3941>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5115, 0.9491, 0.7513>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9710, -2.1045, -0.0547>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2080, 0.0091, 0.2644>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0573, 2.0692, -0.5630>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_68 = union {
  cylinder { <0.8774, 1.5502, 0.0924>, <-0.4495, 1.3076, -0.0557>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8774, 1.5502, 0.0924>, <1.3614, 2.3942, -0.3941>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8774, 1.5502, 0.0924>, <1.5115, 0.9491, 0.7513>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4495, 1.3076, -0.0557>, <-1.1209, 0.0408, 0.1925>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4495, 1.3076, -0.0557>, <-1.0573, 2.0692, -0.5630>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1209, 0.0408, 0.1925>, <-0.4641, -1.1378, 0.0151>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1209, 0.0408, 0.1925>, <-2.2080, 0.0091, 0.2644>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4641, -1.1378, 0.0151>, <0.8550, -1.2504, -0.2097>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4641, -1.1378, 0.0151>, <-0.9710, -2.1045, -0.0547>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8550, -1.2504, -0.2097>, <1.1937, -0.3159, -0.3419>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[68] = union {
    object {AtomsGroup_68}
    object {BondsGroup_68}
}
// End of Section 68
        //Begin of Section 69  
#declare AtomsGroup_69 = union {
  sphere { <0.8759, 1.5496, 0.0895>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4494, 1.3058, -0.0610>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1190, 0.0397, 0.1891>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4659, -1.1384, 0.0140>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8554, -1.2572, -0.2135>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.1992, -0.3239, -0.3440>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3565, 2.4083, -0.3805>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5141, 0.9500, 0.7490>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9765, -2.1030, -0.0461>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2054, 0.0107, 0.2672>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0566, 2.0696, -0.5670>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_69 = union {
  cylinder { <0.8759, 1.5496, 0.0895>, <-0.4494, 1.3058, -0.0610>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8759, 1.5496, 0.0895>, <1.3565, 2.4083, -0.3805>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8759, 1.5496, 0.0895>, <1.5141, 0.9500, 0.7490>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4494, 1.3058, -0.0610>, <-1.1190, 0.0397, 0.1891>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4494, 1.3058, -0.0610>, <-1.0566, 2.0696, -0.5670>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1190, 0.0397, 0.1891>, <-0.4659, -1.1384, 0.0140>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1190, 0.0397, 0.1891>, <-2.2054, 0.0107, 0.2672>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4659, -1.1384, 0.0140>, <0.8554, -1.2572, -0.2135>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4659, -1.1384, 0.0140>, <-0.9765, -2.1030, -0.0461>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8554, -1.2572, -0.2135>, <1.1992, -0.3239, -0.3440>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[69] = union {
    object {AtomsGroup_69}
    object {BondsGroup_69}
}
// End of Section 69
        //Begin of Section 70  
#declare AtomsGroup_70 = union {
  sphere { <0.8749, 1.5508, 0.0876>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4481, 1.3038, -0.0682>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1176, 0.0381, 0.1852>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4657, -1.1385, 0.0127>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8565, -1.2599, -0.2184>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2028, -0.3392, -0.3460>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3507, 2.4245, -0.3651>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5166, 0.9508, 0.7456>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9827, -2.1015, -0.0362>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2033, 0.0125, 0.2706>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0557, 2.0698, -0.5710>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_70 = union {
  cylinder { <0.8749, 1.5508, 0.0876>, <-0.4481, 1.3038, -0.0682>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8749, 1.5508, 0.0876>, <1.3507, 2.4245, -0.3651>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8749, 1.5508, 0.0876>, <1.5166, 0.9508, 0.7456>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4481, 1.3038, -0.0682>, <-1.1176, 0.0381, 0.1852>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4481, 1.3038, -0.0682>, <-1.0557, 2.0698, -0.5710>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1176, 0.0381, 0.1852>, <-0.4657, -1.1385, 0.0127>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1176, 0.0381, 0.1852>, <-2.2033, 0.0125, 0.2706>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4657, -1.1385, 0.0127>, <0.8565, -1.2599, -0.2184>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4657, -1.1385, 0.0127>, <-0.9827, -2.1015, -0.0362>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8565, -1.2599, -0.2184>, <1.2028, -0.3392, -0.3460>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[70] = union {
    object {AtomsGroup_70}
    object {BondsGroup_70}
}
// End of Section 70
        //Begin of Section 71  
#declare AtomsGroup_71 = union {
  sphere { <0.8753, 1.5555, 0.0874>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4453, 1.3021, -0.0777>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1168, 0.0367, 0.1802>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4629, -1.1392, 0.0112>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8563, -1.2619, -0.2241>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2053, -0.3600, -0.3484>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3434, 2.4426, -0.3472>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5183, 0.9516, 0.7405>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9890, -2.0997, -0.0246>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2020, 0.0144, 0.2746>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0544, 2.0694, -0.5751>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_71 = union {
  cylinder { <0.8753, 1.5555, 0.0874>, <-0.4453, 1.3021, -0.0777>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8753, 1.5555, 0.0874>, <1.3434, 2.4426, -0.3472>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8753, 1.5555, 0.0874>, <1.5183, 0.9516, 0.7405>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4453, 1.3021, -0.0777>, <-1.1168, 0.0367, 0.1802>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4453, 1.3021, -0.0777>, <-1.0544, 2.0694, -0.5751>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1168, 0.0367, 0.1802>, <-0.4629, -1.1392, 0.0112>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1168, 0.0367, 0.1802>, <-2.2020, 0.0144, 0.2746>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4629, -1.1392, 0.0112>, <0.8563, -1.2619, -0.2241>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4629, -1.1392, 0.0112>, <-0.9890, -2.0997, -0.0246>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8563, -1.2619, -0.2241>, <1.2053, -0.3600, -0.3484>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[71] = union {
    object {AtomsGroup_71}
    object {BondsGroup_71}
}
// End of Section 71
        //Begin of Section 72  
#declare AtomsGroup_72 = union {
  sphere { <0.8766, 1.5641, 0.0890>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4407, 1.3008, -0.0887>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1153, 0.0358, 0.1736>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4606, -1.1407, 0.0102>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8528, -1.2726, -0.2297>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2099, -0.3780, -0.3524>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3345, 2.4623, -0.3265>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5193, 0.9523, 0.7332>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9940, -2.0974, -0.0112>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2014, 0.0160, 0.2791>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0528, 2.0686, -0.5796>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_72 = union {
  cylinder { <0.8766, 1.5641, 0.0890>, <-0.4407, 1.3008, -0.0887>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8766, 1.5641, 0.0890>, <1.3345, 2.4623, -0.3265>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8766, 1.5641, 0.0890>, <1.5193, 0.9523, 0.7332>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4407, 1.3008, -0.0887>, <-1.1153, 0.0358, 0.1736>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4407, 1.3008, -0.0887>, <-1.0528, 2.0686, -0.5796>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1153, 0.0358, 0.1736>, <-0.4606, -1.1407, 0.0102>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1153, 0.0358, 0.1736>, <-2.2014, 0.0160, 0.2791>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4606, -1.1407, 0.0102>, <0.8528, -1.2726, -0.2297>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4606, -1.1407, 0.0102>, <-0.9940, -2.0974, -0.0112>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8528, -1.2726, -0.2297>, <1.2099, -0.3780, -0.3524>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[72] = union {
    object {AtomsGroup_72}
    object {BondsGroup_72}
}
// End of Section 72
        //Begin of Section 73  
#declare AtomsGroup_73 = union {
  sphere { <0.8767, 1.5753, 0.0925>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4344, 1.3004, -0.0999>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1118, 0.0360, 0.1646>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4633, -1.1430, 0.0106>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8473, -1.2946, -0.2354>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2177, -0.3905, -0.3588>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3248, 2.4845, -0.3038>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5200, 0.9530, 0.7239>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.9974, -2.0950, 0.0040>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.2003, 0.0173, 0.2838>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0508, 2.0678, -0.5848>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_73 = union {
  cylinder { <0.8767, 1.5753, 0.0925>, <-0.4344, 1.3004, -0.0999>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8767, 1.5753, 0.0925>, <1.3248, 2.4845, -0.3038>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8767, 1.5753, 0.0925>, <1.5200, 0.9530, 0.7239>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4344, 1.3004, -0.0999>, <-1.1118, 0.0360, 0.1646>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4344, 1.3004, -0.0999>, <-1.0508, 2.0678, -0.5848>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1118, 0.0360, 0.1646>, <-0.4633, -1.1430, 0.0106>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1118, 0.0360, 0.1646>, <-2.2003, 0.0173, 0.2838>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4633, -1.1430, 0.0106>, <0.8473, -1.2946, -0.2354>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4633, -1.1430, 0.0106>, <-0.9974, -2.0950, 0.0040>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8473, -1.2946, -0.2354>, <1.2177, -0.3905, -0.3588>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[73] = union {
    object {AtomsGroup_73}
    object {BondsGroup_73}
}
// End of Section 73
        //Begin of Section 74  
#declare AtomsGroup_74 = union {
  sphere { <0.8745, 1.5875, 0.0984>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4289, 1.3014, -0.1098>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1073, 0.0386, 0.1528>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4715, -1.1475, 0.0125>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8453, -1.3162, -0.2424>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2250, -0.4091, -0.3671>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3152, 2.5111, -0.2789>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5211, 0.9537, 0.7129>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0007, -2.0949, 0.0218>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1957, 0.0188, 0.2881>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0484, 2.0674, -0.5914>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_74 = union {
  cylinder { <0.8745, 1.5875, 0.0984>, <-0.4289, 1.3014, -0.1098>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8745, 1.5875, 0.0984>, <1.3152, 2.5111, -0.2789>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8745, 1.5875, 0.0984>, <1.5211, 0.9537, 0.7129>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4289, 1.3014, -0.1098>, <-1.1073, 0.0386, 0.1528>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4289, 1.3014, -0.1098>, <-1.0484, 2.0674, -0.5914>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1073, 0.0386, 0.1528>, <-0.4715, -1.1475, 0.0125>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1073, 0.0386, 0.1528>, <-2.1957, 0.0188, 0.2881>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4715, -1.1475, 0.0125>, <0.8453, -1.3162, -0.2424>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4715, -1.1475, 0.0125>, <-1.0007, -2.0949, 0.0218>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8453, -1.3162, -0.2424>, <1.2250, -0.4091, -0.3671>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[74] = union {
    object {AtomsGroup_74}
    object {BondsGroup_74}
}
// End of Section 74
        //Begin of Section 75  
#declare AtomsGroup_75 = union {
  sphere { <0.8724, 1.6017, 0.1052>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4268, 1.3051, -0.1184>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.1040, 0.0403, 0.1395>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4735, -1.1509, 0.0135>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8428, -1.3279, -0.2499>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2286, -0.4410, -0.3766>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.3051, 2.5399, -0.2516>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5233, 0.9547, 0.7011>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0074, -2.1009, 0.0425>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1863, 0.0215, 0.2912>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0457, 2.0681, -0.5996>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_75 = union {
  cylinder { <0.8724, 1.6017, 0.1052>, <-0.4268, 1.3051, -0.1184>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8724, 1.6017, 0.1052>, <1.3051, 2.5399, -0.2516>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8724, 1.6017, 0.1052>, <1.5233, 0.9547, 0.7011>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4268, 1.3051, -0.1184>, <-1.1040, 0.0403, 0.1395>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4268, 1.3051, -0.1184>, <-1.0457, 2.0681, -0.5996>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1040, 0.0403, 0.1395>, <-0.4735, -1.1509, 0.0135>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.1040, 0.0403, 0.1395>, <-2.1863, 0.0215, 0.2912>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4735, -1.1509, 0.0135>, <0.8428, -1.3279, -0.2499>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4735, -1.1509, 0.0135>, <-1.0074, -2.1009, 0.0425>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8428, -1.3279, -0.2499>, <1.2286, -0.4410, -0.3766>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[75] = union {
    object {AtomsGroup_75}
    object {BondsGroup_75}
}
// End of Section 75
        //Begin of Section 76  
#declare AtomsGroup_76 = union {
  sphere { <0.8685, 1.6183, 0.1111>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4252, 1.3121, -0.1240>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0972, 0.0366, 0.1250>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4693, -1.1541, 0.0123>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8316, -1.3448, -0.2548>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2343, -0.4725, -0.3902>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2939, 2.5689, -0.2206>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5272, 0.9568, 0.6889>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0156, -2.1101, 0.0670>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1751, 0.0265, 0.2925>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0443, 2.0726, -0.6104>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_76 = union {
  cylinder { <0.8685, 1.6183, 0.1111>, <-0.4252, 1.3121, -0.1240>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8685, 1.6183, 0.1111>, <1.2939, 2.5689, -0.2206>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8685, 1.6183, 0.1111>, <1.5272, 0.9568, 0.6889>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4252, 1.3121, -0.1240>, <-1.0972, 0.0366, 0.1250>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4252, 1.3121, -0.1240>, <-1.0443, 2.0726, -0.6104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0972, 0.0366, 0.1250>, <-0.4693, -1.1541, 0.0123>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0972, 0.0366, 0.1250>, <-2.1751, 0.0265, 0.2925>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4693, -1.1541, 0.0123>, <0.8316, -1.3448, -0.2548>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4693, -1.1541, 0.0123>, <-1.0156, -2.1101, 0.0670>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8316, -1.3448, -0.2548>, <1.2343, -0.4725, -0.3902>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[76] = union {
    object {AtomsGroup_76}
    object {BondsGroup_76}
}
// End of Section 76
        //Begin of Section 77  
#declare AtomsGroup_77 = union {
  sphere { <0.8655, 1.6333, 0.1182>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4251, 1.3181, -0.1286>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0858, 0.0371, 0.1105>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4805, -1.1696, 0.0123>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8266, -1.3745, -0.2594>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2443, -0.4970, -0.4086>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2817, 2.5975, -0.1867>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9630, 0.6753>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0190, -2.1108, 0.0940>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1646, 0.0323, 0.2919>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0448, 2.0806, -0.6220>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_77 = union {
  cylinder { <0.8655, 1.6333, 0.1182>, <-0.4251, 1.3181, -0.1286>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8655, 1.6333, 0.1182>, <1.2817, 2.5975, -0.1867>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8655, 1.6333, 0.1182>, <1.5304, 0.9630, 0.6753>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4251, 1.3181, -0.1286>, <-1.0858, 0.0371, 0.1105>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4251, 1.3181, -0.1286>, <-1.0448, 2.0806, -0.6220>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0858, 0.0371, 0.1105>, <-0.4805, -1.1696, 0.0123>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0858, 0.0371, 0.1105>, <-2.1646, 0.0323, 0.2919>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4805, -1.1696, 0.0123>, <0.8266, -1.3745, -0.2594>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4805, -1.1696, 0.0123>, <-1.0190, -2.1108, 0.0940>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8266, -1.3745, -0.2594>, <1.2443, -0.4970, -0.4086>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[77] = union {
    object {AtomsGroup_77}
    object {BondsGroup_77}
}
// End of Section 77
        //Begin of Section 78  
#declare AtomsGroup_78 = union {
  sphere { <0.8651, 1.6330, 0.1181>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4244, 1.3185, -0.1292>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0868, 0.0362, 0.1110>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4791, -1.1672, 0.0120>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8280, -1.3735, -0.2596>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2431, -0.4991, -0.4088>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2815, 2.5973, -0.1863>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9636, 0.6754>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0199, -2.1116, 0.0944>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1643, 0.0325, 0.2913>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0448, 2.0805, -0.6214>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_78 = union {
  cylinder { <0.8651, 1.6330, 0.1181>, <-0.4244, 1.3185, -0.1292>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8651, 1.6330, 0.1181>, <1.2815, 2.5973, -0.1863>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8651, 1.6330, 0.1181>, <1.5304, 0.9636, 0.6754>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4244, 1.3185, -0.1292>, <-1.0868, 0.0362, 0.1110>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4244, 1.3185, -0.1292>, <-1.0448, 2.0805, -0.6214>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0868, 0.0362, 0.1110>, <-0.4791, -1.1672, 0.0120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0868, 0.0362, 0.1110>, <-2.1643, 0.0325, 0.2913>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4791, -1.1672, 0.0120>, <0.8280, -1.3735, -0.2596>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4791, -1.1672, 0.0120>, <-1.0199, -2.1116, 0.0944>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8280, -1.3735, -0.2596>, <1.2431, -0.4991, -0.4088>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[78] = union {
    object {AtomsGroup_78}
    object {BondsGroup_78}
}
// End of Section 78
        //Begin of Section 79  
#declare AtomsGroup_79 = union {
  sphere { <0.8646, 1.6326, 0.1179>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4236, 1.3188, -0.1304>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0883, 0.0353, 0.1118>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4767, -1.1641, 0.0115>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8301, -1.3724, -0.2597>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2412, -0.5025, -0.4096>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2811, 2.5968, -0.1856>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9647, 0.6756>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0211, -2.1125, 0.0954>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1641, 0.0329, 0.2903>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0448, 2.0804, -0.6202>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_79 = union {
  cylinder { <0.8646, 1.6326, 0.1179>, <-0.4236, 1.3188, -0.1304>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8646, 1.6326, 0.1179>, <1.2811, 2.5968, -0.1856>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8646, 1.6326, 0.1179>, <1.5304, 0.9647, 0.6756>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4236, 1.3188, -0.1304>, <-1.0883, 0.0353, 0.1118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4236, 1.3188, -0.1304>, <-1.0448, 2.0804, -0.6202>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0883, 0.0353, 0.1118>, <-0.4767, -1.1641, 0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0883, 0.0353, 0.1118>, <-2.1641, 0.0329, 0.2903>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1641, 0.0115>, <0.8301, -1.3724, -0.2597>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1641, 0.0115>, <-1.0211, -2.1125, 0.0954>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8301, -1.3724, -0.2597>, <1.2412, -0.5025, -0.4096>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[79] = union {
    object {AtomsGroup_79}
    object {BondsGroup_79}
}
// End of Section 79
        //Begin of Section 80  
#declare AtomsGroup_80 = union {
  sphere { <0.8646, 1.6326, 0.1179>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4237, 1.3187, -0.1305>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0882, 0.0356, 0.1118>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4767, -1.1643, 0.0115>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8301, -1.3726, -0.2597>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2412, -0.5024, -0.4098>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2810, 2.5967, -0.1855>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9648, 0.6756>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0211, -2.1123, 0.0955>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1642, 0.0329, 0.2902>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0448, 2.0804, -0.6201>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_80 = union {
  cylinder { <0.8646, 1.6326, 0.1179>, <-0.4237, 1.3187, -0.1305>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8646, 1.6326, 0.1179>, <1.2810, 2.5967, -0.1855>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8646, 1.6326, 0.1179>, <1.5304, 0.9648, 0.6756>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3187, -0.1305>, <-1.0882, 0.0356, 0.1118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3187, -0.1305>, <-1.0448, 2.0804, -0.6201>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0882, 0.0356, 0.1118>, <-0.4767, -1.1643, 0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0882, 0.0356, 0.1118>, <-2.1642, 0.0329, 0.2902>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1643, 0.0115>, <0.8301, -1.3726, -0.2597>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1643, 0.0115>, <-1.0211, -2.1123, 0.0955>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8301, -1.3726, -0.2597>, <1.2412, -0.5024, -0.4098>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[80] = union {
    object {AtomsGroup_80}
    object {BondsGroup_80}
}
// End of Section 80
        //Begin of Section 81  
#declare AtomsGroup_81 = union {
  sphere { <0.8647, 1.6326, 0.1179>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4238, 1.3185, -0.1307>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0880, 0.0360, 0.1118>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4767, -1.1648, 0.0115>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8301, -1.3731, -0.2596>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2413, -0.5024, -0.4101>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2810, 2.5966, -0.1853>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9650, 0.6757>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0210, -2.1119, 0.0956>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1644, 0.0330, 0.2900>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0449, 2.0805, -0.6199>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_81 = union {
  cylinder { <0.8647, 1.6326, 0.1179>, <-0.4238, 1.3185, -0.1307>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1179>, <1.2810, 2.5966, -0.1853>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1179>, <1.5304, 0.9650, 0.6757>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4238, 1.3185, -0.1307>, <-1.0880, 0.0360, 0.1118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4238, 1.3185, -0.1307>, <-1.0449, 2.0805, -0.6199>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0880, 0.0360, 0.1118>, <-0.4767, -1.1648, 0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0880, 0.0360, 0.1118>, <-2.1644, 0.0330, 0.2900>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1648, 0.0115>, <0.8301, -1.3731, -0.2596>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4767, -1.1648, 0.0115>, <-1.0210, -2.1119, 0.0956>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8301, -1.3731, -0.2596>, <1.2413, -0.5024, -0.4101>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[81] = union {
    object {AtomsGroup_81}
    object {BondsGroup_81}
}
// End of Section 81
        //Begin of Section 82  
#declare AtomsGroup_82 = union {
  sphere { <0.8647, 1.6326, 0.1178>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4240, 1.3183, -0.1309>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0878, 0.0365, 0.1119>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4766, -1.1653, 0.0115>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8301, -1.3737, -0.2594>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2414, -0.5023, -0.4106>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2808, 2.5965, -0.1851>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9654, 0.6758>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0209, -2.1114, 0.0959>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1646, 0.0330, 0.2897>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0450, 2.0805, -0.6196>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_82 = union {
  cylinder { <0.8647, 1.6326, 0.1178>, <-0.4240, 1.3183, -0.1309>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1178>, <1.2808, 2.5965, -0.1851>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1178>, <1.5304, 0.9654, 0.6758>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4240, 1.3183, -0.1309>, <-1.0878, 0.0365, 0.1119>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4240, 1.3183, -0.1309>, <-1.0450, 2.0805, -0.6196>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0878, 0.0365, 0.1119>, <-0.4766, -1.1653, 0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0878, 0.0365, 0.1119>, <-2.1646, 0.0330, 0.2897>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4766, -1.1653, 0.0115>, <0.8301, -1.3737, -0.2594>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4766, -1.1653, 0.0115>, <-1.0209, -2.1114, 0.0959>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8301, -1.3737, -0.2594>, <1.2414, -0.5023, -0.4106>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[82] = union {
    object {AtomsGroup_82}
    object {BondsGroup_82}
}
// End of Section 82
        //Begin of Section 83  
#declare AtomsGroup_83 = union {
  sphere { <0.8647, 1.6326, 0.1178>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4242, 1.3180, -0.1312>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0876, 0.0369, 0.1119>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4765, -1.1658, 0.0115>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8302, -1.3744, -0.2593>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2415, -0.5023, -0.4113>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2807, 2.5964, -0.1848>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5304, 0.9658, 0.6759>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0207, -2.1108, 0.0962>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1649, 0.0331, 0.2893>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0451, 2.0807, -0.6193>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_83 = union {
  cylinder { <0.8647, 1.6326, 0.1178>, <-0.4242, 1.3180, -0.1312>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1178>, <1.2807, 2.5964, -0.1848>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6326, 0.1178>, <1.5304, 0.9658, 0.6759>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3180, -0.1312>, <-1.0876, 0.0369, 0.1119>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3180, -0.1312>, <-1.0451, 2.0807, -0.6193>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0876, 0.0369, 0.1119>, <-0.4765, -1.1658, 0.0115>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0876, 0.0369, 0.1119>, <-2.1649, 0.0331, 0.2893>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4765, -1.1658, 0.0115>, <0.8302, -1.3744, -0.2593>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4765, -1.1658, 0.0115>, <-1.0207, -2.1108, 0.0962>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8302, -1.3744, -0.2593>, <1.2415, -0.5023, -0.4113>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[83] = union {
    object {AtomsGroup_83}
    object {BondsGroup_83}
}
// End of Section 83
        //Begin of Section 84  
#declare AtomsGroup_84 = union {
  sphere { <0.8647, 1.6325, 0.1177>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4242, 1.3177, -0.1315>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0873, 0.0372, 0.1120>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4764, -1.1661, 0.0116>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8304, -1.3752, -0.2591>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2415, -0.5025, -0.4121>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2805, 2.5963, -0.1844>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5305, 0.9663, 0.6761>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0207, -2.1102, 0.0967>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1652, 0.0333, 0.2888>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0452, 2.0808, -0.6188>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_84 = union {
  cylinder { <0.8647, 1.6325, 0.1177>, <-0.4242, 1.3177, -0.1315>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6325, 0.1177>, <1.2805, 2.5963, -0.1844>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8647, 1.6325, 0.1177>, <1.5305, 0.9663, 0.6761>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3177, -0.1315>, <-1.0873, 0.0372, 0.1120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3177, -0.1315>, <-1.0452, 2.0808, -0.6188>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0873, 0.0372, 0.1120>, <-0.4764, -1.1661, 0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0873, 0.0372, 0.1120>, <-2.1652, 0.0333, 0.2888>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4764, -1.1661, 0.0116>, <0.8304, -1.3752, -0.2591>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4764, -1.1661, 0.0116>, <-1.0207, -2.1102, 0.0967>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8304, -1.3752, -0.2591>, <1.2415, -0.5025, -0.4121>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[84] = union {
    object {AtomsGroup_84}
    object {BondsGroup_84}
}
// End of Section 84
        //Begin of Section 85  
#declare AtomsGroup_85 = union {
  sphere { <0.8645, 1.6324, 0.1176>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4242, 1.3175, -0.1319>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0870, 0.0373, 0.1120>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4762, -1.1661, 0.0116>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8306, -1.3760, -0.2589>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2415, -0.5029, -0.4130>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2802, 2.5962, -0.1840>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5306, 0.9668, 0.6764>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0206, -2.1095, 0.0972>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1654, 0.0334, 0.2882>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0454, 2.0810, -0.6182>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_85 = union {
  cylinder { <0.8645, 1.6324, 0.1176>, <-0.4242, 1.3175, -0.1319>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8645, 1.6324, 0.1176>, <1.2802, 2.5962, -0.1840>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8645, 1.6324, 0.1176>, <1.5306, 0.9668, 0.6764>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3175, -0.1319>, <-1.0870, 0.0373, 0.1120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4242, 1.3175, -0.1319>, <-1.0454, 2.0810, -0.6182>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0870, 0.0373, 0.1120>, <-0.4762, -1.1661, 0.0116>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0870, 0.0373, 0.1120>, <-2.1654, 0.0334, 0.2882>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4762, -1.1661, 0.0116>, <0.8306, -1.3760, -0.2589>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4762, -1.1661, 0.0116>, <-1.0206, -2.1095, 0.0972>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8306, -1.3760, -0.2589>, <1.2415, -0.5029, -0.4130>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[85] = union {
    object {AtomsGroup_85}
    object {BondsGroup_85}
}
// End of Section 85
        //Begin of Section 86  
#declare AtomsGroup_86 = union {
  sphere { <0.8643, 1.6322, 0.1174>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4241, 1.3174, -0.1323>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0867, 0.0372, 0.1121>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4760, -1.1659, 0.0117>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8309, -1.3766, -0.2587>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2414, -0.5035, -0.4141>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2800, 2.5961, -0.1835>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5307, 0.9674, 0.6766>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0207, -2.1090, 0.0978>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1656, 0.0336, 0.2875>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0456, 2.0812, -0.6176>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_86 = union {
  cylinder { <0.8643, 1.6322, 0.1174>, <-0.4241, 1.3174, -0.1323>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8643, 1.6322, 0.1174>, <1.2800, 2.5961, -0.1835>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8643, 1.6322, 0.1174>, <1.5307, 0.9674, 0.6766>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4241, 1.3174, -0.1323>, <-1.0867, 0.0372, 0.1121>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4241, 1.3174, -0.1323>, <-1.0456, 2.0812, -0.6176>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0867, 0.0372, 0.1121>, <-0.4760, -1.1659, 0.0117>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0867, 0.0372, 0.1121>, <-2.1656, 0.0336, 0.2875>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4760, -1.1659, 0.0117>, <0.8309, -1.3766, -0.2587>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4760, -1.1659, 0.0117>, <-1.0207, -2.1090, 0.0978>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8309, -1.3766, -0.2587>, <1.2414, -0.5035, -0.4141>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[86] = union {
    object {AtomsGroup_86}
    object {BondsGroup_86}
}
// End of Section 86
        //Begin of Section 87  
#declare AtomsGroup_87 = union {
  sphere { <0.8641, 1.6320, 0.1172>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4239, 1.3172, -0.1328>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0866, 0.0369, 0.1121>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4757, -1.1654, 0.0118>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8313, -1.3773, -0.2585>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2413, -0.5044, -0.4153>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2797, 2.5961, -0.1830>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5308, 0.9681, 0.6770>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0209, -2.1084, 0.0986>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1658, 0.0338, 0.2865>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0458, 2.0814, -0.6168>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_87 = union {
  cylinder { <0.8641, 1.6320, 0.1172>, <-0.4239, 1.3172, -0.1328>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8641, 1.6320, 0.1172>, <1.2797, 2.5961, -0.1830>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8641, 1.6320, 0.1172>, <1.5308, 0.9681, 0.6770>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3172, -0.1328>, <-1.0866, 0.0369, 0.1121>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3172, -0.1328>, <-1.0458, 2.0814, -0.6168>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0369, 0.1121>, <-0.4757, -1.1654, 0.0118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0369, 0.1121>, <-2.1658, 0.0338, 0.2865>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4757, -1.1654, 0.0118>, <0.8313, -1.3773, -0.2585>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4757, -1.1654, 0.0118>, <-1.0209, -2.1084, 0.0986>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8313, -1.3773, -0.2585>, <1.2413, -0.5044, -0.4153>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[87] = union {
    object {AtomsGroup_87}
    object {BondsGroup_87}
}
// End of Section 87
        //Begin of Section 88  
#declare AtomsGroup_88 = union {
  sphere { <0.8638, 1.6318, 0.1170>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4237, 1.3172, -0.1333>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0865, 0.0364, 0.1121>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4752, -1.1648, 0.0120>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8318, -1.3778, -0.2582>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2410, -0.5057, -0.4169>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2793, 2.5961, -0.1823>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5310, 0.9690, 0.6774>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0211, -2.1079, 0.0996>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1658, 0.0341, 0.2854>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0461, 2.0817, -0.6159>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_88 = union {
  cylinder { <0.8638, 1.6318, 0.1170>, <-0.4237, 1.3172, -0.1333>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8638, 1.6318, 0.1170>, <1.2793, 2.5961, -0.1823>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8638, 1.6318, 0.1170>, <1.5310, 0.9690, 0.6774>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3172, -0.1333>, <-1.0865, 0.0364, 0.1121>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3172, -0.1333>, <-1.0461, 2.0817, -0.6159>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0865, 0.0364, 0.1121>, <-0.4752, -1.1648, 0.0120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0865, 0.0364, 0.1121>, <-2.1658, 0.0341, 0.2854>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4752, -1.1648, 0.0120>, <0.8318, -1.3778, -0.2582>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4752, -1.1648, 0.0120>, <-1.0211, -2.1079, 0.0996>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8318, -1.3778, -0.2582>, <1.2410, -0.5057, -0.4169>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[88] = union {
    object {AtomsGroup_88}
    object {BondsGroup_88}
}
// End of Section 88
        //Begin of Section 89  
#declare AtomsGroup_89 = union {
  sphere { <0.8637, 1.6316, 0.1168>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4236, 1.3170, -0.1338>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0866, 0.0361, 0.1120>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4746, -1.1644, 0.0123>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8323, -1.3783, -0.2579>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2408, -0.5073, -0.4188>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2788, 2.5963, -0.1815>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5312, 0.9699, 0.6780>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0215, -2.1073, 0.1009>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1656, 0.0344, 0.2839>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0465, 2.0821, -0.6148>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_89 = union {
  cylinder { <0.8637, 1.6316, 0.1168>, <-0.4236, 1.3170, -0.1338>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8637, 1.6316, 0.1168>, <1.2788, 2.5963, -0.1815>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8637, 1.6316, 0.1168>, <1.5312, 0.9699, 0.6780>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4236, 1.3170, -0.1338>, <-1.0866, 0.0361, 0.1120>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4236, 1.3170, -0.1338>, <-1.0465, 2.0821, -0.6148>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0361, 0.1120>, <-0.4746, -1.1644, 0.0123>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0361, 0.1120>, <-2.1656, 0.0344, 0.2839>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4746, -1.1644, 0.0123>, <0.8323, -1.3783, -0.2579>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4746, -1.1644, 0.0123>, <-1.0215, -2.1073, 0.1009>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8323, -1.3783, -0.2579>, <1.2408, -0.5073, -0.4188>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[89] = union {
    object {AtomsGroup_89}
    object {BondsGroup_89}
}
// End of Section 89
        //Begin of Section 90  
#declare AtomsGroup_90 = union {
  sphere { <0.8637, 1.6315, 0.1166>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4237, 1.3168, -0.1343>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0866, 0.0361, 0.1118>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4738, -1.1643, 0.0126>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8326, -1.3790, -0.2575>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2408, -0.5091, -0.4213>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2783, 2.5965, -0.1806>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5314, 0.9709, 0.6787>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0219, -2.1065, 0.1024>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1654, 0.0347, 0.2821>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0470, 2.0825, -0.6136>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_90 = union {
  cylinder { <0.8637, 1.6315, 0.1166>, <-0.4237, 1.3168, -0.1343>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8637, 1.6315, 0.1166>, <1.2783, 2.5965, -0.1806>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8637, 1.6315, 0.1166>, <1.5314, 0.9709, 0.6787>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3168, -0.1343>, <-1.0866, 0.0361, 0.1118>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4237, 1.3168, -0.1343>, <-1.0470, 2.0825, -0.6136>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0361, 0.1118>, <-0.4738, -1.1643, 0.0126>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0866, 0.0361, 0.1118>, <-2.1654, 0.0347, 0.2821>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4738, -1.1643, 0.0126>, <0.8326, -1.3790, -0.2575>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4738, -1.1643, 0.0126>, <-1.0219, -2.1065, 0.1024>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8326, -1.3790, -0.2575>, <1.2408, -0.5091, -0.4213>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[90] = union {
    object {AtomsGroup_90}
    object {BondsGroup_90}
}
// End of Section 90
        //Begin of Section 91  
#declare AtomsGroup_91 = union {
  sphere { <0.8639, 1.6315, 0.1165>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4239, 1.3166, -0.1347>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0865, 0.0364, 0.1113>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4732, -1.1647, 0.0132>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8327, -1.3799, -0.2569>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2410, -0.5110, -0.4243>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2776, 2.5970, -0.1796>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5317, 0.9719, 0.6794>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0222, -2.1057, 0.1043>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1651, 0.0349, 0.2800>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0476, 2.0830, -0.6122>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_91 = union {
  cylinder { <0.8639, 1.6315, 0.1165>, <-0.4239, 1.3166, -0.1347>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8639, 1.6315, 0.1165>, <1.2776, 2.5970, -0.1796>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8639, 1.6315, 0.1165>, <1.5317, 0.9719, 0.6794>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3166, -0.1347>, <-1.0865, 0.0364, 0.1113>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3166, -0.1347>, <-1.0476, 2.0830, -0.6122>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0865, 0.0364, 0.1113>, <-0.4732, -1.1647, 0.0132>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0865, 0.0364, 0.1113>, <-2.1651, 0.0349, 0.2800>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1647, 0.0132>, <0.8327, -1.3799, -0.2569>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1647, 0.0132>, <-1.0222, -2.1057, 0.1043>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8327, -1.3799, -0.2569>, <1.2410, -0.5110, -0.4243>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[91] = union {
    object {AtomsGroup_91}
    object {BondsGroup_91}
}
// End of Section 91
        //Begin of Section 92  
#declare AtomsGroup_92 = union {
  sphere { <0.8640, 1.6315, 0.1164>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4240, 1.3163, -0.1350>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0858, 0.0368, 0.1104>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4730, -1.1653, 0.0141>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8326, -1.3813, -0.2562>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2416, -0.5128, -0.4281>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2768, 2.5978, -0.1784>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5320, 0.9729, 0.6803>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0226, -2.1048, 0.1066>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1649, 0.0351, 0.2776>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0483, 2.0836, -0.6107>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_92 = union {
  cylinder { <0.8640, 1.6315, 0.1164>, <-0.4240, 1.3163, -0.1350>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6315, 0.1164>, <1.2768, 2.5978, -0.1784>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6315, 0.1164>, <1.5320, 0.9729, 0.6803>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4240, 1.3163, -0.1350>, <-1.0858, 0.0368, 0.1104>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4240, 1.3163, -0.1350>, <-1.0483, 2.0836, -0.6107>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0858, 0.0368, 0.1104>, <-0.4730, -1.1653, 0.0141>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0858, 0.0368, 0.1104>, <-2.1649, 0.0351, 0.2776>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4730, -1.1653, 0.0141>, <0.8326, -1.3813, -0.2562>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4730, -1.1653, 0.0141>, <-1.0226, -2.1048, 0.1066>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8326, -1.3813, -0.2562>, <1.2416, -0.5128, -0.4281>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[92] = union {
    object {AtomsGroup_92}
    object {BondsGroup_92}
}
// End of Section 92
        //Begin of Section 93  
#declare AtomsGroup_93 = union {
  sphere { <0.8640, 1.6317, 0.1164>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4239, 1.3165, -0.1350>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0848, 0.0368, 0.1091>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4732, -1.1654, 0.0153>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8324, -1.3831, -0.2553>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2425, -0.5148, -0.4325>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2759, 2.5990, -0.1771>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5323, 0.9739, 0.6812>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0231, -2.1042, 0.1094>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1648, 0.0352, 0.2747>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0491, 2.0843, -0.6090>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_93 = union {
  cylinder { <0.8640, 1.6317, 0.1164>, <-0.4239, 1.3165, -0.1350>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6317, 0.1164>, <1.2759, 2.5990, -0.1771>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6317, 0.1164>, <1.5323, 0.9739, 0.6812>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3165, -0.1350>, <-1.0848, 0.0368, 0.1091>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4239, 1.3165, -0.1350>, <-1.0491, 2.0843, -0.6090>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0848, 0.0368, 0.1091>, <-0.4732, -1.1654, 0.0153>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0848, 0.0368, 0.1091>, <-2.1648, 0.0352, 0.2747>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1654, 0.0153>, <0.8324, -1.3831, -0.2553>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1654, 0.0153>, <-1.0231, -2.1042, 0.1094>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8324, -1.3831, -0.2553>, <1.2425, -0.5148, -0.4325>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[93] = union {
    object {AtomsGroup_93}
    object {BondsGroup_93}
}
// End of Section 93
        //Begin of Section 94  
#declare AtomsGroup_94 = union {
  sphere { <0.8640, 1.6322, 0.1165>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4234, 1.3171, -0.1349>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0839, 0.0365, 0.1075>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4735, -1.1652, 0.0166>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8324, -1.3847, -0.2543>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2433, -0.5176, -0.4378>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2748, 2.6006, -0.1756>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5326, 0.9747, 0.6821>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0239, -2.1039, 0.1128>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1645, 0.0353, 0.2713>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0499, 2.0850, -0.6071>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_94 = union {
  cylinder { <0.8640, 1.6322, 0.1165>, <-0.4234, 1.3171, -0.1349>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6322, 0.1165>, <1.2748, 2.6006, -0.1756>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8640, 1.6322, 0.1165>, <1.5326, 0.9747, 0.6821>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4234, 1.3171, -0.1349>, <-1.0839, 0.0365, 0.1075>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4234, 1.3171, -0.1349>, <-1.0499, 2.0850, -0.6071>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0839, 0.0365, 0.1075>, <-0.4735, -1.1652, 0.0166>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0839, 0.0365, 0.1075>, <-2.1645, 0.0353, 0.2713>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4735, -1.1652, 0.0166>, <0.8324, -1.3847, -0.2543>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4735, -1.1652, 0.0166>, <-1.0239, -2.1039, 0.1128>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8324, -1.3847, -0.2543>, <1.2433, -0.5176, -0.4378>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[94] = union {
    object {AtomsGroup_94}
    object {BondsGroup_94}
}
// End of Section 94
        //Begin of Section 95  
#declare AtomsGroup_95 = union {
  sphere { <0.8644, 1.6331, 0.1171>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4229, 1.3179, -0.1348>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0837, 0.0367, 0.1058>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4732, -1.1654, 0.0180>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8325, -1.3860, -0.2533>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2439, -0.5216, -0.4437>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2734, 2.6026, -0.1738>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5328, 0.9755, 0.6828>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0248, -2.1040, 0.1169>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1638, 0.0353, 0.2673>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0507, 2.0858, -0.6051>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_95 = union {
  cylinder { <0.8644, 1.6331, 0.1171>, <-0.4229, 1.3179, -0.1348>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8644, 1.6331, 0.1171>, <1.2734, 2.6026, -0.1738>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8644, 1.6331, 0.1171>, <1.5328, 0.9755, 0.6828>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4229, 1.3179, -0.1348>, <-1.0837, 0.0367, 0.1058>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4229, 1.3179, -0.1348>, <-1.0507, 2.0858, -0.6051>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0837, 0.0367, 0.1058>, <-0.4732, -1.1654, 0.0180>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0837, 0.0367, 0.1058>, <-2.1638, 0.0353, 0.2673>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1654, 0.0180>, <0.8325, -1.3860, -0.2533>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4732, -1.1654, 0.0180>, <-1.0248, -2.1040, 0.1169>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8325, -1.3860, -0.2533>, <1.2439, -0.5216, -0.4437>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[95] = union {
    object {AtomsGroup_95}
    object {BondsGroup_95}
}
// End of Section 95
        //Begin of Section 96  
#declare AtomsGroup_96 = union {
  sphere { <0.8651, 1.6348, 0.1183>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4226, 1.3187, -0.1347>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0837, 0.0375, 0.1038>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4726, -1.1662, 0.0192>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8325, -1.3876, -0.2519>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2443, -0.5264, -0.4505>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2719, 2.6051, -0.1715>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5328, 0.9762, 0.6829>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0257, -2.1042, 0.1218>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1631, 0.0352, 0.2627>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0514, 2.0866, -0.6029>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_96 = union {
  cylinder { <0.8651, 1.6348, 0.1183>, <-0.4226, 1.3187, -0.1347>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8651, 1.6348, 0.1183>, <1.2719, 2.6051, -0.1715>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8651, 1.6348, 0.1183>, <1.5328, 0.9762, 0.6829>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4226, 1.3187, -0.1347>, <-1.0837, 0.0375, 0.1038>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4226, 1.3187, -0.1347>, <-1.0514, 2.0866, -0.6029>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0837, 0.0375, 0.1038>, <-0.4726, -1.1662, 0.0192>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0837, 0.0375, 0.1038>, <-2.1631, 0.0352, 0.2627>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4726, -1.1662, 0.0192>, <0.8325, -1.3876, -0.2519>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4726, -1.1662, 0.0192>, <-1.0257, -2.1042, 0.1218>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8325, -1.3876, -0.2519>, <1.2443, -0.5264, -0.4505>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[96] = union {
    object {AtomsGroup_96}
    object {BondsGroup_96}
}
// End of Section 96
        //Begin of Section 97  
#declare AtomsGroup_97 = union {
  sphere { <0.8650, 1.6369, 0.1199>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4218, 1.3194, -0.1346>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0823, 0.0374, 0.1015>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-0.4727, -1.1662, 0.0204>, atom_rad_2
    pigment { color rgb <.3,.3,.3> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <0.8324, -1.3907, -0.2501>, atom_rad_2
    pigment { color rgb <1,0,0> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2447, -0.5310, -0.4582>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.2704, 2.6084, -0.1687>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <1.5330, 0.9767, 0.6826>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0261, -2.1041, 0.1275>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-2.1632, 0.0351, 0.2577>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
  sphere { <-1.0520, 2.0879, -0.6007>, atom_rad_h
    pigment { color rgb <1,1,1> filter trans_atom }
    finish { AtomFinish }
  }
}
#declare BondsGroup_97 = union {
  cylinder { <0.8650, 1.6369, 0.1199>, <-0.4218, 1.3194, -0.1346>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8650, 1.6369, 0.1199>, <1.2704, 2.6084, -0.1687>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8650, 1.6369, 0.1199>, <1.5330, 0.9767, 0.6826>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4218, 1.3194, -0.1346>, <-1.0823, 0.0374, 0.1015>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4218, 1.3194, -0.1346>, <-1.0520, 2.0879, -0.6007>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0823, 0.0374, 0.1015>, <-0.4727, -1.1662, 0.0204>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-1.0823, 0.0374, 0.1015>, <-2.1632, 0.0351, 0.2577>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4727, -1.1662, 0.0204>, <0.8324, -1.3907, -0.2501>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <-0.4727, -1.1662, 0.0204>, <-1.0261, -2.1041, 0.1275>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
  cylinder { <0.8324, -1.3907, -0.2501>, <1.2447, -0.5310, -0.4582>, bond_rad
    pigment { color rgb <0.7, 0.7, 0.7> filter trans_bd }
    finish { BdFinish }
  }
}
// Combine Atoms, Bonds and Mesh 
#declare molecule[97] = union {
    object {AtomsGroup_97}
    object {BondsGroup_97}
}
// End of Section 97
        