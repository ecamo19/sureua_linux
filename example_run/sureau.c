// SurEau_C
// Version  20-04-2023
// by Hervé Cochard, INRAE, UMR-PIAF, Clermont-Ferrand, France
// still under developement...use it at your own risk !
// This model should not be distributed without prior notice to H Cochard
// thanks for reporting bugs...

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))
  
// Librairies
#include "stdio.h"
#include <stdio.h>
#include "stdlib.h"
#include <math.h>
#include "time.h"    
#include "stddef.h" 
#include <string.h>


#define NPAR 304  //total number of parameters
void print_screen(void);
void print_transient(void);
void print_annual(void);
void print_final(void);
void Compute_ST(void);
void setup(void);
void Phenology(long double);
void Get_DATA(long double);
void Compute_Cavitation(long double);
 
char version[]=" 20-04-2023";
unsigned number=1;
long double COMPET_number=0;

long double debug=0;  // if 1 then print ALL the data in the file out
long double dq_Soil1;
long double DYNAMIC0=1,DYNAMIC=1;
long double PRINT_GRAPH=1;
int WARNING=0;
int New_Year=0;

long double para[NPAR],para0[NPAR];
char filenumber[5]="1";
int year=1,N_days=0,N_days2=0;
unsigned long long indice=1;
long double indice_double=1;
long double FLUID=1;                // to account for the temperature dependance of K
long double SURFACE_TENSION=1;      // to account for the temperature dependance of P50
long double T_g_cuti=1;             // to account for the temperature dependance of g_min
long double T_OSMOTIC=1;            // to accout for the temperature dependance of osmotic potentials
long double TLEAF=1;                // if one then compute radiative Leaf temperature
long double IRRIGATE=1, Irrigation, IRR_DOY_S, IRR_DOY_F;             // for irrigation
long double CONTINUOUS;             // if one the simulation is not reinitialised at the end of the year
long double FRACTAL;                // if one the tree is fractal
long double ACCLIMATE=0,Acc_P1,Acc_P2,Acc_P3;  // for trait acclimation
long double Simul;
long double Test[4];



int Screen_out[NPAR];
int File_out[NPAR];
long double beta=0.5;
char *Label[]=    { 	    "Days    "   , "T_air      ", "T_Leaf     ", "RH_air     ", "PAR        ", "VPD_Leaf   ", "VPD_Cut    ", "VPD_Axil   ", "VPD_Branch ", "VPD_Trunk  ",
						"VPD_Root1  ", "VPD_Root2  ", "VPD_Root3  ", "VPD_Soil   ", "PAR_tot_pot", "P_AIR      ", "E_clim_m2  ", "E_Leaf     ", "E_cuti     ", "E_Branch   ", 
						"E_Trunk    ", "E_Axil     ", "E_Root1_m2 ", "E_Root2_m2 ", "E_Root3_m2 ", "E_Soil_m2  ", "E_plant    ", "E_Plant_m2 ", "E_tot kg   ", "g_canopy   ", 
						"g_s        ", "g_cuti     ", "P_evap_Apo ", "P_Leaf_sym ", "P_Leaf_Apo ", "turgor     ", "Turgor_Leaf", "P_Axil_Symp", "P_Axil_Apo ", "P_Branch_s ", 
						"P_Branch_a ", "P_Trunk_s  ", "P_Trunk_a  ", "Turgor_tr  ", "P_Root1_sy ", "P_endo1    ", "P_Root1_ap ", "P_Root2_sy ", "P_endo2    ", "P_Root2_ap ", 
						"P_Root3_sy ", "P_endo3    ", "P_Root3_a  ", "P_soil1    ", "P_soil2    ", "P_soil3    ", "RWC1       ", "RWC2       ", "RWC3       ", "K_Soil1    ", 
						"K_Soil2    ", "K_Soil3    ", "K_Inter1   ", "K_Inter2   ", "K_Inter3   ", "K_Leaf_s   ", "K_Leaf_a   ", "K_Leaf     ", "K_Root1_s  ", "K_Root1_a  ", 
						"K_Root2_s  ", "K_Root2_a  ", "K_Root3_s  ", "K_Root3_a  ", "K_Root     ", "K_Trunk    ", "K_Plant    ", "K_tot      ", "PLC_Leaf   ", "PLC_Branch ", 
						"PLC_Axil   ", "PLC_Trunk  ", "PLC_Root1  ", "PLC_Root2  ", "PLC_Root3  ", "Q_plant    ", "DQ_soil    ", "Q_evap_a   ", "Q_Leaf_s   ", "Q_Leaf_a   ", 
						"Q_Branch_s ", "Q_Branch_a ", "Q_Axil_s   ", "Q_Axil_a   ", "Q_Trunk_s  ", "Q_Trunk_a  ", "Q_Root_s   ", "Q_Root_a   ", "Q_endo     ", "Q_soil1    ",
						"Q_soil2    ", "Q_soil3    ", "Growth_r   ", "Radius     ", "A_net[0]   ", "A_net_tot_c", "WUE        ", "R_main     ", "Reserve    ", "Leaf_Area  ", 
						"Leaf_wc    ", "Branch_wc  ", "Br_Sy_wc   ", "Br_Apo_wc  ", "RWC_shoot  ", "RWC_Axil   ", "RWC_soil   ", "RWC_min    ", "RWC_int    ", "P_soil     ", 
						"P_soil_min ", "P_soil_int ", "Irrigation ", "Sap_Flow_d ", "Sap_Flow_d2", "VPD_Air    ", "Fruit_Diam ", "ETP_Pen_t  ", "ETP_Pen    ", "Rain_tot   ", 
						"Rain_soil  ", "VPD_air_t  ", "VPD_Leaf_t ", "Rain_Leaf_t", "ETP_Leaf_t ", "T_air_an   ", "T_air_an_l ", "EvapoT     ", "EvapoT mm  ", "ETP_d, mm  ", 
						"Sap_flow   ", "K_Branch   ", "P_min_lf   ", "P_min_br   ", "K_tot2     ", "A_net_c    ", "Intercep   ", "GPP        ", "GPP_day    ", "ETR_d, mm  ", 
						"E_day, mm  ", "T_Soil     ", "PAR_pot    ", "Cloud      ", "F_ls_le    ", "F_la_le    ", "F_ba_la    ", "F_ba_bs    ", "F_ta_ba    ", "F_ta_ts    ", 
						"F_bua_bus  ", "F_ba_bua   ", "F_ra_ta    ", "F_re_ra    ", "F_rs_re    ", "F_s_re     ", "F_s_12     ", "F_s_23     ", "A_net_V    ", "A_net_J    ", 
						"g_Axil     ", "Root-Area  ", "Axil_WC    ", "RU_soil_t  ", "ETR mm     ", "Drainage   ", "R_growth   ", "A_net_tot  ", "A_gros_tot ", "Resp_tot   ", 
						"Export     ", "Export_tot ", "T_Soil1    ", "T_Soil2    ", "T_Soil3    ", "RU_soil_wp ", "REW_t      ", "Q_pet_s    ", "E_Petiole  ", "P_Pet_s    ", 
						"Pet_diam   ", "Q_Plant_a  ", "Q_Plant_s  ", "RWC_Leaf   ", "RWC_Leaf_s ", "RWC_Leaf_a ", "RWC_Branch ", "RWC_Br_s   ", "RWC_Br_a   ", "RWC_Trunk  ", 
						"RWC_Tr_s   ", "RWC_Tr_a   ", "RWC_Root   ", "RWC_Root_s ", "RWC_Root_a ", "RWC_Plant  ", "RWC_Plant_s", "RWC_Plant_a", "E_day mmol ", "P_min_lf_d ", 
						"P_max_lf_d ", "gs_min_d   ", "gs_max_d   ", "SF_min_d   ", "SF_max_d   ", "gs_max_d2  ", "Fruit_gr   ", "Teta_Soil  ", "Teta_Soil1 ", "Teta_Soil2 ", 
						"Teta_Soil3 ", "Ca  ppm    ", "K_Plant_20 ", "PLC_Plant  ", "P50_Leaf   ", "P50_bran   ", "P50_Trunk  ", "P50_Root   ", "PI0_Leaf   ", "Px_gs      ", 
						"PI_Leaf    ", "K_tot_20   ", "PLC_P_soil ", "PAR_tot    ", "P_mean_lf  ", "g_Leaf     ", "E_stomata  ", "gs_Jarvis  ", "g_crown    ", "g_bl       ", 
						"Cm         ", "VcMax      ", "VjMax      ", "Rd         ", "Ci         ", "Cbs        ", "PLC_Leaf1  ", "PLC_Leaf2  ", "PLC_Leaf3  ", "PLC_Branch1", 
						"PLC_Branch2", "PLC_Branch3", "gs_1       ", "gs_2       ", "gs_3       ", "P_Leaf_s1  ", "P_Leaf_s2  ", "P_Leaf_s3  ", "P_Leaf_a1  ", "P_Leaf_a2  ", 
						"P_Leaf_a3  ", "P_Branch_a1", "P_Branch_a2", "P_Branch_a3", "T_Leaf1    ", "T_Leaf2    ", "T_Leaf3    ", "E_Leaf1    ", "E_Leaf2    ", "E_Leaf3    ", 
						"A_net1     ", "A_net2     ", "A_net3     ", "PAR1       ", "PAR2       ", "PAR3       ", "Leaf_area1 ", "Leaf_area2 ", "Leaf_area3 ", "E_layer_t  ",
						"ETR_layer_t", "E_layer_d  ", "ETR_layer_d", "T_Branch1  ", "T_Branch2  ", "T_Branch3  ", "T_Trunk    ", "C_growth_t ", "R_main_t   ", "R_growth_t ",
						"T_Axil     ", "Turgor_Ax  ", "F_s_re1    ", "F_s_re2    ", "F_s_re3    ", "REW_t1     ", "REW_t2     ", "REW_t3     ", "REW_wp     ", "REW_wp1    ", 
						"REW_wp2    ", "REW_wp3    ", "Test1      ","Test2      "};
   
char filename_IN[]="ixxxxxxx_sureau.ini";
char filename_OUT[]="ixxxxxxx_sureau.out1";
char filename_OUT2[]="ixxxxxxx_sureau.out2";
char filename_CLIM[]="cxxxxxxx_climat_day.ini";
char filename_TRANS[]="xxxxxxx_transient.out1";
char filename_SOIL[]="xxxxxxxx_soil.trs";

long double DATA[NPAR];
long TT=0;
int DEAD=0;
int END_CLIMAT=0;
int N=1;
long double CAPILLARITY=1;  // allow flow between soil layers by capillarity
int PREM=1, PREM1=1;
long double RANDOMISE=0;
long double TRANSIENT=0, gs_cst=0,PRINT_SCREEN=0;
long double YEAR1=2000, YEAR2=2001, dt,dt_dyna,dt_stat,T,T0,days_simul,T_max;     // dt time interval in seconds; T time of the day in seconds
long double CUT;                    // to simulate a cut Branch (=1) or  Trunk (=2)
long double CLIMAT;                 // if 1 then climatic data are from a data file; if 0 they are generated
//long j;
long double Lat=45;                     // latitude, between -90 and 90°
long double Day_length;              // length of day light, hours
long double DOY_0=180,DOY;                     // day of year 0-365
long double dPLC_crit=0.01,REW_crit=0.8;

//Atmosphere
long double cloud_cover, POTENTIAL_PAR, CO2_atm, T_air, T_air_min, T_air_max, Cum_T_air,Cum_T_air_l,T_air_an, T_air_an_l,RH_air, RH_air_min, RH_air_max, RH_air_min_0;
long double P_air, VPD,VPD_Air_tot,VPD_Leaf_tot,PAR[5],PAR_max,Rain_1, Rain_2, Rain_tot,Rain_soil, Rain_Leaf_tot, Proba_rain, Wind[4];
long double T_0, T_1, T_air_1, RH_air_1,PAR_1, T_2, T_air_2, RH_air_2,PAR_2,T_air_max_0, T_air_min_1, T_air_max_1, RH_air_min_1, RH_air_max_1,PAR_tot,Pot_PAR_tot;
long double  Max_PAR_tot,PAR_max_1,T_air_min_2, T_air_max_2, RH_air_min_2, RH_air_max_2,PAR_max_2,VPD_Air,PAR_att=1;
long double HW=1;                   //Heat Wave
long double HW_day, HW_duration, HW_T;
FILE *climat_in;
long double HH1=12.00;  // time of Tmax
long double HH2=0.0;    // time of RHmax
long double INTERCEPTION,    Interception_min,    Threshold_rain,    Interception_factor; // for rain interception model
long double ETP_Penman,ETP_Penman_tot,ETP_Penman_day,ETP_Leaf_tot;
//SOIL
long double Rock_f1=0,Rock_f2=0,Rock_f3=0, Teta_Soil,Teta_Soil1,Teta_Soil2,Teta_Soil3, Teta_ini_1, Teta_ini_2, Teta_ini_3, REW_t,REW_t1,REW_t2,REW_t3,  T_REW_Soil, gap, VPD_Soil,  T_Soil,T_Soil1,T_Soil2,T_Soil3;
long double T_Soil_1, T_Soil_2,T_Soil_11, T_Soil_21, T_Soil_12, T_Soil_22,T_Soil_13, T_Soil_23,E_Soil, g_Soil,g_Soil0, dq_Soil,  Q_Soil1,Q_Soil2,Q_Soil3;
long double Q_Soil01,Q_Soil02,Q_Soil03, K_Soil1,K_Soil2,K_Soil3,K_Soil_tot1,K_Soil_tot2,K_Soil_tot3, K_Interface1,K_Interface2,K_Interface3;
long double P_Soil1, P_Soil2, P_Soil3, Teta_s_1,Teta_s_2,Teta_s_3, Teta_fc_1, Teta_fc_2,Teta_fc_3,Teta_r_1,Teta_r_2,Teta_r_3,Teta_wp_1,Teta_wp_2,Teta_wp_3;
long double alpha_1,alpha_2,alpha_3,m_1,m_2,m_3,n_1,n_2,n_3,K_sat_1,K_sat_2,K_sat_3,L,Volume_soil,Volume_soil1,Volume_soil2,Volume_soil3,Soil_Depth, Soil_Width;
long double  K_s,RWC1,RWC2,RWC3,RWC_Soil, RWC_Soil_1,RWC_Soil_2,RWC_Soil_3,RWC_min,RWC_int,REW_int1,REW_int2,REW_wp, REW_wp1,REW_wp2,REW_wp3,RWC_Irr,RWC_fc_1,RWC_fc_2,RWC_fc_3, P_soil,P_soil_min,P_soil_int;
long double Surface_Soil,  dq_Soil_Root_Symp,dq_Soil_12,dq_Soil_23,Fluidity_soil, Fluidity_soil1, Fluidity_soil2, Fluidity_soil3,Osmotic_TSoil,ST_Soil;
long double Layer_1,Layer_2,Layer_3;
long double REHYDRATE=1, PLC_REHYD=99.99, Daily_Irr=1, Drainage ;
long double PENMAN=1, Penman_Coeff=0.75;
long double Teta_soil, teta_c=0.232;
long double PLC_LIMIT=0, T_LIMIT, T_Soil_Crit;
long double Q_soil1_init, Q_soil2_init, Q_soil3_init,COMPET;
long double PI0_Soil=-1; // soil osmotic potential, MPa
long double Psoil_FC=33; // soil potential at field capacity, KPa
 
//PLANT
FILE *File_fractal1,*File_fractal2;
long double crown_diam, Length_term_shoot,side[1001],N_fractal, Nb_segment[1001],Diam_segment[1001],Diam_sapwood[1001],Length_segment[1001],scale,angle,X[1001],Y[1001];
long double LA_term_shoot=100, Length_fine_Root=50, Tapering=1.8,Root_shoot_ratio=0.2,Root_ramif=0.3;  //parameters for fractal
long double SF_min_d,SF_max_d,Q_Plant,Q_Plant_s,Q_Plant_a, Q_Plant0,Q_Plant_s0,Q_Plant_a0,K_Plant,Crown_Area;
long double REFILL,P_REFILL; // to allow xylem refilling above P_REFILL
long double SYMP_CAVIT=0; // then water released by cavitation goes only into the adjacent Symplasm; activated if REFILL=2
long double GPP, EvapoT_day,E_tot_day,GRAVITY, Pg_Leaf[4], Pg_Branch[4], Pg_Trunk,Extinction_Coeff, A_net_c, A_gross, A_gross_tot,A_net_tot_c,A_gross_tot_c,A_net_day_c,Resp_tot_c; // the gravimetric pressures
long double Regul_gs, Regul_gs_para1, Regul_gs_para2, Regul_gs_para3, Regul_gs_para4; // if 1 E regulated by Leaf turgor, if 0 E regulated by P12
long double ST_Leaf=1, ST_Air=1; // Water surface tension at Leaf and air temperatures relative to SF at 20 °C
long double Osmotic_TLeaf[4], Osmotic_TAir;
long double g_Canopy,gs_Jarvis,gs_max[4],gs_night, g_cuti[4], g_cuti_20, g_Branch, g_Trunk,g_Root[4], Jarvis_PAR,g_bl,g_crown,g_crown0,RWC_sens;
long double Reserve,Reserve0,  Rm,Rg, Q_Wood;
long double K_tot,K_tot_20,K_tot_20_0,K_tot0; // Hydraulic conductance soil to Leaf
long double Flow_save, Q_buffer;
long double Fluidity_Leaf[4],Fluidity_Branch[4],Fluidity_Trunk,Fluidity_air;
long double E_max=0, E_max_gs_close=0;
long double Px_gs; // the set xylem pressure for stomatal closure
long double turgor,K_to_Leaf_Apo,E_clim[4], E_t_B;
long double Slope_Leaf_Fall,P50_Leaf_Fall;
long double Tgs_optim=20, Tgs_sens=17; // variable to describe the effect of Leaf temperature on gs_max;
long double Leaf_Fall=0;
long double K_VAR=0, K_VAR_P1,K_VAR_P2,K_VAR_P3; // variable K: 0 = KLeaf CST KRoot CST 1= KLeaf VAR KRoot CST 2: KLeaf CST KRoot VAR 3: KLeaf VAR KRoot VAR
long double END_DEATH=0;
long double PLC_END; // the PLC defining hydraulic failure and computation end
long double gs_tc=0; //  time constant of stomatal response in s
long double P50_Leaf_Apo_0[4], P50_Branch_Apo_0[4], P50_Trunk_Apo_0, P50_Root_Apo_0, Pi0_Leaf_Symp_0,Pi0_Branch_Symp_0,Pi0_Trunk_Symp_0, Px_gs_0, RWC_END;


//Fractal tree
long double Q_Branch_Symp_FR, Q_Branch_Apo_FR, Q_Trunk_Sym_FR, Q_Trunk_Apo_FR, Q_Root_Symp_FR, Q_Root_Apo_FR, Branch_Area_FR, Trunk_Area_FR, Root_Area_FR;
long double Root_Area_fi_0,Root_Area_fi, K_Branch_Apo_FR, K_Trunk_Apo_FR, K_Root_Apo_FR, Length_Root_FR, Diam_Root_FR ;

// variables for Plant
long double Alpha, Leaf_Area[4], LAI_Crown,LAI_Soil, Branch_Area[4], Trunk_Area, Root_Area, Root_Area1,Root_Area2,Root_Area3, E_Leaf[4], dq_cuti[4], dq_stomate[4], dq_Branch[4], dq_Trunk;
long double dq_Root1,dq_Root2,dq_Root3, E_cuti[4], E_Branch[4], E_Trunk, E_Root1,E_Root2,E_Root3, t_out, E_tot, EvapoT, EvapoT_tot, g;
long double P_min_Leaf, P_min_stem, P_min_Leaf_mean_Leafy; // minimal seasonal water potentials
long double K_Plant_20,K_Plant_20_0,PLC_Plant,PLC_Plant_Soil; // Whole plant conductance à 20°C. Use to compute PLC_Plant as 100*(1-K_Plant_20/K_Plant_0)
long double F_Leaf[4], F_Branch[4];
// variabes for legacy of PLC
long double leg_Leaf=0, leg_Branch=0.5, leg_Trunk=0.75, leg_Root=0;

// variables for growth
long double Growth_Trunk, Growth_rate_Trunk, Growth_rate_Root,Growth_Fruit, Growth_rate_Fruit, Extensibility_Trunk, Yield_Trunk, Extensibility_Fruit, Yield_Fruit;
long double Growth_rate_Branch[4],Growth_rate_Root1,Growth_rate_Root2,Growth_rate_Root3;
long double R_main, R_growth, R_main_tot=0, R_growth_tot=0, C_growth_tot=0, Sapwood_V;
long double bark=0.1; // the portion of the Q_symp that counts for diameter variation
long double NSC=0.05; // the concentration of NCS in sapwood or leaves, in % or g/g of dry matter
long double Rm_25=100,Rg_25=0.25 ;  //Rm in umol Co2 /m3/s Rg is a fraction of C_growth
long double Growth_control=1, Growth_para1=180, Growth_para2=1; //sigmoidal growth for trunk

// variables for Leaf
long double Leaf_rain=0, gs_min_d,gs_max_d,gs_max_d2;
long double Pgs_88,Pgs_12,P_min_lf_d ,P_max_lf_d;
long double gs_0, GS_MAX, E_Leaf_night, T_Leaf_max, TP, Q10_1, Q10_2, LMA, Succulence, Leaf_Apo_fraction, PLF[4], g_s[4], Leaf_size, A_net[4], A_net_V, A_net_J;
long double A_net_max, A_net_tot, A_net_day, Turgor_lp, VPD_Leaf[4], VPD_Cuti[4],T_Leaf[4],Px_Leaf_Apo[4];
long double PLCx, P12_Leaf_Apo[4],PLC_Leaf_Apo[4], Q_Leaf_Evap[4], Q_Leaf_Evap0[4], Q_Leaf_Symp[4], Q_Leaf_Apo[4], dQ_Leaf_Apo[4], Q_Leaf_Symp0[4], Q_Leaf_Apo0[4];
long double Q_Leaf_Apo1[4], K_Leaf_Symp[4], K_Leaf_Symp_0[4], K_Leaf_Symp2[4];
long double K_Leaf_Apo[4], K_Leaf_Apo0[4], P_Leaf_Evap[4], P_Leaf_Symp[4], Turgor_Leaf_Symp[4],Turgor_Leaf_Symp_Ref, P_Leaf_Apo[4], C_Leaf_Evap[4], C_Leaf_Apo0[4],C_Leaf_Apo[4];
long double Epsilon_Leaf_Symp, Pi0_Leaf_Symp;
long double P50_Leaf_Symp, Slope_Leaf_Symp, P50_Leaf_Apo[4], Slope_Leaf_Apo[4];
long double LA_Var, LA_max_Pheno,LA_max_Pheno2, LA_min, LA_max, LA_max_init, LA_day1, LA_day2, LA_day3, LA_day4,GDD,S_GDD=200,T_base1=5,T_base2=10,LGE;
long double LA_para1,LA_para2,LA_para3,LA_para4;
long double PhotoS_model, Rd,VcMaxT, VjMaxT, VpMaxT,VcMax=110, VjMax=180, Kc25=40.4, Ko25=24800, Qye=0.12, Rd25=-0.37, Ca=400, Cm, Ci,Cbs;
long double a_Res, Resp_tot;
long double Export, Export_tot=0;   //exportation of µmol of glucose /m2/s
long double Leaf_angle[4];    //angle from horizontal    degrees
long double Gamma=0.8;       //an attenuation factor for E to reach a target potential
long double turgor_ref_factor;
long double gs_CO2_sens=0;
long double E_stomata[4];

// variables for Axil: Axillary bud /Fruit/Flower
long double Type_Axil, T_RWC_Axil, T_PLC_Axil, T_TLP_Axil, RWC_Axil, E_Axil, E_Petiole, T_Axil, N_Axil, dq_Axil, K_Axil_Apo0,K_Axil_Apo, K_Axil_Symp, K_Axil_Symp2, Epsilon_Axil_Symp, Pi0_Axil_Symp;
long double Petiole_diam, P50_Axil_Apo, Slope_Axil_Apo, VPD_Axil, VPD_Petiole, g_Axil, g_Axil_min20, g_Axil_min,g_Axil_max, g_Petiole , C_Axil_Apo, Diam_Axil, Axil_Area, Petiole_area, WC_Axil;
long double PLC_Axil_Apo,Q_Axil_Symp, Q_Axil_Apo, Q_Axil_Symp0, Q_Petiole_Symp, Q_Petiole_Symp0, Q_Axil_Apo0,Q_Axil_Apo1, P_Axil_Symp, P_Petiole_Symp,P_Axil_Apo,Turgor_Axil_Symp, dq_Fruit, Length_Petiole, Diam_Petiole;
long double g_Axil_regul,TP_Axil,Q10_1_Axil,Q10_2_Axil,T_Axil,g_bl_Axil;

// variables for Branch
long double RWC_Branch_s[4],RWC_Branch_s_min[4],RWC_shoot, Density, Branch_Symp_fraction, Branch_Apo_fraction,Length_Branch, Number_Branch, Diam_Branch, VPD_Branch[4],T_Branch,T_Branch1,T_Branch2,T_Branch3,T_Branch1_1,T_Branch2_1,T_Branch3_1,T_Branch1_2,T_Branch2_2,T_Branch3_2,PLC_Branch_Apo[4],Q_Branch_Symp[4], Q_Branch_Apo[4], Q_Branch_Symp0[4];
long double Q_Branch_Apo0[4], Q_Branch_Apo1[4], K_Branch_Symp[4], K_Branch_Symp0[4],K_Branch_Apo[4], K_Branch_Apo0[4],P_Branch_Symp[4], P_Branch_Apo[4], C_Branch_Apo[4], C_Branch_Apo0[4],Epsilon_Branch_Symp, Pi0_Branch_Symp;
long double P50_Branch_Symp, Slope_Branch_Symp, P50_Branch_Apo[4], Slope_Branch_Apo[4],Branch_distri[4];
long double Turgor_Branch_Symp[4];
long double Turgor_Root_Symp1,Turgor_Root_Symp2,Turgor_Root_Symp3;

// variables for Trunk
long double Radius_Trunk,Radius_Trunk_rel,Radius_bark_Trunk,Radius_bark_Trunk_rel,Trunk_sapwood_fraction, Trunk_Symp_fraction, Trunk_Apo_fraction,Length_Trunk, Diam_Trunk,Turgor_Trunk_Symp,VPD_Trunk, T_Trunk,T_Trunk_1,T_Trunk_2,PLC_Trunk_Apo,Q_Trunk_Symp, Q_Trunk_Apo;
long double Q_Trunk_Symp0, Q_Trunk_Apo0,Q_Trunk_Apo1, K_Trunk_Symp, K_Trunk_Symp0,K_Trunk_Apo,K_Trunk_Apo0, P_Trunk_Symp, P_Trunk_Apo, C_Trunk_Apo,C_Trunk_Apo0, Epsilon_Trunk_Symp, Pi0_Trunk_Symp;
long double P50_Trunk_Symp, Slope_Trunk_Symp, P50_Trunk_Apo, Slope_Trunk_Apo;

// variables for Root
long double K_Root_0, Root_Symp_fraction, Root_Apo_fraction,Length_Root, Length_Root_fi, Diam_Root,Q_Root_Endo0, Q_Root_Symp0, Q_Root_Apo_t0, Q_Root_Apo_t1,K_Root_Sympqq,K_Root_Symp1;
long double K_Root_Symp2,K_Root_Symp0,K_Root_Apo,K_Root_Apo0, C_Root_Endo, C_Root_Apo,C_Root_Apo0;
long double Epsilon_Root_Symp, Pi0_Root_Symp, P50_Root_Symp, Slope_Root_Symp, P50_Root_Apo, Slope_Root_Apo, PRF;
long double VPD_Root1,PLC_Root_Apo1,Q_Root_Endo1,Q_Root_Endo01, Q_Root_Symp1, Q_Root_Apo1, Q_Root_Symp01, Q_Root_Apo01, Q_Root_Apo11, K_Root_Symp11,K_Root_Symp21, K_Root_Apo1,K_Root_Apo01;
long double P_Root_Endo1, P_Root_Symp1, P_Root_Apo1, C_Root_Endo1, C_Root_Apo1,C_Root_Apo01, Epsilon_Root_Symp1, Pi0_Root_Symp1, P50_Root_Symp1, Slope_Root_Symp1, P50_Root_Apo1, Slope_Root_Apo1;
long double VPD_Root2,PLC_Root_Apo2,Q_Root_Endo2,Q_Root_Endo02, Q_Root_Symp2, Q_Root_Apo2, Q_Root_Symp02, Q_Root_Apo02, Q_Root_Apo12, K_Root_Symp12,K_Root_Symp22, K_Root_Apo2,K_Root_Apo02;
long double P_Root_Endo2, P_Root_Symp2, P_Root_Apo2, C_Root_Endo2, C_Root_Apo2,C_Root_Apo02, Epsilon_Root_Symp2, Pi0_Root_Symp2, P50_Root_Symp2, Slope_Root_Symp2, P50_Root_Apo2, Slope_Root_Apo2;
long double VPD_Root3,PLC_Root_Apo3,Q_Root_Endo3,Q_Root_Endo03, Q_Root_Symp3, Q_Root_Apo3, Q_Root_Symp03, Q_Root_Apo03, Q_Root_Apo13, K_Root_Symp13,K_Root_Symp23, K_Root_Apo3,K_Root_Apo03;
long double P_Root_Apo, Q_Root_Apo_t, P_Root_Endo3, P_Root_Symp3, P_Root_Apo3, C_Root_Endo3, C_Root_Apo3,C_Root_Apo03, Epsilon_Root_Symp3, Pi0_Root_Symp3, P50_Root_Symp3, Slope_Root_Symp3, P50_Root_Apo3, Slope_Root_Apo3;
long double T_Root,T_Root_1,T_Root_2,T_Root_3,Root_upper,Root_middle,Root_lower,Root_upper0,Root_middle0,Root_lower0;

// variables for dq
long double dq_Leaf_Apo_Leaf_Symp[4],dq_Axil_Apo_Axil_Symp, dq_Leaf_Symp_Leaf_Evap[4],dq_Branch_Apo_Leaf_Apo[4],dq_Branch_Apo_Axil_Apo, dq_Branch_Apo_Branch_Symp[4], dq_Trunk_Apo_Branch_Apo[4];
long double dq_Trunk_Apo_Trunk_Symp, dq_Root_Apo_Trunk_Apo1, dq_Root_Symp_Root_Apo1,dq_Root_Apo_Trunk_Apo2, dq_Root_Symp_Root_Apo2,dq_Root_Apo_Trunk_Apo3, dq_Root_Symp_Root_Apo3;
long double dq_Leaf_Apo_Leaf_Evap[4],dq_Root_Endo_Root_Apo1, dq_Root_Symp_Root_Endo1, dq_Soil_Root_Endo1,dq_Root_Endo_Root_Apo2, dq_Root_Symp_Root_Endo2, dq_Soil_Root_Endo2;
long double dq_Root_Endo_Root_Apo3, dq_Root_Symp_Root_Endo3, dq_Soil_Root_Endo3, dq_Branch_Symp_Axil_Symp,dq_Axil_Apo_Petiole_Symp, dq_Petiole;
long double dq_Root1_Root2,dq_Root1_Root3,dq_Root2_Root3, dq_Root_Apo_Trunk_Apo;

long double T_PLC_Leaf[4],T_TLP_Leaf[4],T_PLC_Branch[4],T_RWC_Branch[4],T_PLC_Trunk,T_PLC_Root,T_PLC_Root1, T_gs_close,T_gs_50mmol,T_gs_regul,T_budbreak, T_max_LAI;

void init(void)  //intilialize a bunch of variables
{
//	FILE *soil;
	int i;
	long double K_Root, K_Root1, K_Root2, K_Root3,K_Rhizo;

	Surface_Soil=3.1416*Soil_Width*Soil_Width/4; //surface of the soil in m2 assuming a circle
	Crown_Area=crown_diam*crown_diam/4*3.1416; //use to compute LAI
   // g_Axil_max= 67.39;
	T_Soil=20;
	T_Soil1=T_Soil;
	T_Soil2=T_Soil;
	T_Soil3=T_Soil;
	Teta_Soil1=Teta_ini_1;
	Teta_Soil2=Teta_ini_2;
	Teta_Soil3=Teta_ini_3;
	T_budbreak=T_max_LAI=0;
	P_min_Leaf_mean_Leafy=0;
	PAR_tot=0;
	Pot_PAR_tot=0;
	Max_PAR_tot=0;
	DOY=DOY_0;
	GDD=0;
	E_max=0, E_max_gs_close=0;
	if (DYNAMIC==0) dt=dt_stat;
	if (DYNAMIC==2) dt=dt_stat;
	if (DYNAMIC==1) dt=dt_dyna;
	PLC_LIMIT=0;
	P_min_lf_d=0;
	P_max_lf_d=-1000;
	gs_min_d=10000;
	gs_max_d=0;
	gs_max_d2=0;
	SF_min_d=10000;
	SF_max_d=0;
	LAI_Soil=LA_max_init/Surface_Soil;   // LAI in m2/m2
	LAI_Crown=LA_max_init/Crown_Area; 
	WARNING=0;
	T_air_an=0;
	T_air_an_l=0;
	Cum_T_air=0;
	Cum_T_air_l=0;
	N_days=0;
	N_days2=0;
	Rain_tot=0;
	Rain_Leaf_tot=0;
	VPD_Air_tot=0;
	VPD_Leaf_tot=0;
	ETP_Leaf_tot=0;
	A_net_day=0;
	A_net_tot=0;
	A_net_tot_c=0;
	E_tot_day=0;
	EvapoT_day=0;
	indice=1;
	indice_double=0;
	DEAD=0;
	END_CLIMAT=0;
	REW_t=1;REW_t1=1;REW_t2=1;REW_t3=1;
	for (i=1;i<4;i++)
	{
		P50_Leaf_Apo[i]=P50_Leaf_Apo_0[i];
		P50_Branch_Apo[i]=P50_Branch_Apo_0[i];
		PLF[i]=0;
		T_PLC_Leaf[i]=8.64e11;
		T_TLP_Leaf[i]=8.64e11;
		T_PLC_Branch[i]=8.64e11;
		T_RWC_Branch[i]=8.64e11;
		RWC_Branch_s_min[i]=8.64e11;
	}
	P50_Trunk_Apo=P50_Trunk_Apo_0;
	P50_Root_Apo=P50_Root_Apo_0;
	Pi0_Leaf_Symp=Pi0_Leaf_Symp_0;
	Pi0_Branch_Symp=Pi0_Branch_Symp_0;
	Pi0_Trunk_Symp=Pi0_Trunk_Symp_0;
	Px_gs=Px_gs_0;
	T_PLC_Trunk=8.64e11;
	T_PLC_Root=8.64e11;
	T_PLC_Root1=8.64e11;
	T_RWC_Axil=8.64e11;
	T_PLC_Axil=8.64e11;
	T_TLP_Axil=8.64e11;
	T_REW_Soil=8.64e11;
	T_gs_close=0;
	T_gs_50mmol=0;
	T_gs_regul=8.64e11;
	g_bl=1000;
	g_bl_Axil=1000;
	gs_0=0;
	for (i=1;i<4;i++) g_s[i]=0;
	ST_Leaf=1.0;
	ST_Air=1.0;
	for (i=1;i<4;i++)
	{
		T_Leaf[i]=T_air;
		g_cuti[i]=g_cuti_20;	
	}
	turgor=1;
	E_max=0;
	E_max_gs_close=0;
	T_Leaf_max=0;
	E_Leaf_night=100;
	if (LA_Var) LA_max_Pheno=LA_min;
	else LA_max_Pheno=LA_max_init;
		
	for (i=1;i<4;i++)
		{
		if (LA_Var)  Leaf_Area[i]=LA_min*Branch_distri[i]; 
		else Leaf_Area[i]=LA_max_init*Branch_distri[i];
		}
	
	Root_Area_fi=Root_Area_fi_0;
	Rain_soil=0;
	ETP_Leaf_tot=0;
	ETP_Penman_tot=0;
	g_Soil=g_Soil0;
	Drainage=0;
	if (CUT==1) // cut between Branch and Trunk
	{
		for (i=1;i<4;i++) K_Branch_Apo0[i]=0;
		Turgor_Leaf_Symp_Ref=-Pi0_Leaf_Symp;
	}
	
	if (CUT==2)  // cut between Trunk and Root
	{
		K_Trunk_Apo0=0;
		Turgor_Leaf_Symp_Ref=-Pi0_Leaf_Symp;
	}

	if (CUT==3)  // cut between leaf and branch
	{
		for (i=1;i<4;i++) K_Leaf_Apo0[i]=0;
		Turgor_Leaf_Symp_Ref=-Pi0_Leaf_Symp;
	}

	for (i=1;i<4;i++) Turgor_Leaf_Symp[i]=-Pi0_Leaf_Symp;
	for (i=1;i<4;i++) Leaf_Area[i]=LA_max_init*Branch_distri[i];
	
	// valeur en attendant
	for (i=1;i<4;i++) 
		{
		C_Leaf_Evap[i]   	= C_Leaf_Apo[i]/10;
		C_Leaf_Apo0[i]   	= C_Leaf_Apo[i];
		C_Branch_Apo0[i] 	= C_Branch_Apo[i];
		Q_Leaf_Evap0[i]  	= Q_Leaf_Apo[i]/1000;
		}
	C_Trunk_Apo0   = C_Trunk_Apo;
	C_Root_Apo01   = C_Root_Apo1;
	C_Root_Apo02   = C_Root_Apo2;
	C_Root_Apo03   = C_Root_Apo3;
	
	

	Q_Root_Endo01  = Q_Root_Apo_t0/1000;
	C_Root_Endo1   = C_Root_Apo/10;
	Q_Root_Endo02  = Q_Root_Apo_t0/1000;
	C_Root_Endo2   = C_Root_Apo/10;
	Q_Root_Endo03  = Q_Root_Apo_t0/1000;
	C_Root_Endo3   = C_Root_Apo/10;
	
   
	
	// initialise Roots
	Root_upper=Root_upper0;
	Root_middle=Root_middle0;
	Root_lower=Root_lower0;
	Root_Area1= Root_Area * Root_upper;
	Root_Area2= Root_Area * Root_middle;
	Root_Area3= Root_Area * Root_lower;
	Q_Root_Symp01= Q_Root_Symp0* Root_upper;
	Q_Root_Symp02= Q_Root_Symp0* Root_middle;
	Q_Root_Symp03= Q_Root_Symp0* Root_lower;
	Q_Root_Apo01= Q_Root_Apo_FR* Root_upper*1000*1000/18;
	Q_Root_Apo02= Q_Root_Apo_FR* Root_middle*1000*1000/18;
	Q_Root_Apo03= Q_Root_Apo_FR* Root_lower*1000*1000/18;
	Q_Root_Apo_t0=Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03;
	g_Root[1]=g_Root[0];
	g_Root[2]=g_Root[0];
	g_Root[3]=g_Root[0];
	C_Root_Apo1=C_Root_Apo;
	C_Root_Apo2=C_Root_Apo;
	C_Root_Apo3=C_Root_Apo;
	Epsilon_Root_Symp1=Epsilon_Root_Symp;
	Epsilon_Root_Symp2=Epsilon_Root_Symp;
	Epsilon_Root_Symp3=Epsilon_Root_Symp;
	Pi0_Root_Symp1=Pi0_Root_Symp;
	Pi0_Root_Symp2=Pi0_Root_Symp;
	Pi0_Root_Symp3=Pi0_Root_Symp;
	K_Root_Apo01=K_Root_Apo0 * Root_upper;
	K_Root_Apo02=K_Root_Apo0 * Root_middle;
	K_Root_Apo03=K_Root_Apo0 * Root_lower;
	K_Root_Symp1    = K_Root_Symp0 * Root_Area_fi;  // only the fine Roots absorb water
	K_Root_Symp2    = K_Root_Symp0 * Root_Area_FR;  // all the Root surface has a Symplasmic compartment
	K_Root_Symp11=K_Root_Symp1 * Root_upper;
	K_Root_Symp12=K_Root_Symp1 * Root_middle;
	K_Root_Symp13=K_Root_Symp1 * Root_lower;
	P50_Root_Apo1=P50_Root_Apo;
	P50_Root_Apo2=P50_Root_Apo;
	P50_Root_Apo3=P50_Root_Apo;
	Slope_Root_Apo1=Slope_Root_Apo;
	Slope_Root_Apo2=Slope_Root_Apo;
	Slope_Root_Apo3=Slope_Root_Apo;
	P_min_Leaf = P_min_stem=0;
	Compute_ST();
	for (i=1;i<4;i++) Px_Leaf_Apo[i]= P50_Leaf_Apo[i]*ST_Leaf+ 25/Slope_Leaf_Apo[i]*log((100-PLCx)/PLCx); //compute the Pressure for x PLC
	A_net_tot=0;
	A_net_max=0;
	// Growth=0;
	Sapwood_V=(Length_Trunk+Length_Branch)*3.1416*Diam_Trunk*Diam_Trunk/4*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction; //whole aboveground tree sapwood volume in m3
	Sapwood_V+=Root_shoot_ratio*Sapwood_V; //includes root sapwood

	Reserve0=NSC*Sapwood_V*1000*1000*Density/1000/180*6;   	// whole NSC content in moles equ CO2 with 180= MW glucose and 6 C par glusose assuming 0.05 g g-1 or 5%
	Reserve0 +=NSC*LA_max_init*LMA/180*6;						// added NSC in the leaves
	Reserve=Reserve0;
	//Test[0]=Sapwood_V;
	
	m_1=(1-1/n_1);
	m_2=(1-1/n_2);
	m_3=(1-1/n_3);
	Flow_save=0;
	//Volume_soil  = Surface_Soil*Soil_Depth*(1-1/3*(Rock_f1+Rock_f2+Rock_f3)); //rock in the soil contain no available water for Roots
	Volume_soil1 = Surface_Soil*Soil_Depth*Layer_1*(1-Rock_f1); //m3
	Volume_soil2 = Surface_Soil*Soil_Depth*Layer_2*(1-Rock_f2);
	Volume_soil3 = Surface_Soil*Soil_Depth*Layer_3*(1-Rock_f3);
	Volume_soil  = Volume_soil1 + Volume_soil2 +Volume_soil3;

	if (CONTINUOUS==0 || PREM==1)
	{
		
		
			RWC1=RWC_fc_1;  // set RWC to field capacity
			RWC2=RWC_fc_2;
			RWC3=RWC_fc_3;
			RWC_Soil_1=RWC_fc_1;
			RWC_Soil_2=RWC_fc_2;
			RWC_Soil_3=RWC_fc_3;
			RWC_Soil=RWC_fc_1*Layer_1+RWC_fc_2*Layer_2+RWC_fc_3*Layer_3;
			RWC_min=RWC_fc_1*Layer_1+RWC_fc_2*Layer_2+RWC_fc_3*Layer_3;
			RWC_int=RWC_fc_1*Layer_1+RWC_fc_2*Layer_2+RWC_fc_3*Layer_3;

			P_soil=-Psoil_FC/1000;  // water potential at field capacity (Pf=2.5=-33kPa)
			P_soil_min=-Psoil_FC/1000;
			P_soil_int=-Psoil_FC/1000;
			P_Soil1=-Psoil_FC/1000;
			P_Soil2=-Psoil_FC/1000;
			P_Soil3=-Psoil_FC/1000;
			if (Teta_fc_1<Teta_ini_1)	Q_Soil01 = Teta_fc_1*Volume_soil1*1000*1000*1000/18;    // water in bulk soil in mmol at field capacity
			else 					Q_Soil01 = Teta_ini_1*Volume_soil1*1000*1000*1000/18;  
			if (Teta_fc_2<Teta_ini_2)	Q_Soil02 = Teta_fc_2*Volume_soil2*1000*1000*1000/18;
			else					Q_Soil02 = Teta_ini_2*Volume_soil2*1000*1000*1000/18;
			if (Teta_fc_3<Teta_ini_3)	Q_Soil03 = Teta_fc_3*Volume_soil3*1000*1000*1000/18;	
			else					Q_Soil03 = Teta_ini_3*Volume_soil3*1000*1000*1000/18;	
		
		
		Q_Soil1        = Q_Soil01;
		Q_Soil2        = Q_Soil02;
		Q_Soil3        = Q_Soil03;
		Q_soil1_init   = Q_Soil01;
		Q_soil2_init   = Q_Soil02;
		Q_soil3_init   = Q_Soil03;
		
		K_Soil1=1E8;
		K_Soil2=1E8;
		K_Soil3=1E8;
		K_Interface1= K_Soil1*10;
		K_Interface2= K_Soil2*10;
		K_Interface3= K_Soil3*10;
		P_Branch_Apo[1]=P_Branch_Apo[2]=P_Branch_Apo[3]=0;
		P_Axil_Apo=0;
		P_Leaf_Apo[1]=P_Leaf_Apo[2]=P_Leaf_Apo[3]=0;
		P_Trunk_Apo=0;
		P_Root_Apo1=0;
		P_Root_Apo2=0;
		P_Root_Apo3=0;
		Drainage=0;
	}
	
	for (i=1;i<4;i++) T_Leaf[i]=T_air_min; //Leaf temperature at least equal to min air temperature to begin with
	E_tot=0;
	EvapoT_tot=0;
	for (i=1;i<4;i++)
		{
		dq_Branch_Apo_Branch_Symp[i]=0;
		dq_Branch_Apo_Leaf_Apo[i]=0;
		dq_Leaf_Apo_Leaf_Symp[i]=0;
		dq_Trunk_Apo_Branch_Apo[i]=0;

		}
	dq_Root_Apo_Trunk_Apo1=0;
	dq_Root_Symp_Root_Apo1=0;
	dq_Root_Apo_Trunk_Apo2=0;
	dq_Root_Symp_Root_Apo2=0;
	dq_Root_Apo_Trunk_Apo3=0;
	dq_Root_Symp_Root_Apo3=0;
	dq_Trunk_Apo_Trunk_Symp=0;
	for (i=1;i<4;i++) dq_Branch[i]=E_Branch[i]*Branch_Area[i]*dt;
	dq_Trunk=E_Trunk*Trunk_Area*dt;
	dq_Root1=E_Root1*Root_Area1*dt;
	dq_Root2=E_Root2*Root_Area2*dt;
	dq_Root3=E_Root3*Root_Area3*dt;
	for (i=1;i<4;i++) K_Leaf_Symp[i] = K_Leaf_Symp_0[i];
	for (i=1;i<4;i++) K_Leaf_Symp2[i]=1*K_Leaf_Symp[i];                   // assume the resitance of the Symplasmic sap pathway is the same as the resitance to the Symplasmis reservoir
	K_Root_Symp21=K_Root_Symp2 * Root_upper;
	K_Root_Symp22=K_Root_Symp2 * Root_middle;
	K_Root_Symp23=K_Root_Symp2 * Root_lower;
	
	for (i=1;i<4;i++)K_Leaf_Apo[i]=K_Leaf_Apo0[i];
	K_Axil_Apo=K_Axil_Apo0;
	for (i=1;i<4;i++) K_Branch_Apo[i]=K_Branch_Apo0[i];
	K_Trunk_Apo=K_Trunk_Apo0;
	K_Root_Apo1=K_Root_Apo01;
	K_Root_Apo2=K_Root_Apo02;
	K_Root_Apo3=K_Root_Apo03;
	
	if (K_Soil1 && K_Interface1 && K_Root_Symp11 && K_Root_Apo1) K_Root1=1/(1/K_Soil1 + 1/K_Interface1 + 1/K_Root_Symp11 + 1/K_Root_Apo1); else K_Root1=0;
	if (K_Soil2 && K_Interface2 && K_Root_Symp12 && K_Root_Apo2) K_Root2=1/(1/K_Soil2 + 1/K_Interface2 + 1/K_Root_Symp12 + 1/K_Root_Apo2); else K_Root2=0;
	if (K_Soil3 && K_Interface3 && K_Root_Symp13 && K_Root_Apo3) K_Root3=1/(1/K_Soil3 + 1/K_Interface3 + 1/K_Root_Symp13 + 1/K_Root_Apo3); else K_Root3=0;
	K_Rhizo=K_Root1+K_Root2+K_Root3;
 
 
   
	K_Leaf_Apo[0] = K_Leaf_Symp[0]= K_Branch_Apo[0]=0;
	for (i=1;i<4;i++) 
	{
	K_Leaf_Apo[0]+=		K_Leaf_Apo[i];
	K_Leaf_Symp[0]+=	K_Leaf_Symp[i];
	K_Branch_Apo[0]+=	K_Branch_Apo[i]; 
	}           
	if (K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo && K_Rhizo) K_tot_20_0=  1/(1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0] + 1/K_Branch_Apo[0]+ 1/K_Trunk_Apo + 1/K_Rhizo);  
	else  K_tot_20_0=0;

	if (K_Root_Symp11 && K_Root_Apo1) K_Root1=1/(1/K_Root_Symp11 + 1/K_Root_Apo1); else K_Root1=0;
	if (K_Root_Symp12 && K_Root_Apo2) K_Root2=1/(1/K_Root_Symp12 + 1/K_Root_Apo2); else K_Root2=0;
	if (K_Root_Symp13 && K_Root_Apo3) K_Root3=1/(1/K_Root_Symp13 + 1/K_Root_Apo3); else K_Root3=0;
	K_Root=K_Root1+K_Root2+K_Root3;
	K_Root_0=K_Root;
	
	if (K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo && K_Root) K_Plant_20_0=  1/(1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0] + 1/K_Branch_Apo[0]+ 1/K_Trunk_Apo + 1/K_Root);  
	else  K_Plant_20_0=0;
	
	for (i=1;i<4;i++) P12_Leaf_Apo[i]=P50_Leaf_Apo[i]+50/Slope_Leaf_Apo[i];
	if (CLIMAT!=2)
	{
		for (i=1;i<4;i++) 
			{
				PLC_Leaf_Apo[i]=0;
				PLC_Branch_Apo[i]=0;
			}
		PLC_Trunk_Apo=0;
		PLC_Root_Apo1=0;
		PLC_Root_Apo2=0;
		PLC_Root_Apo3=0;
	}
	
		
	Irrigation     =0;
	Turgor_lp= Pi0_Leaf_Symp*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp+Epsilon_Leaf_Symp);
	
	for (i=1;i<4;i++)
		{
			if (GRAVITY==1)// gravimetric pressure drop at the top leaves, Trunk and Branches in MPa
			{
				Pg_Leaf[i]		= -9.81*(Length_Branch*(1.5-(double)i/2)+Length_Trunk)/1000;
				Pg_Trunk		= -9.81* Length_Trunk/1000;                     			
				Pg_Branch[i] 	= -9.81*(Length_Branch*(1.5-(double)i/2)+Length_Trunk)/1000;
			}
			else
			{
				Pg_Leaf[i]=0;
				Pg_Trunk =0;
				Pg_Branch[i]=0;
			}
		}
	for (i=1;i<4;i++)
		{
		P_Leaf_Evap[i]=Pg_Branch[i];        // the water potential at the site of evaporation in the Leaf
		P_Leaf_Symp[i]=Pg_Branch[i];         // Leaf Symplasmic water potential
		P_Leaf_Apo[i]=Pg_Branch[i];
		P_Branch_Apo[i]=Pg_Branch[i];
		P_Branch_Symp[i]=Pg_Branch[i];
		}

	P_Axil_Apo=Pg_Branch[1];
	P_Axil_Symp=Pg_Branch[1];
	P_Trunk_Apo=Pg_Trunk;
	P_Trunk_Symp=Pg_Trunk;
	P_Root_Apo1=0;
	P_Root_Symp1=0;
	P_Root_Apo2=0;
	P_Root_Symp2=0;
	P_Root_Apo3=0;
	P_Root_Symp3=0;
	
	if(GRAVITY==0 || GRAVITY==3)
	{
		for (i=1;i<4;i++)
		{
		Q_Leaf_Evap[i]    = Q_Leaf_Evap0[i];
		Q_Leaf_Symp [i]   = Q_Leaf_Symp0[i];
		Q_Leaf_Apo[i]     = Q_Leaf_Apo0[i];
		Q_Leaf_Apo1[i]   = Q_Leaf_Apo0[i];
		Q_Branch_Symp[i]  = Q_Branch_Symp0[i];
		Q_Branch_Apo[i]  = Q_Branch_Apo0[i];
		Q_Branch_Apo1[i]  = Q_Branch_Apo0[i];
		}
		Q_Axil_Symp    = Q_Axil_Symp0;
		Q_Axil_Apo     = Q_Axil_Apo0;
		Q_Axil_Apo1    = Q_Axil_Apo0;
		Q_Trunk_Symp   = Q_Trunk_Symp0;
		Q_Trunk_Apo    = Q_Trunk_Apo0;
		Q_Trunk_Apo1   = Q_Trunk_Apo0;
	}
	
	else // compute steady-state Q at Pg with C (apo) or PV curves (Symp); assume potential is above turgor loss point
	{
		long double discri, tlp;
		for (i=1;i<4;i++)
		{
		Q_Leaf_Evap[i]    = Q_Leaf_Evap0[i]   + C_Leaf_Evap[i]   * P_Leaf_Evap[i] ;
		Q_Leaf_Apo[i]     = Q_Leaf_Apo0[i]    + C_Leaf_Apo[i]    * P_Leaf_Apo[i] ;
		Q_Leaf_Apo1[i]    = Q_Leaf_Apo[i] ;
		Q_Branch_Apo[i]   = Q_Branch_Apo0[i]  + C_Branch_Apo[i]  * P_Branch_Apo[i] ;
		Q_Branch_Apo1[i]  = Q_Branch_Apo[i];
		}
		Q_Trunk_Apo    = Q_Trunk_Apo0   + C_Trunk_Apo   * P_Trunk_Apo;
		Q_Trunk_Apo1   = Q_Trunk_Apo;
		Q_Axil_Apo     = Q_Axil_Apo0   + C_Axil_Apo   * P_Axil_Apo  ;
		Q_Axil_Apo1    = Q_Axil_Apo;
	
		// LEAF
		tlp=Pi0_Leaf_Symp*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp+Epsilon_Leaf_Symp);
		for (i=1;i<4;i++) 
		{
			if (P_Leaf_Symp[i]>tlp)
			{
			discri=pow(Pi0_Leaf_Symp*Osmotic_TLeaf[i] + Epsilon_Leaf_Symp + P_Leaf_Symp[i],2)-4*Pi0_Leaf_Symp*Osmotic_TLeaf[i]*Epsilon_Leaf_Symp;
			Q_Leaf_Symp[i] = ((Pi0_Leaf_Symp*Osmotic_TLeaf[i] + Epsilon_Leaf_Symp + P_Leaf_Symp[i]) + pow(discri,0.5))/(2*Epsilon_Leaf_Symp)*Q_Leaf_Symp0[i];
			}
			else   Q_Leaf_Symp[i] = Q_Leaf_Symp0[i]*Pi0_Leaf_Symp*Osmotic_TLeaf[i]/P_Leaf_Symp[i];
		}
		
		// BUD Flower
		if (Type_Axil)
		{
			tlp=Pi0_Axil_Symp*Epsilon_Axil_Symp/(Pi0_Axil_Symp+Epsilon_Axil_Symp);
			if (P_Axil_Symp>tlp)
			{
				discri=pow(Pi0_Axil_Symp*Osmotic_TLeaf[1] + Epsilon_Axil_Symp + P_Axil_Symp,2)-4*Pi0_Axil_Symp*Osmotic_TLeaf[1]*Epsilon_Axil_Symp;
				Q_Axil_Symp = ((Pi0_Axil_Symp*Osmotic_TLeaf[1] + Epsilon_Axil_Symp + P_Axil_Symp) + pow(discri,0.5))/(2*Epsilon_Axil_Symp/Q_Axil_Symp0);
			}
			else Q_Axil_Symp = Q_Axil_Symp0*Pi0_Axil_Symp*Osmotic_TLeaf[1]/P_Axil_Symp;
			if (Type_Axil>=2)  // petiole
			{
				tlp=Pi0_Axil_Symp*Epsilon_Axil_Symp/(Pi0_Axil_Symp+Epsilon_Axil_Symp);
				if (P_Petiole_Symp>tlp)
				{
					discri=pow(Pi0_Axil_Symp*Osmotic_TLeaf[1] + Epsilon_Axil_Symp + P_Petiole_Symp,2)-4*Pi0_Axil_Symp*Osmotic_TLeaf[1]*Epsilon_Axil_Symp;
					Q_Petiole_Symp = ((Pi0_Axil_Symp*Osmotic_TLeaf[1] + Epsilon_Axil_Symp + P_Petiole_Symp) + pow(discri,0.5))/(2*Epsilon_Axil_Symp/Q_Petiole_Symp0);
				}
				else Q_Petiole_Symp = Q_Petiole_Symp0*Pi0_Axil_Symp*Osmotic_TLeaf[1]/P_Petiole_Symp;
			}
		}
		// BRANCH
		tlp=Pi0_Branch_Symp*Epsilon_Branch_Symp/(Pi0_Branch_Symp+Epsilon_Branch_Symp);
		for (i=1;i<4;i++)
		{
			if (P_Branch_Symp[i]>tlp)
			{
			discri=pow(Pi0_Branch_Symp*Osmotic_TAir + Epsilon_Branch_Symp + P_Branch_Symp[i],2)-4*Pi0_Branch_Symp*Osmotic_TAir*Epsilon_Branch_Symp;
			Q_Branch_Symp[i] = ((Pi0_Branch_Symp*Osmotic_TAir + Epsilon_Branch_Symp + P_Branch_Symp[i]) + pow(discri,0.5))/(2*Epsilon_Branch_Symp/Q_Branch_Symp0[i]);
			}
		else Q_Branch_Symp[i] = Q_Branch_Symp0[i]*Pi0_Branch_Symp*Osmotic_TAir/P_Branch_Symp[i];
		}
		// TRUNK
		tlp=Pi0_Trunk_Symp*Epsilon_Trunk_Symp/(Pi0_Trunk_Symp+Epsilon_Trunk_Symp);
		if (P_Trunk_Symp>tlp)
		{
			discri=pow(Pi0_Trunk_Symp*Osmotic_TAir + Epsilon_Trunk_Symp + P_Trunk_Symp,2)-4*Pi0_Trunk_Symp*Osmotic_TAir*Epsilon_Trunk_Symp;
			Q_Trunk_Symp = ((Pi0_Trunk_Symp*Osmotic_TAir + Epsilon_Trunk_Symp + P_Trunk_Symp) + pow(discri,0.5))/(2*Epsilon_Trunk_Symp/Q_Trunk_Symp0);
		}
		else  Q_Trunk_Symp = Q_Trunk_Symp0*Pi0_Trunk_Symp*Osmotic_TAir/P_Trunk_Symp;
	}
	
	Q_Root_Apo1    = Q_Root_Apo01;
	Q_Root_Apo11   = Q_Root_Apo01;
	Q_Root_Symp1   = Q_Root_Symp01;
	Q_Root_Endo1   = Q_Root_Endo01;
	Q_Root_Apo2    = Q_Root_Apo02;
	Q_Root_Apo12   = Q_Root_Apo02;
	Q_Root_Symp2   = Q_Root_Symp02;
	Q_Root_Endo2   = Q_Root_Endo02;
	Q_Root_Apo3    = Q_Root_Apo03;
	Q_Root_Apo13   = Q_Root_Apo03;
	Q_Root_Symp3   = Q_Root_Symp03;
	Q_Root_Endo3   = Q_Root_Endo03;
	Q_Root_Apo_t	= Q_Root_Apo_t0;
	Q_Root_Apo_t1	= Q_Root_Apo_t0;
	
	Growth_Trunk=Q_Trunk_Symp0;
	
	if (K_VAR==12) //set native embolism in Branches 
		{
		K_Branch_Apo0[1] =	K_Branch_Apo0[1]*(100-K_VAR_P1)/100;
		K_Branch_Apo0[2] =	K_Branch_Apo0[2]*(100-K_VAR_P2)/100;
		K_Branch_Apo0[3] =	K_Branch_Apo0[3]*(100-K_VAR_P3)/100;
		}

}

void legacy(void)
{
	int i;
	for (i=1;i<4;i++) //Legacy effect on PLC
	{
		PLC_Leaf_Apo[i]*=leg_Leaf;
		PLC_Branch_Apo[i]*=leg_Branch;
	}
	PLC_Trunk_Apo*=leg_Trunk;
	PLC_Axil_Apo*=leg_Leaf;  //same as for LeafPLC_Trunk_Apo*=leg_Trunk;
	PLC_Root_Apo1*=leg_Root;
	PLC_Root_Apo2*=leg_Root;
	PLC_Root_Apo3*=leg_Root;
	Compute_Cavitation(dt);
}

void Reset (void) // When it is a new run or a new year. Reset simulation to zero
{
	FILE *transient;
	int i;
	DEAD=0;
	//printf("reset");
	if (DYNAMIC0==0) dt=dt_stat;
	if (DYNAMIC0>=2) dt=dt_stat;
	if (DYNAMIC0==1) dt=dt_dyna;  
	if (TRANSIENT)
	{
		transient = fopen(filename_TRANS,"a");
	//	fprintf(transient, "\n");
		fclose(transient);
	}
	
	if (PRINT_SCREEN) 
	{
		printf(" THF_l: %.1Lf PLC: %.1Lf",T_PLC_Leaf[1]/3600/24/1000,PLC_Leaf_Apo[0]);
		printf("# %.0Lf Yr %.0Lf\n", Simul, YEAR2);		
	}
		else 
			{
				printf(" THF_l: %.1Lf PLC: %.1Lf",T_PLC_Leaf[1]/3600/24/1000,PLC_Leaf_Apo[0]);
				printf("\n# %.0Lf Yr %.0Lf", Simul, YEAR2);		
			}

//  T_budbreak=T_max_LAI=0;
	indice=1;
	indice_double=0;
	gs_0=0;
	T_air_1=T_air_2;
	RH_air_1=RH_air_2;
	PAR_1=PAR_2;
	T_Soil_1=T_Soil_2;
			
	T_1=T_2;
	T_air_min=T_air_min_2;  //set values for next day loaded before
	T_air_max=T_air_max_2;
	RH_air_min=RH_air_min_2;
	RH_air_max=RH_air_max_2;
	PAR_max=PAR_max_2*PAR_att;
	Rain_1=Rain_2;
	Pot_PAR_tot=0;
	PAR_tot=0;
	C_growth_tot=0;
	R_growth_tot=0;
	R_main_tot=0;
	if (CLIMAT==4)
	{
		T_Soil_11=T_Soil_21;
		T_Soil_12=T_Soil_22;
		T_Soil_13=T_Soil_23;
		if (T_Soil_11<T_Soil_Crit || T_Soil_12<T_Soil_Crit || T_Soil_13<T_Soil_Crit) T_LIMIT=1;
		else if (T_Soil_11>T_Soil_Crit+0.1 && T_Soil_12>T_Soil_Crit+0.1 && T_Soil_13>T_Soil_Crit+0.1) T_LIMIT=0;
		else T_LIMIT=1;
	}
	else   T_Soil_1=T_Soil_2;
		
	legacy();
	
	ETP_Penman_tot=0;
	for (i=1;i<4;i++)
	{
		Q_Leaf_Apo[i]=Q_Leaf_Apo0[i];
		Q_Leaf_Apo1[i]=Q_Leaf_Apo0[i];
		Q_Branch_Apo[i]=Q_Branch_Apo0[i];
		Q_Branch_Apo1[i]=Q_Branch_Apo0[i];
	}
	for (i=1;i<4;i++) g_s[i]=0;
	Q_Axil_Apo=Q_Axil_Apo0;
	Q_Axil_Apo1=Q_Axil_Apo0;	
	
	Q_Trunk_Apo=Q_Trunk_Apo0;
	Q_Trunk_Apo1=Q_Trunk_Apo0;
	
	Q_Root_Apo1=Q_Root_Apo01;
	Q_Root_Apo2=Q_Root_Apo02;
	Q_Root_Apo3=Q_Root_Apo03;

	Q_Root_Apo11=Q_Root_Apo01;
	Q_Root_Apo12=Q_Root_Apo02;
	Q_Root_Apo13=Q_Root_Apo03;
	Q_Root_Apo_t=Q_Root_Apo_t0;
	
	E_tot=0;EvapoT_tot=0;A_net_tot=0;A_net_tot_c=0;ETP_Penman_tot=0;Irrigation=0;RWC_int=0;REW_int1=0;REW_int2=0;REW_wp=0;P_soil_int=0;P_soil_min=0;P_min_stem=P_min_Leaf=0;
	RWC_min=1;Rain_tot=0;Rain_soil=0;VPD_Air_tot=0;VPD_Leaf_tot=0;Rain_Leaf_tot=0;ETP_Leaf_tot=0;Cum_T_air_l=0;g_Soil=g_Soil0;A_net_max=0;
	Drainage=0;A_gross_tot=0;A_gross_tot_c=0;Resp_tot=0;Resp_tot_c=0;A_net_day=0;E_tot_day=0;EvapoT_day=0;P_min_Leaf_mean_Leafy=0;P_min_lf_d=0;
	Q_Axil_Symp0=N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil*1000*1000/18;
	Q_Axil_Symp=Q_Axil_Symp0;
	init();
	Phenology(dt);
	//printf(" reset done");
}

void default_para(void) // create a init file with default values when init file is missing; 
{
	FILE *new_ini;
	
	printf("SurEau_ini.txt not found, using default values\n");
	new_ini = fopen("sureau_ini.txt","w");
	fprintf(new_ini,"1	0.0001	0.95	0	11	0	1	0	0.01	60.000	2	0.000	0	0	1	99.9	1	1	1	1	1	1	0	182	45.000	15.400	27.600	32.575	95.000	1532.495	1	12	2	0.000	3.000	3.000	15	0.00E+00	0.28	0.1	0.0005	2	5.00	0.5	0	0	0	1.1195	1.3229	2.5000	0.4000	0.1000	0.0000	0	1	20	2	1	0	0.99	10	0	192	0	90	2	100	5	1.75	1.8	0.2	0.25	6.000	2.500E+00	0.400	0.200	3.000E-03	5.7839E+00	2.4265E+00	4.8529E+00	700.000	10	8.590E-02	0.400	0.200	0.331	2.6912E+00	6.3854E+00	1.2771E+01	1.114E+03	1.000E-03	1.8043E+01	3.5000E+00	3.8358E+03	1.4713E-03	1.9764E+00	3.9528E+00	0.400	0.200	50	45	1	10.5	0	123	140	285	325	0	-0.5	200	100	0.250	45.849	0	0.7	45	0.0000	15	0.4150	-1.50	-1.90	-2.60	200	20.000	0	0.006	7	-5	25.000	17.000	1	3.00	37.500	1.200	4.800	3.000	3.000	1.000	3.00E+01	110.000	180	3.000	-0.370	0.120	41.140	27350.000	-1.815	0.5	10.0000	10.0000	10.0000	10.0000	-2.1000	-1.0000	-1.0000	-1.0000	0	-1.275	181.8	0.169	5.0000	5.00E+03	5.00E+04	2.50E+04	8.2563E+01	3.9442E+02	3.5708E+03	1.8000	1.0000	1.0000	1.4000	-3.4	-3.4	-3.4	-3.4	60	60	60	60	2.000	-1.000	0.000	0.000	0.000	0.000	1.0000E-01	1.0000E-01	1.0000E-01	1.0000E-01	0.0000E+00	9.0000E-01	0	0.0000E+00	0.001	0.010	0.002	0.500	2.0000E+00	5.0000E-03	2.0000E-02	1.000E+01	-2.100	-3.400	60.000	1.000	1.0000E+00	1.0000E+00	5.000E-02	0.000	0.500");
	fclose(new_ini);
}

void default_para_file(void)
{
	FILE *new_para;
	
	printf("SurEau_para.txt not found, using default values\n");
	new_para = fopen("sureau_para.txt","w");
	fprintf(new_para,"1	2	3	4	5	205	6	7	14	15	16	195	10	13	11	12	18	19	20	21	23	63	17	30	31	32	33	34	35	36	37	199	200	38	39	40	41	194	42	43	44	45	46	47	48	196	197	139	140	141	142	143	144	24	52	53	54	55	56	57	58	59	60	61	62	103	206	211	210	207	208	209	88	89	91	92	90	110	104	105	93	94	95	96	97	98	111	106	107	99	100	112	113	114	115	108	109	101	102	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	49	50	51	204	25	26	27	28	29	119	120	8	121	122	198	123	124	22	125	126	127	128	129	130	131	132	80	81	82	83	84	85	86	87	64	145	146	147	148	149	150	151	152	9	201	202	203	153	154	155	156	116	117	118	157	158	159	160	161	162	163	164	165	166	167	168	188	189	190	191	192	193	135	136	137	138	133	134	169	170	182	183	184	185	171	172	173	174	175	176	177	178	179	180	181	186	187");
	fclose(new_para);
}


void Leaf_C_budget(long double dt_long, int i)
{
	long double export0=0.0;             //0.74477;   //exportation of 0.833 µmol of glucose /m2/s
	//long double Reserve0=50000;          // maximum reserve in µmol Assuming 10% mg/g de NSC in leave and a LMA of 90g/m2 then 9g NSC/m2 or Reserve0=50 mmol equ glucose with MW=180
	
	// Reserve in µmol of glucose
	
	Export=-Turgor_Leaf_Symp[i]/Pi0_Leaf_Symp*export0;
	Export_tot+=Export*dt_long*dt;
	//if (g_s) Reserve+=(A_gross*dt_long*dt/6-export*dt_long*dt);   //in µmol of glucose; 6 CO2 for 1 glucose
	//else
	Reserve+=(Rd25*exp(18.72-46390/(8.314*(T_Leaf[i]+273.15)))*dt_long*dt/6);  // if the carbon loss is due only to the photosynthesis respiration if CO2 from Glucose is recycled
	if (Reserve>Reserve0) Reserve=Reserve0;
}

long double Respiration(int i)  //under development;
{
	long double  Tref=25,  RTref, Resp;
	// long double temp, Q10=2, Rm_25=100;
	
	// Q_Wood= (Q_Branch_Symp + Q_Branch_Apo + Q_Trunk_Symp + Q_Trunk_Apo + Q_Root_Symp1 + Q_Root_Apo1 + Q_Root_Endo1+ Q_Root_Symp2 + Q_Root_Apo2 + Q_Root_Endo2+ Q_Root_Symp3 + Q_Root_Apo3 + Q_Root_Endo3)*18/1000/1000/1000;
	// temp=T_air;
	// Rm= Rm_25*pow(Q10,(temp-Tb)/10); //Basic Q10 respiration rate in µmol CO2/s/m3
	// Rg=0;
	
	// Only leaves Rm is included so far; from Heskel et al PNAS 2016
	RTref=exp(a_Res+0.1012*Tref-0.0005*Tref*Tref);
	Resp=RTref*exp(0.1012*(T_Leaf[i]-Tref)-0.0005*(T_Leaf[i]*T_Leaf[i]-Tref*Tref));
	return Resp;
}

void Compute_Growth (void)  //under development. Lockhart model
{
	long double Growth_50,Growth_slope,sigma;
	long double Yield_Root,Yield_Branch;	
	long double Extensibility_Root,Extensibility_Branch;
	long double gamma,etha;
	long double A, Rn=8.314, HA=87.5e3, SD, HD=333e3,Tk,Tcrit=7,Topt=30;
	int i;
	
	Yield_Branch=Yield_Root=Yield_Trunk;
	Extensibility_Root=Extensibility_Branch=Extensibility_Trunk;
	for (i=0;i<3;i++)
		{
			if (Turgor_Branch_Symp[i] > Yield_Branch) 	Growth_rate_Branch[i] = Extensibility_Branch * (Turgor_Branch_Symp[i] - Yield_Branch);  // radial growth rate of the Branch Symplasm in s-1
			else                                 		Growth_rate_Branch[i] = 0;	
		}
	if (Turgor_Root_Symp1 > Yield_Root)  Growth_rate_Root1 = Extensibility_Root * (Turgor_Root_Symp1 - Yield_Root);  // radial growth rate of the Root1 Symplasm in s-1
	else                                 Growth_rate_Root1 = 0;
	if (Turgor_Root_Symp2 > Yield_Root)  Growth_rate_Root2 = Extensibility_Root * (Turgor_Root_Symp2 - Yield_Root);  // radial growth rate of the Root2 Symplasm in s-1
	else                                 Growth_rate_Root2 = 0;
	if (Turgor_Root_Symp3 > Yield_Root)  Growth_rate_Root3 = Extensibility_Root * (Turgor_Root_Symp3 - Yield_Root);  // radial growth rate of the Root3 Symplasm in s-1
	else                                 Growth_rate_Root3 = 0;
	if (Turgor_Trunk_Symp > Yield_Trunk) Growth_rate_Trunk = Extensibility_Trunk * (Turgor_Trunk_Symp - Yield_Trunk);  // radial growth rate of the Trunk Symplasm in s-1
	else                                 Growth_rate_Trunk = 0;
	if (Turgor_Axil_Symp  > Yield_Fruit)  Growth_rate_Fruit = Q_Axil_Symp0*Extensibility_Fruit * (Turgor_Axil_Symp - Yield_Fruit);  // volume growth rate of the Fruit Symplasm in mmol/s
	else                                 Growth_rate_Fruit = 0;	
	
	if (Growth_control==1)	//maximum growth rate according to sigmoidal growth
	{
		Growth_50=Growth_para1;
		Growth_slope=Growth_para2;
		etha= exp(Growth_slope/25*(DOY-Growth_50));
		gamma=4*etha/(1+etha)/(1+etha); 				
		Growth_rate_Trunk=		min(Growth_rate_Trunk,		gamma*Extensibility_Trunk 	* (-Pi0_Trunk_Symp*Osmotic_TAir 	- Yield_Trunk));
		Growth_rate_Branch[0]=	min(Growth_rate_Branch[0], 	gamma*Extensibility_Branch 	* (-Pi0_Branch_Symp*Osmotic_TAir 	- Yield_Branch));
		Growth_rate_Branch[1]=	min(Growth_rate_Branch[1], 	gamma*Extensibility_Branch 	* (-Pi0_Branch_Symp*Osmotic_TAir 	- Yield_Branch));
		Growth_rate_Branch[3]=	min(Growth_rate_Branch[2], 	gamma*Extensibility_Branch 	* (-Pi0_Branch_Symp*Osmotic_TAir 	- Yield_Branch));
		Growth_rate_Root1=		min(Growth_rate_Root1,		gamma*Extensibility_Root 		* (-Pi0_Root_Symp*Osmotic_TSoil 	- Yield_Root));
		Growth_rate_Root2=		min(Growth_rate_Root2,		gamma*Extensibility_Root	 	* (-Pi0_Root_Symp*Osmotic_TSoil 	- Yield_Root));
		Growth_rate_Root3=		min(Growth_rate_Root3,		gamma*Extensibility_Root 		* (-Pi0_Root_Symp*Osmotic_TSoil 	- Yield_Root));	
	}
	
	if (Growth_control==2)	// temperature effect according to Cabon et al NP 2020
	{
		Tcrit=Growth_para1;
		Topt=Growth_para2;
		SD=0.01190448*Topt*Topt-4.343334*Topt+1209.74;
		A=pow(10,0.000167997*Topt*Topt-0.06134963*Topt+14.42771);
		Tk=T_Branch+273.16;
		if (T_Branch<Tcrit) sigma=0; 
		else sigma=Tk*A*exp(-HA/Rn/Tk)/(1+exp(SD/Rn*(1-HD/SD/Tk)));
		for (i=0;i<3;i++) Growth_rate_Branch[i] *= sigma;
	
		Tk=T_Root+273.16;
		if (T_Root<Tcrit) sigma=0;
		else sigma=Tk*A*exp(-HA/Rn/Tk)/(1+exp(SD/Rn*(1-HD/SD/Tk)));
		Growth_rate_Root1 *= sigma;
		Growth_rate_Root2 *= sigma;
		Growth_rate_Root3 *= sigma;
		
		Tk=T_Trunk+273.16;
		if (T_Trunk<Tcrit) sigma=0;
		else sigma=Tk*A*exp(-HA/Rn/Tk)/(1+exp(SD/Rn*(1-HD/SD/Tk)));
		Growth_rate_Trunk *= sigma;
		
		Tk=T_Axil +273.16;		
		if (T_Axil<Tcrit) sigma=0;
		else sigma=Tk*A*exp(-HA/Rn/Tk)/(1+exp(SD/Rn*(1-HD/SD/Tk)));
		Growth_rate_Fruit *= sigma;	
	}
}


void C_budget (long double dt_long) //compute the carbon budget of the whole tree
{
	long double  temp, Tref=25, Q_live;
	long double  Q10=2; 
	long double  C_growth; //growth of carbon for the whole tree

	//Get_DATA(dt_long);
	temp=T_air;
	Q_live=Sapwood_V; 											// whole plant living volume in m3
	R_main= Q_live*Rm_25/1000000*pow(Q10,(temp-Tref)/10); 	// Maintenance respiration  in mol CO2/s 
	C_growth= 	(Growth_rate_Trunk*Q_Trunk_Symp0 +				// growth  in mol CO2/s
				Growth_rate_Branch[0]*Q_Branch_Symp0[0] +
				Growth_rate_Branch[1]*Q_Branch_Symp0[1] +
				Growth_rate_Branch[2]*Q_Branch_Symp0[2] +
				Growth_rate_Root1*Q_Root_Symp01 +
				Growth_rate_Root2*Q_Root_Symp02 +
				Growth_rate_Root3*Q_Root_Symp03)*bark*18/1000/1000/1000*Density*1000/2/12;
	if (Reserve<Reserve0*9/10) C_growth= 0;
	R_growth=C_growth*Rg_25*pow(Q10,(temp-Tref)/10);  // resipation is 25% of growth
	C_growth_tot+=C_growth*dt_long*dt;
	R_main_tot+=R_main*dt_long*dt;
	R_growth_tot+=R_growth*dt_long*dt;
	Reserve+=(A_net[0]*Leaf_Area[0]/1000000-R_growth-R_main-C_growth)*dt_long*dt; //mol
}

long double Temperature_correction (long double temp,long double Ha, long double Hd,long double Hs)
{
	long double Corr;	
	Corr= exp(Ha/(8.32*298.15)*(1-298.15/(temp+273.15)))*(1+exp((Hs*298.15-Hd)/(8.32*298.15)))/(1+exp((Hs*(temp+273.15)-Hd)/(8.32*(temp+273.15))));
	return Corr;
}

long double Net_Photosynthesis_C3 (int i)  // From Osborn & Sack 2012 and Bonan 2019 Chap 11; Von Cammerer 2000, 2021
{
	long double temp, Kc,Ko, KmT, Jv, Cm1, Cm2, gt,gm,Cp,Cp25=46,q;
	long double a,b,c,delta;	
	long double PAR_apparent=PAR[i]*cos(Leaf_angle[i]*3.1416/180);
	long double A_LEAF;
	temp=T_Leaf[i];
	gm=gs_max[i]*1.65;
	
	if (g_bl && g_s[i] && gm) gt=1/(1.4/(g_bl/1000) + 1.6/(g_s[i]/1000) + 1/(gm/1000)) ; // this includes gs_co2 +g_mesophyl and Cm not Ci is computed, converted to mol unit
	else gt=0;
	// gt=1/(1.4/(g_bl/1000) + 1.6/(g_s/1000)); // if only the stomatal conductance is taken into account and Ci is used
	if (gt<0) gt=0;
	
	//  Temperature effects on parameters  	
	VcMaxT = VcMax*Temperature_correction(temp,65330,149250,485); // Ha J mol-1,Hs J mol-1,Hd J K-1 mol-1
	VjMaxT = VjMax*Temperature_correction(temp,43540,152040,495); // Ha J mol-1,Hs J mol-1,Hd J K-1 mol-1
	Rd=Rd25*Temperature_correction(temp,46390,150000,490);  // Ha,Hs,Hd
	
	//  Apparent michaelis constant [umol mol-1] in Rubisco-limited situation from Bernacchi et al 2001
	Kc = Kc25 * exp(32.0365 * (temp - 25) / (temp + 273.15));
	Ko = Ko25*1000 * exp(14.6731 * (temp - 25) / (temp + 273.15));
	Cp = Cp25 * exp(15.258  * (temp - 25) / (temp + 273.15));

	KmT = Kc * (1 + 210000 / Ko);
	
	//  Velocity of RuBP regeneration [umol.s-1.m-2]
	q=PAR_apparent * Qye + VjMaxT;
	Jv= (q - pow((q*q -4*Qye*PAR_apparent*0.9*VjMaxT),0.5))/(2*0.9);
	
	// Quadratic solution for Rubisco-limited rate of photosynthesis ; compute Cm
	a= gt;
	b= VcMaxT-Rd+gt*(KmT-Ca);
	c= -VcMaxT*Cp -gt*KmT*Ca-Rd*KmT;
	delta= b*b - 4*a*c;
	if (a) 
		{
		if (delta>0) Cm1= (-b+pow(delta,0.5))/(2*a); 
		else Cm1=0;
		}
	else Cm1= -c/b;
	
	A_net_V=VcMaxT*(Cm1-Cp)/(Cm1+KmT)-Rd;
	
	//  Quadratic solution for RuBP Light-limited rate of photosynthesis
	a= gt;
	b= Jv/4-Rd+gt*(2*Cp-Ca) ;
	c= -(Jv/4*Cp + (Ca*gt+Rd)*2*Cp);
	delta= b*b - 4*a*c;
	if (a) 
		{
		if (delta>=0) Cm2 = (-b+pow(delta,0.5))/(2*a);
		else          Cm2=0;
		}
	else
		{ 
		if (b) Cm2= -c/b;	
		else   Cm2=0;		
		}
	if (Cm2+2*Cp) A_net_J=Jv*(Cm2-Cp)/4/(Cm2+2*Cp)-Rd;
	else A_net_J=0;
	
	//  Net assimilation rate [umol.s-1.m-2]
	A_LEAF=fminl(A_net_V,A_net_J);
	if (A_net_V<A_net_J)  Cm=Cm1;
	else Cm=Cm2; 
		
	
//	A_LEAF = (A_net_V + A_net_J - sqrt(pow((A_net_V+A_net_J), 2) - 4 * shape * A_net_V * A_net_J)) / (2 * shape); //for a smooth transition
	if (g_bl && g_s[i]) Ci=Ca-(A_LEAF)/(1/(1.4/(g_bl/1000) + 1.6/(g_s[i]/1000)));	
	else Ci=0;
	


	return (A_LEAF); // Net photosynthesis
}

long double Net_Photosynthesis_C4 (int i)  //From Von Cammerer PCE 2021 and Yin et al PCE 2011
{
	long double temp,Kc,Ko,Kp, Kp25=80,Vp,Vpr=80,VpMax=120, gt,gbs=13e-3,gm;
	long double a,b,c,delta;	
	long double I=PAR[i]*cos(Leaf_angle[i]*3.1416/180),I2;
	long double A_LEAF;
	long double ATE, ATT, AEE, AET;
	long double gamma_star=0.5/2590, Oi=210000, alpha=0.1,Jatp,fcyc=0.3,omega=26,Topt=43,VjMax_top;
	long double x1,x2,x3,p,q,r,m,n,o,d,f,k,Q,YY,U;
	
//values from von Caemmerer 2021
/*
gbs=0.003;
Kc25=1210;
Ko25=292;
Kp25=82;
VcMax=40.00;
VpMax=200.00;
VjMax=247.7;
Rd=1;
gm=1000; //mmol
*/
	
	//gm=gs_max*2.65*exp(20.1192 * (temp - 25) / (temp + 273.15)); // if gm depends on T
	gm=gs_max[i]*2.65;
	
	temp=T_Leaf[i];
	if (g_bl && g_s[i] && gs_max[i]) gt=1/(1.4/(g_bl/1000) + 1.6/(g_s[i]/1000) + 1/(gm/1000));// this includes g_boundary layer + gs_co2 + g_mesophyl and Cm not Ci is computed in mol
	else gt=0;
	
	Kc = Kc25 * exp(25.9368 * (temp - 25) / (temp + 273.15));
	Ko = Ko25*1000 * exp(4.242   * (temp - 25) / (temp + 273.15));
	Kp = Kp25 * exp(15.4732 * (temp - 25) / (temp + 273.15));

if (PhotoS_model==4)
{	
	// from RAIA-SILVIA et al PCE 2007 for Maize
	VcMaxT = VcMax*Temperature_correction(temp,67294,144568,472); // T, Ha J mol-1, Hd J mol-1, Hs J K-1 mol-1
	VjMaxT = VjMax*Temperature_correction(temp,77900,191929,627); 
	VpMaxT = VpMax*Temperature_correction(temp,70373,117910,376);
	Rd=Rd25*Temperature_correction(temp,46390,150000,490);  
}
	
if (PhotoS_model==5)
{	
	// from Van Cammerer 2021
	VcMaxT = VcMax* exp(31.512 * (temp - 25) / (temp + 273.15));	
	VpMaxT = VpMax* exp(20.2404 * (temp - 25) / (temp + 273.15));
	VjMax_top=VjMax/exp(-(((25-Topt)/omega))*(((25-Topt)/omega)));
	VjMaxT=VjMax_top*exp(-(((temp-Topt)/omega))*(((temp-Topt)/omega)));
	Rd=Rd25*exp(26.8256 * (temp - 25) / (temp + 273.15));
}	
	Rm=Rd/2;		
	//compute Cm and Vp, the PEP carboxylation rate
	a= gt;
	b= VpMaxT+gt*(Kp-Ca);
	c= -gt*Ca*Kp;
	delta= b*b - 4*a*c;
	if (delta>=0) Cm = (-b+pow(delta,0.5))/(2*a); 	else Cm=0;
	Vp=(Cm*VpMaxT)/(Cm+Kp);
	if (Vp>Vpr) Vp=Vpr;
			//Jatp 
	I2=I*0.85*(1-0.15)*(1-fcyc)/(2-fcyc);
	Jatp=(3-fcyc)/(1-fcyc)/4*(I2+VjMaxT-(pow((I2+VjMaxT)*(I2+VjMaxT)-4*0.7*I2*VjMaxT,0.5)))/(2*0.7);

	// ATE
	x1= VcMaxT;
	x2= Kc/Ko;
	x3= Kc;
	a= x2*gt*alpha/0.047 -gt -gbs;
	b= gt*(Ca*gbs +Vp -Rm) + (x3+x2*Oi)*gt*gbs + (x1*gamma_star + x2*Rd)*gt*alpha/0.047 + (gt+ gbs)*(x1-Rd);
	c= -gt*(Ca*gbs +Vp -Rm)*(x1-Rd)+gt*gbs*(x1*gamma_star*Oi +Rd*(x3+x2*Oi));
	delta= b*b - 4*a*c;
	ATE= (-b +pow(delta,0.5))/(2*a);

	// ATT	
	x1=(1-0.4)*Jatp/3;
	x2=7*gamma_star/3;
	x3=0;	
	a= x2*gt*alpha/0.047 -gt -gbs;
	b= gt*(Ca*gbs +Vp -Rm) + (x3+x2*Oi)*gt*gbs + (x1*gamma_star + x2*Rd)*gt*alpha/0.047 + (gt+ gbs)*(x1-Rd);
	c= -gt*(Ca*gbs +Vp -Rm)*(x1-Rd)+gt*gbs*(x1*gamma_star*Oi +Rd*(x3+x2*Oi));
	delta= b*b - 4*a*c;
	ATT= (-b +pow(delta,0.5))/(2*a);

	//AEE
	x1= VcMaxT;
	x2= Kc/Ko;
	x3= Kc;
	d= gt*(Rm-VpMaxT - Ca*(gt+2*gbs)- Kp*(gt+gbs));
	f= gt*gt*(Ca*VpMaxT+(Ca+Kp)*(gbs*Ca-Rm));
	k= gt*gt*gbs*(Ca+Kp);
	m= d- (x3+x2*Oi)*gt*gbs+(Rd-x1)*(gt+gbs)-(x1*gamma_star*gt + +x2*Rd*gt - x2*k/gbs)*alpha/0.047;
	n= f+(x3+x2*Oi)*k + d*(Rd-x1)-gt*gbs*(x1*gamma_star*Oi + Rd*(x3+x2*Oi))+(x1*gamma_star+x2*Rd)*k*alpha/(0.047*gbs);
	o= Rd*(f+(x3+x2*Oi)*k)-x1*(f-k*gamma_star*Oi);
	p= m/(gt+gbs -x2*gt*alpha/0.047);
	q= n/(gt+gbs -x2*gt*alpha/0.047);
	r= o/(gt+gbs -x2*gt*alpha/0.047);
	Q= (p*p-3*q)/9;
	U= (2*p*p*p -9*p*q +27*r)/54;
	YY= acos(U/pow(Q*Q*Q,0.5));	
	AEE= -2*pow(Q,0.5)*cos(YY/3)-p/3;

	//AET
	x1=(1-0.4)*Jatp/3;
	x2=7*gamma_star/3;
	x3=0;	
	d= gt*(Rm-VpMaxT - Ca*(gt+2*gbs)- Kp*(gt+gbs));
	f= gt*gt*(Ca*VpMaxT+(Ca+Kp)*(gbs*Ca-Rm));
	k= gt*gt*gbs*(Ca+Kp);
	m= d- (x3+x2*Oi)*gt*gbs+(Rd-x1)*(gt+gbs)-(x1*gamma_star*gt + x2*k/gbs)*alpha/(0.047*gbs);
	n= f+(x3+x2*Oi)*k + d*(Rd-x1)-gt*gbs*(x1*gamma_star*Oi + Rd*(x3+x2*Oi))+(x1*gamma_star+x2*Rd)*k*alpha/(0.047*gbs);
	o= Rd*(f+(x3+x2*Oi)*k)-x1*(f-k*gamma_star*Oi);
	p= m/(gt+gbs -x2*gt*alpha/0.047);
	q= n/(gt+gbs -x2*gt*alpha/0.047);
	r= o/(gt+gbs -x2*gt*alpha/0.047);
	Q= (p*p-3*q)/9;
	U= (2*p*p*p -9*p*q +27*r)/54;
	YY= acos(U/pow(Q*Q*Q,0.5));	
	AET= -2*pow(Q,0.5)*cos(YY/3)-p/3;

	A_net_V= fminl(ATE,AEE);
	A_net_J=fminl(ATT,AET);
	A_LEAF=fminl(A_net_V,A_net_J);
	if (g_bl && g_s[i]) Ci=Ca-(A_LEAF)/(1/(1.4/(g_bl/1000) + 1.6/(g_s[i]/1000)));		else Ci=0;
	if (g_bl && g_s[i] && gm) Cm=Ca-(A_LEAF)/(1/(1.4/(g_bl/1000) + 1.6/(g_s[i]/1000)+ 1/(gm/1000)));	else Cm=0;
	Cbs= Cm+(Vp-A_LEAF-Rm)/gbs;
	
	return A_LEAF;
 }



long double Arrhenius(long double T2, long double Ea, long double Hd, long double DS)
{
	long double woo;
	long double Rgas = 8.314;  // forgotten in inital R code... appears to be the ideal gas law coefficient
	woo = exp(Ea * (T2 - 298) / (298 * Rgas * T2)) * (1 + exp((298 * DS - Hd) / (298 * Rgas))) / (1 + exp((T2 * DS - Hd) / (T2 * Rgas)));
	return woo;
}


long double Declination (void)
{
	long double c1 = 0.398749068925246; // =Sin(23.5*pi/180), 23.5 = Earth declination
	long double c2  = 2 * 3.1416 / 365;
	long double c3  = 80; // date of spring
	long double x;
	
	x = c1 * sin((DOY - c3) * c2);
	return atan(x / pow((1 - x*x),0.5));
}

long double Potential_PAR(long double timeOfDay)
{
	long double decl,pn,pz,hRad, se, sn, sz, alt ,azi , pfd ,dpfd;
	long double diffuseFraction = 0.1;
	long double solarConstant = 2084;
	long double attenuationCoef  = -0.174353387144778;
	   
	decl = Declination();
	pn = -cos(Lat * 3.1416 / 180);
	pz = sin(Lat * 3.1416 / 180);
	hRad = (timeOfDay - 6) * 3.1416/12;
	se = cos(hRad) * cos(decl);
	sn = -pz * sin(hRad) * cos(decl) - pn * sin(decl);
	sz = -pn * sin(hRad) * cos(decl) + pz * sin(decl);
	alt = atan(sz / (pow((se*se + sn*sn),0.5)));
	azi = 3.1416 + atan(se / sn);
	if (sn > 0)  azi = azi + 3.1416;
	if (alt > 0)  pfd = solarConstant * exp(attenuationCoef / sin(alt));  else pfd=0;
	if (alt > 0)  dpfd = diffuseFraction * pfd; else dpfd=0;
	return dpfd + pfd * sin(alt);
}

long double BB(long double TTT, long double em)
{
	//Longwave black body irradiation
	//... Input:
	//... T in Celsius
	//... emissivity   none
	//... Output:
	//... LW irradiation for upward and downward faces in W/m2
	return 5.6704e-8*em*pow(TTT+273.15,4);
}

long double gs(long double GsIn)
{
	//Set the stomatal conductance according to microclimate and surface temperature
	//... Input:
	//... Gs Imposed (m/s)
	//... Output:
	//... stomatal conductance gs (m/s)
	
	return  GsIn;
}

long double Transpi(long double TTT,long double gb,long double RH, long double GsIn)
{
	// Computation the evaporative heat flux
	//... Input:
	//... Surface temperature (Temperature in K)
	//... Boundary layer conductance (gb in m/s)
	//... Relative humidity (RH in %)
	//... Stomatal conductance (gs in  m/s)
	//... Output: latent heat flux (w/m2)
	
	long double a = 17.443;
	long double b = 2795.;
	long double c = 3.868;
	long double Po = 100000.;
	long double rhoAir=1.292 ;  // density of dry air    kg/m3
	long double CpAir=1010;   // heat capacity of dry air    J kg-1 K-1
	long double Temp = TTT + 273.15;
	long double esat = Po*exp(log(10.)*(a-b/Temp-c*log10(Temp)));
	long double ea=esat*(RH/100);
	long double VPD1 = esat-ea;
	long double gw = 1/(1/gs(GsIn)+1/gb);
	
	return rhoAir*CpAir*gw*VPD1/66.5;
	
}

long double Convec(long double TTT,long double hh, long double Tair)
{
	// Computation the convective heat flux
	//... Input:
	//... Air Temperature (K): Tair
	//... Heat transfer coefficient (W/K/m2): hh
	//... Surface temperature (K): T
	//... Output: Sensible heat flux (w/m2)
	
	return hh*(TTT-Tair);
}

long double Flux (long double TTT,long double gb,long double hh, long double SWRa,long double GsIn,long double RH,long double Tair,long double em_Leaf,long double em_air)
{
	// Compute the heat flux balance
	//... Input:
	//... Surface Temperature (Temp in K)
	//... Heat transfer conductance (s/m) for Sensible Heat Flux: gb
	//... Heat transfer coefficient (W/K/m2): hh
	//... Absorbed Incoming Shortwave Radiation: SWRa
	//... Parameter to compute the Stomatal Conductance or directly the Stomatal conductance: GsIn
	//... Relative Hulidity of Air: RH
	//... Air Temperature: Tair
	//... Leaf emissivity: em_Leaf
	//... Air emissivity: em_air
	//... Output: Heat balance (W/m2)
	
	return SWRa+BB(Tair,em_air)-BB(TTT,em_Leaf)-Transpi(TTT,gb,RH,GsIn)-Convec(TTT,hh,Tair);
}

long double brentq( long double xa, long double xb, long double xtol, long double rtol, long double gb, long double hh,long double SWRa,long double GsIn,long double RH,long double Tair,long double em_Leaf,long double em_air, int maxiter)
{
	// brentq method taken from scipy.optimize library
	
	double xpre = xa, xcur = xb;
	double xblk = 0.0, fpre, fcur, fblk = 0.0, spre = 0.0, scur = 0.0, sbis, tol;
	double stry, dpre, dblk;
	int i;
	fpre = Flux(xpre,gb,hh,SWRa,GsIn,RH,Tair,em_Leaf,em_air);
	fcur = Flux(xcur,gb,hh,SWRa,GsIn,RH,Tair,em_Leaf,em_air);
	if (fpre == 0) return xpre;
	if (fcur == 0) return xcur;
	for(i = 0; i < maxiter; i++)
	{
		if (fpre*fcur < 0)
		{
			xblk = xpre;
			fblk = fpre;
			spre = scur = xcur - xpre;
		}
		if (fabs(fblk) < fabs(fcur))
		{
			xpre = xcur; xcur = xblk; xblk = xpre;
			fpre = fcur; fcur = fblk; fblk = fpre;
		}
		
		tol = xtol + rtol*fabs(xcur);
		sbis = (xblk - xcur)/2;
		if (fcur == 0 || fabs(sbis) < tol)
			return xcur;
		
		if (fabs(spre) > tol && fabs(fcur) < fabs(fpre))
		{
			if (xpre == xblk)   stry = -fcur*(xcur - xpre)/(fcur - fpre);   /* interpolate */
			else                                                            /* extrapolate */
			{
				dpre = (fpre - fcur)/(xpre - xcur);
				dblk = (fblk - fcur)/(xblk - xcur);
				stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre));
			}
			if (2*fabs(stry) < min(fabs(spre), 3*fabs(sbis) - tol))         /* good short step */
			{
				spre = scur;
				scur = stry;
			}
			else         /* bisect */
			{
				spre = sbis;
				scur = sbis;
			}
		}
		else            /* bisect */
		{
			spre = sbis;
			scur = sbis;
		}
		
		xpre = xcur;
		fpre = fcur;
		if (fabs(scur) > tol)  xcur += scur;
		else                   xcur += (sbis > 0 ? tol : -tol);
		fcur = Flux(xcur,gb,hh,SWRa,GsIn,RH,Tair,em_Leaf,em_air);
	}
	return xcur;
}

void TLeaf(int i)          //Leaf Energy budget from Ecofiz_TLeaf_K2_v3 www.landflux.org/r
{
	long double SWR;   					// short-wave radiation    W m-2
	long double WS;    					// windspeed    m s-1
	long double Tair;   					// air temperature    oC
	long double RH;    					// relative humidity    %
	long double aSWR=0.5;  				// absorptance to SWR     %
	long double em_Leaf=0.97;    			// emissivity    none
	long double d;    						// characteristic dimension    mm
	long double rst;   					// stomatal resistance    s m-1
	long double SB=5.6704e-8;  			// Stefan-Boltzman constant    W m-2 K-4
	long double p=1.292 ;  				// density of dry air    kg/m3
	long double Cp=1010;  					// heat capacity of dry air    J kg-1 K-1
	long double y=0.066  ;  				// psychrometric constant    kPa K-1
	long double a=0.61121  ;  				// coefficient in esat equation    kPa
	long double b=17.502  ;  				// coefficient in esat equation    none
	long double z=240.97  ;  				// coefficient in esat equation    °C
	long double gflat=0.00662;				// coefficient in rbl equation    m
	long double gcyl=0.00403;  			// coefficient in rbl equation    m
	long double jflat=0.5,jcyl=0.6;  		// coefficient in rbl equation    none
	long double esat ;  					// saturation vapor pressure    kPa
	long double ea   ; 					// water vapor pressure of the air    kPa
	long double s   ;						// slope of esat/T curve    kPa oC-1
	long double VPDx  ; 					// water vapor pressure deficit of the air    kPa
	long double SWRabs;   					// absorbed short-wave radiation    W m-2
	long double LWRin  ; 					// incoming long-wave radiation    W m-2
	long double LWRouti;   				// isothermal outgoing long-wave radiation    W m-2
	long double Rni;  						// isothermal net radiation    W m-2
	long double rr ;  						// radiative resistance    s m-1
	long double rblr;   					// boundary-layer + radiative resistance    s m-1
	long double ym ; 						// modified psychrometric constant    kPa K-1
	long double rbl ;  					// Leaf boundary-layer resistance    s m-1
	long double Delta_T  ; 				// Leaf-to-air temperature difference    oC
	long double TLeaf, TLeaf_NonLinear;   // Leaf temperature    oC
	int maxiter = 50;
	long double em_air;
	
	
	if (POTENTIAL_PAR) 
	{
		if (TLEAF==1 || TLEAF==2) cloud_cover=PAR[0]/POTENTIAL_PAR; 
		else cloud_cover=0;
	}
	if (cloud_cover>1)	cloud_cover=1;
	if (CLIMAT==5 || CLIMAT==11)		cloud_cover=1;
	//if (CLIMAT==1)		cloud_cover=1;
	d=Leaf_size;
	Tair=T_air;
	RH= RH_air;
	
	if (g_s[i]+g_cuti[i]) rst=1/(g_s[i]+g_cuti[i])*1000*40;   // conversion fromm mmol/s/m2 to s/m 
	else     rst=9999.99;
	SWR=    PAR[i]*0.5495;    // from µmol/m²/s to Watts/m²
	WS=     Wind[i];
	esat=   a*exp(b*Tair/(Tair+z)); //kPa
	ea=     esat*(RH/100);
	s=      esat*b*z/(pow((Tair+z),2));
	em_air =((1-0.84*cloud_cover)*1.31*pow(10*ea/(Tair+273.15),0.14285714)+0.84*cloud_cover);
	VPDx=   esat-ea;
	SWRabs= aSWR*cos(Leaf_angle[i]*3.1416/180)*SWR;
	LWRin=  em_air*SB*pow(Tair+273.15,4);  // for clear and cloudy sky
	LWRouti=em_Leaf*SB*pow(Tair+273.15,4);
	Rni=    SWRabs+LWRin-LWRouti;
	rr=     p*Cp/(4*em_Leaf*SB*pow(Tair+273.15,3));
	
	if (Leaf_size > 3 ) rbl=1/(1.5*gflat*(pow(WS,jflat)/pow(d/1000,(1-jflat))));     //  a flat Leaf if > 3mm
	else                rbl=1/(1.5*gcyl*(pow(WS,jcyl)/pow(d/1000,(1-jcyl))));       // a needle, formula for a cylinder
	g_bl=1/rbl*1000*40;     //Leaf boundary layer conductance in mmol/s/m2
	rblr=1/(1/rbl+1/rr);
	ym=y*(rst/rblr);
	// compute TLeaf with linear approximation
	if (TLEAF==1 || TLEAF==11)
	{
		Delta_T= (ym*Rni*rblr/(p*Cp)-VPDx)/(s+ym);
		TLeaf=   Tair+Delta_T;
		T_Leaf[i]=  TLeaf;
	}
	// compute non linear TLeaf
	else if (TLEAF==2 || TLEAF==12)
	{
		TLeaf_NonLinear =  brentq(Tair-100.0, Tair+100.0, 1e-16, 1e-16, 1/rbl,p*Cp/rbl,SWRabs,1/rst,RH,Tair,em_Leaf,em_air,maxiter);
		T_Leaf[i]= TLeaf_NonLinear;
	}
	else T_Leaf[i]=T_air;

}

void TAxil(void)          //Flower Energy budget from Ecofiz_TLeaf_K2_v3 www.landflux.org/r  Suppose all flower are on top branch 
{
	long double SWR;   					// short-wave radiation    W m-2
	long double WS;    					// windspeed    m s-1
	long double Tair;   					// air temperature    oC
	long double RH;    					// relative humidity    %
	long double aSWR=0.5;  				// absorptance to SWR     %
	long double em_Leaf=0.97;    			// emissivity    none
	long double d;    					// characteristic dimension    mm
	long double rst;   					// stomatal resistance    s m-1
	long double SB=5.6704e-8;  			// Stefan-Boltzman constant    W m-2 K-4
	long double p=1.292 ;  				// density of dry air    kg/m3
	long double Cp=1010;  					// heat capacity of dry air    J kg-1 K-1
	long double y=0.066  ;  				// psychrometric constant    kPa K-1
	long double a=0.61121  ;  				// coefficient in esat equation    kPa
	long double b=17.502  ;  				// coefficient in esat equation    none
	long double z=240.97  ;  				// coefficient in esat equation    °C
	long double gflat=0.00662;				// coefficient in rbl equation    m
	long double gcyl=0.00403;  			// coefficient in rbl equation    m
	long double jflat=0.5,jcyl=0.6;  		// coefficient in rbl equation    none
	long double esat ;  					// saturation vapor pressure    kPa
	long double ea   ; 					// water vapor pressure of the air    kPa
	long double s   ;						// slope of esat/T curve    kPa oC-1
	long double VPDx  ; 					// water vapor pressure deficit of the air    kPa
	long double SWRabs;   					// absorbed short-wave radiation    W m-2
	long double LWRin  ; 					// incoming long-wave radiation    W m-2
	long double LWRouti;   				// isothermal outgoing long-wave radiation    W m-2
	long double Rni;  						// isothermal net radiation    W m-2
	long double rr ;  						// radiative resistance    s m-1
	long double rblr;   					// boundary-layer + radiative resistance    s m-1
	long double ym ; 						// modified psychrometric constant    kPa K-1
	long double rbl ;  					// Leaf boundary-layer resistance    s m-1
	long double Delta_T  ; 				// Leaf-to-air temperature difference    oC
	long double TLeaf_NonLinear;   	// Leaf temperature    oC
	int maxiter = 50;
	long double em_air;
	
	
	if (POTENTIAL_PAR) 
	{
		if (TLEAF==1 || TLEAF==2) cloud_cover=PAR[0]/POTENTIAL_PAR; 
		else cloud_cover=0;
	}
	if (cloud_cover>1)	cloud_cover=1;
	if (CLIMAT==5 || CLIMAT==11)		cloud_cover=1;
	//if (CLIMAT==1)		cloud_cover=1;
	d=Diam_Axil*1000; // from m to mm
	Tair=T_air;
	RH= RH_air;
	
	if (g_Axil) rst=1/(g_Axil)*1000*40;   // conversion fromm mmol/s/m2 to s/m 
	else     rst=9999.99;
	SWR=    PAR[0]*0.5495;    // from µmol/m²/s to Watts/m²
	WS=     Wind[0];
	esat=   a*exp(b*Tair/(Tair+z)); //kPa
	ea=     esat*(RH/100);
	s=      esat*b*z/(pow((Tair+z),2));
	em_air =((1-0.84*cloud_cover)*1.31*pow(10*ea/(Tair+273.15),0.14285714)+0.84*cloud_cover);
	VPDx=   esat-ea;
	SWRabs= aSWR*1*SWR; // consider flower are flat
	LWRin=  em_air*SB*pow(Tair+273.15,4);  // for clear and cloudy sky
	LWRouti=em_Leaf*SB*pow(Tair+273.15,4);
	Rni=    SWRabs+LWRin-LWRouti;
	rr=     p*Cp/(4*em_Leaf*SB*pow(Tair+273.15,3));
	
	if (d > 3 ) 	rbl=1/(1.5*gflat*(pow(WS,jflat)/pow(d/1000,(1-jflat))));     //  a flat Leaf if > 3mm
	else       	rbl=1/(1.5*gcyl*(pow(WS,jcyl)/pow(d/1000,(1-jcyl))));       // a needle, formula for a cylinder
	g_bl_Axil=1/rbl*1000*40;     //Flower boundary layer conductance in mmol/s/m2
	rblr=1/(1/rbl+1/rr);
	ym=y*(rst/rblr);
	// compute TLeaf with linear approximation
	if (TLEAF==1 || TLEAF==11)
	{
		Delta_T= (ym*Rni*rblr/(p*Cp)-VPDx)/(s+ym);
		T_Axil=   Tair+Delta_T;
	}
	// compute non linear TLeaf
	else if (TLEAF==2 || TLEAF==12)
	{
		TLeaf_NonLinear =  brentq(Tair-100.0, Tair+100.0, 1e-16, 1e-16, 1/rbl,p*Cp/rbl,SWRabs,1/rst,RH,Tair,em_Leaf,em_air,maxiter);
		T_Axil= TLeaf_NonLinear;
	}
	else T_Axil=T_air;
}

void fill_soil_layers(void)
{
	if (Q_Soil1>Q_Soil01)
	{
		Q_Soil2+=(Q_Soil1 - Q_Soil01);
		Q_Soil1=Q_Soil01;
	}
	if (Q_Soil2>Q_Soil02)
	{
		Q_Soil3+=(Q_Soil2 - Q_Soil02);
		Q_Soil2=Q_Soil02;
	}
	if (Q_Soil3>Q_Soil03)
	{
		Drainage+=(Q_Soil3 - Q_Soil03);
		Q_Soil3=Q_Soil03;
	}
}

void Irrigate (void)
{
	if (IRRIGATE==1 && REW_t<RWC_Irr)  //then resature the soil when threshold RWC is reached
	{
		Irrigation += ((Q_Soil01+Q_Soil02+Q_Soil03) - (Q_Soil1+Q_Soil2+Q_Soil3))/(Surface_Soil*1000*1000/18); // irrigation in mm
		Q_Soil1        =Q_Soil01;
		Q_Soil2        =Q_Soil02;
		Q_Soil3        =Q_Soil03;
	}
	else if (IRRIGATE==10 && RWC3<RWC_Irr)  //then resature the soil when threshold RWC3 is reached
	{
		Irrigation += ((Q_Soil03) - (Q_Soil3))/(Surface_Soil*1000*1000/18); // irrigation in mm
		Q_Soil3        =Q_Soil03;
		RWC_Irr=0.99; // to keep the soil saturated
	}
	else if (IRRIGATE==4 && REW_t<RWC_Irr) //then add only the Daily irrigation to top layer when threshold RWC is reached
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18;
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if (IRRIGATE==2 || IRRIGATE==6) //automatic daily irrigation of Daily_irr mm
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if ((IRRIGATE==3 || IRRIGATE==7) && T>IRR_DOY_S*3600*24*1000 && T<IRR_DOY_F*3600*24*1000) //automatic daily irrigation of Daily_irr mm only between DOY_start and DOY_end
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if (IRRIGATE==5 && T>IRR_DOY_S*3600*24*1000 && T<IRR_DOY_F*3600*24*1000 && REW_t<RWC_Irr) //automatic daily irrigation of Daily_irr mm only between DOY_start and DOY_end when below RWC_irr
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if (IRRIGATE==8 && !((int)(T/3600/24/1000-DOY_0+1)%(int)RWC_Irr)) //automatic daily irrigation of Daily_irr mm every Ndays
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if (IRRIGATE==9 && !((int)(T/3600/24/1000-DOY_0+1)%(int)RWC_Irr)) //automatic daily irrigation of Daily_irr mm every Ndays
	{
		Irrigation += ((Q_Soil01+Q_Soil02+Q_Soil03) - (Q_Soil1+Q_Soil2+Q_Soil3))/(Surface_Soil*1000*1000/18); // irrigation in mm
		Q_Soil1        =Q_Soil01;
		Q_Soil2        =Q_Soil02;
		Q_Soil3        =Q_Soil03;
	}
	else if (IRRIGATE==11)
	{
		Q_Soil1+=  RWC_Irr*EvapoT_day; //add rain 
		Irrigation+=RWC_Irr*EvapoT_day/(Surface_Soil*1000*1000/18);
		fill_soil_layers();
	}
}

void Rehydrate(long double dt_long)
{
	if (REHYDRATE==11 || REHYDRATE==21 || REHYDRATE==31) //Leaf PLC
			if (PLC_Leaf_Apo[1]>=PLC_REHYD) //based on the top Branch
				{
					if (REHYDRATE==11 || REHYDRATE==31)
					{
						Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
						Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
						Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
						
					if (REHYDRATE==11)	if (!REFILL) PLC_Leaf_Apo[1]-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration; disable of 31
					}
					if (REHYDRATE==21)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
			
		if (REHYDRATE==12 || REHYDRATE==22 || REHYDRATE==32) //Branch PLC
			if (PLC_Branch_Apo[1]>=PLC_REHYD)
				{
				if (REHYDRATE==12 || REHYDRATE==32)
					{
						
						Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
						Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
						Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
							
					if (REHYDRATE==12)	if (!REFILL) PLC_Branch_Apo[1]-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
				if (REHYDRATE==22)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
			
				}
				
			if (REHYDRATE==13 || REHYDRATE==23|| REHYDRATE==33) //Trunk PLC
				if (PLC_Trunk_Apo>=PLC_REHYD)
				{
					if (REHYDRATE==13 || REHYDRATE==33)
					{
						
								Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
								Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
								Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
							
					if (REHYDRATE==13)	if (!REFILL) PLC_Trunk_Apo-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
					if (REHYDRATE==23)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
			if (REHYDRATE==14 || REHYDRATE==24|| REHYDRATE==34)	//Root1 PLC
				if (PLC_Root_Apo1>=PLC_REHYD)
				{
					if (REHYDRATE==14 || REHYDRATE==34)
					{
						
								Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
								Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
								Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
							
					if (REHYDRATE==14)	if (!REFILL) PLC_Root_Apo1-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
					if (REHYDRATE==24)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
			if (REHYDRATE==15 || REHYDRATE==25|| REHYDRATE==35)	//Whole plant PLC
				if (PLC_Plant>=PLC_REHYD)
				{
					if (REHYDRATE==15 || REHYDRATE==35)
					{
						
								Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
								Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
								Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
							
					if (REHYDRATE==15)	if (!REFILL) PLC_Plant-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
					if (REHYDRATE==25)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
			if (REHYDRATE==16 || REHYDRATE==26 || REHYDRATE==36)	//Whole plant+soil PLC
				if (PLC_Plant_Soil>=PLC_REHYD)
				{
					if (REHYDRATE==16 || REHYDRATE==36)
					{
						
								Q_Soil1        =Teta_fc_1*Volume_soil1*1000*1000*1000/18;
								Q_Soil2        =Teta_fc_2*Volume_soil2*1000*1000*1000/18;
								Q_Soil3        =Teta_fc_3*Volume_soil3*1000*1000*1000/18;
							
					if (REHYDRATE==16)	if (!REFILL) PLC_Plant_Soil-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
					if (REHYDRATE==26)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
}

void Beer_Lambert(double PAR_in) //compute PAR for each Branch with Beer_Lambert's law of attenuation with LAI.
{
    int i;
    double LAI_crown[5];
    LAI_crown[0]=0;
    for (i=1;i<=4;i++)
    {
        PAR[i]= PAR_in*exp(-Extinction_Coeff*LAI_crown[i-1]);
        LAI_crown[i]=LAI_crown[i-1]+Leaf_Area[i]/Crown_Area;  //cumulative LAI in the canopy
    }
	//for (i=1;i<4;i++) Test[i]=PAR[i];
}

void Wind_profile(void) // compute wind profile in the canopy according to Yi (2008); assume canopy height is 2xBranch_length and LAI is evenly distributed
{
	if (Extinction_Coeff) 
		{
		Wind[3]=Wind[0]*exp(-LAI_Crown/3);
		Wind[2]=Wind[0]*exp(-LAI_Crown/6);
		Wind[1]=Wind[0];
		}
		else Wind[3]=Wind[2]=Wind[1]=Wind[0];
}
	
void Climat (void) //copmpute climatic data
{
	long double e_sat, e,  e_sat_air, e_air, e_air_sol1, e_air_sol2, e_air_sol3, slope, tangente,TTTT,day,day_HW;
	int i=0;
	
	if      (CO2_atm==0)  Ca=  400;
	else if (CO2_atm==1)                   // RCP 2.6
	{
		if (YEAR1<2000) 	Ca=  255.1911 + 282.7456/pow(1+exp(-(YEAR1-2058.0442)/11.3900),0.1781); 
		else 				Ca= -0.00000001516648304*pow(YEAR1,5) + 0.0001586265578*pow(YEAR1,4) - 0.663190082*pow(YEAR1,3) + 1385.434635*pow(YEAR1,2) - 1446186.605*YEAR1 + 603458228.5;
	}
	else if (CO2_atm==2)  Ca=  255.1911 + 282.7456/pow(1+exp(-(YEAR1-2058.0442)/11.3900),0.1781);                 // RCP 4.5
	else if (CO2_atm==3)  Ca=  282.8391 + 902.2988/pow(1+exp(-(YEAR1-2104.6814)/17.3348),0.3910);                 // RCP 8.5
	else                  Ca=  CO2_atm;
	Wind_profile();
	
	if (CLIMAT==2 || CLIMAT==3 || CLIMAT==4 || CLIMAT==7) // except for interpolated values
	{
		tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOY-80)/365*2*3.1416)));
		TTTT=(T/3600/1000-DOY*24); // time of day
		if (tangente<-1)     Day_length=0;
		else if (tangente>1) Day_length=24;
		else                 Day_length=    24 -  acos(tangente)*7.6394194;
		PAR[0]=Potential_PAR(TTTT)/Potential_PAR(12)*PAR_max;
		Beer_Lambert(PAR[0]);
	}
	
	if (CLIMAT==0 || CLIMAT==10) //computed from constant data in sureau.init
	{
		if (CLIMAT==0) Day_length=12;
		else
			{
			tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOY-80)/365*2*3.1416)));
			TTTT=(T/3600/1000-DOY*24); // time of day
			if (tangente<-1)     Day_length=0;
			else if (tangente>1) Day_length=24;
			else                 Day_length=    24 -  acos(tangente)*7.6394194;
			}
		
		TTTT=T/3600/1000-24*(float)(int)(T/3600/1000/24);
		PAR[0]=Potential_PAR(TTTT)/Potential_PAR(12)*PAR_max*PAR_att;
		Beer_Lambert(PAR[0]);
	
		T_Soil=20;
		T_air  = (T_air_max+T_air_min)/2+(T_air_max-T_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH1));    //daily cos variation between T_min et T_max et T_max at HH1
	
		RH_air = (RH_air_max+RH_air_min)/2+(RH_air_max-RH_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH2)); //daily cos variation between HR_min et T_max et HR_max at 0h00
		if (HW) // a Heat Wave
		{
			day= T/24/3600/1000;
			//printf("%Lf ",day);
			if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
			{
			 if (HW==1) 
				{
				RH_air*= exp((18.678-T_air/234.5)*T_air/(257.14+T_air))/exp((18.678-(T_air+HW_T)/234.5)*(T_air+HW_T)/(257.14+T_air+HW_T));  //new RH_air assuming constant e_air
				T_air+=HW_T;
				}
			 else if (HW==3) 
				{
				RH_air= 100-(611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air)))/(611.21*exp((18.678-(T_air+HW_T)/234.5)*(T_air+HW_T)/(257.14+(T_air+HW_T))))*(100-RH_air);    
				T_air+=HW_T;
				}
				
			else if (HW==4) // a dry air wave at constant T_air
				{
				e_sat_air= 611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
				e_air =    e_sat_air*RH_air/100;                                       // vapour water pressure at Tair and RHair
				VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
				RH_air= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
				}
			else if (HW==5) //// a Heat wave starting at day HW_day, lasting HW_duration days with a daily temperature increase of HW_T °C/day
				{
					day_HW=day-HW_day; //number of days since HW onset
					T_air+=HW_T*day_HW;
				}
			else if (HW==6) // a wind wave !
				{
					Wind[0]*=HW_T;
					Wind_profile();
				}
			
			 else T_air+=HW_T;
			}
		}
		
	 //   if (IRRIGATE==1 ||IRRIGATE==4 || IRRIGATE==5) Irrigate();  // for cases 2 & 3 irrigate only once at midnight
		
	}
	
	if (CLIMAT==5) //constant values 
	{
		Day_length=12;
		T_Soil=20;
		T_air  = T_air_max;
		RH_air = RH_air_min;
		PAR[0]=PAR_max*PAR_att;
		Beer_Lambert(PAR[0]);
		//if (IRRIGATE==1 ||IRRIGATE==4 || IRRIGATE==5) Irrigate(); // for cases 2 & 3 irrigate only once at midnight
		if (HW) // a Heat Wave
		{
			day= T/24/3600/1000;
			if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
			{
			 if (HW==1) 
				{
				RH_air*= exp((18.678-T_air/234.5)*T_air/(257.14+T_air))/exp((18.678-(T_air+HW_T)/234.5)*(T_air+HW_T)/(257.14+T_air_2+HW_T));  //new RH_air assuming constant e_air
				T_air+=HW_T;
				}
			 else if (HW==3) 
				{
				RH_air= 100-(611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air)))/(611.21*exp((18.678-(T_air+HW_T)/234.5)*(T_air+HW_T)/(257.14+(T_air+HW_T))))*(100-RH_air);    
				T_air+=HW_T;
				}
			 else if (HW==4) // a dry air wave at constant T_air
				{
				e_sat_air= 611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
				e_air =    e_sat_air*RH_air/100;                                       // vapour water pressure at Tair and RHair
				VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
				RH_air= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
				}
			 else if (HW==6) // a wind wave !
				{
				Wind[0]*=HW_T;
				Wind_profile();
				}
			 else T_air+=HW_T;
			}
		}
	}
	
	
	if (CLIMAT==1 || CLIMAT==6 || CLIMAT==8 ||  CLIMAT==11)  // climatic variable interpolated from 2 hourly values
	{
		tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOY-80)/365*2*3.1416)));
		if (tangente<-1)Day_length=0;
		if (tangente>1) Day_length=24;
		else Day_length=    24 -  acos(tangente)*7.6394194;
		T_Soil=20;
		
		if (T_2-T_1)
		{
			T_air  = T_air_1 + (T/1000-T_1) * (T_air_2-T_air_1)/(T_2-T_1);
			RH_air = RH_air_1 + (T/1000-T_1) * (RH_air_2-RH_air_1)/(T_2-T_1);
			PAR[0]    = (PAR_1 + (T/1000-T_1) /(T_2-T_1) * (PAR_2-PAR_1))*PAR_att;
			T_Trunk=T_Trunk_1 + (T/1000-T_1) * (T_Trunk_2-T_Trunk_1)/(T_2-T_1);
			T_Branch1 = T_Branch1_1 + (T/1000-T_1) * (T_Branch1_2-T_Branch1_1)/(T_2-T_1);
			T_Branch2 = T_Branch2_1 + (T/1000-T_1) * (T_Branch2_2-T_Branch2_1)/(T_2-T_1);
			T_Branch3 = T_Branch3_1 + (T/1000-T_1) * (T_Branch3_2-T_Branch3_1)/(T_2-T_1);
			Beer_Lambert(PAR[0]);
			
			
			if (CLIMAT==6)
			{
				T_Soil1= T_Soil_11 + (T/1000-T_1) * (T_Soil_21-T_Soil_11)/(T_2-T_1);
				T_Soil2= T_Soil_12 + (T/1000-T_1) * (T_Soil_22-T_Soil_12)/(T_2-T_1);
				T_Soil3= T_Soil_13 + (T/1000-T_1) * (T_Soil_23-T_Soil_13)/(T_2-T_1);
			}
		}
		else
		{
			T_air  = T_air_1;
			RH_air = RH_air_1;
			PAR[0]    = PAR_1;
			T_Trunk=T_Trunk_1;
			T_Branch1 = T_Branch1_1;
			T_Branch2 = T_Branch2_1;
			T_Branch3 = T_Branch3_1;
			Beer_Lambert(PAR[0]);
			
			
			if (CLIMAT==6)
			{
				T_Soil1= T_Soil_11;
				T_Soil2= T_Soil_12;
				T_Soil3= T_Soil_13;
			}
		
		}
	}
	
	if (CLIMAT==2 || CLIMAT==3 || CLIMAT==4 || CLIMAT==7)  // climatic variable are computed from cos functions with daily values from a file
	{
		//if (IRRIGATE) Irrigate();  //refill soil water content by irrigation
		// T_air
		if (sin(3.14159265359/12*(T/3600/1000-HH1))<0)   // from Tmim_n to Tmax_n starting at 24-HH1 to HH1
			T_air  = (T_air_max+T_air_min)/2+(T_air_max-T_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH1));    //daily cos variation between T_min et T_max et T_max at 14h00
		else   //from T_max_n to T_min_n+1 starting at HH1 to 24-HH1
		{
			if (sin(3.14159265359/12*(T/3600/1000-12))>-0.0001)  
					T_air  = (T_air_max+T_air_min_2)/2+(T_air_max-T_air_min_2)/2*cos(3.14159265359/12*(T/3600/1000-HH1));   // before midnigh from T_max_n to T_min_n+1
			else    T_air  = (T_air_max_0+T_air_min)/2+(T_air_max_0-T_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH1));  // after midnigh from T_max_n-1 to T_min_n
		}
		 
		// T_soil
		if (CLIMAT==4) T_Soil= T_Soil_1 + (T_Soil_2-T_Soil_1)*(T/3600/1000-DOY*24)/24; else T_Soil=20;       // if not measured assumed T_Soil is constant and == 15°C
		
		// RH_air
		if (sin(3.14159265359/12*(T/3600/1000-HH2))>=0)  // from HRmax_n to HRmin_n
			RH_air = (RH_air_max+RH_air_min)/2+(RH_air_max-RH_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH2)); //daily cos variation between HR_min et T_max et HR_max at 4h00
		else //   from HRmin to HRmax
			if (sin(3.14159265359/12*(T/3600/1000-12))>-0.001)  // before midnigh from T_max_n to T_min_n+1
					RH_air = (RH_air_max_2+RH_air_min)/2+(RH_air_max_2-RH_air_min)/2*cos(3.14159265359/12*(T/3600/1000-HH2)); //daily cos variation between HR_min et T_max et HR_max at 4h00
			else    RH_air = (RH_air_max+RH_air_min_0)/2+(RH_air_max-RH_air_min_0)/2*cos(3.14159265359/12*(T/3600/1000-HH2));
	}
	
	if (IRRIGATE==1 || IRRIGATE==4 || IRRIGATE==5 ) Irrigate();
		//if (COMPET && COMPET_number==0) Irrigate();  // for cases 2,3,6,7 irrigate only once at midnight or 19h00
	
	if (RH_air<0.1) RH_air=0.1;
	if (RH_air>100) RH_air=99.999;
	if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
				 e_air=e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
	P_air =     0.4609418*(273.15+T_air)*log(e_air/e_sat_air);              // air water potential, MPa
	VPD_Air=    (e_sat_air - e_air)/1000;                                  // vpd in kPa
	if (VPD_Air<0) VPD_Air=0;
	slope=4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
	ETP_Penman= 0.5625*(0.408*slope*PAR[0]*0.5495*3.6e-3+0.066*37*Wind[0]*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind[0]));  //ETP in mm/h  0.5625 to fit observed ETP
	ETP_Penman/=3600;  // ETP in mm/s
	ETP_Penman=ETP_Penman/18*1000*1000; // ETP in mmol/s  to be multiplied by soil area to obtain the volume
	POTENTIAL_PAR=Potential_PAR((T/3600/1000-DOY*24));
	
	// LEAF
	for (i=1;i<4;i++)
	{
		if (!TLEAF) T_Leaf[i]=T_air;             //Leaf temperature is then assumed to equal the air temperature
		if (T_Leaf[i]>0) e_sat=611.21*exp((18.678-T_Leaf[i]/234.5)*T_Leaf[i]/(257.14+T_Leaf[i]));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else          e_sat=611.15*exp((23.036-T_Leaf[i]/333.7)*T_Leaf[i]/(279.82+T_Leaf[i]));
		e= e_sat*exp(P_Leaf_Evap[i]*2.16947115/(T_Leaf[i]+273.15));
	
		if (Leaf_Area[i])  VPD_Leaf[i]= (e-e_air)/1000; //vpd between Leaf and air in kPa
		//VPD_Leaf= (611.21*exp((18.678-T_Leaf/234.5)*T_Leaf/(257.14+T_Leaf))*exp(P_Leaf_Evap*2.16947115/(T_Leaf+273.15))-611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air))*RH_air/100)/1000
		else            VPD_Leaf[i]=0;
		if (VPD_Leaf[i]<0) VPD_Leaf[i]=0;  // to avoid back flow

		e=        e_sat*exp(P_Leaf_Symp[i]*2.16947115/(T_Leaf[i]+273.15));
		VPD_Cuti[i]= (e-e_air)/1000; //vpd between Leaf and air in kPa
		if (VPD_Cuti[i]<0) VPD_Cuti[i]=0;  // to avoid back flow
	
	}
		// BUD
	if (Type_Axil)
	{
		//T_Axil=T_air;  // could also be T_Leaf, or T_bud!
		if (T_Axil>0) e_sat=611.21*exp((18.678-T_Axil/234.5)*T_Axil/(257.14+T_Axil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else          e_sat=611.15*exp((23.036-T_Axil/333.7)*T_Axil/(279.82+T_Axil));
		e=            e_sat*exp(P_Axil_Symp*2.16947115/(T_Axil+273.15));
		VPD_Axil=(e-e_air)/1000; //vpd between Branch and air in mmol/m3
		if (VPD_Axil<0) VPD_Axil=0;
		if (Type_Axil>=2) //petiole
		{
			//T_Axil=T_air;  // could also be T_Leaf, or T_bud!
			if (T_Axil>0) e_sat=611.21*exp((18.678-T_Axil/234.5)*T_Axil/(257.14+T_Axil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
			else          e_sat=611.15*exp((23.036-T_Axil/333.7)*T_Axil/(279.82+T_Axil));
			e=            e_sat*exp(P_Petiole_Symp*2.16947115/(T_Axil+273.15));
			VPD_Petiole=(e-e_air)/1000; //vpd between Branch and air in mmol/m3
			if (VPD_Petiole<0) VPD_Petiole=0;
		}

	}
	// BRANCH
	if (CLIMAT !=6) T_Branch=T_Branch1=T_Branch2=T_Branch3=T_air;
	else T_Branch=T_Branch1;  																	// until each Branch has it own VPD
	if (T_Branch>0) e_sat=611.21*exp((18.678-T_Branch/234.5)*T_Branch/(257.14+T_Branch));  	// saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Branch/333.7)*T_Branch/(279.82+T_Branch));
	for (i=1;i<4;i++)
	{
		e= e_sat*exp(P_Branch_Symp[i]*2.16947115/(T_Branch+273.15));
		VPD_Branch[i]=(e-e_air)/1000; //vpd between Branch and air in mmol/m3
		if (VPD_Branch[i]<0) VPD_Branch[i]=0;
	}
	// TRUNK
	if (CLIMAT !=6) T_Trunk=T_air;
	if (T_Trunk>0)  e_sat=611.21*exp((18.678-T_Trunk/234.5)*T_Trunk/(257.14+T_Trunk));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Trunk/333.7)*T_Trunk/(279.82+T_Trunk));
	e=              e_sat*exp(P_Trunk_Symp*2.16947115/(T_Trunk+273.15));
	VPD_Trunk=(e-e_air)/1000; //vpd between Trunk and air in mmol/m3
	if (VPD_Trunk<0) VPD_Trunk=0;
	
	// ROOT
	
	if (CLIMAT==4 || CLIMAT==6 ) T_Root_1=T_Soil1;
	else T_Root_1=T_Soil;
	if (T_Root_1>0)		e_sat=611.21*exp((18.678-T_Root_1/234.5)*T_Root_1/(257.14+T_Root_1));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            	e_sat=611.15*exp((23.036-T_Root_1/333.7)*T_Root_1/(279.82+T_Root_1));
	
	e_air_sol1=     e_sat*exp(P_Soil1 *2.16947115/(T_Soil1+273.15));               // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp1*2.16947115/(T_Root_1+273.15));         // Root vapour pressure at Root water potential
	VPD_Root1=     (e-e_air_sol1)/1000; //vpd between Root   and air in mmol/m3
	if (VPD_Root1<0) VPD_Root1=0;
	
	if (CLIMAT==4 || CLIMAT==6 ) T_Root_2=T_Soil2;
	else T_Root_2=T_Soil;
	e_air_sol2=     e_sat*exp(P_Soil2 *2.16947115/(T_Soil2+273.15));           // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp2*2.16947115/(T_Root_2+273.15));
	VPD_Root2=     (e-e_air_sol2)/1000; //vpd between Root and air in mmol/m3
	if (VPD_Root2<0) VPD_Root2=0;
	
	if (CLIMAT==4 || CLIMAT==6 ) T_Root_3=T_Soil3;
	else T_Root_3=T_Soil;
	e_air_sol3=     e_sat*exp(P_Soil3 *2.16947115/(T_Soil3+273.15));           // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp3*2.16947115/(T_Root_3+273.15));
	VPD_Root3=     (e-e_air_sol3)/1000; //vpd between Root and air in mmol/m3
	if (VPD_Root3<0) VPD_Root3=0;
	
	// SOIL
	//if (CLIMAT !=4)T_Root=T_Soil;
	//else T_Root=T_Soil1;
	if (T_Soil>0)   e_sat=611.21*exp((18.678-T_Soil/234.5)*T_Soil/(257.14+T_Soil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Soil/333.7)*T_Soil/(279.82+T_Soil));
	
	e=              e_sat*exp(P_Soil1*2.16947115/(T_Soil+273.15));
	VPD_Soil=      (e-e_air)/1000; //vpd between soil and air in mmol/m3
	if (VPD_Soil<0) VPD_Soil=0;
	
}

void Interception (void) //rain interception from foliage
{
	long double In; //interception in mm
	if (Surface_Soil>Crown_Area) Rain_soil+=Rain_1*(Surface_Soil-Crown_Area)/Surface_Soil; //when canopy is not closed
	Interception_min=10+3*LAI_Crown;
	if (Leaf_Area[1])
	{
		if (Leaf_rain<0) Leaf_rain=0; //Leaf_rain = amount of water on the leaves, mm
		if (Rain_1>Threshold_rain)  In=Rain_1*(Interception_min+(100- Interception_min)/(exp(Interception_factor*(Rain_1 - Threshold_rain))))/100;  // then compute the interception in mm
		else                        In=Rain_1;
		Leaf_rain+=In;                                                                                                  // fill Leaf reservoir with IN, max is Interception_min
		if (Leaf_rain>Threshold_rain)                                                                                   // reservoir is saturated, the excess is not intercepted
		{
			In -= (Leaf_rain-Threshold_rain);
			Leaf_rain=Threshold_rain ;
		}
		Rain_1 = Rain_1-In;
	}
	if (Surface_Soil>Crown_Area) Rain_soil+=Rain_1*Crown_Area/Surface_Soil;
	else Rain_soil+=Rain_1;
}

void acclimate(void)
{
	long double  DOYb;
	long double  P50_rate=Acc_P1;
	long double  PIO_rate=Acc_P2;
	long double  Px_gs_rate=Acc_P3;
	int i;
	
	DOYb=T/3600/24/1000;
	if (DOYb<LA_day2)
	{
	for (i=1;i<4;i++)
		{
		P50_Leaf_Apo[i]=P50_Leaf_Apo_0[i];
		P50_Branch_Apo[i]=P50_Branch_Apo_0[i];
		}
	
	P50_Trunk_Apo=P50_Trunk_Apo_0;
	P50_Root_Apo=P50_Root_Apo_0;
	Pi0_Leaf_Symp=Pi0_Leaf_Symp_0;
	Pi0_Branch_Symp=Pi0_Branch_Symp_0;
	Pi0_Trunk_Symp=Pi0_Trunk_Symp_0;
	Px_gs=Px_gs_0;
	}
	else if (DOYb>LA_day2 && DOYb<LA_day3)  // acclimatation only when LAI is max
	{
		for (i=1;i<4;i++)
		{
		P50_Leaf_Apo[i]=	P50_Leaf_Apo_0[i]		+P50_rate*(DOYb-LA_day2);
		P50_Branch_Apo[i]=	P50_Branch_Apo_0[i]	+P50_rate*(DOYb-LA_day2);
		}
		P50_Trunk_Apo=		P50_Trunk_Apo_0	+P50_rate*(DOYb-LA_day2);
		P50_Root_Apo=		P50_Root_Apo_0		+P50_rate*(DOYb-LA_day2);
		Pi0_Leaf_Symp=		Pi0_Leaf_Symp_0	+PIO_rate*(DOYb-LA_day2);
		Pi0_Branch_Symp=	Pi0_Branch_Symp_0	+PIO_rate*(DOYb-LA_day2);
		Pi0_Trunk_Symp=	Pi0_Trunk_Symp_0	+PIO_rate*(DOYb-LA_day2);
		Px_gs=Px_gs_0+Px_gs_rate*(DOYb-LA_day2);
	}
}

void LA_acclimatation(void)  //Leaf area acclimatation from year to year; work only in continuuous mode!
{
	if (LA_Var==12) LA_max=LA_max_init*(100-PLC_Branch_Apo[1])/100; //based on top Branch
	if (LA_Var==13) LA_max=LA_max_init*(100-PLC_Leaf_Apo[1])/100;
}

void purge(void)
{
//	printf("skip end of year \n");
	while (T_2>T_1 && !END_CLIMAT) 
	{
		++N_days;
	//	if (Leaf_Area && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4))) ++N_days2;
			
		if (!feof(climat_in))
			{
			if (CLIMAT==2)   fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
			else             fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
			}
			else END_CLIMAT=1;
			T_2=T_2*3600*24;
			DOY++;
			Cum_T_air+= (T_air_min_2+T_air_max_2)/2;
			if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
			Rain_tot+=Rain_2;
			if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4))) Rain_Leaf_tot+=Rain_2;
	} 
	
	DEAD=0;
	//if (PRINT_SCREEN) print_screen();
	//if (TRANSIENT) print_transient();
	DOY=T_2/3600/24;
	if (LA_Var>10) LA_acclimatation();
	T0=T_2*1000;
	T=T0;
	T_1=T_0;
	indice=0;
	indice_double=0;
	if (!END_CLIMAT) Reset(); 
	
}

void next_climat (void)  // load new set of climatic data
{
	long double buffer, day , e_sat_air, e_air;
	if (CLIMAT==1 || CLIMAT==6 || CLIMAT==11)  //data on a hour basis
	{
		if (T/1000>T_2)      // reach the end of the time interval; load new values
		{
			Rain_tot+=Rain_1;
			if (Leaf_Area[1] && LA_Var && T/3600/24/1000>=LA_para1 && T/3600/24/1000<LA_para4) Rain_Leaf_tot+=Rain_1;
			if (Surface_Soil>Crown_Area) if ((COMPET && COMPET_number==0) || !COMPET) Q_Soil1+=Rain_1*(Surface_Soil-Crown_Area)*1000*1000/18; //add rain falling directly on the soil from the past interval
			if (INTERCEPTION) Interception();
			else if ((COMPET && COMPET_number==0) || !COMPET) Rain_soil+=Rain_1;
			if ((COMPET && COMPET_number==0) || !COMPET) Q_Soil1+=Rain_1*Crown_Area*1000*1000/18; //add rain falling under the canopy from the past interval to soil; Rain_1 is depleted from foliage interception
			fill_soil_layers();

			T_1=T_2;
			T_air_1=T_air_2;
			RH_air_1=RH_air_2;
			PAR_1=PAR_2;
			T_Soil_1=T_Soil_2;
			if (CLIMAT==6)
			{
				T_Soil_11=T_Soil_21;
				T_Soil_12=T_Soil_22;
				T_Soil_13=T_Soil_23;
				if (T_Soil_11<T_Soil_Crit || T_Soil_12<T_Soil_Crit || T_Soil_13<T_Soil_Crit) T_LIMIT=1;
				else if (T_Soil_11>T_Soil_Crit+0.1 && T_Soil_12>T_Soil_Crit+0.1 && T_Soil_13>T_Soil_Crit+0.1) T_LIMIT=0;
				else T_LIMIT=1;
				T_Trunk_1 = T_Trunk_2;
				T_Branch1_1= T_Branch1_2;
				T_Branch2_1= T_Branch2_2;
				T_Branch3_1= T_Branch3_2;
			}
		
			Rain_1=Rain_2;
			if (!feof(climat_in))
			{
				YEAR1=YEAR2;
				if (CLIMAT==1 || CLIMAT==11) 	fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR2, &T_2, &PAR_2,&T_air_2, &RH_air_2,&Rain_2, &Wind[0]);
				else 						fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &PAR_2, &T_Soil_21, &T_Soil_22, &T_Soil_23, &T_Trunk_2,&T_Branch1_2,&T_Branch2_2,&T_Branch3_2,&T_air_2,&RH_air_2,&Rain_2, &Wind[0]);
				if (PAR_1<0) PAR_1=0; if (PAR_2<0) PAR_2=0;
				T_2=T_2*3600*24;  //time in sec
			   // if (T_air_2 < 0.5) T_air_2=0.5; // otherwise Leaf energy budget doesn't work
				if (HW) // a Heat Wave
				{
					day= T/24/3600/1000;
					if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
					{
						if (HW==1) // spécific humidity is constant
							{
							RH_air_2*= exp((18.678-T_air_2/234.5)*T_air_2/(257.14+T_air_2))/exp((18.678-(T_air_2+HW_T)/234.5)*(T_air_2+HW_T)/(257.14+T_air_2+HW_T));  //new RH_air assuming constant e_air
							T_air_2+=HW_T;
							}
						 else if (HW==3) // vpd is constant
							{
							RH_air_2= 100-(611.21*exp((18.678-T_air_2/234.5)*T_air_2/(257.14+T_air_2)))/(611.21*exp((18.678-(T_air_2+HW_T)/234.5)*(T_air_2+HW_T)/(257.14+(T_air_2+HW_T))))*(100-RH_air_2);   //new RH_air assuming constant VPD_air
							T_air_2+=HW_T;
							}
						else if (HW==4) // a dry air wave at constant T_air
							{
							e_sat_air= 611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
							e_air =    e_sat_air*RH_air_2/100;                                       // vapour water pressure at Tair and RHair
							VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
							RH_air_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
							}
						else if (HW==6) // a wind wave !
							{
							Wind[0]*=HW_T;
							Wind_profile();
							}
						else T_air_2+=HW_T; 
					}
				}
		
				if (T_2<T_1) // then it is a new run or a new year. Reset simulation to zero
				{
					New_Year=1;
					DOY=(long double)(long long)T_2/3600/24;
				//	if (DOY<1) DOY=0;
					if (LA_Var>10) LA_acclimatation();
					 T0=T_2*1000;
					 T=T0;
					 T_1=T_0;
					 indice=0;
					 indice_double=0;
					 print_annual();
					 print_final();
					 if (!CONTINUOUS) Reset(); else legacy();
				}
			}
			else END_CLIMAT=1;
		}
	}
	
	if (CLIMAT==2 || CLIMAT==4) // data on a day basis
	{
		if (T/1000>=T_2) //end of day or DEAD 
		{
			++N_days;
		//	if (Leaf_Area && LA_Var && T/3600/24>=LA_para1 && T/3600/24<LA_para4) ++N_days2;
		
			Rain_tot+=Rain_1;
			if (Leaf_Area[1] && LA_Var && T/3600/24/1000>=LA_para1 && T/3600/24/1000<LA_para4) Rain_Leaf_tot+=Rain_1;
			if (Surface_Soil>Crown_Area) if ((COMPET && COMPET_number==0) || !COMPET) Q_Soil1+=Rain_1*(Surface_Soil-Crown_Area)*1000*1000/18; //add rain falling directly on the soil from the past interval
			if (INTERCEPTION) Interception();
			else if ((COMPET && COMPET_number==0) || !COMPET) Rain_soil+=Rain_1;
			if ((COMPET && COMPET_number==0) || !COMPET) Q_Soil1+=Rain_1*Crown_Area*1000*1000/18; //add rain falling under the canopy from the past interval to soil; Rain_1 is depleted from foliage interception
			fill_soil_layers();
			
			
			
	   //     if (IRRIGATE) Irrigate();   //refill soil water content by irrigation
			T_Soil_1=T_Soil_2;
			if (CLIMAT==4)
			{
				T_Soil_11=T_Soil_21;
				T_Soil_12=T_Soil_22;
				T_Soil_13=T_Soil_23;
				if (T_Soil_11<T_Soil_Crit || T_Soil_12<T_Soil_Crit || T_Soil_13<T_Soil_Crit) T_LIMIT=1;
				else if (T_Soil_11>T_Soil_Crit+0.1 && T_Soil_12>T_Soil_Crit+0.1 && T_Soil_13>T_Soil_Crit+0.1) T_LIMIT=0;
				else T_LIMIT=1;
			}
			
			T_air_max_0=T_air_max;
			T_air_min=T_air_min_2;  //set values for next day loaded before
			T_air_max=T_air_max_2;
			RH_air_min_0=RH_air_min;
			RH_air_min=RH_air_min_2;
			RH_air_max=RH_air_max_2;
			PAR_max=PAR_max_2*PAR_att;
			Rain_1=Rain_2;
			T_Soil=T_Soil_2;
			T_1=T_2;
		   
			if (!feof(climat_in))
			{
				YEAR1=YEAR2;
				DOY=T_2/3600/24;
				if (CLIMAT==2) 
				{
					if (LA_Var!=6) fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",    &YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
					else 		  fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0], &LA_max_Pheno);
				}
				if (CLIMAT==4)   fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_21,&T_Soil_22,&T_Soil_23, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
			
			
				T_2=T_2*3600*24;
				if (HW) // a Heat Wave
				{
					day= T/24/3600/1000;
					if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
					{
						if (HW==1) 
							{
							RH_air_min_2*= exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2))/exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+T_air_max_2+HW_T));  //new RH_air assuming constant e_air
							RH_air_max_2*= exp((18.678-T_air_min_2/234.5)*T_air_min_2/(257.14+T_air_min_2))/exp((18.678-(T_air_min_2+HW_T)/234.5)*(T_air_min_2+HW_T)/(257.14+T_air_min_2+HW_T));  //new RH_air assuming constant e_air
							T_air_max_2+=HW_T;
							T_air_min_2+=HW_T;
							}
					   else if (HW==3) 
						   {
							RH_air_min_2= 100-(611.21*exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2)))/(611.21*exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+(T_air_max_2+HW_T))))*(100-RH_air_min_2);    
							RH_air_max_2= 100-(611.21*exp((18.678-T_air_min_2/234.5)*T_air_min_2/(257.14+T_air_min_2)))/(611.21*exp((18.678-(T_air_min_2+HW_T)/234.5)*(T_air_min_2+HW_T)/(257.14+(T_air_min_2+HW_T))))*(100-RH_air_max_2);    
							T_air_max_2+=HW_T;
							T_air_min_2+=HW_T;
							}
						else if (HW==4) // a dry air wave at constant T_air
							{
							e_sat_air= 611.21*exp((18.678-T_air_max/234.5)*T_air_max/(257.14+T_air_max));    // saturation vapour water pressure at Tair in Pa from Buck's equation
							e_air =    e_sat_air*RH_air_min_2/100;                                       // vapour water pressure at Tair and RHair
							VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
							RH_air_min_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
							
							e_sat_air= 611.21*exp((18.678-T_air_min/234.5)*T_air_min/(257.14+T_air_min));    // saturation vapour water pressure at Tair in Pa from Buck's equation
							e_air =    e_sat_air*RH_air_max_2/100;                                       // vapour water pressure at Tair and RHair
							VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
							RH_air_max_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
							}
						else if (HW==6) // a wind wave !
							{
							Wind[0]*=HW_T;
							Wind_profile();
							}
						else 
						{
						T_air_max_2+=HW_T;
						T_air_min_2+=HW_T;
						}
					}
				}
			  
				Cum_T_air += (T_air_min_2+T_air_max_2)/2;
				if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
				if (T_air_min_2>T_air_max_2) // swap Tmin<>Tmax
				{
					buffer=T_air_min_2;
					T_air_min_2=T_air_max_2;
					T_air_max_2=buffer;
				}
				if (RH_air_min_2>RH_air_max_2) // swap RHmin<>RHmax
				{
					buffer=RH_air_min_2;
					RH_air_min_2=RH_air_max_2;
					RH_air_max_2=buffer;
				}
				if (RH_air_min_2==100) RH_air_min_2=99.9999;  // to avoid a bug if 100%
				if (Wind[0]<0.1) Wind[0]=0.1;
				
				if (T_2<T_1) // then it is a new run or a new year. Reset simulation to zero
				{	
					New_Year=1;
					DOY=T_2/3600/24;
					if (LA_Var>10) LA_acclimatation();
					 T0=T_2*1000;
					 T=T0;
					 T_1=T_0;
					 indice=0;
					 indice_double=0;
					 print_annual();
					 print_final();
					 if (!CONTINUOUS) Reset();
					 else legacy();		
				}
				
			}
			else END_CLIMAT=1; // the end of file
			
		   // T=T_1;
		}
		
	}
	
	if (CLIMAT==3) // data on a month basis
	{
		long double rain,proba;
		srand((unsigned) time(NULL));
		if (!(indice%(unsigned long)(3600*24/dt))) //end of day
		{
			++N_days;
	//		if (Leaf_Area && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4))) ++N_days2;
			Cum_T_air += (T_air_min_1+T_air_max_1)/2;
			if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))   Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
			
			proba=(long double)(rand() % 1000)/1000;                 // probability of rain this day
			if (proba<Proba_rain) rain=Rain_1/Proba_rain;            // rain this day
			else rain=0;
			
			
			//Rain_tot+=Rain_1;
			Rain_tot+=rain;
			if (Leaf_Area[1] && LA_Var && T/3600/24/1000>=LA_para1 && T/3600/24/1000<LA_para4) Rain_Leaf_tot+=Rain_1;
			if (Surface_Soil>Crown_Area) if (COMPET && COMPET_number==0) Q_Soil1+=Rain_1*(Surface_Soil-Crown_Area)*1000*1000/18; //add rain falling directly on the soil from the past interval
			if (INTERCEPTION) Interception();
			else Rain_soil+=Rain_1;
			if (COMPET && COMPET_number==0) Q_Soil1+=Rain_1*Crown_Area*1000*1000/18; //add rain falling under the canopy from the past interval to soil; Rain_1 is depleted from foliage interception
			fill_soil_layers();

			
	   //     if (IRRIGATE) Irrigate();  //refill soil water content by irrigation
		}
		
		if (T>T_2)                // a new month
		{
			T_1=T_2;
			if (!feof(climat_in))
			{
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Proba_rain, &Wind[0]);
				
				if (HW) // a Heat Wave
				{
					day= T/24/3600/1000;
					if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
					{
						if (HW==1) RH_air_min_2*= exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2))/exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+T_air_max_2+HW_T));  //new RH_air assuming constant e_air
						 if (HW==3) RH_air_min_2= 100-(611.21*exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2)))/(611.21*exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+(T_air_max_2+HW_T))))*(100-RH_air_min_2);    
						 //case 4 to be done...
						 T_air_max_2+=HW_T;
						T_air_min_2+=HW_T;
					}
				}
				T_2=(T_2)*3600*24*365/12;
				Cum_T_air+= (T_air_min_2+T_air_max_2)/2;
				if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
				T_air_min=T_air_min_2;  //set values for the first day
				T_air_max=T_air_max_2;
				T_air_max_0=T_air_max_2;
				T_air_min_1=T_air_min_2;
				T_air_max_1=T_air_max_2;
				
				RH_air_min=RH_air_min_2;
				RH_air_max=RH_air_max_2;
				PAR_max=PAR_max_2*PAR_att;
				RH_air_min_1=RH_air_min_2;
				RH_air_max_1=RH_air_max_2;
				PAR_max_1=PAR_max_2*PAR_att;
				Rain_1=Rain_2;
			}
			else END_CLIMAT=1;
			T=T_1;
		}
	}
	if (Wind[0]<=0) Wind[0]=0.1;
	if (CLIMAT!=4 && CLIMAT!=6)
	{
		T_Soil1=T_Soil;
		T_Soil2=T_Soil;
		T_Soil3=T_Soil;
	}
}

void load_climat (void)  // load the first set of climatic data from .txt files
{
	long double day, e_sat_air, e_air;
	
	if (CLIMAT==1 || CLIMAT==6 || CLIMAT==11) // data on a hour basis
	{
		if ((climat_in = fopen(filename_CLIM,"r+"))==NULL) 
		{
			CLIMAT=0;
			printf("\ncan't find climat file\n");
		}
		else
			{
			if (CLIMAT==1 || CLIMAT==11)
				{
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR1,&T_1, &PAR_1,&T_air_1, &RH_air_1,&Rain_1, &Wind[0]);
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR2,&T_2, &PAR_2,&T_air_2, &RH_air_2,&Rain_2, &Wind[0]);
				if (PAR_1<0) PAR_1=0; if (PAR_2<0) PAR_2=0;
				}	
			else   // CLIMAT==6 inludes soil temperature , Trunk and Branch in the climat file
				{
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &PAR_1, &T_Soil_11, &T_Soil_12, &T_Soil_13, &T_Trunk_1,&T_Branch1_1,&T_Branch2_1,&T_Branch3_1,&T_air_1,&RH_air_1,&Rain_1, &Wind[0]);
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &PAR_2, &T_Soil_21, &T_Soil_22, &T_Soil_23, &T_Trunk_2,&T_Branch1_2,&T_Branch2_2,&T_Branch3_2,&T_air_2,&RH_air_2,&Rain_2, &Wind[0]);
				T_Soil1=T_Soil_11;
				T_Soil2=T_Soil_12;
				T_Soil3=T_Soil_13;
				if (T_Soil1<T_Soil_Crit || T_Soil2<T_Soil_Crit || T_Soil3<T_Soil_Crit) T_LIMIT=1;
				else if (T_Soil1>T_Soil_Crit+0.1 && T_Soil2>T_Soil_Crit+0.1 && T_Soil3>T_Soil_Crit+0.1) T_LIMIT=0;
				else T_LIMIT=1;
				T_Trunk = T_Trunk_1;
				T_Branch1 = T_Branch1_1;
				T_Branch2 = T_Branch2_1;
				T_Branch3 = T_Branch3_1;
				}

			if (HW) // a Heat Wave
			{
				day= T/24/3600/1000;
				if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
				{
					if (HW==1)
						{
						RH_air_1*= exp((18.678-T_air_1/234.5)*T_air_1/(257.14+T_air_1))/exp((18.678-(T_air_1+HW_T)/234.5)*(T_air_1+HW_T)/(257.14+T_air_1+HW_T));  //new RH_air assuming constant e_air
						RH_air_2*= exp((18.678-T_air_2/234.5)*T_air_2/(257.14+T_air_2))/exp((18.678-(T_air_2+HW_T)/234.5)*(T_air_2+HW_T)/(257.14+T_air_2+HW_T));  //new RH_air assuming constant e_air
						T_air_1+=HW_T;
						T_air_2+=HW_T;
						}
					else if (HW==3) 
						{
						RH_air_1= 100-(611.21*exp((18.678-T_air_1/234.5)*T_air_1/(257.14+T_air_1)))/(611.21*exp((18.678-(T_air_1+HW_T)/234.5)*(T_air_1+HW_T)/(257.14+(T_air_1+HW_T))))*(100-RH_air_1);   //new RH_air assuming constant VPD_air
						RH_air_2= 100-(611.21*exp((18.678-T_air_2/234.5)*T_air_2/(257.14+T_air_2)))/(611.21*exp((18.678-(T_air_2+HW_T)/234.5)*(T_air_2+HW_T)/(257.14+(T_air_2+HW_T))))*(100-RH_air_2);   //new RH_air assuming constant VPD_air
						T_air_1+=HW_T;
						T_air_2+=HW_T;
						}
					else if (HW==4) // a dry air wave at constant T_air
						{
						e_sat_air= 611.21*exp((18.678-T_air_1/234.5)*T_air_1/(257.14+T_air_1));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_2/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_1= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
					
						e_sat_air= 611.21*exp((18.678-T_air_2/234.5)*T_air_2/(257.14+T_air_2));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_2/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
						}
					else 
						{
						T_air_1+=HW_T;
						T_air_2+=HW_T;
						}
					}
			}
			DOY=(float)(int)(T_1);
			T_1=T_1*24*3600;        // time of the day in secondes
			T_2=T_2*24*3600;
			T0=T_1*1000;
			T=T_1*1000;
			T_air=T_air_1;
			RH_air=RH_air_1;
			PAR[0]=PAR_1*PAR_att;
			//printf("%Lf %Lf",T_1, T);
		  }
	 
	}
	
	if (CLIMAT==2 || CLIMAT==4) // data on a day basis
	{
		if ((climat_in = fopen(filename_CLIM,"r+"))==NULL) 
			{
			CLIMAT=0;
			printf("\ncan't find climat file\n");
			}
		else
		{
			++N_days;
		//	if (Leaf_Area && LA_Var && T/3600/24>=LA_para1 && T/3600/24<LA_para4) ++N_days2;
		   
			if (CLIMAT==2)  //climat_day
			{
				if (LA_Var!=6)
				{
					fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &T_air_min_1, &T_air_max_1, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Wind[0]);
					fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
				}		

				else  //LA is the file
				{
					fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &T_air_min_1, &T_air_max_1, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Wind[0], &LA_max_Pheno);
					fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0], &LA_max_Pheno);
				}
			}
			
			else if (CLIMAT==4)  // CLIMAT==4 inludes soil temperature in the climat file
			{
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &T_air_min_1, &T_air_max_1, &T_Soil_11,&T_Soil_12,&T_Soil_13, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Wind[0]);
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_21, &T_Soil_22, &T_Soil_23, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind[0]);
				T_Soil=T_Soil_1;
				if (CLIMAT==4)
				{
					T_Soil1=T_Soil_11;
					T_Soil2=T_Soil_12;
					T_Soil3=T_Soil_13;
					if (T_Soil1<T_Soil_Crit || T_Soil2<T_Soil_Crit || T_Soil3<T_Soil_Crit) T_LIMIT=1;
					else if (T_Soil1>T_Soil_Crit+0.1 && T_Soil2>T_Soil_Crit+0.1 && T_Soil3>T_Soil_Crit+0.1) T_LIMIT=0;
					else T_LIMIT=1;
				}
				
			}
			if (HW) // a Heat Wave
			{
				day= T_1/24/3600;
				if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
				{
					if (HW==1) 
						{
						RH_air_min_1*= exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1))/exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+T_air_max_1+HW_T));  //new RH_air assuming constant e_air
						RH_air_min_2*= exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2))/exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+T_air_max_2+HW_T));  //new RH_air assuming constant e_air
						RH_air_max_1*= exp((18.678-T_air_min_1/234.5)*T_air_min_1/(257.14+T_air_min_1))/exp((18.678-(T_air_min_1+HW_T)/234.5)*(T_air_min_1+HW_T)/(257.14+T_air_min_1+HW_T));  //new RH_air assuming constant e_air
						RH_air_max_2*= exp((18.678-T_air_min_2/234.5)*T_air_min_2/(257.14+T_air_min_2))/exp((18.678-(T_air_min_2+HW_T)/234.5)*(T_air_min_2+HW_T)/(257.14+T_air_min_2+HW_T));  //new RH_air assuming constant e_air
						T_air_max_1+=HW_T;
						T_air_max_2+=HW_T;
						T_air_min_1+=HW_T;
						T_air_min_2+=HW_T;
						}
					else if (HW==3) 
						{
						RH_air_min_1= 100-(611.21*exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1)))/(611.21*exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+(T_air_max_1+HW_T))))*(100-RH_air_min_1);    
						RH_air_min_2= 100-(611.21*exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2)))/(611.21*exp((18.678-(T_air_max_2+HW_T)/234.5)*(T_air_max_2+HW_T)/(257.14+(T_air_max_2+HW_T))))*(100-RH_air_min_2);    
						RH_air_max_1= 100-(611.21*exp((18.678-T_air_min_1/234.5)*T_air_min_1/(257.14+T_air_min_1)))/(611.21*exp((18.678-(T_air_min_1+HW_T)/234.5)*(T_air_min_1+HW_T)/(257.14+(T_air_min_1+HW_T))))*(100-RH_air_max_1);    
						RH_air_max_2= 100-(611.21*exp((18.678-T_air_min_2/234.5)*T_air_min_2/(257.14+T_air_min_2)))/(611.21*exp((18.678-(T_air_min_2+HW_T)/234.5)*(T_air_min_2+HW_T)/(257.14+(T_air_min_2+HW_T))))*(100-RH_air_max_2);    
						T_air_max_1+=HW_T;
						T_air_max_2+=HW_T;
						T_air_min_1+=HW_T;
						T_air_min_2+=HW_T;
						}
					else if (HW==4) // a dry air wave at constant T_air
						{
						e_sat_air= 611.21*exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_min_1/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_min_1= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
					
						e_sat_air= 611.21*exp((18.678-T_air_max_2/234.5)*T_air_max_2/(257.14+T_air_max_2));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_min_2/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_min_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
						
						e_sat_air= 611.21*exp((18.678-T_air_min_1/234.5)*T_air_min_1/(257.14+T_air_min_1));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_max_1/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_max_1= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
					
						e_sat_air= 611.21*exp((18.678-T_air_min_2/234.5)*T_air_min_2/(257.14+T_air_min_2));    // saturation vapour water pressure at Tair in Pa from Buck's equation
						e_air =    e_sat_air*RH_air_max_2/100;                                       // vapour water pressure at Tair and RHair
						VPD_Air=  (e_sat_air - e_air)/1000;                                   // vpd in kPa
						RH_air_max_2= (e_sat_air- VPD_Air*HW_T*1000)/e_sat_air*100;
						}
					
					else 
						{
						T_air_max_1+=HW_T;
						T_air_max_2+=HW_T;
						T_air_min_1+=HW_T;
						T_air_min_2+=HW_T;
						}
				}
			 }
			
			Cum_T_air +=   (T_air_min_1+T_air_max_1)/2;
			if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))  Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
			T_air_min=      T_air_min_1;  //set values for the first day
			T_air_max=      T_air_max_1;
			T_air_max_0=    T_air_max_1;
			RH_air_min=     RH_air_min_1;
			RH_air_max=     RH_air_max_1;
			PAR_max=        PAR_max_1*PAR_att;
			DOY=T_1;
			T_1=T_1*24*3600;            // time of the day in secondes
			T_2=T_2*24*3600;
			T=T_1*1000;
		}
	}
	
	if (CLIMAT==3) // monthly average of daily values (including for precipitation)
	{
		if ((climat_in = fopen("climat_month_in.txt","r+"))==NULL) 
			{
			CLIMAT=0;
			printf("\ncan't find climat file\n");
			}
		else 
		{
		++N_days;
	//	if (Leaf_Area && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4))) ++N_days2;
	   
		fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&T_1, &T_air_min_1, &T_air_max_1, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Proba_rain, &Wind[0]);
		if (HW) // a Heat Wave
		{
			day= T/24/3600/1000;
			if ((day>HW_day && day<(HW_day+HW_duration)) || HW_duration==0)  // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
			{
				if (HW==1) RH_air_min_1*= exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1))/exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+T_air_max_1+HW_T));  //new RH_air assuming constant e_air
				 if (HW==3) RH_air_min_1= 100-(611.21*exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1)))/(611.21*exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+(T_air_max_1+HW_T))))*(100-RH_air_min_1);    
				   
				 T_air_max_1+=HW_T;
				T_air_min_1+=HW_T;
			}
		}
		Cum_T_air+= (T_air_min_1+T_air_max_1)/2;
		if (Leaf_Area[1] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))   Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
		T_air_min=T_air_min_1;  //set values for the first day
		T_air_max=T_air_max_1;
		T_air_max_0=T_air_max;
		RH_air_min=RH_air_min_1;
		RH_air_max=RH_air_max_1;
		PAR_max=PAR_max_1*PAR_att;
		
		T_air_min_2=T_air_min_1;
		T_air_max_2=T_air_max_1;
		RH_air_min_2=RH_air_min_1;
		RH_air_max_2=RH_air_max_1;
		PAR_max_2=PAR_max_1*PAR_att;
		Rain_2=Rain_1;
		T_2=(T_1)*24*3600*365/12;
		T_1=0 ;          // time of the day in secondes
		T=T_1*1000;
		}
	}
	if (Wind[0]<=0) Wind[0]=0.1;
	if (CLIMAT!=4 && CLIMAT!=6)
	{
		T_Soil=20;
		T_Soil1=T_Soil;
		T_Soil2=T_Soil;
		T_Soil3=T_Soil;
	}
	if (PRINT_SCREEN) printf("\n# %.0Lf Yr %.0Lf\n", Simul, YEAR1);
	else printf("\n# %.0Lf Yr %.0Lf", Simul, YEAR1);
}

void E_day(long double dt_court, long double dt_long, int i) //compute plant transpiration based on ETP, gs and VPD on mmol/m2/s
{
	long double E_Soil1, E_Soil2;
	//FILE *out;
	long double K_Root, K_Root1, K_Root2, K_Root3, K_to_Branch, K_to_Leaf, P_Soil;
	long double gs_max2, gs_night2, PL_gs, slope_gs, P50_gs,PTrunk; //PBranch,
	long double SLmax;
	int k;
	
	gs_0=g_s[i];
	g_crown=g_crown0*pow(Wind[0],0.6);             // wind effect on canopy conductance
	
	if (GS_MAX==1)      gs_max2 = gs_max[i]/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));                	// temperature effect on gs_max
	else if (GS_MAX==2) gs_max2=  gs_max[i]*300/Ca;                                                		// CO2 effect on gs_max; assume 300ppm in 1950
	else if (GS_MAX==3) gs_max2 = gs_max[i]*300/Ca/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));         	// CO2+T effects on gs_max
	else if (GS_MAX==4) gs_max2=  gs_max[i]*pow(300/Ca,0.5);  
	else if (GS_MAX==5) gs_max2 = gs_max[i]*pow(300/Ca,0.5)/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2)); 
	else if (GS_MAX==6) gs_max2 = gs_max[i]*(1 + gs_CO2_sens/100*(Ca-300)/100);
	else if (GS_MAX==7) gs_max2 =(gs_max[i]*(1 + gs_CO2_sens/100*(Ca-300)/100))/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));
	else                gs_max2=  gs_max[i];  
	if (gs_max2<0) gs_max2=0;                                                     // gs_max does not vary with temperature and CO2
	
	if (GS_MAX==1)      gs_night2 = gs_night/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));                // temperature effect on gs_min
	else if (GS_MAX==2) gs_night2=  gs_night*300/Ca;                                                // CO2 effect on gs_night; assume 300ppm in 1950
	else if (GS_MAX==3) gs_night2 = gs_night*300/Ca/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));         // CO2+T effects on gs_night
	else if (GS_MAX==4) gs_night2=  gs_night*pow(300/Ca,0.5);                                                // CO2 effect on gs_night; assume 300ppm in 1950
	else if (GS_MAX==5) gs_night2 = gs_night*pow(300/Ca,0.5)/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));         // CO2+T effects on gs_night
	else if (GS_MAX==6) gs_night2 = gs_night* (1 + gs_CO2_sens/100*(Ca-300)/100);
	else if (GS_MAX==7) gs_night2 = (gs_night* (1 + gs_CO2_sens/100*(Ca-300)/100))/(1+pow((T_Leaf[i]-Tgs_optim)/Tgs_sens,2));
	else                gs_night2=  gs_night;                                                       // gs_night does not vary with temperature and CO2
	if (gs_night2<0) gs_night2=0;
	
	if (!gs_cst) gs_Jarvis=   (gs_night2 + (gs_max2-gs_night2)*(1-exp(-Jarvis_PAR*PAR[i])));           // limitation of gs by light following Jarvis
	else         //gs_Jarvis=   (gs_night2 + (gs_max2-gs_night2)*(1-exp(-Jarvis_PAR*PAR_max)));       // gs does not respond to PAR
				  gs_Jarvis=   gs_max2;
					
	if (Leaf_Area[i] && gs_Jarvis)  E_clim[i] =    1/ (1/g_bl + 1/(gs_Jarvis+g_cuti[i])  + 1/g_crown) * VPD_Leaf[i]/101.3;                //Leaf  evaporation rate due to climatic demand in mmol s-1 m-2
	else E_clim[i]=0;
	
	if (PENMAN)
	{
		if (Penman_Coeff==0) Penman_Coeff= -0.006*LAI_Crown*LAI_Crown +0.134*LAI_Crown + 0.036;           // Empirical formula from Granier
		if (Leaf_Area[i]) 
		{
			if ((E_clim[i]*Leaf_Area[i]) >(ETP_Penman*Penman_Coeff*Surface_Soil)) 
				E_clim[i]=ETP_Penman*Penman_Coeff*Surface_Soil/Leaf_Area[i];        //Evaporation limited to ETP_Penman*coeff following Granier. soil evaporation is not included!
		}
		else E_clim[i]=0;
	}
	
	if (Leaf_Area[i] && g_cuti[i]) E_cuti[i] =    1/( 1/g_bl + 1/g_cuti[i] + 1/g_crown)  * VPD_Cuti[i]/101.3;                //evaporation through Leaf cuticle due to climatic demand
	else                     	    E_cuti[i] = 0;
	
	if (Type_Axil)  // evaporation from bud or flower or fruit
	{
		//g_Axil=(g_Axil_max-g_Axil_min)*(P_Axil_Symp + 2) +g_Axil_min;
		if (g_Axil_regul==1) g_Axil=g_Axil_max*(0.76*(P_Axil_Symp-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
		else if (g_Axil_regul==2) // g_Axil regulation by turgor_max=-PI0*turgor_ref_factor passed as Regul_gs_para1
		{			 
			if (Turgor_Axil_Symp<(-Pi0_Axil_Symp*Regul_gs_para1))  turgor=-Turgor_Axil_Symp/(Pi0_Axil_Symp*Regul_gs_para1);
			else  turgor=1;
			//Test[0]=turgor;
			g_Axil=g_Axil_max*turgor*(1-exp(-Jarvis_PAR*PAR[0]));  // including light effect
		}
		else if (g_Axil_regul==3) g_Axil=g_Axil_min;
		else g_Axil=g_Axil_max;
		
		if (g_Axil<g_Axil_min) g_Axil=g_Axil_min;
		if (g_Axil>g_Axil_max) if (g_Axil_regul!=3) g_Axil=g_Axil_max;
		if (g_Axil) E_Axil =  (1/ (1/g_bl_Axil + 1/g_Axil + 1/g_crown)) * VPD_Axil/101.3;           //for all the Axillary organs
		else E_Axil=0;
	}
	else 
		{
		g_Axil=0;
		E_Axil=0;
		}
	

	E_Branch[i]=   g_Branch*VPD_Branch[i]/101.3;        //evaporation through Branch cuticle due to climatic demand
	E_Trunk=    g_Trunk*VPD_Trunk/101.3;          //evaporation through Trunk cuticle due to climatic demand
	if (Type_Axil >= 2) E_Petiole=  g_Petiole*VPD_Petiole/101.3;  else E_Petiole=0;
	E_Root1=    g_Root[1]*VPD_Root1/101.3;          //evaporation through Root cuticle due to climatic demand
	E_Root2=    g_Root[2]*VPD_Root2/101.3;
	E_Root3=    g_Root[3]*VPD_Root3/101.3;
	
	if (RWC1<0.4) g_Soil=	g_Soil0*RWC1/0.4; // a kind of Biljou relationship
	else g_Soil=	g_Soil0;
	E_Soil1=	g_Soil*VPD_Soil/101.3;                      //VPD effect
	
	if (g_Soil0) E_Soil2=    g_Soil/g_Soil0*ETP_Penman*exp(-0.5*LAI_Soil);    // limitation by ETP depending on radiation reaching the soil
	else         E_Soil2=0;
	E_Soil=     min(E_Soil1,E_Soil2);
	if (T_Soil<0 || T_Soil1<0) E_Soil=0;                       //no evaporation from frozen soil

	dq_Branch[i]=   E_Branch[i]*Branch_Area[i]*dt_court;           // flows during dt in mmol
	dq_Axil=     E_Axil*Axil_Area*dt_court;               // for all the Axillary organs
	if (Type_Axil >= 2) dq_Petiole=  E_Petiole*Petiole_area*dt_court;   else dq_Petiole=0;     // for all petioles
	dq_cuti[i]=     Leaf_Area[i]*E_cuti[i]*dt_court;              // Plant evaporation through cuticle only during dt in mmol on the time interval dt; use one LA !
	if (Type_Axil == 2) dq_Fruit=    Growth_rate_Fruit*dt_court;  else dq_Fruit=0;           // water flow due to Fruit growth
	dq_Trunk=    E_Trunk*Trunk_Area*dt_court;
	dq_Root1=    E_Root1*Root_Area1*dt_court;
	dq_Root2=    E_Root2*Root_Area2*dt_court;
	dq_Root3=    E_Root3*Root_Area3*dt_court;
	dq_Soil=     E_Soil*Surface_Soil*dt_court;
 
	if (K_Soil1 && K_Interface1 &&K_Root_Symp11 && K_Root_Apo1) K_Root1=1/(1/K_Soil1 + 1/K_Interface1 + 1/K_Root_Symp11 + 1/K_Root_Apo1); else K_Root1=0;
	if (K_Soil2 && K_Interface2 &&K_Root_Symp12 && K_Root_Apo2) K_Root2=1/(1/K_Soil2 + 1/K_Interface2 + 1/K_Root_Symp12 + 1/K_Root_Apo2); else K_Root2=0;
	if (K_Soil3 && K_Interface3 &&K_Root_Symp13 && K_Root_Apo3) K_Root3=1/(1/K_Soil3 + 1/K_Interface3 + 1/K_Root_Symp13 + 1/K_Root_Apo3); else K_Root3=0;
	K_Root=K_Root1+K_Root2+K_Root3;
	K_Leaf_Apo[0]=K_Leaf_Symp[0]=K_Branch_Apo[0]=Leaf_Area[0]=0;
	for (k=1;k<4;k++) // conductances in parallel
		{
			K_Leaf_Apo[0]+=K_Leaf_Apo[k];
			K_Leaf_Symp[0]+=K_Leaf_Symp[k];
			K_Branch_Apo[0]+=K_Branch_Apo[k];
			Leaf_Area[0]+=Leaf_Area[k];
		}

	if (Leaf_Area[0] && K_Root && K_Trunk_Apo && K_Branch_Apo[0] && K_Leaf_Apo[0] && K_Leaf_Symp[0])  K_tot=1/(1/K_Root+1/K_Trunk_Apo+1/K_Branch_Apo[0]+1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0])/Leaf_Area[0];
	else K_tot=0;
	
	
	switch((int)Regul_gs) // if gs is computed then E is derived and vice versa
	{
		case 0: // gs is not regulated by water stress but still respond to PAR
			g_s[i]=gs_Jarvis;
			break;
			
			
		case 1:            // gs is decreased proportionnaly to the decrease in turgor pressure; ref is midday turgor of control
			if (Turgor_Leaf_Symp[i]<Turgor_Leaf_Symp_Ref) turgor=Turgor_Leaf_Symp[i]/Turgor_Leaf_Symp_Ref;
			else  turgor=1;
			g_s[i]=gs_max2*turgor;
			break;
			
		case 2:  // gs regulation by turgor_max=-PI0*turgor_ref_factor passed as Regul_gs_para1
			if (Turgor_Leaf_Symp[i]<(-Pi0_Leaf_Symp*Regul_gs_para1))  turgor=-Turgor_Leaf_Symp[i]/(Pi0_Leaf_Symp*Regul_gs_para1);
			else  turgor=1;
			g_s[i]=gs_max2*turgor;
			break;
		
		case 3:     // E respond to para2 = P_Branch_Apo; gamma attenuation factor is para1 (0.8 is typical)
			K_to_Branch=1/(1/K_Root+1/K_Trunk_Apo);
			if (K_Root) P_Soil= (P_Soil1*K_Root1 + P_Soil2*K_Root2 + P_Soil3*K_Root3)/K_Root;
				else    P_Soil=Regul_gs_para2+Pg_Leaf[1];
			PTrunk= P_Soil + Pg_Trunk - (F_Leaf[0]+F_Branch[0])/K_to_Branch;
			if (Leaf_Area[i]) E_Leaf[i]=Regul_gs_para1*(K_Branch_Apo[i]*(PTrunk-Regul_gs_para2+Pg_Leaf[i]))/Leaf_Area[i];
				else          E_Leaf[i]=0;
			if (E_Leaf[i]<0)E_Leaf[i]=0;
			E_stomata[i]=E_Leaf[i]-Gamma*E_cuti[i];
			if (VPD_Leaf[i]) g_s[i]=E_stomata[i]/VPD_Leaf[i]*101.3;  else g_s[i] = 0;
		
			break;
		
		case 4:     // E respond to  para2= P_Leaf_Apo; gamma attenuation factor is para1 (0.8 is typical)
			K_to_Leaf_Apo=1/(1/K_Root+1/K_Trunk_Apo+1/K_Branch_Apo[i]+1/K_Leaf_Apo[i]);
			if (K_Root) P_Soil= (P_Soil1*K_Root1 + P_Soil2*K_Root2 + P_Soil3*K_Root3)/K_Root;
			else        P_Soil=Regul_gs_para2+Pg_Leaf[1];
			if (Leaf_Area[i]) 	E_Leaf[i]=Regul_gs_para1*(K_to_Leaf_Apo*(P_Soil-Regul_gs_para2+Pg_Leaf[i]))/Leaf_Area[i];
			else           	E_Leaf[i]=0;
			if (E_Leaf[i]<0)E_Leaf[i]=0;
			E_stomata[i]= E_Leaf[i] - Gamma*E_cuti[i];
			if (VPD_Leaf[i]) g_s[i]=E_stomata[i]/VPD_Leaf[i]*101.3;  else g_s[i] = 0;
			break;
			
		case 5:     // E respond to para2=P_Leaf_evap; gamma attenuation factor is para1 (0.8 is typical)
			if (ACCLIMATE) Regul_gs_para2=Px_gs;
			K_to_Leaf=1/(1/K_Root+1/K_Trunk_Apo+1/K_Branch_Apo[i]+1/K_Leaf_Apo[i] +1/K_Leaf_Symp[i]);
			if (K_Root) P_Soil= (P_Soil1*K_Root1 + P_Soil2*K_Root2 + P_Soil3*K_Root3)/K_Root;
			else P_Soil=Regul_gs_para2+Pg_Leaf[1];
			if (Leaf_Area[i]) E_Leaf[i]=Regul_gs_para1*(K_to_Leaf*(P_Soil-Regul_gs_para2+Pg_Leaf[i]))/Leaf_Area[i];
			else           E_Leaf[i]=0;
			if (E_Leaf[i]<0)E_Leaf[i]=0;
			E_stomata[i]= E_Leaf[i] - Gamma*E_cuti[i];
			if (VPD_Leaf[i]) g_s[i]=E_stomata[i]/VPD_Leaf[i]*101.3;  else g_s[i] = 0;
			break;
			
		case 6:   //  E respond to P_Leaf_Apo set to para2 = Leaf_PLC; gamma attenuation factor is para1 (0.8 is typical)
			Px_Leaf_Apo[i]= P50_Leaf_Apo[i]*ST_Leaf+ 25/Slope_Leaf_Apo[i]*log((100-Regul_gs_para2)/Regul_gs_para2); //compute the Pressure for x PLC	Surface tension removed
			K_to_Leaf_Apo=1/(1/K_Root+1/K_Trunk_Apo+1/K_Branch_Apo[i]+1/K_Leaf_Apo[i]);
			if (K_Root) P_Soil= (P_Soil1*K_Root1 + P_Soil2*K_Root2 + P_Soil3*K_Root3)/K_Root;
				else    P_Soil= Px_Leaf_Apo[i]+Pg_Leaf[i];
			if (Leaf_Area[i]) E_Leaf[i]=Regul_gs_para1*(K_to_Leaf_Apo*(P_Soil-Px_Leaf_Apo[i]+Pg_Leaf[i]))/Leaf_Area[i];
				else          E_Leaf[i]=0;
			if (E_Leaf[i]<0)E_Leaf[i]=0;
	
			E_stomata[i]= E_Leaf[i] - Gamma*E_cuti[i];
			if (VPD_Leaf[i]) g_s[i]=E_stomata[i]/VPD_Leaf[i]*101.3;  else g_s[i] = 0;
			break;
			
		case 8: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Soil
			g_s[i]=gs_max2*(0.76*(P_soil-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
			
		case 9: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Leaf_Symp
			g_s[i]=gs_max2*(0.76*(P_Leaf_Symp[i]-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
		
		case 10: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Leaf_Symp
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Symp[i]-P50_gs)));
			g_s[i] = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 11: //   a power function of  g_s = Regul_gs_para1*pow(-P_Leaf_Apo,Regul_gs_para2)
			if (P_Leaf_Apo[i]<0) g_s[i] = Regul_gs_para1*pow(-P_Leaf_Apo[i],Regul_gs_para2); else g_s[i]=gs_max2;
			break;
			
		case 12: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Leaf_Apo
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Apo[i]-P50_gs)));
			g_s[i] = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 13: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Branch_Apo
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Branch_Apo[i]-P50_gs)));
			g_s[i] = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 14: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Branch_Apo
			g_s[i]=gs_max2*(0.76*(P_Branch_Apo[i]-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
			
		case 15: //like 2 but with a sigmoid fit between P_Leaf_Symp at 0.12 *max_turgor*para1 and 0.88*max_turgor*para1
			if (ACCLIMATE) Regul_gs_para2=Px_gs;
			Pgs_12=-Pi0_Leaf_Symp*Regul_gs_para1*Osmotic_TLeaf[i]*0.88+Epsilon_Leaf_Symp*Pi0_Leaf_Symp*Osmotic_TLeaf[i]/(Epsilon_Leaf_Symp+Pi0_Leaf_Symp*Osmotic_TLeaf[i]-Pi0_Leaf_Symp*Regul_gs_para1*Osmotic_TLeaf[i]*0.88);
			Pgs_88=-Pi0_Leaf_Symp*Regul_gs_para1*Osmotic_TLeaf[i]*0.12+Epsilon_Leaf_Symp*Pi0_Leaf_Symp*Osmotic_TLeaf[i]/(Epsilon_Leaf_Symp+Pi0_Leaf_Symp*Osmotic_TLeaf[i]-Pi0_Leaf_Symp*Regul_gs_para1*Osmotic_TLeaf[i]*0.12);
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Symp[i]-P50_gs)));
			g_s[i] = gs_max2 *(100-PL_gs)/100;
			break;
		
		case 16:
			g_s[i]= 1000*exp(P_Branch_Apo[i]*Regul_gs_para1);
			break;
			
		case 17: //gs decreases wit soil REW_wp below a threshold at Reguls_ge_para1
			REW_t=   (Teta_Soil1- Teta_r_1)/(Teta_fc_1-Teta_r_1)*Layer_1 + (Teta_Soil2- Teta_r_2)/(Teta_fc_2-Teta_r_2)*Layer_2 + (Teta_Soil3- Teta_r_3)/(Teta_fc_3-Teta_r_3)*Layer_3;  //REW based on Teta_r
			g_s[i]= gs_max2 *(REW_t/(Regul_gs_para1-Regul_gs_para2)-Regul_gs_para2/(Regul_gs_para1-Regul_gs_para2));
			break;
		
		case 18: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Root_Apo
			g_s[i]=gs_max2*(0.76*(P_Root_Apo-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;

		case 19: // gs regulated by P_Leaf_symp between Regul_gs_para1 and Regul_gs_para2
			g_s[i]=gs_max2/(Regul_gs_para1 - Regul_gs_para2)*(P_Leaf_Symp[i] - Regul_gs_para2);
			break;
			
		default: g_s[i]=gs_Jarvis; //E_Leaf=E_clim;
			break;
	}
	
	if (g_s[i]<0) g_s[i]=0;
	if (!gs_cst) 
	{
		if (g_s[i] > gs_Jarvis) g_s[i] = gs_Jarvis;
		else if (T<T_gs_regul) T_gs_regul=(T-T0);    // define the timing of onset of stomatal regulation
	}
	if (g_s[i] > gs_max2) g_s[i] = gs_max2;
	if (g_s[i]) g_Canopy=1/ (1/g_bl + 1/g_s[i] + 1/g_crown); else g_Canopy=0;
	E_stomata[i]=g_Canopy*VPD_Leaf[i]/101.3;
	E_Leaf[i] = E_stomata[i] + E_cuti[i];

   if (gs_tc) // if the time constant gs_tc of stomatal response is not null 
	{
	
		if ((g_s[i]-gs_0)>0) 	SLmax= (g_s[i]-gs_0)/(Regul_gs_para3*2.718);  	//slope of the gs increase with time
		else 				SLmax= (g_s[i]-gs_0)/(Regul_gs_para4*2.718);	//slope of the gs decrease with time
		g_s[i]=gs_0+SLmax*dt_long*dt;			//new stomatal conductance for canopy layer i
		if (g_s[i]) g_Canopy=1/ (1/g_bl + 1/g_s[i] + 1/g_crown); else g_Canopy=0; 
		E_stomata[i]=g_Canopy*VPD_Leaf[i]/101.3;
		E_Leaf[i] = E_stomata[i]+ E_cuti[i];
	} 
	 Test[0]=gs_0;
	 Test[1]=g_s[1];
	 Test[2]=SLmax;
	 Test[3]=gs_0+SLmax*dt_long*dt;	
	 
	dq_stomate[i]=  Leaf_Area[i]*E_stomata[i]*dt_court;     // Plant evaporation through stomata only during dt in mmol
	if (!Leaf_Area[i]) g_s[i]=0;
	if ((g_s[i] + g_cuti[i])>gs_max_d) gs_max_d = g_s[i] + g_cuti[i];
	if ((g_s[i] + g_cuti[i])<gs_min_d) gs_min_d = g_s[i] + g_cuti[i];
  
 // STATS  
	if (g_s[i] == gs_night2) if (E_Leaf[i] < E_Leaf_night) E_Leaf_night = E_Leaf[i];
	if (E_Leaf[i] <= (E_cuti[i]*1.01))  //when sigmoid fit E never reaches E_cuti so 1.01 is necessary
	{
		E_Leaf[i] = E_cuti[i];                                  	 // E cannot be less than E_cuti
		if (E_Leaf[i]>E_max_gs_close) E_max_gs_close=E_Leaf[i];     // Max E_Leaf when stomata are closed
		if (END_DEATH==13) DEAD=1;
	}
	if (Leaf_Area[1]) if (g_s[1] > (gs_max[1]/100)) T_gs_close=(T-T0); // define the timing of total of stomatal regulation when sigmoid fit g_s never reaches 0 so 1% of gs_max is the limit 
	if (Leaf_Area[1]) if (g_s[1] > 50) T_gs_50mmol=(T-T0);
	if (E_Leaf[i] > E_max) E_max=E_Leaf[i];  // grand Max E_Leaf
 
 if (!Leaf_Area[i])
	{
	gs_max_d  = 0;
	gs_min_d  = 0;
	gs_max_d2 = 0;
	}
	
 if (Leaf_rain)                                 // evaporate first the water on the leaves until dry
	{
	// set all E to 0
	E_Branch[i]=0; E_Trunk=0; E_Axil=0; E_cuti[i]=0; E_Leaf[i]=0; E_clim[i]=0; E_Petiole=0;
	// set all dq to 0
	dq_Branch[i]=0; dq_Axil=0; dq_Petiole=0;dq_cuti[i]=0; dq_stomate[i]=0; dq_Fruit=0; dq_Trunk=0;
	if (DYNAMIC==0) Leaf_rain-= ETP_Penman*18/1000/1000*dt;  //water loss by ETP in mm during dt_stat;
	else Leaf_rain-= ETP_Penman*18/1000/1000*dt_long*dt;
	if (Leaf_rain<0) Leaf_rain=0;
	}
}

void soil(long double time)   //pedotransfer functions following van Genuchten
{	
	Teta_Soil1=Q_Soil1/1000/1000/1000*18/Volume_soil1;
	Teta_Soil2=Q_Soil2/1000/1000/1000*18/Volume_soil2;
	Teta_Soil3=Q_Soil3/1000/1000/1000*18/Volume_soil3;

	RWC1=(Teta_Soil1- Teta_r_1)/(Teta_s_1 - Teta_r_1); 		// Relative volumetric water content
	if(RWC1>RWC_fc_1) RWC1=RWC_fc_1;		// cannot contain more water than field capacity 
	if(RWC1>0.9999) RWC1=0.9999;			// cannot contain more water than saturation  
	if (RWC1) P_Soil1=-(pow(((pow((1/RWC1),(1/m_1)))-1),(1/n_1)))/alpha_1/10000  + PI0_Soil/RWC1;  else P_Soil1=-9999;
	K_s=1000*K_sat_1*2*3.1416*Length_Root_fi*Root_upper/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_upper/Volume_soil1,0.5)/(Diam_Root/2))*Fluidity_soil1; //based on total Root system
	K_Soil1=K_s*pow(RWC1,L)*pow(1-pow(1-pow(RWC1,1/m_1),m_1),2);
	if (K_Soil1) K_Interface1=K_Soil1*10*pow(Q_Root_Symp1/Q_Root_Symp01,gap);    // Root-soil interface Conductance linked to Root diameter
	else         K_Interface1=0;
	K_Soil_tot1=Surface_Soil*1000*K_sat_1/(Soil_Depth*Layer_1)*pow(RWC1,L)*pow(1-pow(1-pow(RWC1,1/m_1),m_1),2)*Fluidity_soil1;

	RWC2=(Teta_Soil2- Teta_r_2)/(Teta_s_2 - Teta_r_2);
	if(RWC2>RWC_fc_2) RWC2=RWC_fc_2;
	if(RWC2>0.9999) RWC2=0.9999;
	if (RWC2) P_Soil2=-(pow(((pow((1/RWC2),(1/m_2)))-1),(1/n_2)))/alpha_2/10000 + PI0_Soil/RWC2; else P_Soil2=-9999;
	K_s=1000*K_sat_2*2*3.1416*Length_Root_fi*Root_middle/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_middle/Volume_soil2,0.5)/(Diam_Root/2))*Fluidity_soil1;
	K_Soil2=K_s*pow(RWC2,L)*pow(1-pow(1-pow(RWC2,1/m_2),m_2),2);
	if (K_Soil2)  K_Interface2=K_Soil2*10*pow(Q_Root_Symp2/Q_Root_Symp02,gap);
	else          K_Interface2 =0;
	K_Soil_tot2=Surface_Soil*1000*K_sat_2/(Soil_Depth*Layer_2)*pow(RWC2,L)*pow(1-pow(1-pow(RWC2,1/m_2),m_2),2)*Fluidity_soil2;

	RWC3=(Teta_Soil3- Teta_r_3)/(Teta_s_3 - Teta_r_3);
	if(RWC3>RWC_fc_3) RWC3=RWC_fc_3;
	if(RWC3>0.9999) RWC3=0.9999;
	if (RWC3) P_Soil3=-(pow(((pow((1/RWC3),(1/m_3)))-1),(1/n_3)))/alpha_3/10000 + PI0_Soil/RWC3; else P_Soil3=-9999;
	K_s=1000*K_sat_3*2*3.1416*Length_Root_fi*Root_lower/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_lower/Volume_soil3,0.5)/(Diam_Root/2))*Fluidity_soil2;
	K_Soil3=K_s*pow(RWC3,L)*pow(1-pow(1-pow(RWC3,1/m_3),m_3),2);
	if (K_Soil3)  K_Interface3=K_Soil3*10*pow(Q_Root_Symp3/Q_Root_Symp03,gap);
	else          K_Interface3 =0;
	K_Soil_tot3=Surface_Soil*1000*K_sat_3/(Soil_Depth*Layer_3)*pow(RWC3,L)*pow(1-pow(1-pow(RWC3,1/m_3),m_3),2)*Fluidity_soil3;
			
	Teta_Soil=Teta_Soil1*Layer_1 + Teta_Soil2*Layer_2 + Teta_Soil3*Layer_3;
	RWC_Soil_1= (Teta_Soil1- Teta_r_1)/(Teta_s_1 - Teta_r_1);
	RWC_Soil_2= (Teta_Soil2- Teta_r_2)/(Teta_s_2 - Teta_r_2);
	RWC_Soil_3= (Teta_Soil3- Teta_r_3)/(Teta_s_3 - Teta_r_3);
	RWC_Soil=   RWC_Soil_1*Layer_1 + RWC_Soil_2*Layer_2 + RWC_Soil_3*Layer_3; //mean soil RWC

	REW_t1=   (Teta_Soil1- Teta_r_1)/(Teta_fc_1-Teta_r_1);  //REW based on Teta_r
	REW_t2=   (Teta_Soil2- Teta_r_2)/(Teta_fc_2-Teta_r_2);  //REW based on Teta_r
	REW_t3=   (Teta_Soil3- Teta_r_3)/(Teta_fc_3-Teta_r_3);  //REW based on Teta_r
	REW_t=    REW_t1*Layer_1 + REW_t2*Layer_2 + REW_t3*Layer_3;  //REW based on Teta_r

	REW_wp1=  (Teta_Soil1- Teta_wp_1)/(Teta_fc_1-Teta_wp_1);  //REW based on Teta_wp pf4.2
	REW_wp2=  (Teta_Soil2- Teta_wp_2)/(Teta_fc_2-Teta_wp_2);  //REW based on Teta_wp pf4.2
	REW_wp3=  (Teta_Soil3- Teta_wp_3)/(Teta_fc_3-Teta_wp_3);  //REW based on Teta_wp pf4.2
	REW_wp=   REW_wp1*Layer_1 + REW_wp2*Layer_2 + REW_wp3*Layer_3;  //REW based on Teta_wp pf4.2
	
	if (RWC_Soil<RWC_min) RWC_min=RWC_Soil;  // min soil RWC
	if (RWC_Soil_1>RWC_fc_1) RWC_Soil_1=RWC_fc_1;
	if (RWC_Soil_2>RWC_fc_2) RWC_Soil_2=RWC_fc_2;
	if (RWC_Soil_3>RWC_fc_3) RWC_Soil_3=RWC_fc_3;
	
	RWC_int+=(1-RWC_Soil)*time/3600/24; 					// cumulated RWC deficit
	if (REW_t<=0.4) REW_int2+=(0.4-REW_t)*time/3600/24;   	// cumulated REW deficit based on Teta_r		
	if (REW_wp<=0.4)	REW_int1+=(0.4-REW_wp)*time/3600/24; 	// cumulated REW deficit, like in BilJou based on Teta_wp at Pf=4.2
	if (RWC_Soil>0) P_soil=P_Soil1*Layer_1 + P_Soil2*Layer_2 + P_Soil3*Layer_3; 
	else P_soil=-9999;//mean soil water potential
	P_soil_int +=P_soil*time/3600/24;
	if (P_soil<P_soil_min) P_soil_min=P_soil; // min soil WP
}


void  Compute_g_cuti(int i)  //temperature and RWC dependance of g_cuti
{
	long double g_cuti_tp, g_Axil_tp,RWC_Leaf,RWC_tlp, PLC, P_tlp;
	
	
	if (T_g_cuti==2) // computed with whole leaf RWC
	{
		P_tlp=(Epsilon_Leaf_Symp * Pi0_Leaf_Symp*Osmotic_TLeaf[i])/(Epsilon_Leaf_Symp + Pi0_Leaf_Symp*Osmotic_TLeaf[i]);
		PLC=100/(1+exp(Slope_Leaf_Apo[i]/25*(P_tlp-P50_Leaf_Apo[i]*ST_Leaf)));
		RWC_tlp=(Epsilon_Leaf_Symp + Pi0_Leaf_Symp*Osmotic_TLeaf[i])/Epsilon_Leaf_Symp*Leaf_Apo_fraction +  (100-PLC)/100*(1-Leaf_Apo_fraction);
		RWC_Leaf= (Q_Leaf_Apo0[i]*(100-PLC_Leaf_Apo[i])/100+Q_Leaf_Symp[i])/(Q_Leaf_Apo0[i] + Q_Leaf_Symp0[i]);
	}
	
	if (T_g_cuti==3) // computed with leaf symplasmic RWC only
	{
		RWC_tlp=(Epsilon_Leaf_Symp + Pi0_Leaf_Symp*Osmotic_TLeaf[i])/Epsilon_Leaf_Symp;
		RWC_Leaf= Q_Leaf_Symp[i]/Q_Leaf_Symp0[i];	
	}
	
	if (T_g_cuti==1)  // T effect
	{
		g_cuti_tp= g_cuti_20*pow(Q10_1,(TP-20)/10);
		g_Axil_tp= g_Axil_min20*pow(Q10_1_Axil,(TP_Axil-20)/10);
		
		if (T_Leaf[i]<TP) g_cuti[i]= g_cuti_20*pow(Q10_1,(T_Leaf[i]-20)/10);
		else 			 g_cuti[i]= g_cuti_tp*pow(Q10_2,(T_Leaf[i]-TP)/10);
			
		if (T_Axil<TP_Axil) g_Axil_min= g_Axil_min20*pow(Q10_1_Axil,(T_Axil-20)/10);
		else 			   g_Axil_min= g_Axil_tp*pow(Q10_2_Axil,(T_Axil-TP_Axil)/10);
	}
	
	else if (T_g_cuti==2 || T_g_cuti==3)  // linear RWC effect
	{
		if (RWC_Leaf<RWC_tlp) g_cuti[i]= 100*RWC_sens*(RWC_Leaf-RWC_tlp)+g_cuti_20;
		else g_cuti[i]= g_cuti_20;
	}
	
	
	else if (T_g_cuti==12 || T_g_cuti==13) // T + RWC linear
	{
		if (RWC_Leaf<RWC_tlp) g_cuti[i] = 100*RWC_sens*(RWC_Leaf-RWC_tlp)+g_cuti_20;
		else g_cuti[i]= g_cuti_20;
		
		g_cuti_tp= g_cuti[i]*pow(Q10_1,(TP-20)/10);
		
		if (T_Leaf[i]<TP) g_cuti[i]= g_cuti[i]*pow(Q10_1,(T_Leaf[i]-20)/10);
		else 			 g_cuti[i]= g_cuti_tp*pow(Q10_2,(T_Leaf[i]-TP)/10);
		
	}

	
}

void Compute_C (void)   //The apoplasmic capacitance if function of the apoplaspic water content;
{
	long double a,b;
	int i=1;
	a=0.0 ;  // fraction to prevent null C at 100 PLC
	b=1;
	
	for (i=1;i<4;i++) C_Leaf_Apo[i]    = b*C_Leaf_Apo0[i]   *(a+(1-a)*(100-PLC_Leaf_Apo[i])   /100);
	for (i=1;i<4;i++) C_Branch_Apo[i]  = b*C_Branch_Apo0[i] *(a+(1-a)*(100-PLC_Branch_Apo[i]) /100);
	C_Trunk_Apo   = b*C_Trunk_Apo0  *(a+(1-a)*(100-PLC_Trunk_Apo)  /100);
	C_Root_Apo1   = b*C_Root_Apo01  *(a+(1-a)*(100-PLC_Root_Apo1)  /100);
	C_Root_Apo2   = b*C_Root_Apo02  *(a+(1-a)*(100-PLC_Root_Apo2)  /100);
	C_Root_Apo3   = b*C_Root_Apo03  *(a+(1-a)*(100-PLC_Root_Apo3)  /100);
}

void Compute_Q_steady(void)
{
	long double Rs, Qold;
	int i;
	
	for (i=1;i<4;i++)
	{
	//Leaf
	Qold=Q_Leaf_Symp[i];
	Rs=(-(P_Leaf_Symp[i]+Pi0_Leaf_Symp-Epsilon_Leaf_Symp)-pow((pow((P_Leaf_Symp[i]+Pi0_Leaf_Symp-Epsilon_Leaf_Symp),2)+ 4*Epsilon_Leaf_Symp*P_Leaf_Symp[i]),0.5))/(2*Epsilon_Leaf_Symp);
	if (Rs<(1-Pi0_Leaf_Symp/P_Leaf_Symp[i])) Rs=1-Pi0_Leaf_Symp/P_Leaf_Symp[i];
	Turgor_Leaf_Symp[i]=-Pi0_Leaf_Symp*Osmotic_TLeaf[i] - Epsilon_Leaf_Symp*Rs;
	if (Turgor_Leaf_Symp[i]<0) Turgor_Leaf_Symp[i]=0;
	Q_Leaf_Symp[i]=Q_Leaf_Symp0[i]*(1-Rs);
	Q_Soil3+= Qold-Q_Leaf_Symp[i];
	
	//Branch
	Qold=Q_Branch_Symp[i];
	Rs=(-(P_Branch_Symp[i]+Pi0_Branch_Symp-Epsilon_Branch_Symp)-pow((pow((P_Branch_Symp[i]+Pi0_Branch_Symp-Epsilon_Branch_Symp),2)+ 4*Epsilon_Branch_Symp*P_Branch_Symp[i]),0.5))/(2*Epsilon_Branch_Symp);
	if (Rs<(1-Pi0_Branch_Symp/P_Branch_Symp[i])) Rs=1-Pi0_Branch_Symp/P_Branch_Symp[i];
	Q_Branch_Symp[i]=Q_Branch_Symp0[i]*(1-Rs);
	Q_Soil3+= Qold-Q_Branch_Symp[i];
	}
	//Axil
	if (Type_Axil)
	{
		Qold=Q_Axil_Symp;
		Rs=(-(P_Axil_Symp+Pi0_Axil_Symp-Epsilon_Axil_Symp)-pow((pow((P_Axil_Symp+Pi0_Axil_Symp-Epsilon_Axil_Symp),2)+ 4*Epsilon_Axil_Symp*P_Axil_Symp),0.5))/(2*Epsilon_Axil_Symp);
		if (Rs<(1-Pi0_Axil_Symp/P_Axil_Symp)) Rs=1-Pi0_Axil_Symp/P_Axil_Symp;
		Turgor_Axil_Symp=-Pi0_Axil_Symp*Osmotic_TAir - Epsilon_Axil_Symp*Rs;
		if (Turgor_Axil_Symp<0) Turgor_Axil_Symp=0;
		Q_Axil_Symp=Q_Axil_Symp0*(1-Rs);
		Q_Soil3+= Qold-Q_Axil_Symp;
		if (Type_Axil>=2) // a flower or a Fruit
		{
			Qold=Q_Petiole_Symp;
			Rs=(-(P_Petiole_Symp+Pi0_Axil_Symp-Epsilon_Axil_Symp)-pow((pow((P_Petiole_Symp+Pi0_Axil_Symp-Epsilon_Axil_Symp),2)+ 4*Epsilon_Axil_Symp*P_Petiole_Symp),0.5))/(2*Epsilon_Axil_Symp);
			if (Rs<(1-Pi0_Axil_Symp/P_Petiole_Symp)) Rs=1-Pi0_Axil_Symp/P_Petiole_Symp;
			Q_Petiole_Symp=Q_Petiole_Symp0*(1-Rs);
			Q_Soil3+= Qold-Q_Petiole_Symp;
		}
	}
	
	//Trunk
	Qold=Q_Trunk_Symp;
	Rs=(-(P_Trunk_Symp+Pi0_Trunk_Symp-Epsilon_Trunk_Symp)-pow((pow((P_Trunk_Symp+Pi0_Trunk_Symp-Epsilon_Trunk_Symp),2)+ 4*Epsilon_Trunk_Symp*P_Trunk_Symp),0.5))/(2*Epsilon_Trunk_Symp);
	if (Rs<(1-Pi0_Trunk_Symp/P_Trunk_Symp)) Rs=1-Pi0_Trunk_Symp/P_Trunk_Symp;
	Turgor_Trunk_Symp=-Pi0_Trunk_Symp*Osmotic_TAir - Epsilon_Trunk_Symp*Rs;
	if (Turgor_Trunk_Symp<0) Turgor_Trunk_Symp=0;
	Q_Trunk_Symp=Q_Trunk_Symp0*(1-Rs);
	Q_Soil3+= Qold-Q_Trunk_Symp;
	
	//Roots
	Qold=Q_Root_Symp1;
	Rs=(-(P_Root_Symp1+Pi0_Root_Symp-Epsilon_Root_Symp)-pow((pow((P_Root_Symp1+Pi0_Root_Symp-Epsilon_Root_Symp),2)+ 4*Epsilon_Root_Symp*P_Root_Symp1),0.5))/(2*Epsilon_Root_Symp);
	if (Rs<(1-Pi0_Root_Symp/P_Root_Symp1)) Rs=1-Pi0_Root_Symp/P_Root_Symp1;
	Q_Root_Symp1= Q_Root_Symp01 *(1-Rs);
	Q_Soil3+= Qold-Q_Root_Symp1;
	
	Qold=Q_Root_Symp2;
	Rs=(-(P_Root_Symp2+Pi0_Root_Symp-Epsilon_Root_Symp)-pow((pow((P_Root_Symp2+Pi0_Root_Symp-Epsilon_Root_Symp),2)+ 4*Epsilon_Root_Symp*P_Root_Symp2),0.5))/(2*Epsilon_Root_Symp);
	if (Rs<(1-Pi0_Root_Symp/P_Root_Symp2)) Rs=1-Pi0_Root_Symp/P_Root_Symp2;
	Q_Root_Symp2= Q_Root_Symp02 *(1-Rs);
	Q_Soil3+= Qold-Q_Root_Symp2;
	
	Qold=Q_Root_Symp3;
	Rs=(-(P_Root_Symp3+Pi0_Root_Symp-Epsilon_Root_Symp)-pow((pow((P_Root_Symp3+Pi0_Root_Symp-Epsilon_Root_Symp),2)+ 4*Epsilon_Root_Symp*P_Root_Symp3),0.5))/(2*Epsilon_Root_Symp);
	if (Rs<(1-Pi0_Root_Symp/P_Root_Symp3)) Rs=1-Pi0_Root_Symp/P_Root_Symp3;
	Q_Root_Symp3= Q_Root_Symp03 *(1-Rs);
	Q_Soil3+= Qold-Q_Root_Symp3;
}

void Compute_P_steady(long double dt_court)
{
	long double K_Root1, K_Root2, K_Root3, K_Root1b, K_Root2b, K_Root3b,FE_Root1, FE_Root2, FE_Root3;
	long double F_Trunk, F_Axil, F_Petiole, F_Cuti[4], F_Root1,F_Root2,F_Root3;
	long double LIMIT=0;
	int i;
	
	printf("PRINT 1");

	if (K_Soil1 && K_Interface1 && K_Root_Apo1 && K_Root_Symp11) K_Root1=1/(1/K_Soil1+1/K_Interface1+1/K_Root_Symp11+1/K_Root_Apo1);  else K_Root1=LIMIT;
	if (K_Soil2 && K_Interface2 && K_Root_Apo2 && K_Root_Symp12) K_Root2=1/(1/K_Soil2+1/K_Interface2+1/K_Root_Symp12+1/K_Root_Apo2);  else K_Root2=LIMIT;
	if (K_Soil3 && K_Interface3 && K_Root_Apo3 && K_Root_Symp13) K_Root3=1/(1/K_Soil3+1/K_Interface3+1/K_Root_Symp13+1/K_Root_Apo3);  else K_Root3=LIMIT;
	
	if (K_Soil1 && K_Interface1 && K_Root_Symp11) K_Root1b=1/(1/K_Soil1+1/K_Interface1+1/K_Root_Symp11); else K_Root1b=LIMIT;
	if (K_Soil2 && K_Interface2 && K_Root_Symp12) K_Root2b=1/(1/K_Soil2+1/K_Interface2+1/K_Root_Symp12); else K_Root2b=LIMIT;
	if (K_Soil3 && K_Interface3 && K_Root_Symp13) K_Root3b=1/(1/K_Soil3+1/K_Interface3+1/K_Root_Symp13); else K_Root3b=LIMIT;
	F_Leaf[0]= 0;
	F_Cuti[0]= 0;  
	F_Branch[0]=0; 

	for (i=1;i<4;i++)
		
		{
		if (Branch_distri[i])
		{
		F_Leaf[i]=  (E_Leaf[i]-E_cuti[i])*Leaf_Area[i]; //flow only through stomata
		F_Cuti[i]=   E_cuti[i]*Leaf_Area[i];
		F_Branch[i]= E_Branch[i]*Branch_Area[i]; // for all Branches
		}
		else
			{
			F_Leaf[i]=  0; //flow only through stomata
			F_Cuti[i]=   0;
			F_Branch[i]= 0; // for all Branches
		}

		F_Leaf[0]+= F_Leaf[i];
		F_Cuti[0]+= F_Cuti[i];  
		F_Branch[0]+=F_Branch[i]; 
		}
	F_Trunk=  E_Trunk*Trunk_Area;

	if (Type_Axil) F_Axil=   E_Axil*Axil_Area;  else F_Axil=0;   // for all flowers
	if (Type_Axil>=2) F_Petiole=E_Petiole*Petiole_area; else F_Petiole=0; //for one petioles
	
	FE_Root1= E_Root1*Root_Area1;
	FE_Root2= E_Root2*Root_Area2;
	FE_Root3= E_Root3*Root_Area3;
	
	P_Root_Apo1= -((F_Leaf[0]+F_Branch[0]+F_Trunk+F_Axil+F_Petiole+F_Cuti[0]) - K_Root1*P_Soil1  - K_Root2*P_Soil2 - K_Root3*P_Soil3)/(K_Root1+K_Root2+K_Root3);
	P_Root_Apo2=P_Root_Apo1;
	P_Root_Apo3=P_Root_Apo1;
	F_Root1= (P_Soil1-P_Root_Apo1)*K_Root1;
	F_Root2= (P_Soil2-P_Root_Apo2)*K_Root2;
	F_Root3= (P_Soil3-P_Root_Apo3)*K_Root3;
	//printf(""); // needed for PelleC to compile ???
	printf("PRINT 1");
	if (COMPET)
	{
		Q_Soil1-=F_Root1*dt_court*COMPET;
		Q_Soil2-=F_Root2*dt_court*COMPET;
		Q_Soil3-=F_Root3*dt_court*COMPET;
	}
	else
	{
		Q_Soil1-=F_Root1*dt_court;
		Q_Soil2-=F_Root2*dt_court;
		Q_Soil3-=F_Root3*dt_court;
	}
	
	if (K_Root1b) P_Root_Endo1= P_Soil1 - F_Root1/K_Root1b; else P_Root_Endo1=P_Root_Apo1;
	if (K_Root2b) P_Root_Endo2= P_Soil2 - F_Root2/K_Root2b; else P_Root_Endo2=P_Root_Apo2;
	if (K_Root3b) P_Root_Endo3= P_Soil3 - F_Root3/K_Root3b; else P_Root_Endo3=P_Root_Apo3;
	
	if (K_Root_Symp21)  P_Root_Symp1= P_Root_Endo1 -FE_Root1/(K_Root_Symp21); else P_Root_Symp1= P_Root_Endo1;
	if (K_Root_Symp22)  P_Root_Symp2= P_Root_Endo2 -FE_Root2/(K_Root_Symp22); else P_Root_Symp2= P_Root_Endo2;
	if (K_Root_Symp23)  P_Root_Symp3= P_Root_Endo3 -FE_Root3/(K_Root_Symp23); else P_Root_Symp3= P_Root_Endo3;
	
	if (K_Trunk_Apo) P_Trunk_Apo= Pg_Trunk+ P_Root_Apo1-(F_Leaf[0]+F_Branch[0]+F_Trunk+F_Axil+F_Petiole+F_Cuti[0])/K_Trunk_Apo; else P_Trunk_Apo=P_Root_Apo1;
	if (K_Trunk_Symp)P_Trunk_Symp=P_Trunk_Apo-F_Trunk/K_Trunk_Symp; else P_Trunk_Symp=P_Trunk_Apo;

	for (i=1;i<4;i++)
		{
			if (Branch_distri[i])
			{
			P_Branch_Apo[i]=Pg_Branch[i]- Pg_Trunk + P_Trunk_Apo-(F_Leaf[i]+F_Cuti[i]+F_Branch[i]+F_Axil+F_Petiole)/K_Branch_Apo[i]; 
			P_Branch_Symp[i]=P_Branch_Apo[i]-F_Branch[i]/K_Branch_Symp[i]; 
			}
			else 
			{
			P_Branch_Apo[i]=P_Trunk_Apo;
			P_Branch_Symp[i]=P_Trunk_Apo;	
			}
		}
	if (Type_Axil)
	{
		if (K_Axil_Apo) P_Axil_Apo=P_Branch_Apo[i]-(F_Axil+F_Petiole)/(K_Axil_Apo); else P_Axil_Apo=P_air; // no Symplasmic connexion here !
		P_Axil_Symp=P_Axil_Apo-F_Axil/(K_Axil_Symp);
		P_Petiole_Symp=P_Axil_Apo-F_Petiole/(K_Axil_Symp2);
	}
	for (i=1;i<4;i++)
		{
		if (Branch_distri[i])
			{
			if (Leaf_Area[i])   P_Leaf_Apo [i]=Pg_Leaf[i]-Pg_Branch[i] + P_Branch_Apo[i]-(F_Leaf[i]+F_Cuti[i])/(K_Leaf_Apo[i]); 	else P_Leaf_Apo[i] =P_Branch_Apo[i];
			if (Leaf_Area[i])   P_Leaf_Evap[i]=P_Leaf_Apo[i]-(F_Leaf[i]+F_Cuti[i])/(K_Leaf_Symp[i]); 								else P_Leaf_Evap[i]=P_Branch_Apo[i];
			if (Leaf_Area[i])   P_Leaf_Symp[i]=P_Leaf_Evap[i]-F_Cuti[i]/(K_Leaf_Symp2[i]); 										else P_Leaf_Symp[i]=P_Branch_Apo[i];
			}	
		else P_Leaf_Apo [i]=P_Leaf_Evap[i]=P_Leaf_Symp[i]=P_Trunk_Apo;
		}	
}

void Fluidity (void)   // compute water fluidity variation with temperature from Dehaoui et al PNAS 2015; works between -40°C to +60°C
{
	int i;
	if (FLUID)
	{
		for (i=1;i<4;i++) 
			Fluidity_Leaf[i] 	= 1.22783E-04 * T_Leaf[i]	* T_Leaf[i]  	+ 1.91158E-02 * T_Leaf[i] 	+ 5.64601E-01;
			Fluidity_Trunk 	= 1.22783E-04 * T_Trunk    * T_Trunk   	+ 1.91158E-02 * T_Trunk  	+ 5.64601E-01;
			Fluidity_Branch[1]	= 1.22783E-04 * T_Branch1  * T_Branch1   	+ 1.91158E-02 * T_Branch1  + 5.64601E-01;
			Fluidity_Branch[2]	= 1.22783E-04 * T_Branch2  * T_Branch2   	+ 1.91158E-02 * T_Branch2  + 5.64601E-01;
			Fluidity_Branch[3]	= 1.22783E-04 * T_Branch3  * T_Branch3   	+ 1.91158E-02 * T_Branch3  + 5.64601E-01;
			Fluidity_air  		= 1.22783E-04 * T_air  	* T_air   		+ 1.91158E-02 * T_air  	+ 5.64601E-01;
			Fluidity_soil  	= 1.22783E-04 * T_Soil 	* T_Soil  		+ 1.91158E-02 * T_Soil 	+ 5.64601E-01;
			Fluidity_soil1  	= 1.22783E-04 * T_Soil1	* T_Soil1 		+ 1.91158E-02 * T_Soil1	+ 5.64601E-01;
			Fluidity_soil2  	= 1.22783E-04 * T_Soil2	* T_Soil2 		+ 1.91158E-02 * T_Soil2	+ 5.64601E-01;
			Fluidity_soil3  	= 1.22783E-04 * T_Soil3	* T_Soil3		+ 1.91158E-02 * T_Soil3	+ 5.64601E-01;
	}
	else
	{
		for (i=1;i<4;i++) 
		{
			Fluidity_Leaf[i]  = 1;
			Fluidity_Branch[i]= 1;
		}
		Fluidity_air   = 1;
		Fluidity_Trunk = 1;
		Fluidity_soil  = 1;
		Fluidity_soil1 = 1;
		Fluidity_soil2 = 1;
		Fluidity_soil3 = 1;
	}
}

void Compute_K(void) // changes of K due to cavitation and phenology
{
	long double K,Cortex_Gap,R;
	int i;

	for (i=1;i<4;i++)
		{
		K_Leaf_Apo[i]    = K_Leaf_Apo0[i]    	*(100-PLC_Leaf_Apo[i])/100   * (100-PLF[i])/100 * Fluidity_Leaf[i];
		K_Branch_Apo[i]  = K_Branch_Apo0[i]  	*(100-PLC_Branch_Apo[i])/100 * Fluidity_Branch[i];
		K_Branch_Symp[i] = K_Branch_Symp0[i] 	* Fluidity_Branch[i];
		}
	
	K_Axil_Apo    = K_Axil_Apo0       *(100-PLC_Axil_Apo)/100   * Fluidity_air;
	K_Trunk_Apo   = K_Trunk_Apo0      *(100-PLC_Trunk_Apo)/100  * Fluidity_Trunk;
	K_Root_Apo1   = K_Root_Apo01      *(100-PLC_Root_Apo1)/100  * Fluidity_soil1 ;  // T_soil assumed constant
	K_Root_Apo2   = K_Root_Apo02      *(100-PLC_Root_Apo2)/100  * Fluidity_soil2 ;
	K_Root_Apo3   = K_Root_Apo03      *(100-PLC_Root_Apo3)/100  * Fluidity_soil3 ;

	K_Trunk_Symp  = K_Trunk_Symp0  * Fluidity_Trunk;
	
	// Leaf
	if (K_VAR==1)  // K_Leaf is variable
	{
		//K= 6.83 + 81.4* exp (7.56 * P_Leaf_Symp);
		for (i=1;i<4;i++)
		{
			K= K_VAR_P1 + K_VAR_P2* exp (K_VAR_P3 * P_Leaf_Symp[i]) ;
			if (Leaf_Area[i]) K_Leaf_Symp[i]= 1/(1/K -1/(K_Leaf_Apo0[i]/Leaf_Area[i])) * (100-PLF[i])/100* Fluidity_Leaf[i]; else K_Leaf_Symp[i]=99999;
			if (K_Leaf_Symp[i] > 106) K_Leaf_Symp[i]=106;
			K_Leaf_Symp[i]*=Leaf_Area[i];
		}
	}
	else if (K_VAR==11) // KLeaf varies with Y Leaf; for Lawren C3 C4
	{
		for (i=1;i<4;i++)
		{
			K= K_VAR_P1 + K_VAR_P2* P_Leaf_Symp[i] ;
			if (K<=0) K=K_VAR_P1/1000; //whole leak K
			if (Leaf_Area[i]) K_Leaf_Symp[i]= 1/(1/K -1/(K_Leaf_Apo0[i]/Leaf_Area[i])) * (100-PLF[i])/100* Fluidity_Leaf[i]; else K_Leaf_Symp[i]=99999;
			//if (K_Leaf_Symp[i] > K_VAR_P1) K_Leaf_Symp[i]=K_VAR_P1;
			K_Leaf_Symp[i]*=Leaf_Area[i];
		}
	}
	else  
	{
		if (!LA_Var) //Leaf area is constant
		{
			for (i=1;i<4;i++)
			{
			K_Leaf_Symp[i]  = K_Leaf_Symp_0[i]/LA_max_init*LA_max_Pheno*(100-PLF[i])/100 * Fluidity_Leaf[i];
			K_Leaf_Symp2[i] = K_Leaf_Symp[i];
			}
		}
		else 		// Leaf area varies with phenology
		{
			for (i=1;i<4;i++)
			{
			if (Leaf_Area[i])	K_Leaf_Symp[i]=K_Leaf_Symp_0[i]/LA_max_init*Leaf_Area[i]*(100-PLF[i])/100 * Fluidity_Leaf[i]; else K_Leaf_Symp[i]=0;//K_Leaf_Symp_0[i]/100;
			if (Leaf_Area[i])	K_Leaf_Apo[i]=K_Leaf_Apo0[i] *(100-PLC_Leaf_Apo[i])/100   * (100-PLF[i])/100 * Fluidity_Leaf[i]/LA_max_init*Leaf_Area[i];    else K_Leaf_Apo[i]=0;//K_Leaf_Apo0[i]/100;
			if (Leaf_Area[i])	K_Leaf_Symp2[i]=K_Leaf_Symp[i]; else K_Leaf_Symp2[i]=K_Leaf_Symp[i]; //A VERIFIER ?
			}
		}
		
	}
	
	// Root
	if (K_VAR>1)  //K is variable
	{
		switch((int)K_VAR)
		{
		case 2: // Arabidopsis Scoffoni
		
			//K=  4.25+495*exp(5.85*P_Leaf_Symp);
			K=K_VAR_P1+K_VAR_P2*exp(K_VAR_P3*P_Leaf_Symp[1]);
			K*=LA_max;
			K_Root_Symp1= 1/(1/K -1/(K_Root_Apo0)) ;
			if (K_Root_Symp1>1.6) K_Root_Symp1=1.6;
			break;
		
		case 3: // Effet température sol pour Doux-Glace
		
			if (T_Soil1 >20) K_Root_Symp11 = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi * Root_upper;                  //  follows fluidity
			else if (T_Soil1<T_Soil_Crit) K_Root_Symp11 = 0;                                                             //  zero when below T_Soil_Crit
			else K_Root_Symp11 = K_Root_Symp0 * Root_Area_fi * Root_upper * T_Soil1/(20-T_Soil_Crit);                    //  between T_Soil_Crit and 20 linear with T_soil
			
			if (T_Soil2 >20) K_Root_Symp12 = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi * Root_middle;
			else if (T_Soil2<T_Soil_Crit) K_Root_Symp12 = 0;
			else K_Root_Symp12 = K_Root_Symp0 * Root_Area_fi * Root_middle * T_Soil2/(20-T_Soil_Crit);
			
			if (T_Soil3 >20) K_Root_Symp13 = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi * Root_lower;
			else if (T_Soil3<T_Soil_Crit) K_Root_Symp13 = 0;
			else K_Root_Symp13 = K_Root_Symp0 * Root_Area_fi * Root_lower * T_Soil3/(20-T_Soil_Crit);
			
			if (T_Root_1 < K_VAR_P1) K_Root_Apo1=0; 		// frozen !
			if (T_Root_2 < K_VAR_P1) K_Root_Apo2=0;
			if (T_Root_3 < K_VAR_P1) K_Root_Apo3=0;
			if (T_Trunk  < K_VAR_P2) K_Trunk_Apo=0;
			if (T_Branch1< K_VAR_P3) K_Branch_Apo[1]=0;
			if (T_Branch2< K_VAR_P3) K_Branch_Apo[2]=0;
			if (T_Branch3< K_VAR_P3) K_Branch_Apo[3]=0;
			if (T_Leaf[1]< K_VAR_P3) K_Leaf_Apo[1]=0;
			if (T_Leaf[2]< K_VAR_P3) K_Leaf_Apo[2]=0;
			if (T_Leaf[3]< K_VAR_P3) K_Leaf_Apo[3]=0;
			break;
		
	  case 4: // air gap in the cortex following a sigmoidal curve K_VAR_P1=P50 K_VAR_P2=slope with P_Root_Symp; non reversible loss of K
		
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp1-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi*Root_upper*(100-Cortex_Gap)/100;
			if (K<K_Root_Symp11)  K_Root_Symp11=K;
			//K_Root_Symp21= 1 * K_Root_Symp11;
			//Test[1]=Cortex_Gap;
			//Test[2]=T_Soil1;
			//Test[3]=Fluidity_soil1;
			
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp2-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi* Root_middle *(100-Cortex_Gap)/100;
			if (K<K_Root_Symp12)  K_Root_Symp12=K;
			//K_Root_Symp22= 1 * K_Root_Symp12;
			
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp3-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi*Root_lower*(100-Cortex_Gap)/100;
			if (K<K_Root_Symp13)  K_Root_Symp13=K;
			//K_Root_Symp23= 1 * K_Root_Symp13;
			break;
		
	  case 5: // For Maddy and Tim
	  
		  R=pow(K_Root_Symp0/K_VAR_P1,-1/K_VAR_P2);
		  K_Root_Symp11=K_VAR_P1*pow(R-P_Soil1,K_VAR_P2)*Fluidity_soil1 * Root_Area_fi* Root_upper;
		  K_Root_Symp12=K_VAR_P1*pow(R-P_Soil2,K_VAR_P2)*Fluidity_soil2 * Root_Area_fi* Root_middle;
		  K_Root_Symp13=K_VAR_P1*pow(R-P_Soil3,K_VAR_P2)*Fluidity_soil3 * Root_Area_fi* Root_lower;
		break;
	  
	  case 6: // air gap in the cortex following a sigmoidal curve with P_soil; non reversible loss of K
	  
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil1-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi*Root_upper*(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp11)  K_Root_Symp11=K;
		 // K_Root_Symp21= 1 * K_Root_Symp11;
		  
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil2-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi* Root_middle *(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp12)  K_Root_Symp12=K;
		//  K_Root_Symp22= 1 * K_Root_Symp12;
		  
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil3-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi*Root_lower*(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp13)  K_Root_Symp13=K;
	   //   K_Root_Symp23= 1 * K_Root_Symp13;
		break;
	
	case 7: //K_Branch_Symp variable with PLC
	
		for (i=1;i<4;i++) K_Branch_Symp[i]=K_Branch_Symp0[i]*(100-PLC_Branch_Apo[i])/100 * Fluidity_Branch[i];
		K_Trunk_Symp=K_Trunk_Symp0*(100-PLC_Trunk_Apo)/100 * Fluidity_Trunk;
		break;
	
	case 8: //For Carola and Tim with K_VAR_P1=P50  K_VAR_P2=slope K_VAR_P3=Kmax fraction to have Kmin
	
		K_Root_Symp11=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp1-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil1 * Root_Area_fi* Root_upper;
		K_Root_Symp12=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp2-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil2 * Root_Area_fi* Root_middle;
		K_Root_Symp13=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp3-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil3 * Root_Area_fi* Root_lower;
		break;
	
	case 9 ://For URI
	
		Cortex_Gap= 0.5/(-10 * pow(P_Soil1,3) + 0.5);
		K = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi*Root_upper*Cortex_Gap;
	//  if (K<K_Root_Symp11)  //activate if irriversible
		K_Root_Symp11=K;
		
		Cortex_Gap= 0.5/(-10 * pow(P_Soil2,3) + 0.5);
		K = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi*Root_middle*Cortex_Gap;
	//	if (K<K_Root_Symp12)  
		K_Root_Symp12=K;
	
		Cortex_Gap= 0.5/(-10 * pow(P_Soil3,3) + 0.5);
		K = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi*Root_lower*Cortex_Gap;
	//	if (K<K_Root_Symp13) 
		K_Root_Symp13=K;	
		break;
	
	case 10: //For URI
	
		if ((T/3600/24/1000-DOY_0) > K_VAR_P1) 
		{
			K_Interface1=0;
			K_Interface2=0;
			K_Interface3=0;	
		}
		break;
	
	case 13: // air gap in the cortex following a lineear curve K_VAR_P1=REW init K_VAR_P2=REW end; non reversible loss of K
			Teta_Soil1=Q_Soil1/1000/1000/1000*18/Volume_soil1;
			R=(Teta_Soil1- Teta_r_1)/(Teta_fc_1-Teta_r_1);//REW based on Teta_r			
			Cortex_Gap= (R/(K_VAR_P1-K_VAR_P2)-K_VAR_P2/(K_VAR_P1-K_VAR_P2));			
			if (Cortex_Gap < 0) Cortex_Gap= 0;
			if (Cortex_Gap>1) Cortex_Gap=1;
			K = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi*Root_upper*Cortex_Gap;			
			if (K<K_Root_Symp11) K_Root_Symp11=K;
			
			Teta_Soil2=Q_Soil2/1000/1000/1000*18/Volume_soil2;
			R=(Teta_Soil2- Teta_r_2)/(Teta_fc_2-Teta_r_2);			
			Cortex_Gap= (R/(K_VAR_P1-K_VAR_P2)-K_VAR_P2/(K_VAR_P1-K_VAR_P2));			
			if (Cortex_Gap < 0) Cortex_Gap= 0;
			if (Cortex_Gap>1) Cortex_Gap=1;
			K = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi* Root_middle *Cortex_Gap;
			if (K<K_Root_Symp12) K_Root_Symp12=K;
			
			Teta_Soil3=Q_Soil3/1000/1000/1000*18/Volume_soil3;			
			R=(Teta_Soil3- Teta_r_3)/(Teta_fc_3-Teta_r_3);
			Cortex_Gap= (R/(K_VAR_P1-K_VAR_P2)-K_VAR_P2/(K_VAR_P1-K_VAR_P2));
			if (Cortex_Gap < 0) Cortex_Gap= 0;
			if (Cortex_Gap>1) Cortex_Gap=1;
			K = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi*Root_lower*Cortex_Gap;
			if (K<K_Root_Symp13)  K_Root_Symp13=K;
			break;
		
	default : // K_Root is constant
		K_Root_Symp1 =	K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi; //follows fluidity
		K_Root_Symp11=	K_Root_Symp1 * Root_upper;
		K_Root_Symp12=	K_Root_Symp1 * Root_middle;
		K_Root_Symp13=	K_Root_Symp1 * Root_lower;
		K_Root_Symp2 = 	K_Root_Symp0 * Fluidity_soil1 * Root_Area_FR;
		K_Root_Symp21=	K_Root_Symp2 * Root_upper;
		K_Root_Symp22=	K_Root_Symp2 * Root_middle;
		K_Root_Symp23=	K_Root_Symp2 * Root_lower;
		break;
	}
   }
}

void Compute_Leaf_Fall(int i)   //force leaves to fall according to Leaf water potential
{
	long double LA, PLF_i;
	
	if		(Leaf_Fall==1) PLF_i= 100/(1+exp(Slope_Leaf_Fall/25*(P_Leaf_Symp[i]-P50_Leaf_Fall)));
	else if	(Leaf_Fall==2) PLF_i= PLC_Leaf_Apo[i] ; //PLF_i=100/(1+exp(Slope_Leaf_Fall/25*(P_Leaf_Apo[i] -P50_Leaf_Fall)));
	else PLF_i=0;
	
	LA = Branch_distri[i]*LA_max_Pheno *(100-PLF_i)/100;
	if (LA==0) LA=0.001;
//  if (LA<Leaf_Area[i])
	if (PLF_i>PLF[i])
	{
//		if (Leaf_Area[i])  printf("%d %Lf %Lf\n",i, PLF[i],Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Apo1[i] *= (LA/Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Apo[i]  *= (LA/Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Symp0[i]*= (LA/Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Symp[i] *= (LA/Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Evap[i] *= (LA/Leaf_Area[i]);
		if (Leaf_Area[i])  Q_Leaf_Evap0[i]*= (LA/Leaf_Area[i]);
		Leaf_Area[i]=LA;
		PLF[i]=PLF_i;		
	}	
}

void Compute_Root_Fall(void)   //force Root to die according to Root water potential
{
	long double RA, PRF_i;
	
	PRF_i=100/(1+exp(Slope_Leaf_Fall/25*(P_Root_Symp1-P50_Leaf_Fall)));
	RA = Root_Area_fi_0 *(100-PRF_i)/100;
	if (RA==0) RA=0.001;
	if (RA<Root_Area_fi)
	{
		Root_Area_fi=RA;
		PRF=PRF_i;
	}
}

void Phenology (long double dtt)  //computes seasonal changes in Leaf area
{
	long double dGDD=0,Q_old_s,Root_min,LA0;
	int i;
	
	Root_min=0.1;  //min fraction of active Roots when LA=0
//	if (LA_Var==0) //Leaf area and Root area are constant
/*	{
		LA_max_Pheno=LA_max;
		Root_upper	= Root_upper0;
		Root_middle	= Root_middle0;
		Root_lower	= Root_lower0;	
	}
*/ 
	LA0=LA_max_Pheno;
	if (LA_Var==1 || LA_Var==3 || LA_Var==4) // Leaf phenology based on fixed dates; Root area constant to max
	{
		if (T/24/3600/1000<LA_day1) LA_max_Pheno=LA_min;
		else if (T/24/3600/1000>LA_day1 && T/24/3600/1000<LA_day2) LA_max_Pheno=LA_min + (LA_max - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1);
		else if (T/24/3600/1000>LA_day2 && T/24/3600/1000<LA_day3) LA_max_Pheno=LA_max;
		else if (T/24/3600/1000>LA_day3 && T/24/3600/1000<LA_day4) 
			{
			if (LA_Var==4) 	LA_max_Pheno=LA_max - LA_max*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3);
			else LA_max_Pheno=LA_max- (LA_max - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3);
			}
		else if (T/24/3600/1000>LA_day4) 
		{
			if (LA_Var==4) LA_max_Pheno=0;		// case for wheat
			else LA_max_Pheno=LA_min;
		}
	}
	
	if (LA_Var==2 || LA_Var==3|| LA_Var==4) // Root phenology based on fixe dates; assume that Root_upper_middle_lower are proportional to Leaf Area
	{
		if (LA_Var==2) LA_max_Pheno=LA_max;
		if (T/24/3600/1000<=LA_day1)
		{	
			if (LA_min) 
				{
				Root_upper	= LA_min/LA_max*Root_upper0 ;  // at leat 1% of Root area
				Root_middle	= LA_min/LA_max*Root_middle0;
				Root_lower	= LA_min/LA_max*Root_lower0;
				}
			else
				{
				Root_upper	= (LA_min+ LA_max*Root_min)/LA_max*Root_upper0 ;  // at leat 1% of Root area
				Root_middle	= (LA_min+ LA_max*Root_min)/LA_max*Root_middle0;
				Root_lower	= (LA_min+ LA_max*Root_min)/LA_max*Root_lower0;
				}
		}
		else if (T/24/3600/1000>LA_day1 && T/24/3600/1000<LA_day2) 
		{
			if (LA_min) 
				{
				Root_upper  = (LA_min+ (LA_max - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_upper0;
				Root_middle = (LA_min+ (LA_max - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_middle0;
				Root_lower  = (LA_min+ (LA_max - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_lower0;
				}
			else
				{
				Root_upper  = (LA_min+ LA_max*Root_min + (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_upper0;
				Root_middle = (LA_min+ LA_max*Root_min + (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_middle0;
				Root_lower  = (LA_min+ LA_max*Root_min + (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day1)/(LA_day2-LA_day1))/LA_max*Root_lower0;
				}				
		}
		else if (T/24/3600/1000>LA_day2 && T/24/3600/1000<LA_day3) 
			{
			Root_upper	= Root_upper0;
			Root_middle	= Root_middle0;
			Root_lower	= Root_lower0;
			}
		else if (T/24/3600/1000>LA_day3 && T/24/3600/1000<LA_day4) 
			{
			if (LA_Var==4) // case for wheat back to zero
				{
				Root_upper	=	(1-(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))*Root_upper0;
				Root_middle	=	(1-(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))*Root_middle0;
				Root_lower 	=	(1-(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))*Root_lower0;	
				}
			else
				{		
				if (LA_min) 
					{
					Root_upper	=	(LA_max- (LA_max - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_upper0;
					Root_middle	=	(LA_max- (LA_max - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_middle0;
					Root_lower 	=	(LA_max- (LA_max - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_lower0;	
					}
				else 
					{
					Root_upper	=	(LA_max- (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_upper0;
					Root_middle	=	(LA_max- (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_middle0;
					Root_lower 	=	(LA_max- (LA_max*(1-Root_min) - LA_min)*(T/24/3600/1000-LA_day3)/(LA_day4-LA_day3))/LA_max*Root_lower0;	
					}
				}
			}
			
		else if (T/24/3600/1000>=LA_day4) 
			{
			if (LA_min) 
				{
				Root_upper	= LA_min/LA_max*Root_upper0 ;
				Root_middle	= LA_min/LA_max*Root_middle0;
				Root_lower	= LA_min/LA_max*Root_lower0;
				}
			else
				{
				Root_upper	= (LA_min+ LA_max*Root_min)/LA_max*Root_upper0 ;
				Root_middle	= (LA_min+ LA_max*Root_min)/LA_max*Root_middle0;
				Root_lower	= (LA_min+ LA_max*Root_min)/LA_max*Root_lower0;
				}
			//if (LA_Var==4) Root_upper=Root_middle=Root_lower=0.1;		// case for wheat
			}
			
	
		Root_Area1= Root_Area * Root_upper;
		Root_Area2= Root_Area * Root_middle;
		Root_Area3= Root_Area * Root_lower;
		
		Q_old_s=Q_Root_Symp01;
		Q_Root_Symp01= Q_Root_Symp0* Root_upper;
		Q_Root_Symp1+=(Q_Root_Symp01-Q_old_s);
		
		Q_old_s=Q_Root_Symp02;
		Q_Root_Symp02= Q_Root_Symp0* Root_middle;
		Q_Root_Symp2+=(Q_Root_Symp02-Q_old_s);
		
		Q_old_s=Q_Root_Symp03;
		Q_Root_Symp03= Q_Root_Symp0* Root_lower;
		Q_Root_Symp3+=(Q_Root_Symp03-Q_old_s);
		
		Q_Root_Apo01= Q_Root_Apo_FR* Root_upper*1000*1000/18;
		Q_Root_Apo02= Q_Root_Apo_FR* Root_middle*1000*1000/18;
		Q_Root_Apo03= Q_Root_Apo_FR* Root_lower*1000*1000/18;
		Q_Root_Apo_t0=Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03;
		K_Root_Apo01=K_Root_Apo0 * Root_upper;
		K_Root_Apo02=K_Root_Apo0 * Root_middle;
		K_Root_Apo03=K_Root_Apo0 * Root_lower;
		K_Root_Symp11=K_Root_Symp1 * Root_upper;
		K_Root_Symp12=K_Root_Symp1 * Root_middle;
		K_Root_Symp13=K_Root_Symp1 * Root_lower;
	}
	
	if (LA_Var==5)	//phenology based on a temperature model
	{
		if (T_air>T_base1) // minimum Tair for T accumulation
			{
				dGDD=(T_air-T_base1)*dtt/24/3600; 
				GDD+=dGDD; //GDD in °C*day
			}
			if (GDD<S_GDD)  //before budbreak accumulating T above T_base until S_GDD
			{
				LA_max_Pheno=LA_min;
				for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno*Branch_distri[i]*(100-PLF[i])/100;
				T_budbreak=T/24/3600/1000;
				//Root_upper	= LA_max_Pheno/LA_max_init*Root_upper0 ;
				//Root_middle	= LA_max_Pheno/LA_max_init*Root_middle0;
				//Root_lower	= LA_max_Pheno/LA_max_init*Root_lower0;
			}
		else 									// after budbreak= Leaf growth proportional (LGE) to T-Tbase
			{
				if (T_air>T_base1) dGDD=(T_air-T_base1)*dtt/24/3600; else dGDD=0;
				if(LA_max_Pheno<=LA_max_init)// && T/24/3600<T_max_LAI+10) 
				{
					LA_max_Pheno+=dGDD*LGE*LA_max_init/100;
					T_max_LAI=T/24/3600/1000;
					LA_max=LA_max_Pheno;
					for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno*Branch_distri[i]*(100-PLF[i])/100;
				//	Root_upper	= LA_max_Pheno/LA_max_init*Root_upper0 ;
				//	Root_middle	= LA_max_Pheno/LA_max_init*Root_middle0;
				//	Root_lower	= LA_max_Pheno/LA_max_init*Root_lower0;
				}
				if (!LA_para4) // if LA_para=0 start to scenesce 10 days after T_max_LAI
				{
					if (T/24/3600/1000>(T_max_LAI+10) && T/24/3600/1000<=(T_max_LAI+60)) 
					{
					LA_max_Pheno2=LA_max - LA_max*(T/24/3600/1000-T_max_LAI-10)/(50);
					for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno2*Branch_distri[i]*(100-PLF[i])/100;
					}
					else if (T/24/3600/1000>(T_max_LAI+60)) 
						{
							LA_max_Pheno2=0;	
							for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno2*Branch_distri[i]*(100-PLF[i])/100;
						}
				}
				
				if (LA_para4) //start to scenesce at DOY=LA_para4 for 40 days
				{
					if (T/24/3600/1000>(LA_para4) && T/24/3600/1000<=(LA_para4+40)) 
					{
					LA_max_Pheno2=LA_max - LA_max*(T/24/3600/1000-LA_para4)/(40);
					for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno2*Branch_distri[i]*(100-PLF[i])/100;
					}
					else if (T/24/3600/1000>(LA_para4+40)) 
						{
							LA_max_Pheno2=0;	
							for (i=1;i<4;i++) Leaf_Area[i]=LA_max_Pheno2*Branch_distri[i]*(100-PLF[i])/100;
						}
				}
			}
			
		
		// for Roots assume Root_Area follow Leaf_Area		
	/*	Root_Area1= Root_Area * Root_upper;
		Root_Area2= Root_Area * Root_middle;
		Root_Area3= Root_Area * Root_lower;
	
		Q_old_s=Q_Root_Symp01;
		Q_Root_Symp01= Q_Root_Symp0* Root_upper;
		Q_Root_Symp1+=(Q_Root_Symp01-Q_old_s);
		
		Q_old_s=Q_Root_Symp02;
		Q_Root_Symp02= Q_Root_Symp0* Root_middle;
		Q_Root_Symp2+=(Q_Root_Symp02-Q_old_s);
	
		Q_old_s=Q_Root_Symp03;
		Q_Root_Symp03= Q_Root_Symp0* Root_lower;
		Q_Root_Symp3+=(Q_Root_Symp03-Q_old_s);
	
		Q_Root_Apo01= Q_Root_Apo_FR* Root_upper*1000*1000/18;
		Q_Root_Apo02= Q_Root_Apo_FR* Root_middle*1000*1000/18;
		Q_Root_Apo03= Q_Root_Apo_FR* Root_lower*1000*1000/18;
		Q_Root_Apo_t0=Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03;
		K_Root_Apo01=K_Root_Apo0 * Root_upper;
		K_Root_Apo02=K_Root_Apo0 * Root_middle;
		K_Root_Apo03=K_Root_Apo0 * Root_lower;
		K_Root_Symp11=K_Root_Symp1 * Root_upper;
		K_Root_Symp12=K_Root_Symp1 * Root_middle;
		K_Root_Symp13=K_Root_Symp1 * Root_lower;
		 */ 
	}
	

	if (LA_Var)
	{
		for (i=1;i<4;i++)
		{  
		if (LA_Var !=5)  Leaf_Area[i]=LA_max_Pheno*Branch_distri[i]*(100-PLF[i])/100;
		Q_old_s=Q_Leaf_Symp0[i];
		//if (Leaf_Area[i]<(LA_max_init*Branch_distri[i]/100)) 
		Q_Leaf_Symp0[i] = Succulence/1000 * (1- Leaf_Apo_fraction) * Leaf_Area[i]*1000*1000/18; 
		//else Q_Leaf_Symp0[i] = Succulence/1000 * (1- Leaf_Apo_fraction) * LA_max_init*Branch_distri[i]/100 *1000*1000/18; 
		Q_Leaf_Symp[i]+=(Q_Leaf_Symp0[i]- Q_old_s);
		
		LAI_Soil=Leaf_Area[i]/Surface_Soil;
		LAI_Crown=Leaf_Area[i]/Crown_Area;
		}
	if (LA_max_Pheno>LA0) Reserve-=(LA_max_Pheno-LA0)*LMA/180*6;;
	}
  
}

void Compute_T_Osmotic(void) // temperature dependance of osmotic potentials
{
	int i;
	if (T_OSMOTIC)
	{
		for (i=1;i<4;i++) Osmotic_TLeaf[i] =  (T_Leaf[i] + 273.16)/293.16;
		Osmotic_TAir  =  (T_air  + 273.16)/293.16;
		Osmotic_TSoil=   (T_Soil + 273.16)/293.16;
	}
	else
	{
		for (i=1;i<4;i++) Osmotic_TLeaf[i] = 1.0;
		Osmotic_TAir  = 1.0;
		Osmotic_TSoil = 1.0;
	}
}

void Compute_ST(void)  // temperature dependance of water surface tension
{
	if (SURFACE_TENSION)
	{
		ST_Leaf =  (75.6986 - (2.64569E-4) * T_Leaf[1] * T_Leaf[1] - 1.42361E-1 * T_Leaf[1] )/72.7455;
		ST_Air  =  (75.6986 - (2.64569E-4) * T_air  * T_air  - 1.42361E-1 * T_air  )/72.7455 ;
		ST_Soil =  (75.6986 - (2.64569E-4) * T_Soil * T_Soil - 1.42361E-1 * T_Soil )/72.7455 ;
	}
	else
	{
		ST_Leaf =1.0;
		ST_Air = 1.0;
		ST_Soil= 1.0;
	}
}
void Compute_THF(void)
{
	int i;
	
	for (i=1;i<4;i++)
		{
		if (PLC_Leaf_Apo[i]>PLC_END     && T_PLC_Leaf[i]>(T-T0))      T_PLC_Leaf[i]= T-T0;             // time to Leaf hydraulic failure
		if (PLC_Branch_Apo[i]>PLC_END   && T_PLC_Branch[i]>(T-T0))    T_PLC_Branch[i]= T-T0;           // time to Branch hydraulic failure
		if (Turgor_Leaf_Symp[i] > 0) T_TLP_Leaf[i]= T-T0;    										// time to leaf turgor loss
		if (RWC_Branch_s[i]>RWC_END && T_RWC_Branch[i]>(T-T0)) T_RWC_Branch[i]= T-T0;    										// time to reach critical branch symp value
		}
	//Test[0]=T_PLC_Leaf[1];
	//Test[1]=T;
	//Test[2]=T0;
	
	if (PLC_Trunk_Apo>PLC_END    	&& T_PLC_Trunk>(T-T0))     	T_PLC_Trunk= T-T0;                                          // time to Trunk hydraulic failure
	if ((PLC_Root_Apo1>PLC_END	&& PLC_Root_Apo2>PLC_END && PLC_Root_Apo3>PLC_END) && T_PLC_Root>(T-T0)) T_PLC_Root= T-T0;   // time to ALL Roots hydraulic failure
	if (PLC_Root_Apo1>PLC_END  	&& T_PLC_Root>(T-T0)) 		T_PLC_Root1= T-T0;
	if (PLC_Axil_Apo>PLC_END     	&& T_PLC_Axil>(T-T0))     	T_PLC_Axil= T-T0;
	if (Turgor_Axil_Symp > 0) T_TLP_Axil= T-T0;    	
}

void Compute_Cavitation(long double dtt)
{
	long double PLC=0.0,dQ,KP;
	int i;
   	
	Compute_ST();
	//LEAF
	for (i=1;i<4;i++)
	{
		if (Leaf_Area[i]) PLC=100/(1+exp(Slope_Leaf_Apo[i]/25*(P_Leaf_Apo[i]-P50_Leaf_Apo[i]*ST_Leaf)));
		else PLC=0;
		if (Regul_gs==6) Px_Leaf_Apo[i]= P50_Leaf_Apo[i]*ST_Leaf + 25/Slope_Leaf_Apo[i]*log((100-PLCx)/PLCx); 
		if (DYNAMIC0==3) if (PAR[i]==0) if (REW_t>REW_crit) PLC_LIMIT=0;
		if (DYNAMIC0==4) if (((PLC-PLC_Leaf_Apo[i])/dtt*1000000/10)<dPLC_crit) PLC_LIMIT=0;
		if (Leaf_Area[0]==0)	PLC_LIMIT=0;
		if (DYNAMIC0>=2) if (((PLC-PLC_Leaf_Apo[i])/dtt*1000000)>dPLC_crit) //hybrid mode; look for limit
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Leaf_Apo[i];
		}
		
		if ((REFILL && P_Leaf_Apo[i]>P_REFILL) || (PLC>PLC_Leaf_Apo[i]))
		{
			if (!New_Year) dQ=(PLC-PLC_Leaf_Apo[i])/100*Q_Leaf_Apo0[i]; else dQ=0;

			if (SYMP_CAVIT) // water goes only to the adjacent Symplasm
			{
				if (DYNAMIC && dQ>0) Q_Leaf_Symp[i]+=dQ; 
			}
			else if (DYNAMIC && dQ>0)     // water released is distributed to the connected Symplasmic reservoirs except for Refilling (dQ<0)
			{
			KP=K_Leaf_Symp[i]*P_Leaf_Symp[i] + K_Leaf_Apo[i]*P_Branch_Apo[i];
			if (KP)
				{
				Q_Leaf_Symp[i]+=     dQ*(K_Leaf_Symp[i]*P_Leaf_Symp[i])/KP;  //water released is distributed to the connected reservoirs according to the F=K*P ratio
				Q_Branch_Apo[i]+=    dQ*(K_Leaf_Apo[i]*P_Branch_Apo[i])/KP;
				}
			}
			else Q_Soil3+=  	dQ;                         // water released is distributed the deeper soil layer
			Q_Leaf_Apo1[i]-=   dQ;
			Q_Leaf_Apo[i]-=    dQ;
			PLC_Leaf_Apo[i]=   PLC;
		}
	}
	
	//Axil
	if (Type_Axil)
	{
	PLC=100/(1+exp(Slope_Axil_Apo/25*(P_Axil_Apo-P50_Axil_Apo*ST_Air)));
	if ((PLC-PLC_Axil_Apo)/dtt*1000000>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Axil_Apo>P_REFILL) || (PLC>PLC_Axil_Apo))
	{
		if (!New_Year)  dQ=(PLC-PLC_Axil_Apo)/100*Q_Axil_Apo0; else dQ=0;
		if (DYNAMIC && dQ>0)
		{
			KP=K_Axil_Symp*P_Axil_Symp  + N_Axil *K_Axil_Apo*P_Branch_Apo[1];                  // the water released by cavitation is distributed
			if (KP)
			{
				Q_Branch_Apo[1]+=      N_Axil * dQ*(K_Axil_Apo*P_Branch_Apo[1])/KP;                                            // according to the relative potential flows = K*pressure gradient
				if (Type_Axil>=2) Q_Petiole_Symp+=dQ*(K_Axil_Symp*P_Axil_Symp)/KP;                                     // into the connected apo and Symplasmic compartements
				else Q_Axil_Symp+=dQ*(K_Axil_Symp*P_Axil_Symp)/KP ; ;
			}  
		}
		else Q_Soil3+=dQ;
		Q_Axil_Apo-=dQ;
		Q_Axil_Apo1-=dQ;
		PLC_Axil_Apo=PLC;
	}
	}
	
	//BRANCH
	for (i=1;i<4;i++)
	{
		PLC=100/(1+exp(Slope_Branch_Apo[i]/25*(P_Branch_Apo[i]-P50_Branch_Apo[i]*ST_Air)));
		if (DYNAMIC0==4) if (((PLC-PLC_Branch_Apo[i])/dtt*1000000/10)<dPLC_crit) PLC_LIMIT=0;
	
		if (DYNAMIC0>=2) if ((PLC-PLC_Branch_Apo[i])/dtt*1000000>dPLC_crit)
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Branch_Apo[i];
		}
		if ((REFILL && P_Branch_Apo[i]>P_REFILL) || (PLC>PLC_Branch_Apo[i]))
		{
			if (!New_Year) dQ=(PLC-PLC_Branch_Apo[i])/100*Q_Branch_Apo0[i]; else dQ=0;
			if (SYMP_CAVIT && DYNAMIC && dQ>0)  Q_Branch_Symp[i]+=dQ;  
			else if (DYNAMIC && dQ>0)  // into the connected apo and Symplasmic compartements
			{                             
				if (Type_Axil) KP=K_Branch_Symp[i]*P_Branch_Symp[i] + K_Leaf_Apo[i]*P_Leaf_Apo[i]+ K_Axil_Apo*P_Axil_Apo+ K_Branch_Apo[i]*P_Trunk_Apo;   // the water released by cavitation is distributed
				else KP=K_Branch_Symp[i]*P_Branch_Symp[i] + K_Leaf_Apo[i]*P_Leaf_Apo[i]+K_Branch_Apo[i]*P_Trunk_Apo;   // the water released by cavitation is distributed
				if (KP )
				{
					Q_Branch_Symp[i]+=   dQ*K_Branch_Symp[i]*P_Branch_Symp[i]/KP ;                                     // into the connected apo and Symplasmic compartements
					Q_Leaf_Apo[i]+=      dQ*K_Leaf_Apo[i]*P_Leaf_Apo[i]/KP;                                            // according to the relative potential flows = K*pressure gradient
					Q_Trunk_Apo+=     dQ*K_Branch_Apo[i]*P_Trunk_Apo/KP;
					if (Type_Axil) 	Q_Axil_Apo+=      dQ*K_Axil_Apo*P_Axil_Apo/KP;
				}
			}
			else Q_Soil3+=dQ;
			Q_Branch_Apo[i]-=dQ;
			Q_Branch_Apo1[i]-=dQ;
			PLC_Branch_Apo[i]=PLC;
		}
	}
	//TRUNK
	PLC=100/(1+exp(Slope_Trunk_Apo/25*(P_Trunk_Apo-P50_Trunk_Apo*ST_Air)));
	if (DYNAMIC0==4) if (((PLC-PLC_Trunk_Apo)/dtt*1000000/10)<dPLC_crit) PLC_LIMIT=0;
	
	if (DYNAMIC0>=2) if ((PLC-PLC_Trunk_Apo)/dtt*1000000>dPLC_crit) 
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Trunk_Apo;
		}
	if ((REFILL && P_Trunk_Apo>P_REFILL) || (PLC>PLC_Trunk_Apo))
	{
		if (!New_Year) dQ=(PLC-PLC_Trunk_Apo)/100*Q_Trunk_Apo0; else dQ=0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0) Q_Trunk_Symp+=dQ;
		else if (DYNAMIC && dQ>0) 
		{
			KP=K_Trunk_Symp*P_Trunk_Symp + K_Branch_Apo[1]*P_Branch_Apo[1] +K_Branch_Apo[2]*P_Branch_Apo[2]+K_Branch_Apo[3]*P_Branch_Apo[3]+ K_Trunk_Apo*P_Root_Apo;
			if (KP)
			{
			Q_Trunk_Symp+=  	dQ*K_Trunk_Symp*P_Trunk_Symp/KP;
			Q_Branch_Apo[1]+= 	dQ*K_Branch_Apo[1]*P_Branch_Apo[1]/KP;
			Q_Branch_Apo[2]+= 	dQ*K_Branch_Apo[2]*P_Branch_Apo[1]/KP;
			Q_Branch_Apo[3]+= 	dQ*K_Branch_Apo[3]*P_Branch_Apo[1]/KP;
			Q_Root_Apo_t+= 	dQ*K_Trunk_Apo*P_Root_Apo/KP;
			Q_Root_Apo1+=     	Root_upper * dQ*K_Trunk_Apo*P_Root_Apo1/KP;
			Q_Root_Apo2+=		Root_middle * dQ*K_Trunk_Apo*P_Root_Apo2/KP;
			Q_Root_Apo3+=  	Root_lower * dQ*K_Trunk_Apo*P_Root_Apo3/KP;
			}
		}
		else Q_Soil3+=dQ;
		Q_Trunk_Apo-=dQ;
		Q_Trunk_Apo1-=dQ;
		PLC_Trunk_Apo=PLC;
	}
	
	//ROOTS
	PLC=100/(1+exp(Slope_Root_Apo1/25*(P_Root_Apo1-P50_Root_Apo1*ST_Soil)));
	if (DYNAMIC0>=2) if ((PLC-PLC_Root_Apo1)/dtt*1000000>dPLC_crit) 		
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Root_Apo1;
		}
	if ((REFILL && P_Root_Apo1>P_REFILL) || (PLC>PLC_Root_Apo1))
	{
		if (!New_Year) dQ=(PLC-PLC_Root_Apo1)/100*Q_Root_Apo01; else dQ=0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0) Q_Root_Symp1+=  dQ;
		else if (DYNAMIC && dQ>0)
		{
			KP=K_Root_Symp11*P_Root_Symp1 + K_Trunk_Apo*P_Trunk_Apo;
			if (KP) Q_Root_Symp1+=  dQ*K_Root_Symp11*P_Root_Symp1/KP;
			if (KP) Q_Trunk_Apo+=   dQ*K_Trunk_Apo*P_Trunk_Apo/KP;
		}
		else Q_Soil3+=dQ;
		Q_Root_Apo_t1-=dQ;
		Q_Root_Apo_t-=dQ;
		Q_Root_Apo1-=dQ;
		Q_Root_Apo11-=dQ;
		PLC_Root_Apo1=PLC;
	}
	
	PLC=100/(1+exp(Slope_Root_Apo2/25*(P_Root_Apo2-P50_Root_Apo2*ST_Soil)));
	if (DYNAMIC0>=2) if ((PLC-PLC_Root_Apo2)/dtt*1000000>dPLC_crit) 		
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Root_Apo2;
		}
	if ((REFILL && P_Root_Apo2>P_REFILL) || (PLC>PLC_Root_Apo2))
	{
		if (!New_Year) dQ=(PLC-PLC_Root_Apo2)/100*Q_Root_Apo02; else dQ=0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0) Q_Root_Symp2+=  dQ;
		else if (DYNAMIC && dQ>0)		
		{
		KP= K_Root_Symp12*P_Root_Symp2 + K_Trunk_Apo*P_Trunk_Apo;
		if (KP)
			{
			Q_Root_Symp2+=  dQ*K_Root_Symp12*P_Root_Symp2/KP;
			Q_Trunk_Apo+=  dQ*K_Trunk_Apo*P_Trunk_Apo/KP;
			}
		}
		else Q_Soil3+=dQ;
		Q_Root_Apo_t1-=dQ;
		Q_Root_Apo_t-=dQ;
		Q_Root_Apo2-=dQ;
		Q_Root_Apo12-=dQ;
		PLC_Root_Apo2=PLC;
	}
	
	PLC=100/(1+exp(Slope_Root_Apo3/25*(P_Root_Apo3-P50_Root_Apo3*ST_Soil)));
	if (DYNAMIC0>=2) if ((PLC-PLC_Root_Apo3)/dtt*1000000>dPLC_crit) 		
		{
			PLC_LIMIT=1;
			if (!DYNAMIC) PLC=PLC_Root_Apo3;
		}
	if ((REFILL && P_Root_Apo3>P_REFILL) || (PLC>PLC_Root_Apo3))
	{
		if (!New_Year) dQ=(PLC-PLC_Root_Apo3)/100*Q_Root_Apo03; else dQ=0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0) Q_Root_Symp3+=  dQ;
		else if (DYNAMIC && dQ>0)
		{
			KP= K_Root_Symp13*P_Root_Symp3 + K_Trunk_Apo*P_Trunk_Apo;
			if (KP)
			{
			Q_Root_Symp3+=  dQ*K_Root_Symp13*P_Root_Symp3/KP;
			Q_Trunk_Apo+=  dQ*K_Trunk_Apo*P_Trunk_Apo/KP;
			}
		}
		else Q_Soil3+=dQ;
		Q_Root_Apo_t1-=dQ;
		Q_Root_Apo_t-=dQ;
		Q_Root_Apo3-=dQ;
		Q_Root_Apo13-=dQ;
		PLC_Root_Apo3=PLC;
	}
	New_Year=0;	

}

void init_Cavitation(void)
{
	int i;
	
	Compute_ST();
	for (i=1;i<4;i++)
		{
		if (Leaf_Area[i]) PLC_Leaf_Apo[i]   = 100/(1+exp(Slope_Leaf_Apo[i]  /25*(P_Leaf_Apo[i]  -P50_Leaf_Apo[i]  * ST_Leaf)));  
		else PLC_Leaf_Apo[i]   =0;
		PLC_Branch_Apo[i] = 100/(1+exp(Slope_Branch_Apo[i]/25*(P_Branch_Apo[i]-P50_Branch_Apo[i]* ST_Air )));  
		}
	PLC_Axil_Apo = 100/(1+exp(Slope_Axil_Apo  /25*(P_Axil_Apo  -P50_Axil_Apo  * ST_Air ))); 
	PLC_Trunk_Apo = 100/(1+exp(Slope_Trunk_Apo /25*(P_Trunk_Apo -P50_Trunk_Apo * ST_Air )));  
	PLC_Root_Apo1 = 100/(1+exp(Slope_Root_Apo1 /25*(P_Root_Apo1 -P50_Root_Apo1 * ST_Soil))); 
	PLC_Root_Apo2 = 100/(1+exp(Slope_Root_Apo2 /25*(P_Root_Apo2 -P50_Root_Apo2 * ST_Soil))); 
	PLC_Root_Apo3 = 100/(1+exp(Slope_Root_Apo3 /25*(P_Root_Apo3 -P50_Root_Apo3 * ST_Soil)));
}

void Compute_Turgor_Ref(void)
{
	// compute the midday Leaf turgor for a well watered plant to use this value to compute the regulation of E by Leaf turgor;
	// assume no drougth no cavitation and an ETP=4mm and 25°C
	long double ETP_ref0 = 2;  // 2mm per day
	long double ETP_ref;     // ETP in mmol/s/m2 of Leaf area
	long double K_tott,K_Root;           // whole plant hydraulic LS conductance in mmol/s/m2/MPA
	long double LWP;         //Leaf Water Potential in MPa
	long double tlp, discri, Rs,Tp,Q=0;
	long double T1,T2,T3;
	
	
	T1=T_air;
	T2=T_Leaf[1];
	T3=T_Soil;
	T_air=25;
	T_Leaf[1]=25;
	T_Soil=20;
	//LAI=LA_max_init/Surface_Soil;   
	LAI_Soil=Leaf_Area[0]/Surface_Soil;// LAI in m2/m2
	LAI_Crown=Leaf_Area[0]/Crown_Area;
	if (Penman_Coeff==0) Penman_Coeff= -0.006*LAI_Crown*LAI_Crown +0.134*LAI_Crown + 0.036;
	ETP_ref= (ETP_ref0*1.6106-0.5616)*Surface_Soil*Penman_Coeff;  //ETP at midday in mmol/s/m2
	Fluidity();
	Compute_K();
	Compute_T_Osmotic();
	K_Root=1/(1/K_Root_Symp11+1/K_Root_Apo1) + 1/(1/K_Root_Symp12+1/K_Root_Apo2) + 1/(1/K_Root_Symp13+1/K_Root_Apo3);	
	K_tott= Fluidity_air*(1/(1/(K_Leaf_Apo[1]+K_Leaf_Apo[2]+K_Leaf_Apo[3]) +1/(K_Leaf_Symp[1]+K_Leaf_Symp[2]+K_Leaf_Symp[3]) + 1/(K_Branch_Apo[1]+ K_Branch_Apo[2]+ K_Branch_Apo[3])+ 1/K_Trunk_Apo + 1/K_Root));

	LWP=-ETP_ref/(K_tott);
	if (GRAVITY==1) LWP-=(Length_Trunk+Length_Branch)*0.01;
	
	
	tlp=Pi0_Leaf_Symp*Osmotic_TLeaf[1]*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp*Osmotic_TLeaf[1]+Epsilon_Leaf_Symp);
	if (LWP>tlp)
	{
		discri=pow(Pi0_Leaf_Symp*Osmotic_TLeaf[1] + Epsilon_Leaf_Symp + LWP,2)-4*Pi0_Leaf_Symp*Osmotic_TLeaf[1]*Epsilon_Leaf_Symp;
		Q = ((Pi0_Leaf_Symp*Osmotic_TLeaf[1] + Epsilon_Leaf_Symp + LWP) + pow(discri,0.5))/(2*Epsilon_Leaf_Symp/Q_Leaf_Symp0[1]);
	}
	else   Turgor_Leaf_Symp_Ref=0;
	
	Rs=(Q_Leaf_Symp0[1]-Q)/Q_Leaf_Symp0[1];
	Tp=-Pi0_Leaf_Symp*Osmotic_TLeaf[1] - Epsilon_Leaf_Symp*Rs;
	if (Tp<0)Tp=0;
	Turgor_Leaf_Symp_Ref =Tp;
	if (PRINT_SCREEN) printf("LWP=%Lf TLP=%Lf Turgor_ref=%Lf", LWP, tlp, Tp);
	T_air=T1;           //restore previous values
	T_Leaf[1]=T2;
	T_Soil=T3;
		
}

void Get_DATA(long double dt_court) 
{
	long double Radius_Trunk0,Radius_bark_Trunk0, E_Plant,  K_Root, K_tot2, Diam_Fruit, Specific_Leaf_WC[4],Specific_Br_WC[4], Specific_Br_Symp_WC[4], Specifc_Br_Apo_WC[4], Volume_Branch[4],Sap_Flow_d,Sap_Flow_d2,Axil_WC;
	long double K_Rhizo,K_Rhizo_20,K_Root1, K_Root2, K_Root3, P_Soil;
	long double RWC_Leaf[4], RWC_Leaf_s[4],RWC_Leaf_a[4],RWC_Branch[4],RWC_Branch_a[4],RWC_Trunk,RWC_Trunk_s,RWC_Trunk_a,RWC_Root,RWC_Root_s,RWC_Root_a,RWC_Plant,RWC_Plant_s,RWC_Plant_a;
	long double E_tot_Leaf=0, E_tot_Branch=0,Q_Leaf_Symp_tot=0,Q_Branch_Symp_tot=0,Q_Leaf_Apo_tot=0,Q_Branch_Apo_tot=0,Q_Leaf_Apo00=0,Q_Leaf_Symp00=0,Q_Branch_Apo00=0,Q_Branch_Symp00=0;
	
	int i;  
	
	//for (i=1;i<4;i++)Test[i]=K_Branch_Apo[i];
	for (i=1;i<4;i++) 
	{
		if (FRACTAL==1)     Volume_Branch[i]= Branch_distri[i]*Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4;      //morphometric                                           // Branch volume in m3
		else                Volume_Branch[i]=  Branch_distri[i]*Q_Branch_Apo_FR/Branch_Apo_fraction;                                                                                           // if Fractal then sapwood volume in m3
		if (Leaf_Area[i])      Specific_Leaf_WC[i]= (Q_Leaf_Symp[i]   + Q_Leaf_Apo0[i]   *(100-PLC_Leaf_Apo[i])/100) /1000*18 / (LMA * Leaf_Area[i]); else  Specific_Leaf_WC[i]=0;    // g H20 per g Leaf dry mass
		if (Volume_Branch[i])
		{
			Specific_Br_WC[i]=     (Q_Branch_Symp[i] + Q_Branch_Apo0[i] *(100-PLC_Branch_Apo[i])/100) /1000/1000*18 / (Density * Volume_Branch[i]);                            // kg H20 per kg total Branch dry mass
			Specific_Br_Symp_WC[i]= Q_Branch_Symp[i]/1000/1000*18 / (Density * Volume_Branch[i]* Branch_Symp_fraction);                                                  // kg H20 per kg Symplasm Branch dry mass
			Specifc_Br_Apo_WC[i]=  (Q_Branch_Apo0[i] *(100-PLC_Branch_Apo[i])/100) /1000/1000*18 / (Density * Volume_Branch[i]* Branch_Apo_fraction);                       // kg H20 per kg apoplasm Branch dry mass
		}
		else Specific_Br_WC[i]= Specific_Br_Symp_WC[i]=Specifc_Br_Apo_WC[i]=0;
	}
	
	Radius_bark_Trunk0= (Diam_Trunk -sqrt((Diam_Trunk*Diam_Trunk -4*bark*Q_Trunk_Symp0*18/1000/1000/1000/Length_Trunk/3.1416)))*1000/2; //equivalent initial bark thickness in mm
	Radius_bark_Trunk = (Diam_Trunk -sqrt((Diam_Trunk*Diam_Trunk -4*bark*Growth_Trunk*18/1000/1000/1000/Length_Trunk/3.1416)))*1000/2; //equivalent  bark thickness in mm
	Radius_bark_Trunk_rel= Radius_bark_Trunk-Radius_bark_Trunk0;
	
	Radius_Trunk0=  Diam_Trunk/2*1000;  //equivalent initial bark thickness in mm
	Radius_Trunk =  Radius_Trunk0 + Radius_bark_Trunk_rel;
	Radius_Trunk_rel= Radius_bark_Trunk_rel ;  // will be different when temperature dilatation is added.
	
	
	Petiole_diam =      1000*1000*pow(4/3.1416*Q_Petiole_Symp*18/1000/1000/Length_Petiole/Branch_Symp_fraction/1000,0.5);                                                                    //in µm
	if (N_Axil) Diam_Fruit   =      pow(Q_Axil_Symp/1000/1000/1000*18*6/N_Axil/(WC_Axil*3.1416),0.333333333333333); else Diam_Fruit =0;
	if (N_Axil) Growth_Fruit =      pow(Q_Axil_Symp0/1000/1000/1000*18*6/N_Axil/(WC_Axil*3.1416),0.333333333333333); else  Growth_Fruit =0;
	for (i=1;i<4;i++)
	{
		E_tot_Leaf+=E_Leaf[i]*Leaf_Area[i];
		E_tot_Branch+=E_Branch[i]*Branch_Area[i];
		Q_Leaf_Symp_tot+=Q_Leaf_Symp[i];
		Q_Branch_Symp_tot+=Q_Branch_Symp[i];
		Q_Leaf_Apo_tot+=Q_Leaf_Apo0[i]*(100-PLC_Leaf_Apo[i])/100;
		Q_Branch_Apo_tot+=Q_Branch_Apo[i]*(100-PLC_Branch_Apo[i])/100;
	}
	if      (CUT==1)    E_Plant= (E_tot_Leaf + E_tot_Branch + E_Axil*Axil_Area + E_Petiole*Petiole_area);                                     //cut base of Branch
	else if (CUT==2)    E_Plant= (E_tot_Leaf + E_tot_Branch + E_Axil*Axil_Area + E_Petiole*Petiole_area + E_Trunk*Trunk_Area     );           //cut base of Trunk
	else                E_Plant= (E_tot_Leaf + E_tot_Branch + E_Axil*Axil_Area + E_Petiole*Petiole_area + E_Trunk*Trunk_Area + E_Root1*Root_Area1 + E_Root2*Root_Area2 + E_Root3*Root_Area3 ); //total plant water loss in mmol/s
	EvapoT =            E_Plant + E_Soil*Surface_Soil;                                                                                             //total plant + soil water loss in mmol/s
	if      (CUT==1)
	{
		Q_Plant_s=(Q_Axil_Symp + Q_Petiole_Symp + Q_Leaf_Symp_tot  + Q_Branch_Symp_tot)*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo_tot + Q_Axil_Apo0*(100-PLC_Axil_Apo)/100 + Q_Branch_Apo_tot)*18/1000/1000 ;
	}
	else if (CUT==2)
	{
		Q_Plant_s=(Q_Axil_Symp +Q_Petiole_Symp + Q_Leaf_Symp_tot  + Q_Branch_Symp_tot + Q_Trunk_Symp )*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo_tot + Q_Axil_Apo0*(100-PLC_Axil_Apo)/100 + Q_Branch_Apo_tot + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100)*18/1000/1000 ;
	}
	else
	{
		Q_Plant_s=(Q_Axil_Symp +Q_Petiole_Symp + Q_Leaf_Symp_tot  + Q_Branch_Symp_tot + Q_Trunk_Symp +Q_Root_Symp1 +Q_Root_Symp2 +Q_Root_Symp3 )*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo_tot + Q_Axil_Apo0*(100-PLC_Axil_Apo)/100 + Q_Branch_Apo_tot + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100)*18/1000/1000 ;
	}
	Q_Plant= Q_Plant_a + Q_Plant_s;
	
	if (K_Soil1 && K_Interface1 && K_Root_Symp11 && K_Root_Apo1) K_Root1=1/(1/K_Soil1 + 1/K_Interface1 + 1/K_Root_Symp11 + 1/K_Root_Apo1); else K_Root1=0;
	if (K_Soil2 && K_Interface2 && K_Root_Symp12 && K_Root_Apo2) K_Root2=1/(1/K_Soil2 + 1/K_Interface2 + 1/K_Root_Symp12 + 1/K_Root_Apo2); else K_Root2=0;
	if (K_Soil3 && K_Interface3 && K_Root_Symp13 && K_Root_Apo3) K_Root3=1/(1/K_Soil3 + 1/K_Interface3 + 1/K_Root_Symp13 + 1/K_Root_Apo3); else K_Root3=0;
	K_Rhizo=K_Root1+K_Root2+K_Root3;
	K_Rhizo_20=(K_Root1/Fluidity_soil1+K_Root2/Fluidity_soil2+K_Root3/Fluidity_soil3);
	
	if (K_Root_Symp11 && K_Root_Apo1) K_Root1=1/(1/K_Root_Symp11 + 1/K_Root_Apo1); else K_Root1=0;
	if (K_Root_Symp12 && K_Root_Apo2) K_Root2=1/(1/K_Root_Symp12 + 1/K_Root_Apo2); else K_Root2=0;
	if (K_Root_Symp13 && K_Root_Apo3) K_Root3=1/(1/K_Root_Symp13 + 1/K_Root_Apo3); else K_Root3=0;
	K_Root=K_Root1+K_Root2+K_Root3;

	K_Leaf_Apo[0]=K_Leaf_Symp[0]=K_Branch_Apo[0]=0;	
	for (i=1;i<4;i++) // conductances in parallel
	{
		K_Leaf_Apo[0]+=K_Leaf_Apo[i]*Fluidity_Leaf[i];
		K_Leaf_Symp[0]+=K_Leaf_Symp[i]*Fluidity_Leaf[i];
		K_Branch_Apo[0]+=K_Branch_Apo[i]*Fluidity_Branch[i];
	}
	if (Leaf_Area[0] && K_Rhizo && K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo)  K_tot_20  =  1/(1/K_Leaf_Apo[0]+1/K_Leaf_Symp[0]+ 1/K_Branch_Apo[0]+ 1/K_Trunk_Apo*Fluidity_air + 1/K_Rhizo_20);  
		else K_tot_20=0;
	if (K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo && K_Root)            		K_Plant_20=  1/(1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0] + 1/K_Branch_Apo[0] + 1/K_Trunk_Apo*Fluidity_air + 1/K_Root*Fluidity_soil);  
		else	K_Plant_20=0;
		
	if (K_tot_20_0) PLC_Plant_Soil= 100*(1-K_tot_20/K_tot_20_0); else PLC_Plant_Soil=0;
	if (K_Plant_20_0) PLC_Plant=100*(1-K_Plant_20/K_Plant_20_0); else PLC_Plant=0;	
	
	K_Leaf_Apo[0]=K_Leaf_Symp[0]=K_Branch_Apo[0]=0;
	Leaf_Area[0]=P_Leaf_Evap[0]=0;
	for (i=1;i<4;i++) // conductances in parallel
	{
		K_Leaf_Apo[0]+=K_Leaf_Apo[i];
		K_Leaf_Symp[0]+=K_Leaf_Symp[i];
		K_Branch_Apo[0]+=K_Branch_Apo[i];
		Leaf_Area[0]+=Leaf_Area[i];
	}	
	P_Leaf_Evap[0]= 		P_Leaf_Evap[1]*		Branch_distri[1] +P_Leaf_Evap[2]*		Branch_distri[2]+	P_Leaf_Evap[3]*		Branch_distri[3];
	P_Leaf_Symp[0]= 		P_Leaf_Symp[1]*		Branch_distri[1] +P_Leaf_Symp[2]*		Branch_distri[2]+	P_Leaf_Symp[3]*		Branch_distri[3];
	P_Leaf_Apo[0] = 		P_Leaf_Apo[1]*		Branch_distri[1] +P_Leaf_Apo[2]*		Branch_distri[2]+	P_Leaf_Apo[3]*		Branch_distri[3];
	VPD_Leaf[0]=			VPD_Leaf[1]*			Branch_distri[1] +VPD_Leaf[2]*		Branch_distri[2]+	VPD_Leaf[3]*			Branch_distri[3];
	VPD_Cuti[0]=			VPD_Cuti[1]*			Branch_distri[1] +VPD_Cuti[2]*		Branch_distri[2]+	VPD_Cuti[3]*			Branch_distri[3];
	T_Leaf[0]=			T_Leaf[1]*			Branch_distri[1] +T_Leaf[2]*			Branch_distri[2]+	T_Leaf[3]*			Branch_distri[3];
	VPD_Branch[0]=		VPD_Branch[1]*		Branch_distri[1] +VPD_Branch[2]*		Branch_distri[2]+	VPD_Branch[3]*		Branch_distri[3];
	E_Leaf[0]=			E_Leaf[1]*			Branch_distri[1] +E_Leaf[2]*			Branch_distri[2]+	E_Leaf[3]*			Branch_distri[3];
	E_cuti[0]=			E_cuti[1]*			Branch_distri[1] +E_cuti[2]*			Branch_distri[2]+	E_cuti[3]*			Branch_distri[3];
	E_Branch[0]=			E_Branch[1]*			Branch_distri[1] +E_Branch[2]*		Branch_distri[2]+	E_Branch[3]*			Branch_distri[3];
	g_s[0]=				g_s[1]*				Branch_distri[1] +g_s[2]*				Branch_distri[2]+	g_s[3]*				Branch_distri[3];
	g_cuti[0]=			g_cuti[1]*			Branch_distri[1] +g_cuti[2]*			Branch_distri[2]+	g_cuti[3]*			Branch_distri[3];
	Turgor_Leaf_Symp[0]=	Turgor_Leaf_Symp[1]*	Branch_distri[1] +Turgor_Leaf_Symp[2]*	Branch_distri[2]+	Turgor_Leaf_Symp[3]*	Branch_distri[3];
	P_Branch_Symp[0]=	P_Branch_Symp[1]*	Branch_distri[1] +P_Branch_Symp[2]*	Branch_distri[2]+	P_Branch_Symp[3]*	Branch_distri[3];
	P_Branch_Apo[0]=		P_Branch_Apo[1]*		Branch_distri[1] +P_Branch_Apo[2]*		Branch_distri[2]+	P_Branch_Apo[3]*		Branch_distri[3];
	PLC_Leaf_Apo[0]=		PLC_Leaf_Apo[1]*		Branch_distri[1] +PLC_Leaf_Apo[2]*		Branch_distri[2]+	PLC_Leaf_Apo[3]*		Branch_distri[3];
	PLC_Branch_Apo[0]=	PLC_Branch_Apo[1]*	Branch_distri[1] +PLC_Branch_Apo[2]*	Branch_distri[2]+	PLC_Branch_Apo[3]*	Branch_distri[3];
	E_clim[0]=			E_clim[1]*			Branch_distri[1] +E_clim[2]*			Branch_distri[2]+	E_clim[3]*			Branch_distri[3];
	A_net[0] = 			A_net[1]*			Branch_distri[1] +A_net[2]*			Branch_distri[2]+	A_net[3]*			Branch_distri[3];
	
	if (Leaf_Area[0] && K_Rhizo && K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo)  K_tot  =  1/(1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0] + 1/K_Branch_Apo[0] + 1/K_Trunk_Apo + 1/K_Rhizo);  
		else K_tot=0;
	if (K_Leaf_Apo[0] && K_Leaf_Symp[0] && K_Branch_Apo[0] && K_Trunk_Apo && K_Root)            		K_Plant=  1/(1/K_Leaf_Apo[0] +1/K_Leaf_Symp[0] + 1/K_Branch_Apo[0] + 1/K_Trunk_Apo + 1/K_Root);             
		else K_Plant=0;
	if (K_Rhizo)   P_Soil= (P_Soil1*K_Root1 + P_Soil2*K_Root2 + P_Soil3*K_Root3)/K_Rhizo; else P_Soil=-10;
	if (Leaf_Area[0]) K_tot2= E_Plant/(P_Soil - P_Leaf_Evap[0] + Pg_Leaf[1] ); else K_tot2=0;
	
	for (i=1;i<4;i++) // conductances in parallel
		{
		if (Leaf_Area[i])
			{
			RWC_Leaf[i]= (Q_Leaf_Apo0[i]*(100-PLC_Leaf_Apo[i])/100+Q_Leaf_Symp[i])/(Q_Leaf_Apo0[i] + Q_Leaf_Symp0[i]);
			RWC_Leaf_s[i]= Q_Leaf_Symp[i]/Q_Leaf_Symp0[i];
			RWC_Leaf_a[i]= Q_Leaf_Apo0[i]*(100-PLC_Leaf_Apo[i])/Q_Leaf_Apo0[i]/100;
			}
		else RWC_Leaf[i]=RWC_Leaf_a[i]=RWC_Leaf_s[i]=1;
		if (Branch_distri[i])
			{
			RWC_Branch[i]= (Q_Branch_Apo0[i]*(100-PLC_Branch_Apo[i])/100 + Q_Branch_Symp[i])/(Q_Branch_Apo0[i] + Q_Branch_Symp0[i]);
			RWC_Branch_s[i]= Q_Branch_Symp[i] / Q_Branch_Symp0[i];
			RWC_Branch_a[i]= Q_Branch_Apo0[i]*(100-PLC_Branch_Apo[i])/Q_Branch_Apo0[i]/100;
			}
			else RWC_Branch[i]=RWC_Branch_s[i]=RWC_Branch_a[i]=1;
		}
		
	
	RWC_Trunk= (Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Trunk_Symp)/(Q_Trunk_Apo0 + Q_Trunk_Symp0);
	RWC_Trunk_s= Q_Trunk_Symp/ Q_Trunk_Symp0;
	RWC_Trunk_a= Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/Q_Trunk_Apo0/100;
	RWC_Root=(Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 )/(Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03+Q_Root_Symp01+Q_Root_Symp02+Q_Root_Symp03);
	RWC_Root_s= (Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3)/(Q_Root_Symp01+Q_Root_Symp02+Q_Root_Symp03);
	RWC_Root_a= (Q_Root_Apo01*(100-PLC_Root_Apo1) + Q_Root_Apo02*(100-PLC_Root_Apo2) + Q_Root_Apo03*(100-PLC_Root_Apo3))/(Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03)/100;
	
	Q_Leaf_Apo0[0]=Q_Leaf_Symp[0]=Q_Branch_Apo0[0]=Q_Branch_Symp[0]=0;
	for (i=1;i<4;i++) // conductances in parallel
		{
		Q_Leaf_Apo0[0]+=Q_Leaf_Apo0[i]*(100-PLC_Leaf_Apo[i])/100;
		Q_Leaf_Apo00+=Q_Leaf_Apo0[i];
		Q_Leaf_Symp[0]+=Q_Leaf_Symp[i];
		Q_Leaf_Symp00+=Q_Leaf_Symp0[i];
		Q_Branch_Apo0[0]+=Q_Branch_Apo0[i]*(100-PLC_Branch_Apo[i])/100;
		Q_Branch_Apo00+=Q_Branch_Apo0[i];
		Q_Branch_Symp[0]+=Q_Branch_Symp[i];
		Q_Branch_Symp00+=Q_Branch_Symp0[i];
		}	
	RWC_Plant= (Q_Leaf_Apo0[0]+ Q_Leaf_Symp[0] + Q_Branch_Apo0[0] + Q_Branch_Symp[0] + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Trunk_Symp + Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp + Q_Axil_Apo*(100-PLC_Axil_Apo)/100)/(Q_Leaf_Apo00 +Q_Leaf_Symp00 +Q_Branch_Apo00 + Q_Branch_Symp00 + Q_Trunk_Apo0 + Q_Trunk_Symp0 + Q_Root_Apo01 + Q_Root_Apo02 + Q_Root_Apo03 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp + Q_Axil_Apo);
	RWC_Plant_s=(Q_Leaf_Symp[0]  + Q_Branch_Symp[0] + Q_Trunk_Symp + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp)/(Q_Leaf_Symp00 + Q_Branch_Symp00 +Q_Trunk_Symp0 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp);
	RWC_Plant_a=(Q_Leaf_Apo0[0]+Q_Branch_Apo0[0] + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Axil_Apo*(100-PLC_Axil_Apo)/100)/(Q_Leaf_Apo00 +Q_Branch_Apo00 +  Q_Trunk_Apo0 +  Q_Root_Apo01 + Q_Root_Apo02 + Q_Root_Apo03 + Q_Axil_Apo);
	RWC_shoot=  (Q_Leaf_Symp[0] + Q_Leaf_Apo0[0] + Q_Branch_Symp[0] + Q_Branch_Apo0[0]) / (Q_Leaf_Symp00 + Q_Leaf_Apo00   + Q_Branch_Symp00 + Q_Branch_Apo00);
	Teta_Soil=(Q_Soil1+Q_Soil2+Q_Soil3)/1000/1000/1000*18/Volume_soil;
	Teta_Soil1=(Q_Soil1)/1000/1000/1000*18/Volume_soil1;
	Teta_Soil2=(Q_Soil2)/1000/1000/1000*18/Volume_soil2;
	Teta_Soil3=(Q_Soil3)/1000/1000/1000*18/Volume_soil3;
	GPP=A_net_tot/Surface_Soil*12 ; // gC/m2/year
	if (Type_Axil)
	{
		RWC_Axil=(Q_Axil_Symp) / (Q_Axil_Symp0);
		Axil_WC=(Q_Axil_Symp*18/1000)/(Q_Axil_Symp*18/1000+0.125*3.1416/6*pow(Diam_Axil,3)*1000*1000); //organ water content per dry mass with 0.125g dry mass per cm3
	}
	else
	{
		RWC_Axil=0;
		Axil_WC=0;
	}
	dq_Trunk_Apo_Branch_Apo[0]=0;
	for (i=1;i<4;i++) dq_Trunk_Apo_Branch_Apo[0]+=dq_Trunk_Apo_Branch_Apo[i];
	if (DYNAMIC)    Sap_Flow_d=(dq_Trunk_Apo_Branch_Apo[0]+dq_Trunk_Apo_Trunk_Symp)/1000*18/dt_court/(Diam_Trunk * Diam_Trunk * 3.1416 / 4 *(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction); //sap flow density in the Trunk, in g/s/m2 of sapwood
	else            Sap_Flow_d=(E_Plant)/1000*18/(Diam_Trunk * Diam_Trunk * 3.1416 / 4 *(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction);
					 Sap_Flow_d2=Sap_Flow_d/100000*3600; //sap flow density in the Trunk, in dm3/dm2/h of sapwood
	if (Sap_Flow_d2<SF_min_d) SF_min_d=Sap_Flow_d2;
	if (Sap_Flow_d2>SF_max_d) SF_max_d=Sap_Flow_d2;
	//if (CLIMAT==0 || CLIMAT==1 || CLIMAT==3 || CLIMAT==5) DATA[0]= T/3600/24; else 
		
		
	



	if (t_out==24 && TRANSIENT==2) DATA[0]= T/1000/3600/24-1;
	else DATA[0]= T/3600/24/1000;
	if (File_out[0]==2 || Screen_out[0]==2) DATA[0]=T/1000/3600/24-DOY_0;
	DATA[1]=  T_air;
	DATA[2]=  T_Leaf[0];
	DATA[3]=  RH_air;
	DATA[4]=  PAR[0];
	DATA[5]=  VPD_Leaf[0];
	DATA[6]=  VPD_Cuti[0];
	DATA[7]=  VPD_Axil;
	DATA[8]=  VPD_Branch[0];
	DATA[9]=  VPD_Trunk;
	DATA[10]= VPD_Root1;
	DATA[11]= VPD_Root2;
	DATA[12]= VPD_Root3;
	DATA[13]= VPD_Soil;
	DATA[14]= Pot_PAR_tot;
	DATA[15]= P_air;
	DATA[16]= E_clim[0];
	if (File_out[17]==1 || Screen_out[17]==1) DATA[17]= E_Leaf[0];
	if (File_out[17]==2 || Screen_out[17]==2) DATA[17]= E_Leaf[0]*Leaf_Area[0];
	if (File_out[18]==1 || Screen_out[18]==1) DATA[18]= E_cuti[0];
	if (File_out[18]==2 || Screen_out[18]==2) DATA[18]= E_cuti[0]*Leaf_Area[0];
	if (File_out[19]==1 || Screen_out[19]==1) DATA[19]= E_Branch[0];
	if (File_out[19]==2 || Screen_out[19]==2) DATA[19]= E_Branch[0]*Branch_Area[0];
	if (File_out[20]==1 || Screen_out[20]==1) DATA[20]= E_Trunk;
	if (File_out[20]==2 || Screen_out[20]==2) DATA[20]= E_Trunk*Trunk_Area;
	if (File_out[21]==1 || Screen_out[21]==1) DATA[21]= E_Axil;
	if (File_out[21]==2 || Screen_out[21]==2) DATA[21]= E_Axil*Axil_Area;
	DATA[22]= E_Root1;
	DATA[23]= E_Root2;
	DATA[24]= E_Root3;
	DATA[25]= E_Soil;
	if (File_out[26]==1 || Screen_out[26]==1) DATA[26]= E_Plant;					// in mmol/s
	if (File_out[26]==2 || Screen_out[26]==2) DATA[26]= E_Plant*18/1000*3600; 	// in g/h
	if (Leaf_Area[0]) DATA[27]= E_Plant/Leaf_Area[0]; else DATA[27]=0;			
	DATA[28]= E_tot*18/1000/1000;
	DATA[29]= g_Canopy;	
	DATA[30]= g_s[0];
	DATA[31]= g_cuti[0];
	DATA[32]= P_Leaf_Evap[0];
	DATA[33]=P_Leaf_Symp[0];
	DATA[34]= P_Leaf_Apo[0];
	DATA[35]= turgor;
	DATA[36]= Turgor_Leaf_Symp[0];
	DATA[37]= P_Axil_Symp;
	DATA[38]= P_Axil_Apo;
	DATA[39]= P_Branch_Symp[0];
	DATA[40]= P_Branch_Apo[0];
	DATA[41]= P_Trunk_Symp;
	DATA[42]= P_Trunk_Apo;
	DATA[43]= Turgor_Trunk_Symp;
	DATA[44]= P_Root_Symp1;
	DATA[45]= P_Root_Endo1;
	DATA[46]= P_Root_Apo1;
	DATA[47]= P_Root_Symp2;
	DATA[48]= P_Root_Endo2;
	DATA[49]= P_Root_Apo2;
	DATA[50]= P_Root_Symp3;
	DATA[51]= P_Root_Endo3;
	DATA[52]= P_Root_Apo3;
	DATA[53]= P_Soil1;
	DATA[54]= P_Soil2;
	DATA[55]= P_Soil3;
	DATA[56]= RWC1;
	DATA[57]= RWC2;
	DATA[58]= RWC3;
	DATA[59]= K_Soil1;
	DATA[60]= K_Soil2;
	DATA[61]= K_Soil3;
	DATA[62]= K_Interface1;
	DATA[63]= K_Interface2;
	DATA[64]= K_Interface3;
	if (File_out[65]==2 || Screen_out[65]==2) DATA[65]= K_Leaf_Symp[1] +K_Leaf_Symp[2]+K_Leaf_Symp[3];
	if (File_out[65]==1 || Screen_out[65]==1) 
	{
		if (Leaf_Area[0]) DATA[65]= (K_Leaf_Symp[1]+K_Leaf_Symp[2]+K_Leaf_Symp[3])/Leaf_Area[0]; 
		else DATA[65]=0;
	}
	if (File_out[66]==2 || Screen_out[66]==2) DATA[66]= K_Leaf_Apo[1] +K_Leaf_Apo[2]+K_Leaf_Apo[3];
	if (File_out[66]==1 || Screen_out[66]==1) 
	{
		if (Leaf_Area[0]) DATA[66]= (K_Leaf_Apo[1]+K_Leaf_Apo[2]+K_Leaf_Apo[3])/Leaf_Area[0]; 
		else DATA[66]=0;
	}
	if (File_out[67]==2 || Screen_out[67]==2) DATA[67]= 1/(1/(K_Leaf_Symp[1] +K_Leaf_Symp[2]+K_Leaf_Symp[3])+1/(K_Leaf_Apo[1]+K_Leaf_Apo[2]+K_Leaf_Apo[3]));
	if (File_out[67]==1 || Screen_out[67]==1)
	{
		if (Leaf_Area[0]) DATA[67]= 1/(1/(K_Leaf_Apo[1]+K_Leaf_Apo[2]+K_Leaf_Apo[3])+1/(K_Leaf_Symp[1]+K_Leaf_Symp[2]+K_Leaf_Symp[3]))/Leaf_Area[0]; 
		else DATA[67]=0;
	}
	if (File_out[68]==2 || Screen_out[68]==2) DATA[68]= K_Root_Symp11;
	if (File_out[68]==1 || Screen_out[68]==1)
	{
		if (Leaf_Area[0]) DATA[68]=K_Root_Symp11/Leaf_Area[0];
		else DATA[68]=0;
	}
	
	if (File_out[69]==2 || Screen_out[69]==2) DATA[69]= K_Root_Apo1;
	if (File_out[69]==1 || Screen_out[69]==1)
	{
		if (Leaf_Area[0]) DATA[69]=K_Root_Apo1/Leaf_Area[0];
		else DATA[69]=0;
	}
	if (File_out[70]==2 || Screen_out[70]==2) DATA[70]= K_Root_Symp12;
	if (File_out[70]==1 || Screen_out[70]==1)
	{
		if (Leaf_Area[0]) DATA[70]=K_Root_Symp12/Leaf_Area[0];
		else DATA[70]=0;
	}
	if (File_out[71]==2 || Screen_out[71]==2) DATA[71]= K_Root_Apo2;
	if (File_out[71]==1 || Screen_out[71]==1)
	{
		if (Leaf_Area[0]) DATA[71]=K_Root_Apo2/Leaf_Area[0];
		else DATA[71]=0;
	}
	if (File_out[72]==2 || Screen_out[72]==2) DATA[72]= K_Root_Symp13;
	if (File_out[72]==1 || Screen_out[72]==1)
	{
		if (Leaf_Area[0]) DATA[72]=K_Root_Symp13/Leaf_Area[0];
		else DATA[72]=0;
	}
	if (File_out[73]==2 || Screen_out[73]==2) DATA[73]= K_Root_Apo3;
	if (File_out[73]==1 || Screen_out[73]==1)
	{
		if (Leaf_Area[0]) DATA[73]=K_Root_Apo3/Leaf_Area[0];
		else DATA[73]=0;
	}
	if (File_out[74]==2 || Screen_out[74]==2) DATA[74]= K_Root;
	if (File_out[74]==1 || Screen_out[74]==1)
	{
		if (Leaf_Area[0]) DATA[74]= K_Root/Leaf_Area[0]; 
		else DATA[74]=0;
	}
	if (File_out[75]==2 || Screen_out[75]==2) DATA[75]= K_Trunk_Apo;
	if (File_out[75]==1 || Screen_out[75]==1)
	{
		if (Leaf_Area[0]) DATA[75]= K_Trunk_Apo/Leaf_Area[0]; 
		else DATA[75]=0;
	}
	if (File_out[76]==2 || Screen_out[76]==2) DATA[76]= K_Plant; 	// plant conductance
	if (File_out[76]==1 || Screen_out[76]==1)
	{
		if (Leaf_Area[0]) DATA[76]= K_Plant/Leaf_Area[0]; 
		else DATA[76]= 0;
	}
	if (File_out[77]==2 || Screen_out[77]==2) DATA[77]= K_tot; 	// plant + soil conductance
	if (File_out[77]==1 || Screen_out[77]==1)
	{
		if (Leaf_Area[0]) DATA[77]= K_tot/Leaf_Area[0]; 
		else DATA[77]= 0;
	}
	
	DATA[78]= PLC_Leaf_Apo[0];	
	DATA[79]= PLC_Branch_Apo[0];	
	DATA[80]= PLC_Axil_Apo;
	DATA[81]= PLC_Trunk_Apo;
	DATA[82]= PLC_Root_Apo1;
	DATA[83]= PLC_Root_Apo2;
	DATA[84]= PLC_Root_Apo3;
	DATA[85]= Q_Plant;
	if (File_out[86]==2 || Screen_out[86]==2) DATA[86]= (Q_Soil01-Q_Soil1+Q_Soil02-Q_Soil2+Q_Soil03-Q_Soil3)*18/1000/1000; //soil water decreased
	if (File_out[86]==1 || Screen_out[86]==1) DATA[86]= (Q_Soil1+Q_Soil2+Q_Soil3)*18/1000/1000;
	DATA[87]= (Q_Leaf_Evap[1] +Q_Leaf_Evap[2]+Q_Leaf_Evap[3])*18/1000/1000;
	DATA[88]= (Q_Leaf_Symp[1] +Q_Leaf_Symp[2]+Q_Leaf_Symp[3])*18/1000/1000;
	DATA[89]= (Q_Leaf_Apo[1] +Q_Leaf_Apo[2]+Q_Leaf_Apo[3])*18/1000/1000;
	DATA[90]= (Q_Branch_Symp[1] +Q_Branch_Symp[2]+Q_Branch_Symp[3])*18/1000/1000;
	DATA[91]= (Q_Branch_Apo[1] +Q_Branch_Apo[2]+Q_Branch_Apo[3])*18/1000/1000;
	DATA[92]= Q_Axil_Symp*18/1000/1000;
	DATA[93]= Q_Axil_Apo*18/1000/1000;
	DATA[94]= Q_Trunk_Symp*18/1000/1000;
	DATA[95]= Q_Trunk_Apo*18/1000/1000;
	DATA[96]= (Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp2)*18/1000/1000;
	DATA[97]= Q_Root_Apo_t*18/1000/1000;
	DATA[98]= (Q_Root_Endo1+Q_Root_Endo2+Q_Root_Endo3)*18/1000/1000;
	DATA[99]=  Q_Soil1*18/1000/1000;
	DATA[100]= Q_Soil2*18/1000/1000;
	DATA[101]= Q_Soil3*18/1000/1000;
	DATA[102]= Growth_rate_Trunk;
	if (File_out[103]==2 || Screen_out[103]==2) DATA[103]= Radius_Trunk;   // absolute Trunk radius, mm
	else DATA[103]= Radius_Trunk_rel;												// relative Trunk radius, mm
	DATA[104]= A_net[0]; 				//Mean Leaf Photosynthesis per leaf area, µmol s-1 m-2
	DATA[105]= A_net_tot_c; 			 //Annual plant canopy A_net in mol of CO2
	DATA[106]= A_net_tot/E_tot/1000;		//WUE
	DATA[107]= R_main; 					//maintenance respiration
	DATA[108]= Reserve;
	DATA[109]= Leaf_Area[0];
	DATA[110]= Specific_Leaf_WC[1]*Branch_distri[1] +Specific_Leaf_WC[2]*Branch_distri[2]+Specific_Leaf_WC[3]*Branch_distri[3];
	DATA[111]= Specific_Br_WC[1]*Branch_distri[1] +Specific_Br_WC[2]*Branch_distri[2]+Specific_Br_WC[3]*Branch_distri[3];
	DATA[112]= Specific_Br_Symp_WC[1]*Branch_distri[1] +Specific_Br_Symp_WC[2]*Branch_distri[2]+Specific_Br_Symp_WC[3]*Branch_distri[3];
	DATA[113]= Specifc_Br_Apo_WC[1]*Branch_distri[1] +Specifc_Br_Apo_WC[2]*Branch_distri[2]+Specifc_Br_Apo_WC[3]*Branch_distri[3];
	DATA[114]= RWC_shoot;
	DATA[115]= RWC_Axil;
	DATA[116]= RWC_Soil;
	DATA[117]= RWC_min;
	if (File_out[118]==2 || Screen_out[118]==2) DATA[118]= REW_int1;  // cumulated REW deficit, like in BilJou based on Teta_wp at Pf=4.2		
	else if (File_out[118]==3 || Screen_out[118]==3) DATA[118]= REW_int2;  // cumulated REW deficit based on Teta_r	
	else DATA[118]= RWC_int;
	DATA[119]= P_soil;
	DATA[120]= P_soil_min;
	DATA[121]= P_soil_int;
	DATA[122]= Irrigation;
	DATA[123]= Sap_Flow_d;
	DATA[124]= Sap_Flow_d2;
	DATA[125]= VPD_Air;
	DATA[126]= Diam_Fruit*1000; // in mm
	DATA[127]= ETP_Penman_tot;  // cumul EPT in mm
	DATA[128]= ETP_Penman;      // ETP Penman in mmol/s
	DATA[129]= Rain_tot;
	DATA[130]= Rain_soil;
	DATA[131]= VPD_Air_tot;
	DATA[132]= VPD_Leaf_tot;
	DATA[133]= Rain_Leaf_tot;
	DATA[134]= ETP_Leaf_tot;
	if (N_days) DATA[135]= Cum_T_air/(long double)N_days; else DATA[135]=0;
	if (N_days2) DATA[136]= Cum_T_air_l/(long double)N_days2; else DATA[136]=0;
	DATA[137]= EvapoT;
	DATA[138]= EvapoT_tot*18/1000/1000/Surface_Soil;  // cumul plant + soil evaporation in mm
	DATA[139]= ETP_Penman_day;           // cumul ETP per day in mm
	DATA[140]= (dq_Trunk_Apo_Branch_Apo[0]+dq_Trunk_Apo_Trunk_Symp)/dt_court;  //sap flow in the Trunk in mmol/s; includes Trunk bark transpiration
	if (Leaf_Area[0]) DATA[141]=K_Branch_Apo[0]/Leaf_Area[0]; else DATA[141]=0;
	DATA[142]= P_min_Leaf;
	DATA[143]= P_min_stem;
	if (File_out[144]==2 || Screen_out[144]==2) DATA[144]= K_tot2; 	// plant + soil conductance
	if (File_out[144]==1 || Screen_out[144]==1)
	{
		if (Leaf_Area[0]) DATA[144]= K_tot2/Leaf_Area[0]; 
		else DATA[144]= 0;
	}
	DATA[145]= A_net[0]*Leaf_Area[0];                              // Canopy A_net µmol/s
	DATA[146]= Rain_tot-Rain_soil ;					// Canopy interception, mm
	DATA[147]= GPP ; //GPP in g of C per soil surface per year
	DATA[148]= A_net_day/Surface_Soil*12 ; //GPP in g of C per soil surface per day
	DATA[149]= EvapoT_day*18/1000/1000/Surface_Soil;
	if (File_out[150]==2 || Screen_out[150]==2) DATA[150]= E_tot_day*18/1000/1000; 	//E in kg per day
	else DATA[150]= E_tot_day*18/1000/1000/Crown_Area;		//E in mm per day
	DATA[151]= T_Soil;
	DATA[152]= POTENTIAL_PAR;
	DATA[153]= cloud_cover;
	DATA[154]= -dq_Leaf_Symp_Leaf_Evap[1]/dt_court;
	DATA[155]= dq_Leaf_Apo_Leaf_Evap[1]/dt_court;
	DATA[156]= dq_Branch_Apo_Leaf_Apo[1]/dt_court;
	DATA[157]= dq_Branch_Apo_Branch_Symp[1]/dt_court;
	DATA[158]= dq_Trunk_Apo_Branch_Apo[1]/dt_court;
	DATA[159]= dq_Trunk_Apo_Trunk_Symp/dt_court;
	DATA[160]= dq_Axil_Apo_Axil_Symp/dt_court;
	DATA[161]= dq_Branch_Apo_Axil_Apo/dt_court;
	DATA[162]= (dq_Root_Apo_Trunk_Apo1+dq_Root_Apo_Trunk_Apo2+dq_Root_Apo_Trunk_Apo3)/dt_court;
	DATA[163]= (dq_Root_Endo_Root_Apo1+dq_Root_Endo_Root_Apo2+dq_Root_Endo_Root_Apo3)/dt_court;
	DATA[164]= (dq_Root_Symp_Root_Endo1+dq_Root_Symp_Root_Endo2+dq_Root_Symp_Root_Endo3)/dt_court;
	DATA[165]= (dq_Soil_Root_Endo1+dq_Soil_Root_Endo2+dq_Soil_Root_Endo3)/dt_court;
	if (File_out[166]==2 || Screen_out[166]==2) DATA[166]= -dq_Soil_12/dt_court*18/1000/1000/Surface_Soil*3600*24;  //capillary rise in mm/day
	else DATA[166]= -dq_Soil_12/dt_court;																				//capillary rise in mmol/s
	if (File_out[167]==2 || Screen_out[167]==2) DATA[167]= -dq_Soil_23/dt_court*18/1000/1000/Surface_Soil*3600*24;
	else DATA[167]= -dq_Soil_23/dt_court;
	DATA[168]= A_net_V;
	DATA[169]= A_net_J;
	DATA[170]= g_Axil;
	DATA[171]= Root_Area_fi*(Root_upper+Root_lower+Root_middle)/(Root_upper0+Root_lower0+Root_middle0);
	DATA[172]= Axil_WC;
	DATA[173]= 		(Q_Soil1*18/1000/1000/1000/Volume_soil1 -Teta_r_1)*Soil_Depth*Layer_1*1000*(1-Rock_f1) 
				+ 	(Q_Soil2*18/1000/1000/1000/Volume_soil2 -Teta_r_2)*Soil_Depth*Layer_2*1000*(1-Rock_f2) 
				+ 	(Q_Soil3*18/1000/1000/1000/Volume_soil3 -Teta_r_3)*Soil_Depth*Layer_3*1000*(1-Rock_f3) ; //RU tot
	DATA[174]= E_tot*18/1000/1000/Crown_Area;
	DATA[175]= Drainage*18/1000/1000/Surface_Soil;
	DATA[176]= R_growth;          		 //Growth respiration in µmol
	DATA[177]= A_net_tot;          	 //Annual A_net in mol
	DATA[178]= A_gross_tot;         	//Annual A_gross in mol
	DATA[179]= Resp_tot;            	//Annual Leaf Respiration in mol
	DATA[180]= Export;
	DATA[181]= Export_tot/1000;
	DATA[182]= T_Soil1;
	DATA[183]= T_Soil2;
	DATA[184]= T_Soil3;
	DATA[185]=		(Q_Soil1*18/1000/1000/1000/Volume_soil1 -Teta_wp_1)*Soil_Depth*Layer_1*1000*(1-Rock_f1)  
				+ 	(Q_Soil2*18/1000/1000/1000/Volume_soil2 -Teta_wp_2)*Soil_Depth*Layer_2*1000*(1-Rock_f2) 
				+ 	(Q_Soil3*18/1000/1000/1000/Volume_soil3 -Teta_wp_3)*Soil_Depth*Layer_3*1000*(1-Rock_f3) ;
	
	DATA[186]= REW_t; //RWC based on Teta_r
	DATA[187]= Q_Petiole_Symp*18/1000/1000;   // in Kg
	if (File_out[188]==1 || Screen_out[188]==1) DATA[188]= E_Petiole;
	if (File_out[188]==2 || Screen_out[188]==2) DATA[188]= E_Petiole*Petiole_area;

	DATA[189]= P_Petiole_Symp;
	DATA[190]= Petiole_diam;
	DATA[191]= Q_Plant_a;
	DATA[192]= Q_Plant_s;
	DATA[193]= RWC_Leaf[1]*Branch_distri[1] +RWC_Leaf[2]*Branch_distri[2]+RWC_Leaf[3]*Branch_distri[3];
	DATA[194]= RWC_Leaf_s[1]*Branch_distri[1] +RWC_Leaf_s[2]*Branch_distri[2]+RWC_Leaf_s[3]*Branch_distri[3];
	DATA[195]= RWC_Leaf_a[1]*Branch_distri[1] +RWC_Leaf_a[2]*Branch_distri[2]+RWC_Leaf_a[3]*Branch_distri[3];
	DATA[196]= RWC_Branch[1]*Branch_distri[1] +RWC_Branch[2]*Branch_distri[2]+RWC_Branch[3]*Branch_distri[3];
	DATA[197]= RWC_Branch_s[1]*Branch_distri[1] +RWC_Branch_s[2]*Branch_distri[2]+RWC_Branch_s[3]*Branch_distri[3];
	DATA[198]= RWC_Branch_a[1]*Branch_distri[1] +RWC_Branch_a[2]*Branch_distri[2]+RWC_Branch_a[3]*Branch_distri[3];
	DATA[199]= RWC_Trunk;
	DATA[200]= RWC_Trunk_s;
	DATA[201]= RWC_Trunk_a;
	DATA[202]= RWC_Root;
	DATA[203]= RWC_Root_s;
	DATA[204]= RWC_Root_a;
	DATA[205]= RWC_Plant;
	DATA[206]= RWC_Plant_s;
	DATA[207]= RWC_Plant_a;
	if (Leaf_Area[0]) DATA[208]= E_tot_day/3600/24/Leaf_Area[0]; else DATA[208]=0;  // in mmol/s/m2
	DATA[209]= P_min_lf_d;
	DATA[210]= P_max_lf_d;
	DATA[211]= gs_min_d;
	DATA[212]= gs_max_d;
	DATA[213]= SF_min_d;
	DATA[214]= SF_max_d;
	DATA[215]= gs_max_d2;   //gs_max at Pmin
	DATA[216]= Growth_Fruit*1000; //Fruit diameter in mm without elastic variation
	DATA[217]= Teta_Soil; 
	DATA[218]= Teta_Soil1; 
	DATA[219]= Teta_Soil2; 
	DATA[220]= Teta_Soil3; 
	DATA[221]= Ca; 
	if (File_out[222]==2 || Screen_out[222]==2) DATA[222]= K_Plant_20; 	// plant + soil conductance
	if (File_out[222]==1 || Screen_out[222]==1)
	{
		if (Leaf_Area[0]) DATA[222]= K_Plant_20/Leaf_Area[0]; 
		else DATA[222]= 0;
	}
	DATA[223]= PLC_Plant;
	DATA[224]= P50_Leaf_Apo[1];
	DATA[225]= P50_Branch_Apo[1];
	DATA[226]= P50_Trunk_Apo;
	DATA[227]= P50_Root_Apo;
	DATA[228]= Pi0_Leaf_Symp;
	DATA[229]= Px_gs;
	DATA[230]= P_Leaf_Symp[0]-Turgor_Leaf_Symp[0] ;
	if (File_out[231]==2 || Screen_out[231]==2) DATA[231]= K_tot_20; 	// plant + soil conductance
	if (File_out[231]==1 || Screen_out[231]==1)
	{
		if (Leaf_Area[0]) DATA[231]= K_tot_20/Leaf_Area[0]; 
		else DATA[231]= 0;
	}
	DATA[232]= PLC_Plant_Soil;
	DATA[233]=PAR_tot;
	if (N_days2) DATA[234]=P_min_Leaf_mean_Leafy/N_days2; else DATA[234]=0;
	DATA[235]=(g_s[1]+g_cuti[1])*Branch_distri[1]+ (g_s[2]+g_cuti[2])*Branch_distri[2]+(g_s[3]+g_cuti[3])*Branch_distri[3]; 
	DATA[236]=E_stomata[1]*Branch_distri[1] +E_stomata[2]*Branch_distri[2]+E_stomata[3]*Branch_distri[3];
	DATA[237]=gs_Jarvis;
	DATA[238]=g_crown;
	DATA[239]=g_bl;
	DATA[240]=Cm;
	DATA[241]=VcMaxT;
	DATA[242]=VjMaxT;
	DATA[243]=Rd;
	DATA[244]=Ci;
	DATA[245]=Cbs;
	DATA[246]=PLC_Leaf_Apo[1];
	DATA[247]=PLC_Leaf_Apo[2];
	DATA[248]=PLC_Leaf_Apo[3];
	DATA[249]=PLC_Branch_Apo[1];
	DATA[250]=PLC_Branch_Apo[2];
	DATA[251]=PLC_Branch_Apo[3];
	DATA[252]=g_s[1];
	DATA[253]=g_s[2];
	DATA[254]=g_s[3];
	DATA[255]=P_Leaf_Symp[1];
	DATA[256]=P_Leaf_Symp[2];
	DATA[257]=P_Leaf_Symp[3];
	DATA[258]=P_Leaf_Apo[1];
	DATA[259]=P_Leaf_Apo[2];
	DATA[260]=P_Leaf_Apo[3];
	DATA[261]=P_Branch_Apo[1];
	DATA[262]=P_Branch_Apo[2];
	DATA[263]=P_Branch_Apo[3];
	DATA[264]=T_Leaf[1];
	DATA[265]=T_Leaf[2];
	DATA[266]=T_Leaf[3];
	DATA[267]=E_Leaf[1];
	DATA[268]=E_Leaf[2];
	DATA[269]=E_Leaf[3];
	DATA[270]=A_net[1];
	DATA[271]=A_net[2];
	DATA[272]=A_net[3];
	DATA[273]=PAR[2]; //the PAR at the base of branch 1. At the top it is PAR[1]=PAR[0]
	DATA[274]=PAR[3];
	DATA[275]=PAR[4]; // the PAR reaching the soil
	DATA[276]=Leaf_Area[1];
	DATA[277]=Leaf_Area[2];	
	DATA[278]=Leaf_Area[3];
	if (COMPET)
	{
	DATA[279]=E_tot*18/1000/1000*COMPET;
	DATA[280]=E_tot*18/1000/1000/Crown_Area*COMPET;
	DATA[281]=E_tot_day*18/1000/1000*COMPET;
	DATA[282]=E_tot_day*18/1000/1000/Crown_Area*COMPET;
	}
	else{
	DATA[279]=E_tot*18/1000/1000*1;
	DATA[280]=E_tot*18/1000/1000/Crown_Area*1;
	DATA[281]=E_tot_day*18/1000/1000*1;
	DATA[282]=E_tot_day*18/1000/1000/Crown_Area*1;
	} 
	DATA[283]=T_Branch1;
	DATA[284]=T_Branch2;
	DATA[285]=T_Branch3;
	DATA[286]=T_Trunk;
	DATA[287]=C_growth_tot;
	DATA[288]=R_main_tot;
	DATA[289]=R_growth_tot;
	DATA[290]=T_Axil;
	DATA[291]=Turgor_Axil_Symp;
	DATA[292]= (dq_Soil_Root_Endo1)/dt_court;
	DATA[293]= (dq_Soil_Root_Endo2)/dt_court;
	DATA[294]= (dq_Soil_Root_Endo3)/dt_court;
	DATA[295]=REW_t1;
	DATA[296]=REW_t2;
	DATA[297]=REW_t3;
	DATA[298]=REW_wp;
	DATA[299]=REW_wp1;
	DATA[300]=REW_wp2;
	DATA[301]=REW_wp3;
	DATA[302]=Test[1];
	DATA[303]=Test[2];
}

void compute_climatic_stats(void) //computes some stats on the annual climatic data when given on a daily basis
{
	long double YEAR, YEAR0=2000, DOI,DOI0,Tmin,Tmax,Tsoil,RHmin,RHmax,PAR_day,Rain,Wind_day;
	long double T_annual_mean=0,T_summer_mean=0,T_1_mean=0,T_2_mean=0,T_3_mean=0,T_4_mean=0;
	long double T_annual_min=0,T_summer_min=0,T_1_min=0,T_2_min=0,T_3_min=0,T_4_min=0;
	long double T_annual_max=0,T_summer_max=0,T_1_max=0,T_2_max=0,T_3_max=0,T_4_max=0;
	long double Wind_annual_mean=0,Wind_summer_mean=0,Wind_1_mean=0,Wind_2_mean=0,Wind_3_mean=0,Wind_4_mean=0;
	long double HR_annual_mean=0,HR_summer_mean=0,HR_1_mean=0,HR_2_mean=0,HR_3_mean=0,HR_4_mean=0;
	long double HR_annual_min=0,HR_summer_min=0,HR_1_min=0,HR_2_min=0,HR_3_min=0,HR_4_min=0;
	long double PAR_annual_max=0,PAR_summer_max=0,PAR_1_max=0,PAR_2_max=0,PAR_3_max=0,PAR_4_max=0;
	long double Rain_annual=0,Rain_summer=0,Rain_11=0,Rain_22=0,Rain_3=0,Rain_4=0;
	long double ETP_annual=0,ETP_summer=0,ETP_1=0,ETP_2=0,ETP_3=0,ETP_4=0;
	long double VPD_annual=0,VPD_summer=0,VPD_1=0,VPD_2=0,VPD_3=0,VPD_4=0;
	long double Days=0,days_summer=0,days_1=0,days_2=0,days_3=0,days_4=0,tangente,TTTT;
	long double e_sat_air,e_air,slope,vpd,ETP;
	char filename_Stat[]="                                      ";
	FILE *out;
	int i;
	
	printf(".");
	sprintf(filename_Stat, "Stat_%.0Lf.csv",Simul);	
	out = fopen(filename_Stat,"w");
	fprintf(out,"Grid;YEAR;T_annual_mean;T_annual_min;T_annual_max;HR_annual_mean;HR_annual_min;Wind_annual_mean;PAR_annual_max;Rain_annual;VPD_annual;ETP_annual;P-ETP_annual;;");
	fprintf(out,"T_summer_mean;T_summer_min;T_summer_max;HR_summer_mean;HR_summer_min;Wind_summer_mean;PAR_summer_max;Rain_summer;VPD_summer;ETP_summer;P-ETP_summer;\n");
//	fprintf(out,"T_1_mean;T_1_min;T_1_max;HR_1_mean;HR_1_min;Wind_1_mean;PAR_1_max;Rain_11;VPD_1;ETP_1;P-ETP_1;;");
//	fprintf(out,"T_2_mean;T_2_min;T_2_max;HR_2_mean;HR_2_min;Wind_2_mean;PAR_2_max;Rain_22;VPD_2;ETP_2;P-ETP_2;;");
//	fprintf(out,"T_3_mean;T_3_min;T_3_max;HR_3_mean;HR_3_min;Wind_3_mean;PAR_3_max;Rain_3;VPD_3;ETP_3;P-ETP_3;;");
//	fprintf(out,"T_4_mean;T_4_min;T_4_max;HR_4_mean;HR_4_min;Wind_4_mean;PAR_4_max;Rain_4;VPD_4;ETP_4;P-ETP_4;\n");

	climat_in = fopen(filename_CLIM,"r+");
	if (CLIMAT==4) 		 fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI0, &Tmin, &Tmax, &Tsoil, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
	else           	 fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",    &YEAR, &DOI0, &Tmin, &Tmax, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
	Days++;
	vpd=0;
	ETP=0;
	
	for(i=0;i<24;i++)
	{
		TTTT=(float)i; // time of day
		T_air  = (Tmax+Tmin)/2+(Tmax-Tmin)/2*cos(3.14159265359/12*(TTTT-HH1));    //daily cos variation between T_min et T_max et T_max at 12h00
		RH_air = (RHmax+RHmin)/2+(RHmax-RHmin)/2*cos(3.14159265359/12*(TTTT-HH2)); //daily cos variation between HR_min et T_max et HR_max at 0h00
		tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOI0-80)/365*2*3.1416)));
		if (tangente<-1)     Day_length= 0;
		else if (tangente>1) Day_length= 24;
		else                 Day_length= 24 -  acos(tangente)*7.6394194;
		if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR[0]=0;
		else PAR[0] = PAR_day*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
		
		if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
		
		e_air=e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
		vpd+= (e_sat_air - e_air)/1000;
		slope=4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
		ETP+= 0.5625*(0.408*slope*PAR[0]*0.5495*3.6e-3+0.066*37*Wind_day*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind_day));  //ETP in mm/h  0.5625 to fit observed ETP
		//  printf("%Lf T %Lf RH %Lf PAR %Lf vpd %Lf ETP %Lf \n", TTTT, T_air,RH_air, PAR,vpd,ETP);
	}
	
	T_annual_mean=  (Tmin+Tmax)/2;
	T_annual_min=   Tmin;
	T_annual_max=   Tmax;
	HR_annual_mean= (RHmin+RHmax)/2;
	HR_annual_min=  RHmin;
	Wind_annual_mean=Wind_day;
	PAR_annual_max= PAR_day;
	Rain_annual=    Rain;
	VPD_annual=     vpd;
	ETP_annual=     ETP;
	
	
	while (!feof(climat_in))
	{
		YEAR0=YEAR;
		if (CLIMAT==4) 	fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI, &Tmin, &Tmax, &Tsoil, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
		else 			fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI, &Tmin, &Tmax, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
		if (DOI<DOI0) // a new year
		{
			fprintf(out,"%.0Lf;%.0Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;; ",Simul, YEAR0,T_annual_mean,T_annual_min,T_annual_max,HR_annual_mean,HR_annual_min,Wind_annual_mean,PAR_annual_max,Rain_annual,VPD_annual,ETP_annual,Rain_annual-ETP_annual);
			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;\n", T_summer_mean,T_summer_min,T_summer_max,HR_summer_mean,HR_summer_min,Wind_summer_mean,PAR_summer_max,Rain_summer,VPD_summer,ETP_summer,Rain_summer-ETP_summer);
//			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;;", T_1_mean,T_1_min,T_1_max,HR_1_mean,HR_1_min,Wind_1_mean,PAR_1_max,Rain_11,VPD_1,ETP_1,Rain_1-ETP_1);
//			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;;", T_2_mean,T_2_min,T_2_max,HR_2_mean,HR_2_min,Wind_2_mean,PAR_2_max,Rain_22,VPD_2,ETP_2,Rain_2-ETP_2);
//			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;;", T_3_mean,T_3_min,T_3_max,HR_3_mean,HR_3_min,Wind_3_mean,PAR_3_max,Rain_3,VPD_3,ETP_3,Rain_3-ETP_3);
//			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;\n", T_4_mean,T_4_min,T_4_max,HR_4_mean,HR_4_min,Wind_4_mean,PAR_4_max,Rain_4,VPD_4,ETP_4,Rain_4-ETP_4);
		//	 printf("%.0Lf %Lf T %Lf RH %Lf vpd %Lf ETP %Lf \n",YEAR0, TTTT, T_air,RH_air,vpd,ETP);			
	
			T_annual_mean=T_annual_min=T_annual_max=HR_annual_mean=HR_annual_min=Wind_annual_mean=PAR_annual_max=Rain_annual=VPD_annual=ETP_annual=0;
			T_summer_mean=T_summer_min=T_summer_max=HR_summer_mean=HR_summer_min=Wind_summer_mean=PAR_summer_max=Rain_summer=VPD_summer=ETP_summer=0;
			T_1_mean=T_1_min=T_1_max=HR_1_mean=HR_1_min=Wind_1_mean=PAR_1_max=Rain_11=VPD_1=ETP_1=0;
			T_2_mean=T_2_min=T_2_max=HR_2_mean=HR_2_min=Wind_2_mean=PAR_2_max=Rain_22=VPD_2=ETP_2=0;
			T_3_mean=T_3_min=T_3_max=HR_3_mean=HR_3_min=Wind_3_mean=PAR_3_max=Rain_3=VPD_3=ETP_3=0;
			T_4_mean=T_4_min=T_4_max=HR_4_mean=HR_4_min=Wind_4_mean=PAR_4_max=Rain_4=VPD_4=ETP_4=0;
				
			Days=0;
			days_summer=0;
			days_1=days_2=days_3=days_4=0;
		}
		DOI0=DOI;
		if (INTERCEPTION) Interception();
		Days++;
		ETP=0;
		vpd=0;
		for(i=0;i<24;i++)
		{
			TTTT=(float)i; // time of day
			T_air  = (Tmax+Tmin)/2+(Tmax-Tmin)/2*cos(3.14159265359/12*(TTTT-HH1));    //daily cos variation between T_min et T_max et T_max at 12h00
			RH_air = (RHmax+RHmin)/2+(RHmax-RHmin)/2*cos(3.14159265359/12*(TTTT-HH2)); //daily cos variation between HR_min et T_max et HR_max at 0h00
			tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOI-80)/365*2*3.1416)));
			if      (tangente<-1) Day_length= 0;
			else if (tangente>1)  Day_length= 24;
			else                  Day_length= 24 -  acos(tangente)*7.6394194;
			if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR[0]=0;
			else PAR[0]= PAR_day*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
			
			if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
			else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
			
			e_air=   e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
			vpd+=   (e_sat_air - e_air)/1000;
			slope=  4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
			ETP+=   0.5625*(0.408*slope*PAR[0]*0.5495*3.6e-3+0.066*37*Wind_day*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind_day));  //ETP in mm/h  0.5625 to fit observed ETP
		}
		T_annual_mean=  (T_annual_mean*(Days-1)+(Tmin+Tmax)/2)/Days;
		T_annual_min=   (T_annual_min*(Days-1)+Tmin)/Days;
		T_annual_max=   (T_annual_max*(Days-1)+Tmax)/Days;
		HR_annual_mean= (HR_annual_mean*(Days-1)+(RHmin+RHmax)/2)/Days;
		HR_annual_min=  (HR_annual_min*(Days-1) + RHmin)/Days;
		Wind_annual_mean=(Wind_annual_mean*(Days-1)+Wind_day)/Days;
		PAR_annual_max=(PAR_annual_max*(Days-1)+PAR_day)/Days;
		//VPD_annual=(VPD_annual*(Days-1)+vpd)/Days;
		//ETP_annual=(ETP_annual*(Days-1)+ETP)/Days; // ce n'est pas l'ETP journalière donc !
		VPD_annual+=vpd;
		ETP_annual+=ETP;
		Rain_annual+=Rain;
		
		if (DOI>171 && DOI<265)  //summer period June 20 to Sept 22
		{
			days_summer++;
			T_summer_mean=  (T_summer_mean*(days_summer-1)+(Tmin+Tmax)/2)/days_summer;
			T_summer_min=   (T_summer_min*(days_summer-1)+Tmin)/days_summer;
			T_summer_max=   (T_summer_max*(days_summer-1)+Tmax)/days_summer;
			HR_summer_mean= (HR_summer_mean*(days_summer-1)+(RHmin+RHmax)/2)/days_summer;
			HR_summer_min=  (HR_summer_min*(days_summer-1) + RHmin)/days_summer;
			Wind_summer_mean=(Wind_summer_mean*(days_summer-1)+Wind_day)/days_summer;
			PAR_summer_max= (PAR_summer_max*(days_summer-1)+PAR_day)/days_summer;
			ETP_summer+=ETP;
			VPD_summer+=vpd;
			Rain_summer+=Rain;
		}
		
		if (DOI>=0 && DOI<91) // Jan 1 to April 1
			{
			days_1++;
			T_1_mean=  (T_1_mean*(days_1-1)+(Tmin+Tmax)/2)/days_1;
			T_1_min=   (T_1_min*(days_1-1)+Tmin)/days_1;
			T_1_max=   (T_1_max*(days_1-1)+Tmax)/days_1;
			HR_1_mean= (HR_1_mean*(days_1-1)+(RHmin+RHmax)/2)/days_1;
			HR_1_min=  (HR_1_min*(days_1-1) + RHmin)/days_1;
			Wind_1_mean=(Wind_1_mean*(days_1-1)+Wind_day)/days_1;
			PAR_1_max= (PAR_1_max*(days_1-1)+PAR_day)/days_1;
			ETP_1+=ETP;
			VPD_1+=vpd;
			Rain_11+=Rain;
			}
			
		if (DOI>=91 && DOI<182) // April 1 to July 1
			{
			days_2++;
			T_2_mean=  (T_2_mean*(days_2-1)+(Tmin+Tmax)/2)/days_2;
			T_2_min=   (T_2_min*(days_2-1)+Tmin)/days_2;
			T_2_max=   (T_2_max*(days_2-1)+Tmax)/days_2;
			HR_2_mean= (HR_2_mean*(days_2-1)+(RHmin+RHmax)/2)/days_2;
			HR_2_min=  (HR_2_min*(days_2-1) + RHmin)/days_2;
			Wind_2_mean=(Wind_2_mean*(days_2-1)+Wind_day)/days_2;
			PAR_2_max= (PAR_2_max*(days_2-1)+PAR_day)/days_2;
			ETP_2+=ETP;
			VPD_2+=vpd;
			Rain_22+=Rain;
			}
			
		if (DOI>=182 && DOI<273) // July 1 to Oct 1
			{
			days_3++;
			T_3_mean=  (T_3_mean*(days_3-1)+(Tmin+Tmax)/2)/days_3;
			T_3_min=   (T_3_min*(days_3-1)+Tmin)/days_3;
			T_3_max=   (T_3_max*(days_3-1)+Tmax)/days_3;
			HR_3_mean= (HR_3_mean*(days_3-1)+(RHmin+RHmax)/2)/days_3;
			HR_3_min=  (HR_3_min*(days_3-1) + RHmin)/days_3;
			Wind_3_mean=(Wind_3_mean*(days_3-1)+Wind_day)/days_3;
			PAR_3_max= (PAR_3_max*(days_3-1)+PAR_day)/days_3;
			ETP_3+=ETP;
			VPD_3+=vpd;
			Rain_3+=Rain;
			}
			
		if (DOI>=273 && DOI<366) // Oct 1 to Jan 1
			{
			days_4++;
			T_4_mean=  (T_4_mean*(days_4-1)+(Tmin+Tmax)/2)/days_4;
			T_4_min=   (T_4_min*(days_4-1)+Tmin)/days_4;
			T_4_max=   (T_4_max*(days_4-1)+Tmax)/days_4;
			HR_4_mean= (HR_4_mean*(days_4-1)+(RHmin+RHmax)/2)/days_4;
			HR_4_min=  (HR_4_min*(days_4-1) + RHmin)/days_4;
			Wind_4_mean=(Wind_4_mean*(days_4-1)+Wind_day)/days_4;
			PAR_4_max= (PAR_4_max*(days_4-1)+PAR_day)/days_4;
			ETP_4+=ETP;
			VPD_4+=vpd;
			Rain_4+=Rain;		
			}
			
	}
	fprintf(out,"%.0Lf;%.0Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;; ",Simul, YEAR0,T_annual_mean,T_annual_min,T_annual_max,HR_annual_mean,HR_annual_min,Wind_annual_mean,PAR_annual_max,Rain_annual,VPD_annual,ETP_annual,Rain_annual-ETP_annual);
	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf;\n", T_summer_mean,T_summer_min,T_summer_max,HR_summer_mean,HR_summer_min,Wind_summer_mean,PAR_summer_max,Rain_summer,VPD_summer,ETP_summer,Rain_summer-ETP_summer);
//	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf;;", T_1_mean,T_1_min,T_1_max,HR_1_mean,HR_1_min,Wind_1_mean,PAR_1_max,Rain_11,VPD_1,ETP_1,Rain_1-ETP_1);
//	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf;;", T_2_mean,T_2_min,T_2_max,HR_2_mean,HR_2_min,Wind_2_mean,PAR_2_max,Rain_22,VPD_2,ETP_2,Rain_2-ETP_2);
//	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf;;", T_3_mean,T_3_min,T_3_max,HR_3_mean,HR_3_min,Wind_3_mean,PAR_3_max,Rain_3,VPD_3,ETP_3,Rain_3-ETP_3);
//	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf;\n", T_4_mean,T_4_min,T_4_max,HR_4_mean,HR_4_min,Wind_4_mean,PAR_4_max,Rain_4,VPD_4,ETP_4,Rain_4-ETP_4);
	
	fclose(climat_in);
	fclose(out);
	//exit(1);
}


void Synchronise_Soil_compet (void)  //when plants compet for water in the same soil volume
{
	long double  Q_soil1_saved, Q_soil2_saved, Q_soil3_saved, T_init;
	FILE *soil;
	T_init=0;			
	do
	{
		soil = fopen(filename_SOIL,"r+");
		fscanf(soil,"%Lf %Lf %Lf %Lf",&T_init, &Q_soil1_saved,&Q_soil2_saved,&Q_soil3_saved);
		fclose(soil);
	}
	while (T/1000>T_init); // wait for the other tree

}

void Soil_compet (void)  //when plants compet for water in the same soil volume
{
	long double  Q_soil1_saved, Q_soil2_saved, Q_soil3_saved,T_init;
	long double  dQ_soil1, dQ_soil2, dQ_soil3;
	FILE *soil;
	
	dQ_soil1=Q_soil1_init-Q_Soil1;
	dQ_soil2=Q_soil2_init-Q_Soil2;
	dQ_soil3=Q_soil3_init-Q_Soil3;
	
	soil = fopen(filename_SOIL,"r+");  
	fscanf(soil,"%Lf %Lf %Lf %Lf",&T_init,&Q_soil1_saved,&Q_soil2_saved,&Q_soil3_saved);
	
	if ((T/1000-T_init)>6000) 
		{
		fclose(soil);
		Synchronise_Soil_compet();
		soil = fopen(filename_SOIL,"r+");  
		fscanf(soil,"%Lf %Lf %Lf %Lf",&T_init, &Q_soil1_saved,&Q_soil2_saved,&Q_soil3_saved);
		}
		  
	rewind (soil);
	Q_soil1_init= Q_soil1_saved - dQ_soil1;
	Q_soil2_init= Q_soil2_saved - dQ_soil2;
	Q_soil3_init= Q_soil3_saved - dQ_soil3;
	fprintf(soil,"%Lf %Lf %Lf %Lf",T/1000, Q_soil1_init,Q_soil2_init,Q_soil3_init);
	fclose(soil);

	Q_Soil1=Q_soil1_init;
	Q_Soil2=Q_soil2_init;
	Q_Soil3=Q_soil3_init;
}

void Gravity_test(void)
{
	if(T/1000>(DOY+6.25)*3600*24)
		{
		Pg_Leaf[1]  = -9.81*(Length_Branch+Length_Trunk*4)/1000;
		Pg_Trunk = -9.81*Length_Trunk*4/1000;                     // gravimetric pressure drop at the top leaves, Trunk and Branches in MPa
		Pg_Branch[1]= -9.81*(Length_Branch+Length_Trunk*4)/1000;
		P_Trunk_Apo=Pg_Trunk;
		}
	else if(T/1000>(DOY+4.25)*3600*24)
		{
		Pg_Leaf[1]  = -9.81*(Length_Branch+Length_Trunk*3)/1000;
		Pg_Trunk = -9.81*Length_Trunk*3/1000;                     // gravimetric pressure drop at the top leaves, Trunk and Branches in MPa
		Pg_Branch[1]= -9.81*(Length_Branch+Length_Trunk*3)/1000;
		P_Trunk_Apo=Pg_Trunk;
		}
	else if(T/1000>(DOY+2.25)*3600*24)
		{
		Pg_Leaf[1]  = -9.81*(Length_Branch+Length_Trunk*2)/1000;
		Pg_Trunk = -9.81*Length_Trunk*2/1000;                     // gravimetric pressure drop at the top leaves, Trunk and Branches in MPa
		Pg_Branch[1]= -9.81*(Length_Branch+Length_Trunk*2)/1000;
		P_Trunk_Apo=Pg_Trunk;
		}
	else if(T/1000>(DOY+0.25)*3600*24)
		{
		Pg_Leaf[1]  = -9.81*(Length_Branch+Length_Trunk)/1000;
		Pg_Trunk = -9.81*Length_Trunk/1000;                     // gravimetric pressure drop at the top leaves, Trunk and Branches in MPa
		Pg_Branch[1]= -9.81*(Length_Branch+Length_Trunk)/1000;
		P_Trunk_Apo=Pg_Trunk;
		}
}

void Fractal_loop_shoot (int ordre,int i)
{
	int j;
	long double x,y;	
	x=X[ordre-1];
	y=Y[ordre-1];
	if ((rand()%2)) side[ordre]=side[ordre-1]+1; else side[ordre]=side[ordre-1]-1;
	for(j=i-1;j>=0;j--) 
		{
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_segment[j]*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[j]*scale*sin(side[ordre-1]*angle),y+Length_segment[j]*scale*cos(side[ordre-1]*angle));	
			x+=Length_segment[j]*scale*sin(side[ordre-1]*angle);
			y+=Length_segment[j]*scale*cos(side[ordre-1]*angle);
			X[ordre]=x;
			Y[ordre]=y;			
			if ((rand()%1000)<((N_fractal-1)*1000)) Fractal_loop_shoot(ordre+1,j);
		}
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",pow(LA_term_shoot,0.5)/100*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x,y+pow(LA_term_shoot,0.5)/100*scale);	
	
}

void Fractal_loop_Root (int ordre,int i)
{
	int j;
	long double x,y;	
	x=X[ordre-1];
	y=Y[ordre-1];
	if ((rand()%2)) side[ordre]=side[ordre-1]+1; else side[ordre]=side[ordre-1]-1;
	for(j=i-1;j>=0;j--) 
		{
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_segment[j]*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[j]*scale*cos(side[ordre-1]*angle+3.1416/2),y-Length_segment[j]*scale*sin(side[ordre-1]*angle+3.1416/2));	
			x+=Length_segment[j]*scale*cos(side[ordre-1]*angle+3.1416/2);
			y-=Length_segment[j]*scale*sin(side[ordre-1]*angle+3.1416/2);
			X[ordre]=x;
			Y[ordre]=y;			
			if ((rand()%1000)<((N_fractal-1)*1000)) Fractal_loop_Root(ordre+1,j);
		}
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",Diam_Root*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_fine_Root/100*scale*cos(side[ordre-1]*angle+3.1416/2),y-Length_fine_Root/100*scale*sin(side[ordre-1]*angle+3.1416/2));	
	
	//fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_fine_Root/100*scale,y-Length_fine_Root/100*scale);	
	
}

void Fractal(void)  //compute a fractal tree for the whole plant when FRACTAL=2 or for the Root only when FRACTAL=1;
{
	long double SEGMENTS=100; //discretisation of the shoot and Root systems 100 is enough
	long double Diam_term_sapwood, Diam_Trunk_sapwood;
	long double A_term_Root, Length_term_Root,Diam_tapRoot, Length_tapRoot, Diam_tapRoot_sapwood,N_Root, Mean_diam_Root;
	long double Nb_term_shoot,Nb_term_Root,FineRoot_area, TapRoot_area, Volume_fineRoot, Volume_tapRoot, Volume_fineRoot_sapwood, Volume_tapRoot_sapwood, R_fineRoot, R_tapRoot;
	long double F,Delta1,Delta2;
	long double Discriminant, Node_base_crown=0;
	long double Volume_Branch ,Volume_Branch_sapwood, Volume_Trunk ,Volume_Trunk_sapwood, R_Branch, R_Trunk; 
	long double l_Branch,Leaf; //side1,side2,side3,side4,side5,side6,side7,side8;
	int i,j;
	long double X0=108,Y0=100,XMAX=216,YMAX=279,Y_crown,x,y;
	
	
	
	for (i=0;i<=1000;i++) side[i]=1;
	
// for shoots
if(FRACTAL==2) //shoot is fractal
{
	Diam_term_sapwood=Diam_Branch;
	Diam_Trunk_sapwood=pow(Diam_Trunk*Trunk_sapwood_fraction*(2*Diam_Trunk-Diam_Trunk*Trunk_sapwood_fraction),0.5);
	if (!Length_term_shoot) Length_term_shoot=(Length_Branch+Length_Trunk)/SEGMENTS; else Length_term_shoot=Length_term_shoot/100;
	Nb_term_shoot= LA_max_init / LA_term_shoot*10000;
	Nb_segment[0]=Nb_term_shoot;
	Diam_segment[0]=Diam_Branch;
	Length_segment[0]=Length_term_shoot;
	Branch_Area_FR=Length_term_shoot*3.1416*Diam_segment[0]*Nb_segment[0];
	Trunk_Area_FR=0;
	Volume_Branch=Nb_segment[0]*Diam_segment[0]*Diam_segment[0]*3.1416/4*Length_segment[0];
	Volume_Branch_sapwood=Nb_segment[0]*Diam_term_sapwood*Diam_term_sapwood*3.1416/4*Length_segment[0];
	Volume_Trunk=0;
	Volume_Trunk_sapwood=0;

	Q_Branch_Symp_FR=Volume_Branch_sapwood*Branch_Symp_fraction*1000; //in Kg
	Q_Branch_Apo_FR=Volume_Branch_sapwood*Branch_Apo_fraction*1000;
	Q_Trunk_Sym_FR=Volume_Trunk_sapwood*Trunk_Symp_fraction*1000;
	Q_Trunk_Apo_FR=Volume_Trunk_sapwood*Trunk_Apo_fraction*1000;
	R_Branch=Length_segment[0]/(Diam_term_sapwood*Diam_term_sapwood*3.1416/4)*1/K_Branch_Apo0[0]/Nb_term_shoot*pow((Diam_term_sapwood/Diam_term_sapwood),Tapering);
	R_Trunk=0;
	
	F=2*(Length_Branch+Length_Trunk-Length_term_shoot*(SEGMENTS+1))/((SEGMENTS+1)*SEGMENTS);
	Discriminant=pow((F+2*Length_term_shoot),2)+8*F*Length_Branch;
	Node_base_crown=(-(F+2*Length_term_shoot)+pow(Discriminant,0.5))/(2*F);
	N_fractal=pow(Nb_term_shoot,1/Node_base_crown);
	Delta1=Node_base_crown*log(N_fractal)/log(Diam_Trunk/Diam_Branch);
	Delta2=Node_base_crown*log(N_fractal)/log(Diam_Trunk_sapwood/Diam_term_sapwood);
	
	for(i=1;i<(int)SEGMENTS;i++)
	{		
		Nb_segment[i]=max(Nb_segment[i-1]/N_fractal,1);
		Length_segment[i]=Length_term_shoot+i*F;
		if ((long double)i<Node_base_crown) // for the Branches
		{
				Diam_segment[i]=(Diam_Branch)*pow(N_fractal,((long double)i/Delta1));
				Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,((long double)i/Delta2));
				Volume_Branch+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
				Volume_Branch_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
				Branch_Area_FR+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
				R_Branch+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Branch_Apo0[0]/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);	
		}
		else // for the Trunk
		{
			Diam_segment[i]=(Diam_Branch)*pow(N_fractal,(Node_base_crown/Delta1));
			Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(Node_base_crown/Delta2));
			Trunk_Area_FR+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
			Volume_Trunk+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
			Volume_Trunk_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
			R_Trunk+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Branch_Apo0[0]/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
		}		
	}
	
	Q_Branch_Symp_FR=Volume_Branch_sapwood*Branch_Symp_fraction*1000;  //in Kg
	Q_Branch_Apo_FR=Volume_Branch_sapwood*Branch_Apo_fraction*1000;
	Q_Trunk_Sym_FR=Volume_Trunk_sapwood*Trunk_Symp_fraction*1000;
	Q_Trunk_Apo_FR=Volume_Trunk_sapwood*Trunk_Apo_fraction*1000;
	K_Branch_Apo_FR=1/R_Branch;
	K_Trunk_Apo_FR=1/R_Trunk;	
}

// Fractal map1
    if (PRINT_GRAPH==4 && FRACTAL==2) //fractal shoot
	{
	srand((unsigned)time(NULL));
	YMAX=200;
	angle=3.1416/8;
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	Leaf=pow(LA_term_shoot,0.5)/100;
	File_fractal1 = fopen("fractal_map1.ps", "w+" );   
	fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal1);
	fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal1);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal1);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal1);

	//soil limit
	fprintf(File_fractal1,"stroke 0.5 M setlinewidth 1 0 0 setrgbcolor  newpath ");
	fprintf(File_fractal1,"0 M %Le M %Le M %Le M T\n",Y[0],XMAX,Y[0]);

	//Trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal1);
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //Trunk
	{
		fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_segment[i]*scale);
		Y[0]+=Length_segment[i]*scale;
	}
	Y_crown=Y[0];

	//Branch wood
	fputs("%Branches\n",File_fractal1);
	fprintf(File_fractal1,"stroke 0.5 0 0 setrgbcolor newpath ");
	for(i=(int)(Node_base_crown+1);i>=0;i--)  //ordre 1
		{
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_segment[i]*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_segment[i]*scale);
			Y[0]+=Length_segment[i]*scale;
			//	if ((rand()%1000)<((N_fractal-1)*1000)) 
			//	if ((i)%(int)((N_fractal-1)*100)==(int)(Node_base_crown)%(int)((N_fractal-1)*100))
			{
				side[0]=-side[0];
				Fractal_loop_shoot(1,i);
			}
		}
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",pow(LA_term_shoot,0.5)/100*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0] ,X[0],Y[0]+ pow(LA_term_shoot,0.5)/100*scale);	
	
		
	//heartwood
	
	fputs("\n%%Heartwood\n",File_fractal1);
	X[1]=X0;
	Y[1]=Y0;		
	for(i=(int)(SEGMENTS-1);i>=0;i--)
		{
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.67 0 0 setrgbcolor newpath ",pow(Diam_segment[i]*Diam_segment[i]-Diam_sapwood[i]*Diam_sapwood[i],0.5) *scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[1],Y[1], X[1],Y[1]+Length_segment[i]*scale);
			Y[1]+=Length_segment[i]*scale;
		}	
	
// Fractal map2
//shoot 
	File_fractal2 = fopen("fractal_map2.ps", "w+" );   
	fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal2);
	fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal2);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal2);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal2);

	X[0]=X0;Y[0]=Y0;
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal2);
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //Trunk
	{
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_segment[i]*scale);
		Y[0]+=Length_segment[i]*scale;
	}
	Y_crown=Y[0];
	fputs("%Branches\n",File_fractal2);
	for(i=(int)(Node_base_crown);i>=0;i--)  //ordre 1
		{
			fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_segment[i]*scale);
			fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_segment[i]*scale);
			Y[0]+=Length_segment[i]*scale;
		}
	
	l_Branch=0;
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(Node_base_crown);i>=0;i--) //Branches
	{
		for (j=1;j<=Nb_segment[i];j++)
		{		
		x=X[0]+l_Branch*scale*cos((long double)j*3.1416/Nb_segment[i]);
		y=Y_crown+l_Branch*scale*sin((long double)j*3.1416/Nb_segment[i]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[i]*scale*cos((long double)j*3.1416/Nb_segment[i]),y+Length_segment[i]*scale*sin((long double)j*3.1416/Nb_segment[i]));
		}
		l_Branch+=Length_segment[i];
	}
	
	for (j=0;j<Nb_segment[0];j++)//leaves
	{	
		x=X[0]+(l_Branch)*scale*cos((long double)j*3.1416/Nb_segment[0]);
		y=Y_crown+(l_Branch)*scale*sin((long double)j*3.1416/Nb_segment[0]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",Leaf*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Leaf*scale*cos((long double)j*3.1416/Nb_segment[0]),y+Leaf*scale*sin((long double)j*3.1416/Nb_segment[0]));	
	}
			
}

if (PRINT_GRAPH==4 && FRACTAL==1) //shoot is morphometric
{
	YMAX=200;	
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	Leaf=pow(LA_max_init/Number_Branch,0.5);
	File_fractal1 = fopen("fractal_map1.ps", "w+" );   
	fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal1);
	fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal1);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal1);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal1);
	
	//soil limit

	fprintf(File_fractal1,"stroke 0.5 M setlinewidth 1 0 0 setrgbcolor  newpath ");
	fprintf(File_fractal1,"0 M %Le M %Le M %Le M T\n",Y[0],XMAX,Y[0]);
	//Trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal1);
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_Trunk*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
	Y[0]+=Length_Trunk*scale;
	

	//heartwood
	fputs("\n%%Heartwood\n",File_fractal1);
	X[0]=X0;Y[0]=Y0;	
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.67 0 0 setrgbcolor newpath ",Diam_Trunk*Trunk_sapwood_fraction*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
			
	//Branch wood
	Y[0]+=Length_Trunk*scale;
	fputs("%Branches\n",File_fractal1);
	for(i=1;i<=Number_Branch;i++)  //ordre 1
		{
			angle=3.1416/(Number_Branch+1);
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_Branch*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0]+Length_Branch*scale*cos(3.1416-(long double)i*angle),Y[0]+Length_Branch*scale*sin(3.1416-(long double)i*angle));
			//Leaf
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",Leaf*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0]+Length_Branch*scale*cos(3.1416-(long double)i*angle),Y[0]+Length_Branch*scale*sin(3.1416-(long double)i*angle),X[0]+(Length_Branch+Leaf)*scale*cos(3.1416-(long double)i*angle),Y[0]+(Length_Branch+Leaf)*scale*sin(3.1416-(long double)i*angle));	

		}
	fputs("\nstroke\nshowpage\n%%Trailer\n",File_fractal1);
	fclose(File_fractal1);
	}	
		
// for one Root
	//A_term_Root = LA_max_init/3; //assume the area of fine Roots is equal to the Leaf area cm2
	A_term_Root =Root_Area_fi_0;
	Length_Root=(Length_Branch+Length_Trunk)*Root_shoot_ratio;
	Length_tapRoot= Length_Root*(1-Root_ramif);
	
	Diam_term_sapwood=Diam_Root ;
	Diam_tapRoot=Diam_Trunk; //assume the cross section of the 3 tapRoots=dbh
	Diam_tapRoot_sapwood=pow(Diam_tapRoot*Trunk_sapwood_fraction*(2*Diam_tapRoot-Diam_tapRoot*Trunk_sapwood_fraction),0.5);  //assume same sapwood fraction in the tapRoot that in the Trunk
	//Length_term_Root=Length_Root/SEGMENTS; 
	Length_term_Root=Length_fine_Root/100;
	Nb_term_Root= A_term_Root / (Length_fine_Root/100*Diam_Root*3.1416);
	Nb_segment[0]=Nb_term_Root;
	N_Root=Nb_segment[0];
	
	Diam_segment[0]=Diam_Root;
	Mean_diam_Root=Diam_segment[0]*Nb_segment[0];
	Length_segment[0]=Length_term_Root;
	FineRoot_area=Length_term_Root*3.1416*Diam_segment[0]*Nb_segment[0];
	TapRoot_area=0;
	Volume_fineRoot=Nb_segment[0]*Diam_segment[0]*Diam_segment[0]*3.1416/4*Length_segment[0];
	Volume_fineRoot_sapwood=Nb_segment[0]*Diam_term_sapwood*Diam_term_sapwood*3.1416/4*Length_segment[0];
	Volume_tapRoot=0;
	Volume_tapRoot_sapwood=0;
	Q_Root_Symp_FR=Volume_fineRoot_sapwood*Root_Symp_fraction*1000+Volume_tapRoot_sapwood*Trunk_Symp_fraction*1000;  //in Kg
	Q_Root_Apo_FR=Volume_fineRoot_sapwood*Root_Apo_fraction*1000 + Volume_tapRoot_sapwood*Trunk_Apo_fraction*1000;
	R_fineRoot=Length_segment[0]/(Diam_term_sapwood*Diam_term_sapwood*3.1416/4)*1/K_Root_Apo0/Nb_term_Root*pow((Diam_term_sapwood/Diam_term_sapwood),Tapering);
	R_tapRoot=0;
	
	F=2*(Length_Root-Length_term_Root*(SEGMENTS+1))/((SEGMENTS+1)*SEGMENTS);
	Discriminant=pow((F+2*Length_term_Root),2)+8*F*(Length_Root-Length_tapRoot);
	Node_base_crown=(-(F+2*Length_term_Root)+pow(Discriminant,0.5))/(2*F);
	N_fractal=pow(Nb_term_Root,1/Node_base_crown);
	Delta1=Node_base_crown*log(N_fractal)/log(Diam_tapRoot/Diam_Root);
	Delta2=Node_base_crown*log(N_fractal)/log(Diam_tapRoot_sapwood/Diam_term_sapwood);
	Length_Root_FR=Nb_segment[0]*Length_segment[0];
	
	for(i=1;i<(int)SEGMENTS;i++)
	{
		Nb_segment[i]=max(Nb_segment[i-1]/N_fractal,1);
		N_Root+=Nb_segment[i];
		Length_segment[i]=Length_term_Root+i*F;
		if ((long double)i<Node_base_crown) // for the fineRoot
		{
				Diam_segment[i]=(Diam_Root)*pow(N_fractal,((long double)i/Delta1));
				Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,((long double)i/Delta2));
				Volume_fineRoot+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
				Volume_fineRoot_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
				FineRoot_area+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
				R_fineRoot+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Root_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
				Length_Root_FR+=Nb_segment[i]*Length_segment[i];
				Mean_diam_Root+=Diam_segment[i]*Nb_segment[i];
		}
		else // for the tapRoot
		{
			Diam_segment[i]=(Diam_Root)*pow(N_fractal,(Node_base_crown/Delta1));
			Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(Node_base_crown/Delta2));
			TapRoot_area+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
			Volume_tapRoot+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
			Volume_tapRoot_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
			R_tapRoot+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Root_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
			Length_Root_FR+=Nb_segment[i]*Length_segment[i];
			Mean_diam_Root+=Diam_segment[i]*Nb_segment[i];
		}	
	}
	Q_Root_Symp_FR=Volume_fineRoot_sapwood*Root_Symp_fraction*1000+Volume_tapRoot_sapwood*Trunk_Symp_fraction*1000;  //in Kg
	Q_Root_Apo_FR=Volume_fineRoot_sapwood*Root_Apo_fraction*1000 + Volume_tapRoot_sapwood*Trunk_Apo_fraction*1000;
	K_Root_Apo_FR=1/(R_fineRoot+R_tapRoot);
	Root_Area_FR=FineRoot_area+TapRoot_area;
	Diam_Root_FR=Mean_diam_Root/N_Root;
	//Root_Area_fi_0=A_term_Root;
	Length_Root_fi=Root_Area_fi_0/(3.1416*Diam_Root);

if (PRINT_GRAPH==4) // Fractal map2
{
	if (FRACTAL==2) //print fractal Root in File_fractal1
	{
	angle=3.1416/8;
	X[0]=X0;Y[0]=Y0;
	fputs("%TapRoot\n",File_fractal1);
	fprintf(File_fractal1,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //tap Root
	{
		fprintf(File_fractal1,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]-Length_segment[i]*scale);
		Y[0]-=Length_segment[i]*scale;
	}
	Y_crown=Y[0];
	fputs("%Fine Roots\n",File_fractal1);
	l_Branch=Length_segment[0];
	fprintf(File_fractal1,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(Node_base_crown);i>=0;i--)  
		{
			fprintf(File_fractal1,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]-Length_segment[i]*scale);
			Y[0]-=Length_segment[i]*scale;
			//if ((rand()%1000)<((N_fractal-1)*1000)) 
			//if ((i)%(int)((N_fractal-1)*100)==(int)(Node_base_crown)%(int)((N_fractal-1)*100))
			{
				side[0]=-side[0];
				Fractal_loop_Root(1,i);
			}
		}	
	fputs("\nstroke\nshowpage\n%%Trailer\n",File_fractal1);
	fclose(File_fractal1);	

	}
	if (FRACTAL==1)
	{
		File_fractal2 = fopen("fractal_map2.ps", "w+" );   
		fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal2);
		fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal2);
		fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal2);
		fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal2);
	}

	//soil limit
	X[0]=X0;Y[0]=Y0;
	fprintf(File_fractal2,"stroke 0.5 M setlinewidth 1 0 0 setrgbcolor  newpath ");
	fprintf(File_fractal2,"0 M %Le M %Le M %Le M T\n",Y[0],XMAX,Y[0]);
//Root	
	fputs("%TapRoot\n",File_fractal2);
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //tap Root
	{
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]-Length_segment[i]*scale);
		Y[0]-=Length_segment[i]*scale;
	}
	Y_crown=Y[0];
	fputs("%Fine Roots\n",File_fractal2);
	l_Branch=Length_segment[0];
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
		for(i=(int)(Node_base_crown);i>=0;i--) //fine Roots
		{
			for (j=1;j<=Nb_segment[i];j++)
			{		
			x=X[0]+l_Branch*scale*cos((long double)j*3.1416/Nb_segment[i]);
			y=Y_crown-l_Branch*scale*sin((long double)j*3.1416/Nb_segment[i]);
			fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
			fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[i]*scale*cos((long double)j*3.1416/Nb_segment[i]),y-Length_segment[i]*scale*sin((long double)j*3.1416/Nb_segment[i]));
			}
			l_Branch+=Length_segment[i];
		}
	fprintf(File_fractal2,"stroke  0 0.5 0 setrgbcolor newpath ");	
	for (j=0;j<Nb_term_Root;j++)//terminal Roots
	{	
		x=X[0]+(l_Branch)*scale*cos((long double)j*3.1416/Nb_segment[0]);
		y=Y_crown-(l_Branch)*scale*sin((long double)j*3.1416/Nb_segment[0]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_Root*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_term_Root*scale*cos((long double)j*3.1416/Nb_segment[0]),y-Length_term_Root*scale*sin((long double)j*3.1416/Nb_segment[0]));	
	}
	

if (FRACTAL==1)
	{
	YMAX=200;	
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	Leaf=pow(LA_max_init/Number_Branch,0.5);

	//Trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal2);
	fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_Trunk*scale);
	fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
	Y[0]+=Length_Trunk*scale;
	

	//heartwood
	fputs("\n%%Heartwood\n",File_fractal1);
	X[0]=X0;Y[0]=Y0;	
	fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.67 0 0 setrgbcolor newpath ",Diam_Trunk*Trunk_sapwood_fraction*scale);
	fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
			
	//Branch wood
	Y[0]+=Length_Trunk*scale;
	fputs("%Branches\n",File_fractal2);
	for(i=1;i<=Number_Branch;i++)  //ordre 1
	{
		angle=3.1416/(Number_Branch+1);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_Branch*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0]+Length_Branch*scale*cos(3.1416-(long double)i*angle),Y[0]+Length_Branch*scale*sin(3.1416-(long double)i*angle));
		//Leaf
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",Leaf*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0]+Length_Branch*scale*cos(3.1416-(long double)i*angle),Y[0]+Length_Branch*scale*sin(3.1416-(long double)i*angle),X[0]+(Length_Branch+Leaf)*scale*cos(3.1416-(long double)i*angle),Y[0]+(Length_Branch+Leaf)*scale*sin(3.1416-(long double)i*angle));	
	}
	}	

	fputs("\nstroke\nshowpage\n%%Trailer\n",File_fractal2);
	fclose(File_fractal2);
	}
}

//***************************//
// SHORT TIME STEP FUNCTIONS //
//***************************//
void Compute_dq(long double dt_court)  // in mmol H20 exchanged between compartments during dt
{
	int i;
	
	for(i=1;i<4;i++)
	{
									dq_Leaf_Symp_Leaf_Evap[i]       	=   dt_court * K_Leaf_Symp2[i]   		*   (P_Leaf_Symp[i]    -  P_Leaf_Evap[i]); 
		if (PLC_Leaf_Apo[i]<100)   	dq_Leaf_Apo_Leaf_Evap[i]      	=   dt_court * K_Leaf_Symp[i]   		*   (P_Leaf_Apo[i]     -  P_Leaf_Evap[i]); 							else dq_Leaf_Apo_Leaf_Evap[i]=0;
		if (PLC_Branch_Apo[i]<100) 	dq_Branch_Apo_Leaf_Apo[i]        	=   dt_court * K_Leaf_Apo[i]     		*   (P_Branch_Apo[i]   -  Pg_Branch[i] - P_Leaf_Apo[i] + Pg_Leaf[i]); 	else dq_Branch_Apo_Leaf_Apo[i]=0;
		if (PLC_Branch_Apo[i]<100) 	dq_Branch_Apo_Branch_Symp[i]    	=   dt_court * K_Branch_Symp[i]  		*   (P_Branch_Apo[i]   -  P_Branch_Symp[i]); 							else dq_Branch_Apo_Branch_Symp[i]=0;		
		if (PLC_Trunk_Apo<100)  		dq_Trunk_Apo_Branch_Apo[i]       	=   dt_court * K_Branch_Apo[i]   		*   (P_Trunk_Apo    	 -  Pg_Trunk - P_Branch_Apo[i] + Pg_Branch[i]); 	else dq_Trunk_Apo_Branch_Apo[i]=0;
	}
	
	if (PLC_Trunk_Apo<100)  	dq_Trunk_Apo_Trunk_Symp    	=   dt_court * K_Trunk_Symp   *   (P_Trunk_Apo     -  P_Trunk_Symp); else dq_Trunk_Apo_Trunk_Symp=0;
							dq_Root_Apo_Trunk_Apo      	=   dt_court*K_Trunk_Apo      *   (P_Root_Apo      -  P_Trunk_Apo + Pg_Trunk);
			
	if (Type_Axil) // flow for one Axillary organ
	{
		if (PLC_Axil_Apo<100)  		dq_Axil_Apo_Axil_Symp       =   dt_court*K_Axil_Symp   *   (P_Axil_Apo      	-  P_Axil_Symp);      else dq_Axil_Apo_Axil_Symp=0;
									dq_Branch_Symp_Axil_Symp    =   dt_court*K_Axil_Symp2  *   (P_Branch_Symp[1]  	-  P_Axil_Symp);
		if (PLC_Branch_Apo[1]<100) 	dq_Branch_Apo_Axil_Apo    	 =   dt_court*K_Axil_Apo    *   (P_Branch_Apo[1] 	-  P_Axil_Apo);       else dq_Branch_Apo_Axil_Apo =0;
		if (PLC_Axil_Apo<100)  		dq_Axil_Apo_Petiole_Symp    =   dt_court*K_Axil_Symp2  *   (P_Axil_Apo      	-  P_Petiole_Symp);   else dq_Axil_Apo_Petiole_Symp=0;
	}

	if (Root_upper)
	{
		if (PLC_Root_Apo1<100) 	dq_Root_Apo_Trunk_Apo1   =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_Trunk); else dq_Root_Apo_Trunk_Apo1 =0;
		if (PLC_Root_Apo1<100) 	dq_Root_Endo_Root_Apo1   =   dt_court*K_Root_Apo1    *   (P_Root_Endo1     -  P_Root_Apo); else  dq_Root_Endo_Root_Apo1=0;
								dq_Root_Symp_Root_Endo1  =   dt_court*K_Root_Symp21  *   (P_Root_Symp1     -  P_Root_Endo1);
								dq_Soil_Root_Endo1       =   dt_court*(1/(1/K_Soil1 +1/K_Root_Symp11 +1/K_Interface1))       *  (P_Soil1 -  P_Root_Endo1);
	}
	else 
	{
		dq_Root_Apo_Trunk_Apo1 	= 0;
		dq_Root_Endo_Root_Apo1	= 0;
		dq_Root_Symp_Root_Endo1  	= 0;
		dq_Soil_Root_Endo1       	= 0;
	}
	
	if (Root_middle)
	{
		if (PLC_Root_Apo2<100) 	dq_Root_Apo_Trunk_Apo2      =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_Trunk); else dq_Root_Apo_Trunk_Apo2 =0;
		if (PLC_Root_Apo2<100) 	dq_Root_Endo_Root_Apo2      =   dt_court*K_Root_Apo2    *   (P_Root_Endo2     -  P_Root_Apo); else dq_Root_Endo_Root_Apo2=0;
								dq_Root_Symp_Root_Endo2     =   dt_court*K_Root_Symp22  *   (P_Root_Symp2     -  P_Root_Endo2);
								dq_Soil_Root_Endo2          =   dt_court*(1/(1/K_Soil2 +1/K_Root_Symp12 +1/K_Interface2))       *  (P_Soil2 -  P_Root_Endo2);
	}
	else 
	{
		dq_Root_Apo_Trunk_Apo2 	= 0;
		dq_Root_Endo_Root_Apo2	= 0;
		dq_Root_Symp_Root_Endo2  	= 0;
		dq_Soil_Root_Endo2       	= 0;
	}
	
	if (Root_lower)
	{
	    if (PLC_Root_Apo3<100) 	dq_Root_Apo_Trunk_Apo3      =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_Trunk); else dq_Root_Apo_Trunk_Apo3 =0;
		if (PLC_Root_Apo3<100) 	dq_Root_Endo_Root_Apo3      =   dt_court*K_Root_Apo3    *   (P_Root_Endo3     -  P_Root_Apo); else dq_Root_Endo_Root_Apo3=0;
								dq_Root_Symp_Root_Endo3     =   dt_court*K_Root_Symp23  *   (P_Root_Symp3     -  P_Root_Endo3);
								dq_Soil_Root_Endo3          =   dt_court*(1/(1/K_Soil3 +1/K_Root_Symp13 +1/K_Interface3))       *  (P_Soil3 -  P_Root_Endo3);
	}
	else 
	{
		dq_Root_Apo_Trunk_Apo3 	= 0;
		dq_Root_Endo_Root_Apo3	= 0;
		dq_Root_Symp_Root_Endo3  	= 0;
		dq_Soil_Root_Endo3       	= 0;
	}
	
	if (CAPILLARITY==1 || CAPILLARITY==11) // if salt in the soil the equation should be modified !
	{
		dq_Soil_12                  =   dt_court*K_Soil_tot1*K_Soil_tot2/(K_Soil_tot1+K_Soil_tot2) *(P_Soil1-P_Soil2);
		dq_Soil_23                  =   dt_court*K_Soil_tot2*K_Soil_tot3/(K_Soil_tot2+K_Soil_tot3) *(P_Soil2-P_Soil3);
	}
	else dq_Soil_12 = dq_Soil_23 =0;
}

void Compute_Q(void)  // changes in water content of the different compartments
{
	int i;
	for (i=1;i<4;i++) if (Branch_distri[i])
	{
		Q_Leaf_Evap[i]+=           dq_Leaf_Apo_Leaf_Evap[i]       +   dq_Leaf_Symp_Leaf_Evap[i]          - dq_stomate[i];
		Q_Leaf_Symp[i]+=        -  dq_Leaf_Symp_Leaf_Evap[i]      -   dq_cuti[i]  ;
		Q_Leaf_Apo[i]+=         -  dq_Leaf_Apo_Leaf_Evap[i]       +   dq_Branch_Apo_Leaf_Apo[i];
	}
		
	if (Type_Axil)				// if there is an Axillary organ on branch 1 only
	{
		if (Type_Axil==3) 		// a flower
		{
			Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp       -   dq_Axil;
			Q_Petiole_Symp+=       dq_Axil_Apo_Petiole_Symp    -   dq_Petiole;
			Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo      -   dq_Axil_Apo_Axil_Symp  -  dq_Axil_Apo_Petiole_Symp;
			Q_Branch_Symp[1]+=     dq_Branch_Apo_Branch_Symp[1]   -   dq_Branch[1];
			Q_Branch_Apo[1]+=    - dq_Branch_Apo_Leaf_Apo[1]      -   dq_Branch_Apo_Branch_Symp[1]       + dq_Trunk_Apo_Branch_Apo[1] -  dq_Branch_Apo_Axil_Apo;
		}
		else if (Type_Axil==2) // a Fruit
		{
			Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp      	 -   dq_Axil   + dq_Fruit;
			Q_Petiole_Symp+=       dq_Axil_Apo_Petiole_Symp  	 -   dq_Petiole;
			Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo        -   dq_Axil_Apo_Axil_Symp         - dq_Fruit;
			Q_Branch_Symp[1]+=     dq_Branch_Apo_Branch_Symp[1]  -   dq_Branch[1];
			Q_Branch_Apo[1]+=    - dq_Branch_Apo_Leaf_Apo[1]     -   dq_Branch_Apo_Branch_Symp[1]     + dq_Trunk_Apo_Branch_Apo[1] -  dq_Branch_Apo_Axil_Apo;
		}
		else if (Type_Axil==1) // a bud
		{
			Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp      +   dq_Branch_Symp_Axil_Symp       - dq_Axil  ;
			Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo     -   dq_Axil_Apo_Axil_Symp        ;
			Q_Branch_Symp[1]+=        dq_Branch_Apo_Branch_Symp[1]  -   dq_Branch_Symp_Axil_Symp       - dq_Branch[1];
			Q_Branch_Apo[1]+=       - dq_Branch_Apo_Leaf_Apo[1]     -   dq_Branch_Apo_Branch_Symp[1]       + dq_Trunk_Apo_Branch_Apo[1] -   dq_Branch_Apo_Axil_Apo;
		}
	}
	
	else //no organ
	{
		for (i=1;i<4;i++)
		{
		Q_Branch_Apo[i]+=       -  dq_Branch_Apo_Leaf_Apo[i]      -   dq_Branch_Apo_Branch_Symp[i]       + dq_Trunk_Apo_Branch_Apo[i];
		Q_Branch_Symp[i]+=         dq_Branch_Apo_Branch_Symp[i]   -   dq_Branch[i];
		}
	}
	
	Q_Trunk_Symp+=          dq_Trunk_Apo_Trunk_Symp     -   dq_Trunk;
	Q_Trunk_Apo+=        -  dq_Trunk_Apo_Branch_Apo[1]  -   dq_Trunk_Apo_Branch_Apo[2] -  dq_Trunk_Apo_Branch_Apo[3]     -   dq_Trunk_Apo_Trunk_Symp      + dq_Root_Apo_Trunk_Apo;
	Q_Root_Apo1+=        -  dq_Root_Apo_Trunk_Apo1      +   dq_Root_Endo_Root_Apo1 		;
	Q_Root_Apo2+=        -  dq_Root_Apo_Trunk_Apo2      +   dq_Root_Endo_Root_Apo2 		;
	Q_Root_Apo3+=        -  dq_Root_Apo_Trunk_Apo3      +   dq_Root_Endo_Root_Apo3 		;
	Q_Root_Apo_t+=		-  dq_Root_Apo_Trunk_Apo	   	+   dq_Root_Endo_Root_Apo1	    +  dq_Root_Endo_Root_Apo2 +   dq_Root_Endo_Root_Apo3; 
	Q_Root_Endo1+=       -  dq_Root_Endo_Root_Apo1      +   dq_Root_Symp_Root_Endo1       + dq_Soil_Root_Endo1;
	Q_Root_Endo2+=       -  dq_Root_Endo_Root_Apo2      +   dq_Root_Symp_Root_Endo2       + dq_Soil_Root_Endo2;
	Q_Root_Endo3+=       -  dq_Root_Endo_Root_Apo3      +   dq_Root_Symp_Root_Endo3       + dq_Soil_Root_Endo3;
	Q_Root_Symp1+=       -  dq_Root_Symp_Root_Endo1     -   dq_Root1;
	Q_Root_Symp2+=       -  dq_Root_Symp_Root_Endo2     -   dq_Root2;
	Q_Root_Symp3+=       -  dq_Root_Symp_Root_Endo3     -   dq_Root3;
	
	if (COMPET)
	{
		Q_Soil1+=            -  dq_Soil_Root_Endo1*COMPET          -   dq_Soil                         - dq_Soil_12;
		Q_Soil2+=            -  dq_Soil_Root_Endo2*COMPET          +   dq_Soil_12                      - dq_Soil_23;
		Q_Soil3+=            -  dq_Soil_Root_Endo3*COMPET          +   dq_Soil_23;
	}
	else
	{
		Q_Soil1+=            -  dq_Soil_Root_Endo1         -   dq_Soil                         - dq_Soil_12;
		Q_Soil2+=            -  dq_Soil_Root_Endo2         +   dq_Soil_12                      - dq_Soil_23;
		Q_Soil3+=            -  dq_Soil_Root_Endo3         +   dq_Soil_23;
	}
}

void Compute_P_dynamic(void)   // compute water potentials
{
	long double Tp, Pi, Rs;
	int i;
	

	for (i=1;i<4;i++) 
		if (Branch_distri[i]) 
		{
			//Branch
			Rs=(Q_Branch_Symp0[i]-Q_Branch_Symp[i])/Q_Branch_Symp0[i];
			Tp=-Pi0_Branch_Symp*Osmotic_TAir - Epsilon_Branch_Symp*Rs;
			if (Tp<0)Tp=0;
			if (1-Rs) Pi=Pi0_Branch_Symp*Osmotic_TAir/(1-Rs); else Pi=-1000;
			P_Branch_Symp[i]=Tp+Pi;
			P_Branch_Apo[i]=(Q_Branch_Apo[i]-Q_Branch_Apo1[i])/C_Branch_Apo[i]+Pg_Branch[i];
			Turgor_Branch_Symp[i]=Tp;
			
			if (Leaf_Area[i])//Leaf
			{
				P_Leaf_Evap[i]=(Q_Leaf_Evap[i]-Q_Leaf_Evap0[i])/C_Leaf_Evap[i] + Pg_Leaf[i];
				Rs=(Q_Leaf_Symp0[i]-Q_Leaf_Symp[i])/Q_Leaf_Symp0[i];
				Tp=-Pi0_Leaf_Symp*Osmotic_TLeaf[i] - Epsilon_Leaf_Symp*Rs;
				if (Tp<0)Tp=0;
				if (1-Rs) Pi=Pi0_Leaf_Symp*Osmotic_TLeaf[i]/(1-Rs); else Pi=-1000;
				P_Leaf_Symp[i]=Tp+Pi;
				P_Leaf_Apo[i]=(Q_Leaf_Apo[i]-Q_Leaf_Apo1[i])/C_Leaf_Apo[i]+Pg_Leaf[i];
				Turgor_Leaf_Symp[i]=Tp;
				}
			else 
			{
				P_Leaf_Apo[i]=P_Leaf_Symp[i]=P_Leaf_Evap[i]=P_Branch_Apo[i];
			}
			
		}
	
	//Axil
	if (Type_Axil) // not implemented yet with 3 Branches
	{
		Rs=(Q_Axil_Symp0-Q_Axil_Symp)/Q_Axil_Symp0;
		Tp=-Pi0_Axil_Symp*Osmotic_TAir - Epsilon_Axil_Symp*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Axil_Symp* Osmotic_TAir /(1-Rs); else Pi=-1000;
		P_Axil_Symp=Tp+Pi;
		P_Axil_Apo=(Q_Axil_Apo-Q_Axil_Apo1)/C_Axil_Apo+Pg_Branch[1];
		Turgor_Axil_Symp=Tp;
		
		//Axil Petiole
		if (Type_Axil>=2) // for flower or Fruit with a petiole only
		{
			Rs=(Q_Petiole_Symp0-Q_Petiole_Symp)/Q_Petiole_Symp0;
			Tp=-Pi0_Axil_Symp*Osmotic_TAir - Epsilon_Axil_Symp*Rs;
			if (Tp<0)Tp=0;
			if (1-Rs) Pi=Pi0_Axil_Symp* Osmotic_TAir /(1-Rs); else Pi=-1000;
			P_Petiole_Symp=Tp+Pi;
		}
	}
	else P_Axil_Symp=P_Axil_Apo=P_Petiole_Symp=P_Branch_Apo[1];
	
	//Trunk
	Rs=(Q_Trunk_Symp0-Q_Trunk_Symp)/Q_Trunk_Symp0;
	Tp=-Pi0_Trunk_Symp*Osmotic_TAir - Epsilon_Trunk_Symp*Rs;
	if (Tp<0)Tp=0;
	if (1-Rs) Pi=Pi0_Trunk_Symp*Osmotic_TAir/(1-Rs); else Pi=-1000;
	P_Trunk_Symp=Tp+Pi;
	P_Trunk_Apo=(Q_Trunk_Apo-Q_Trunk_Apo1)/C_Trunk_Apo+Pg_Trunk;
	Turgor_Trunk_Symp=Tp;
	
	//Roots
	P_Root_Apo=(Q_Root_Apo_t-Q_Root_Apo_t1)/C_Root_Apo;
	if (Root_upper)
	{
		Rs=(Q_Root_Symp01-Q_Root_Symp1)/Q_Root_Symp01;
		Tp=-Pi0_Root_Symp1*Osmotic_TSoil - Epsilon_Root_Symp1*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Root_Symp1*Osmotic_TSoil/(1-Rs); else Pi=-1000;
		P_Root_Symp1=Tp+Pi;
		P_Root_Apo1=(Q_Root_Apo1-Q_Root_Apo11)/C_Root_Apo1;
		P_Root_Endo1=(Q_Root_Endo1-Q_Root_Endo01)/C_Root_Endo1;
		Turgor_Root_Symp1=Tp;

	}
	else  P_Root_Symp1=P_Root_Endo1=P_Root_Apo1;
	
	if (Root_middle)
	{
		Rs=(Q_Root_Symp02-Q_Root_Symp2)/Q_Root_Symp02;
		Tp=-Pi0_Root_Symp2*Osmotic_TSoil - Epsilon_Root_Symp2*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Root_Symp2*Osmotic_TSoil/(1-Rs); else Pi=-1000;
		P_Root_Symp2=Tp+Pi;
		P_Root_Apo2=(Q_Root_Apo2-Q_Root_Apo12)/C_Root_Apo2;
		P_Root_Endo2=(Q_Root_Endo2-Q_Root_Endo02)/C_Root_Endo2;
		Turgor_Root_Symp3=Tp;
	}
	else  P_Root_Symp2=P_Root_Endo2=P_Root_Apo2;
	
	if (Root_lower)
	{
		Rs=(Q_Root_Symp03-Q_Root_Symp3)/Q_Root_Symp03;
		Tp=-Pi0_Root_Symp3*Osmotic_TSoil - Epsilon_Root_Symp3*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Root_Symp3*Osmotic_TSoil/(1-Rs); else Pi=-1000;
		P_Root_Symp3=Tp+Pi;
		P_Root_Apo3=(Q_Root_Apo3-Q_Root_Apo13)/C_Root_Apo3;
		P_Root_Endo3=(Q_Root_Endo3-Q_Root_Endo03)/C_Root_Endo3;
		Turgor_Root_Symp3=Tp;
	}
	else  P_Root_Symp3=P_Root_Endo3=P_Root_Apo3;
	P_Root_Apo=(Q_Root_Apo_t-Q_Root_Apo_t1)/C_Root_Apo;
	P_Root_Apo1=P_Root_Apo2=P_Root_Apo3=P_Root_Apo;			 	 
}

void print_screen(void)
{
	int i;
	
	Get_DATA(dt);
	printf("  %5.0Lf",YEAR2);
	if (PRINT_SCREEN==2) printf("%7.0Lf %7.2Lf", DATA[0],PLC_Leaf_Apo[0]);
	else
		for (i=0;i<NPAR;i++)
			if (Screen_out[i])
			{
				if (i==0) printf("%7.2Lf ", DATA[0]);
				else printf("%9.2Le  ", DATA[i]);
			}
	if (WARNING)
	{
		printf("  WARNING %Lf", dt);
		WARNING=0;
	}
	//printf (" %Lf ",DYNAMIC);
	printf("\n"); 
}

void print_transient(void)
{
	FILE *transient;
	int i;
	
	transient = fopen(filename_TRANS,"a");
	Get_DATA(dt);
	fprintf(transient,"%.0Lf; ",YEAR1);
	for (i=0;i<NPAR;i++)
		if (debug || File_out[i] )
		{
			if (i==0) fprintf(transient, "%.5Lf; ", DATA[i]);
			else fprintf(transient, "%.8Le; ", DATA[i]);
		}
	fprintf(transient,"\n");
	fclose(transient);
}

void print_annual(void)
{
	FILE *out;
	int i;
	
	if ((out = fopen(filename_OUT2,"a"))==NULL) printf("\nCan't create file annual.out!!");
	else
	{
		Get_DATA(dt);
		fprintf(out, "%s; %d;  %.0Lf; %.0Lf;", filenumber, N, Simul, YEAR1);
		for (i=0;i<NPAR;i++)  if (File_out[i])  fprintf(out, "%Lf; ", DATA[i]);
		fprintf(out,"\n");
	}
	fclose(out);
	
}

void print_final(void)
{
	int  random1; //para number
	long double random2, random3, random4; //CV, min% max%
	FILE *out, *random_file;
	
	if ((out = fopen(filename_OUT,"a"))==NULL) printf("\nCan't create file sureau.out!!");
	else 
	{
		Get_DATA(dt);
		fprintf(out,"%s; %d; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;", 
		filenumber,  N, Simul, YEAR1, T/3600/24/1000, T_PLC_Leaf[1]/3600/24/1000,T_PLC_Leaf[2]/3600/24/1000,T_PLC_Leaf[3]/3600/24/1000,T_PLC_Axil/3600/24/1000,T_PLC_Branch[1]/3600/24/1000,T_PLC_Branch[2]/3600/24/1000,T_PLC_Branch[3]/3600/24/1000,T_PLC_Trunk/3600/24/1000,T_PLC_Root/3600/24/1000,T_PLC_Root1/3600/24/1000,T_RWC_Axil/3600/24/1000,T_REW_Soil/3600/24/1000, T_gs_regul/3600/24/1000,T_gs_50mmol/3600/24/1000, T_gs_close/3600/24/1000, T_budbreak,T_max_LAI,E_tot*18/1000/1000,A_net_tot, GPP, Radius_Trunk_rel,  T_Leaf_max, PLC_Leaf_Apo[0],PLC_Leaf_Apo[1],PLC_Leaf_Apo[2],PLC_Leaf_Apo[3], PLC_Branch_Apo[0],PLC_Branch_Apo[1],PLC_Branch_Apo[2],PLC_Branch_Apo[3],REW_int1,REW_int2,RWC_int,RWC_min,P_soil_min,P_min_Leaf,P_min_stem,RWC_Branch_s_min[1],RWC_Branch_s_min[2],RWC_Branch_s_min[3],T_TLP_Leaf[1]/3600/24/1000,T_TLP_Axil/3600/24/1000,T_RWC_Branch[1]/3600/24/1000);
		if (RANDOMISE)
		{
		if ((random_file = fopen("random_para.txt","r"))==NULL) printf("\a\nCan't find file random_para.txt!!");
			else while (!feof(random_file)) 
			{
			fscanf(random_file,"%d %Le %Le %Le\n",&random1,&random2,&random3,&random4);
			fprintf(out," %Lf;", para[random1]);
			}
			fprintf(out,"\n");
			fclose(random_file);
		}
		else fprintf(out,"\n");
		fclose(out);
	 
	}	
}
void test_END(void)
	{
	// test for END
				if (END_DEATH== 1)  // dead when all Leaf xylem are fully embolised
					{
					if ((PLC_Leaf_Apo[1]>=PLC_END || !Branch_distri[1]) && (PLC_Leaf_Apo[2]>=PLC_END || !Branch_distri[2]) && (PLC_Leaf_Apo[3]>=PLC_END || !Branch_distri[3])) 
					DEAD=1;   
					}
				else if (END_DEATH== 2)	// dead when all Branches xylem are fully embolised
					{
					if ((PLC_Branch_Apo[1]>=PLC_END || !Branch_distri[1]) && (PLC_Branch_Apo[2]>=PLC_END || !Branch_distri[2]) && (PLC_Branch_Apo[3]>=PLC_END || !Branch_distri[3])) 
					DEAD=1;  
					}
				else if (END_DEATH== 14)	// dead when all leaves and Branches xylem are fully embolised
					{	
					if ((PLC_Leaf_Apo[1]>=PLC_END || !Branch_distri[1]) && (PLC_Leaf_Apo[2]>=PLC_END || !Branch_distri[2]) && (PLC_Leaf_Apo[3]>=PLC_END || !Branch_distri[3])) 
					if ((PLC_Branch_Apo[1]>=PLC_END || !Branch_distri[1]) && (PLC_Branch_Apo[2]>=PLC_END || !Branch_distri[2]) && (PLC_Branch_Apo[3]>=PLC_END || !Branch_distri[3])) 
					DEAD=1;   // printf("dead ici");
					}  	
				else if (END_DEATH== 3 && PLC_Trunk_Apo>=PLC_END) DEAD=1;   // dead when Trunk xylem is fully embolised
				else if (END_DEATH== 4 && (PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END)) DEAD=1;        // dead when ALL Root xylem is fully embolised
				else if (END_DEATH== 5 && (PLC_Leaf_Apo[1]>=PLC_END  || PLC_Leaf_Apo[2]>=PLC_END  || PLC_Leaf_Apo[3]>=PLC_END  || PLC_Branch_Apo[1]>=PLC_END || PLC_Branch_Apo[2]>=PLC_END || PLC_Branch_Apo[3]>=PLC_END || PLC_Axil_Apo>=PLC_END || PLC_Trunk_Apo>=PLC_END || PLC_Root_Apo1>=PLC_END || PLC_Root_Apo2>=PLC_END || PLC_Root_Apo3>=PLC_END)) DEAD=1; // dead when ONE organ is fully embolised
				else if (END_DEATH== 6 && (PLC_Leaf_Apo[1]>=PLC_END  && PLC_Leaf_Apo[2]>=PLC_END  && PLC_Leaf_Apo[3]>=PLC_END  && PLC_Branch_Apo[1]>=PLC_END && PLC_Branch_Apo[2]>=PLC_END && PLC_Branch_Apo[3]>=PLC_END && PLC_Axil_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
				else if (END_DEATH== 7 && (REW_t<=0.2)) DEAD=1;  // Stop when RWC<0.1
				else if (END_DEATH== 8 && PLC_Axil_Apo>=PLC_END) DEAD=1;
				else if (END_DEATH== 9 && (PLC_Axil_Apo>=PLC_END && PLC_Leaf_Apo[1]>=PLC_END)) DEAD=1;
				else if (END_DEATH==10 && (RWC_Axil<=0.02)) DEAD=1;  // Stop when 2% RWC
				else if (END_DEATH==11 && Type_Axil  && (RWC_Axil<=0.02) && (PLC_Leaf_Apo[1]>=PLC_END  &&  PLC_Leaf_Apo[2]>=PLC_END  &&PLC_Leaf_Apo[3]>=PLC_END  && PLC_Branch_Apo[1]>=PLC_END && PLC_Branch_Apo[2]>=PLC_END && PLC_Branch_Apo[3]>=PLC_END && PLC_Trunk_Apo>=PLC_END  &&  PLC_Axil_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END) ) DEAD=1;  // Stop when 2% RWC         
				else if (END_DEATH==11 && !Type_Axil && (PLC_Leaf_Apo[1]>=PLC_END  &&  PLC_Leaf_Apo[2]>=PLC_END  && PLC_Leaf_Apo[3]>=PLC_END  && PLC_Branch_Apo[1]>=PLC_END && PLC_Branch_Apo[2]>=PLC_END && PLC_Branch_Apo[3]>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END) ) DEAD=1;  // Stop when 2% RWC
				else if (END_DEATH==12 && Type_Axil  && (PLC_Leaf_Apo[1]>=PLC_END  &&  PLC_Leaf_Apo[2]>=PLC_END  && PLC_Leaf_Apo[3]>=PLC_END  && PLC_Branch_Apo[1]>=PLC_END && PLC_Branch_Apo[2]>=PLC_END && PLC_Branch_Apo[3]>=PLC_END && PLC_Axil_Apo>=PLC_END  && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
				else if (END_DEATH==12 && !Type_Axil && (PLC_Leaf_Apo[1]>=PLC_END  &&  PLC_Leaf_Apo[2]>=PLC_END  && PLC_Leaf_Apo[3]>=PLC_END  && PLC_Branch_Apo[1]>=PLC_END && PLC_Branch_Apo[2]>=PLC_END && PLC_Branch_Apo[3]>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
				else if (END_DEATH==15 && Type_Axil  && (Turgor_Leaf_Symp[1]==0 && Turgor_Axil_Symp==0)) DEAD=1; 
				else if (END_DEATH==15 && !Type_Axil && Turgor_Leaf_Symp[1]==0) DEAD=1; 
				else if (END_DEATH==16 && RWC_Branch_s[1]<=RWC_END) DEAD=1; // stop when branch symplasmic water content is less thant PLC_END. Typical value is 0.35
				
	}
	
void Test_Time(void)
{
	if (t_out==24 && (TRANSIENT==1 || TRANSIENT==11 || TRANSIENT==21) ) // if T_out is 24h then print only midday max values at the time of T_max HH1
		{
		if (!((unsigned long long)(T)%(unsigned long long)(3600*t_out*1000)-(unsigned long long)(3600*HH1*1000)))
			{
				if (PRINT_SCREEN) print_screen();
				if (TRANSIENT) print_transient();
			}
		}
			
	else if(!((unsigned long long)(T)%(unsigned long long)(3600*t_out*1000)))
		{
			if (PRINT_SCREEN) print_screen();
			if (TRANSIENT) print_transient();		
		}
	
}

//***************************//
//        MAIN FUNCTION      //
//***************************//

void compute (void) //the main function where the job is done
{
	long i=0;
	FILE *transient, *_soil;
	long double Q1,Q2, dt_long, threshold;
	char buffer[40];
	year++;
	Q1=Q_Trunk_Symp0;
	Q2=Q_Axil_Symp0;
	if (CLIMAT!=0 && CLIMAT!=5) load_climat();
	if (CLIMAT !=1) T0=DOY*24*3600*1000;
	T=T0;
	//printf("T0= %Lf %Lf ",T0,T);
	if (COMPET) 
		{
		sprintf(buffer, "%.0Lf", Simul/10);		// the soil.trs file is named with de decade of Simul							
		strcat(buffer,"_SOIL.trs");
		strcpy(filename_SOIL,buffer);
		COMPET_number=Simul-((int)Simul/10)*10;  // decimal part; if =0 then main tree
		if (COMPET_number==0) // only for the first simulation in case of competition
			{
			if((_soil = fopen(filename_SOIL,"r"))==NULL) //create the file if it doesn't exist
				{
				_soil = fopen(filename_SOIL,"w+");
				rewind (_soil);
				fprintf(_soil,"%Lf %Lf %Lf %Lf",T/1000, Q_Soil1,Q_Soil2,Q_Soil3);
				//printf("%Lf %Lf %Lf %Lf",(Q_Soil01-Q_Soil1 + Q_Soil02-Q_Soil2+Q_Soil03-Q_Soil3)*18/1000/1000, Q_Soil1*18/1000/1000,Q_Soil2*18/1000/1000,Q_Soil3*18/1000/1000);
				fclose(_soil);
				}
			else fclose(_soil);
			}
		}
			
	
	if (CLIMAT==0 || CLIMAT==5) printf("#%.0Lf Yr %.0Lf", Simul, YEAR1);
	
	if (DYNAMIC==1) 
	{
		if (!gs_tc) dt_long=60/dt;  // 1 min 
		else dt_long=0.1/dt;  	   // 0.1s si gs dynamique     		
	}
	else   dt_long=600/dt;           // 10 min si statique
	
	Get_DATA(dt_long); //initialise some variable

	Climat();
	if (LA_Var) Phenology(dt);
	if (Leaf_Fall==1 || Leaf_Fall==2) for (i=1;i<4;i++) if(Leaf_Area[i]) Compute_Leaf_Fall(i);
	if (Leaf_Fall==3 && Root_Area_fi) Compute_Root_Fall();
	for (i=1;i<4;i++) TLeaf(i); 
	if (Type_Axil) TAxil();
	Fluidity();
	Compute_ST();
	Compute_T_Osmotic();
	if (T_g_cuti) for (i=1;i<4;i++) Compute_g_cuti(i);
	A_net[0]=0;
	for (i=1;i<4;i++)
	{
		if (Leaf_Area[i])
		{
			if (PhotoS_model==3) A_net[i]=Net_Photosynthesis_C3(i);
			if (PhotoS_model>=4) A_net[i]=Net_Photosynthesis_C4(i);
		}
		else        A_net[i]=0;
	A_net[0]+=	A_net[i]*Branch_distri[i]; //mean photosynthesis rate in umol/s/m2
	}
	Rm=Respiration(i);
	Compute_Cavitation(dt);
	Compute_THF();
	Respiration(i);
	soil(dt);
	init_Cavitation();
	Compute_K();
	
	
	if (PRINT_SCREEN)
	{
		printf("   \n   Year ");
		if (PRINT_SCREEN==2)  printf("   Days  PLC_Leaf");
		else for (i=0;i<NPAR;i++)  if (Screen_out[i])  printf("%s", Label[i]);
		printf("\n");
	}

	while(!DEAD && !END_CLIMAT)
	{		
			indice++;
			indice_double+=1;
			//T=T0+indice_double*dt;
			
			if (DYNAMIC==1)     // Dynamic model
			{
				Compute_dq(dt);
				Compute_Q();
				Compute_P_dynamic();
				//if (GRAVITY==3)   Gravity_test();// simulates a suden change in Trunk xylem potential
				if (DYNAMIC0>=2)         // mix model
				{
					if (!PLC_LIMIT && !T_LIMIT)
					{
						DYNAMIC=0;
						dt=dt_stat;
						dt_long=600/dt;
						indice=indice*dt_dyna/dt_stat+1;
						indice_double=indice_double*dt_dyna/dt_stat+1;
					}
					//if (REW_t>=REW_crit || !Leaf_Area[0]) PLC_LIMIT=0;
				}
			}
			
			else if (DYNAMIC==0)  // quasi-STEADY
			{
				if (LA_Var) Phenology(dt);
				if (Leaf_Fall==1 || Leaf_Fall==2) for (i=1;i<4;i++) if (Leaf_Area[i]) Compute_Leaf_Fall(i);
				if (Leaf_Fall==3 && Root_Area_fi) Compute_Root_Fall();
				Compute_Cavitation(dt);
				Compute_THF();
				Compute_K();
				for (i=1;i<4;i++) if (Branch_distri[i]) E_day(dt,dt_long,i);
				Compute_P_steady(dt);
				Compute_Q_steady();
				Climat();
				soil(dt);
				
				if (DYNAMIC0>=2)         // mix model
				{
					if (PLC_LIMIT || T_LIMIT)
					{
						DYNAMIC=1;
						dt=dt_dyna;
						dt_long=60/dt;
						indice=indice*dt_stat/dt_dyna+1;
						indice_double=indice_double*dt_stat/dt_dyna+1;
						for (i=1;i<4;i++) if (Branch_distri[i]) E_day(dt, dt_long,i);
					}
				}			
			}
			
	
			//T=floorl(T*1000)/1000; // round up Time to 1/100 of sec
			
			//  compute these data only once per mn for Dynamic and once per 10min for Steady
			 
			if(!(indice%(unsigned long)(dt_long)))
			{
				//printf ("loop \n");
				Fluidity();
				Compute_T_Osmotic();
				
				for (i=1;i<4;i++) if (Branch_distri[i]) TLeaf(i);
				if (Type_Axil) TAxil();
			//	if (Leaf_Fall==1 || Leaf_Fall==2) for (i=1;i<4;i++) if(Leaf_Area[i]) Compute_Leaf_Fall(i);
			//	if (Leaf_Fall==3 && Root_Area_fi) Compute_Root_Fall();

				if (ACCLIMATE) acclimate();
				if (T_g_cuti) for (i=1;i<4;i++) Compute_g_cuti(i);
				
				if (COMPET) if (!(indice%(unsigned long)(600/dt))) Soil_compet(); //soil synchronisation every 6 min
 	
				if (DYNAMIC==1)     //long loop for Dynamic model
				{
					if (LA_Var) Phenology(dt_long*dt);
					if (Leaf_Fall==1 || Leaf_Fall==2) for (i=1;i<4;i++) if (Leaf_Area[i]) Compute_Leaf_Fall(i);
					if (Leaf_Fall==3 && Root_Area_fi) Compute_Root_Fall();
					soil(dt_long*dt);
					Compute_Cavitation(dt_long*dt);
					Compute_THF();
					Compute_K();
					for (i=1;i<4;i++) if (Branch_distri[i]) E_day(dt, dt_long,i);
					Climat();   // compute climatic variables
				}
	
		
				/*if (Leaf_Area[0]) Compute_Growth();        // compute Trunk and Fruit growth only when Leafy
					else
				{
					Growth_rate_Fruit=0;
					Growth_rate_Trunk=0;
				}*/
				if (Growth_control) Compute_Growth();        // compute Trunk and Fruit growth only when Leafy
				if (T_Leaf[1] > T_Leaf_max) T_Leaf_max=T_Leaf[1];
				if (P_Leaf_Symp[1] < P_min_Leaf) P_min_Leaf=P_Leaf_Symp[1];
				if (P_Leaf_Symp[1] < P_min_lf_d)
				{
					P_min_lf_d=P_Leaf_Symp[1];
					gs_max_d2= g_cuti[1]+g_s[1];
				}
				if (P_Leaf_Symp[1]  > P_max_lf_d) P_max_lf_d=P_Leaf_Symp[1];
				if (P_Branch_Apo[1] < P_min_stem) P_min_stem=P_Branch_Apo[1];
				
				F_Leaf[0]=F_Branch[0]=0;
				for (i=1;i<4;i++) 
					{
						F_Leaf[0]+=E_Leaf[i]*Leaf_Area[i];
						F_Branch[0]+=E_Branch[i]*Branch_Area[i];
					}
				
				if (CUT==1)      E_tot+=dt_long*dt*(F_Leaf[0]+F_Branch[0]+E_Axil*Axil_Area + E_Petiole*Petiole_area);  // in mmol
				else if (CUT==2) E_tot+=dt_long*dt*(F_Leaf[0]+F_Branch[0]+E_Trunk*Trunk_Area+E_Axil*Axil_Area + E_Petiole*Petiole_area);  // in mmol
				else
				{
					E_tot+=         dt_long*dt*(E_Axil*Axil_Area + E_Petiole*Petiole_area+F_Leaf[0]+F_Branch[0]+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3);  // in mmol
					E_tot_day+=     dt_long*dt*(E_Axil*Axil_Area + E_Petiole*Petiole_area+F_Leaf[0]+F_Branch[0]+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3);  // in mmol
					EvapoT_tot+=    dt_long*dt*(E_Axil*Axil_Area + E_Petiole*Petiole_area+F_Leaf[0]+F_Branch[0]+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3 + E_Soil*Surface_Soil);
					EvapoT_day+=    dt_long*dt*(E_Axil*Axil_Area + E_Petiole*Petiole_area+F_Leaf[0]+F_Branch[0]+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3 + E_Soil*Surface_Soil);
				}
				
				// Leaf Photosynthesis rates
				A_net[0]=0;	
				for (i=1;i<4;i++)
					{
						if (Leaf_Area[i])
						{
							if (PhotoS_model==3) A_net[i]=Net_Photosynthesis_C3(i);
							if (PhotoS_model>=4) A_net[i]=Net_Photosynthesis_C4(i);
						}
						else        A_net[i]=0;
					A_net[0]+=	A_net[i]*Branch_distri[i]; //mean photosynthesis rate in umol/s/m2
					}
				Rm=Respiration(1);
				if (A_net[0]>A_net_max) A_net_max=A_net[0];
			
				A_net_tot+=     A_net[0]*dt_long*dt*Leaf_Area[0]/1000000 ;        //total tree  annual A_net in mol of C02
				A_net_day+=     A_net[0]*dt_long*dt*Leaf_Area[0] /1000000;        //total tree daily  A_net in mol of C02
				Resp_tot+=      Rm*dt_long*dt*Leaf_Area[0] /1000000;              //total tree annual respiration in mol of C02
				
				A_gross=        A_net[0] - Rm;                			 // gross A
				A_gross_tot+=   A_gross*dt_long*dt*Leaf_Area[0] /1000000;    		 // total tree annual A_gross in mol of C02
				
				//Canopy Photosynthesis data
				if (Leaf_Area[0]) 
					{
					if (Extinction_Coeff) A_net_c=A_net[1]/Extinction_Coeff*(1-exp(-Extinction_Coeff*LAI_Crown)); //Canopy photosynthesis per ground area
					else A_net_c=A_net[1];
					}
				else A_net_c=0;
			
				A_net_tot_c +=A_net_c*Crown_Area*dt_long*dt/1000000;  //total plant canopy annual A_net in mol/year
				A_net_day_c +=A_net_c*Crown_Area*dt_long*dt/1000000;  //total plant canopy daily  A_net in mol/day
				Resp_tot_c  +=Rm*Leaf_Area[0]*dt_long*dt/1000000;          //total plant canopy annual  Respiration in mol/year
				A_gross_tot_c=A_net_tot_c-Resp_tot_c;
				
				Pot_PAR_tot+=POTENTIAL_PAR*dt_long*dt/24/3600;
				PAR_tot+=PAR[0]*dt_long*dt/24/3600;
				if (POTENTIAL_PAR>Max_PAR_tot) Max_PAR_tot=POTENTIAL_PAR;
				VPD_Air_tot+=VPD_Air*dt_long*dt/3600/24;
				VPD_Leaf_tot+=  VPD_Leaf[1]*dt_long*dt/24/3600;

				C_budget(dt_long);
				ETP_Penman_tot+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP totale in mm or Kg
				ETP_Penman_day+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP per day in mm or Kg
				
				if (Leaf_Area[0] && (LA_Var && (DOY>=LA_para1 && DOY<LA_para4)))  ETP_Leaf_tot+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP totale in mm or Kg
				if (Type_Axil==2)
				{
					Q_Axil_Symp0+= Growth_rate_Fruit*dt_long*dt;                            //plastic growth for a Fruit
				  //  Growth_Fruit+= (Q_Axil_Symp-Q2);                                        //elastic growth
				}
				Growth_Trunk+= (Q_Trunk_Symp-Q1+ Growth_rate_Trunk*Q_Trunk_Symp0*dt_long*dt);  // cumulative whole Trunk elastic growth in mmol of water
				Q1=Q_Trunk_Symp;
				Q2=Q_Axil_Symp;
				
				
				//refill soil at XXX%
				if (REHYDRATE) if (!COMPET || (COMPET && COMPET_number==0)) Rehydrate(dt_long);
				if (CLIMAT!=0 && CLIMAT!=5)  next_climat();
				if (REW_t<0.2 && T_REW_Soil>(T-T0))      T_REW_Soil= T-T0;
				if (RWC_Axil<0.02 && T_RWC_Axil>(T-T0))  T_RWC_Axil= T-T0;
				for (i=1;i<4;i++) 
				{
					RWC_Branch_s[i]= Q_Branch_Symp[i] / Q_Branch_Symp0[i];
					if (RWC_Branch_s[i]<RWC_Branch_s_min[i]) RWC_Branch_s_min[i]=RWC_Branch_s[i];
				}

				
				
			if (!(indice%(unsigned long)(3600*24/dt))) if (Leaf_Area[0] && LA_Var && T/3600/1000/24>=LA_para1 && T/3600/1000/24<LA_para4) //mean PLeaf_min when Leafy
			{
				P_min_Leaf_mean_Leafy+=P_min_lf_d;
				N_days2++;
			}
			
					
			
		/*	if (!(indice%(unsigned long)(3600*24/dt))) */	if (IRRIGATE == 1 || IRRIGATE == 2 || IRRIGATE == 3 || IRRIGATE==8 ||IRRIGATE==9 || IRRIGATE==10 ) /*if (COMPET && COMPET_number==0)*/ Irrigate(); // irrigate at 24h00
			if (!(indice%(unsigned long)(3600*24/dt)))		if (IRRIGATE == 6 || IRRIGATE == 7 || IRRIGATE==11) /*if (COMPET && COMPET_number==0) */ Irrigate(); // irrigate at 24h00
			if (!(indice%(unsigned long)(3600*24/dt)))    indice=0;  // a new day
			
			if (END_DEATH) test_END();
			if (DEAD) if (COMPET)
			{
				_soil = fopen(filename_SOIL,"r+");
				rewind (_soil);
				fprintf(_soil,"%Lf %Lf %Lf %Lf",T*2/1000, Q_Soil1,Q_Soil2,Q_Soil3);
				fclose(_soil);
			}
			
			if (days_simul && T>(T0+(days_simul*3600*1000*24))) //if days=0 never stops, otherwise compute for x days
				{
					DEAD=1; //printf("dead ici");
					print_final(); 
					print_annual(); 
					//purge(); 
					END_CLIMAT=1;
					//Reset();
				} 
			else if (END_DEATH && DEAD && (CLIMAT==2 || CLIMAT==1 || CLIMAT==11)) // if DEAD skip the end of the year if END_DEATH not zero
				{
					print_final(); 
					print_annual();  
					purge(); 					
				}
			else if (DEAD || END_CLIMAT)
				{
					print_final(); 
					print_annual();  
				}
			
			
		fill_soil_layers();
	
	   }//of long loop
	
	  T=T+dt*1000;	// time in milliseconds
	  Test_Time();
	 // if (!(indice%(unsigned long)(3600*24/dt))) //if (CLIMAT==0 || CLIMAT==10 || CLIMAT==5)    // a new day reset daily values
		if(!((unsigned long long)(T)%(unsigned long long)(3600*t_out*1000)))		
			{
				ETP_Penman_day=0;
				A_net_day=0;
				A_net_day_c=0;
				E_tot_day=0;
				EvapoT_day=0;
				P_min_lf_d=0;
				P_max_lf_d=-1000;
				gs_min_d=10000;
				gs_max_d=0;
				gs_max_d2=0;
				SF_min_d=10000;
				SF_max_d=0;
				if (CLIMAT==10) DOY=DOY+1;
				Pot_PAR_tot=0;
			}	
		
		if (isnan(E_Leaf[0]) || isnan(g_s[0]) || isnan(P_Leaf_Apo[0]) || isnan(P_Trunk_Apo) || isnan(P_Branch_Apo[0]) || isnan(P_Root_Apo1) ||isnan(P_Root_Apo2) || isnan(P_Root_Apo3))
			{print_transient(); DEAD=1;}
		
		threshold=10;
		if(P_Leaf_Apo[1]>threshold || P_Trunk_Apo>threshold || P_Branch_Apo[1]>threshold || P_Root_Apo1>threshold || P_Root_Apo2>threshold || P_Root_Apo3>threshold)
			{
			printf(" Overflow ");
			printf("%Lf %Lf %Lf",P_Leaf_Apo[1],P_Trunk_Apo,P_Branch_Apo[1]);
			DEAD=1;
			PLC_Leaf_Apo[1]=9999.0;
			if (RANDOMISE)setup(); 
			else exit(1);
			}			
	} //of while not dead loop
	//if (END_CLIMAT) printf("end climat");
	if (TRANSIENT)
	{ 	
		//print_transient();
		transient = fopen(filename_TRANS,"a");
		fprintf(transient,"\n");
		fclose(transient);
	}
	
	printf(" THF_l: %.1Lf PLC: %.1Lf",T_PLC_Leaf[1]/3600/24/1000,PLC_Leaf_Apo[0]);

	if (END_CLIMAT) fclose (climat_in);
	for (i=1;i<4;i++)
	{
	PLC_Leaf_Apo[i]=  0;
	PLC_Branch_Apo[i]=0;
	}
	PLC_Axil_Apo=  0;
	PLC_Trunk_Apo= 0;
	PLC_Root_Apo1= 0;
	PLC_Root_Apo2= 0;
	PLC_Root_Apo3= 0;
	
	year=1;
}

void convert(void)  //compute different parameters and convert Kg to mmol for capacitance Q
{
	int i;
	//	Q_Leaf_Apo0   = Succulence/1000 * Leaf_Apo_fraction * LA_max_Pheno;  // in Kg

	for (i=1;i<4;i++)
	{
	Q_Leaf_Apo0[i]   = Succulence/1000 * Leaf_Apo_fraction * LA_max_init*Branch_distri[i];
	Q_Leaf_Symp0[i]  = Succulence/1000 * (1- Leaf_Apo_fraction) * LA_max_Pheno*Branch_distri[i]; // in Kg
	}
	// compute the total Q for all organs
	if (FRACTAL==2 || FRACTAL==3) //then the total number of Axillary organs is N_Axil
	{
		if (Type_Axil==3) //a flower, surface of a disc
		{
			Q_Axil_Symp0    = N_Axil * WC_Axil * 3.1416 /4 * Diam_Axil * Diam_Axil;  // organ surface times surfacic water content  in Kg, WC in kg/m2
			Q_Petiole_Symp0 = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Symp_fraction * 1000; //total petiole Symp volume in Kg
			Q_Axil_Apo0     = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Apo_fraction * 1000;
			Petiole_area    = N_Axil *Length_Petiole * Diam_Petiole * 3.1416;
		}
		else if (Type_Axil==2) //  a Fruit, volume of a sphere
		{
			Q_Axil_Symp0    = N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
			Q_Petiole_Symp0 = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Symp_fraction * 1000; //total Branch Symp volume in Kg
			Q_Axil_Apo0     = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Apo_fraction * 1000;
			Petiole_area    = N_Axil *Length_Petiole * Diam_Petiole * 3.1416;
		}
		
		else if (Type_Axil==1)  // a bud  volume of a sphere
		{
			Q_Axil_Symp0   = N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
			Q_Axil_Apo0    = N_Axil *Q_Axil_Symp0/20;   // assume 5%  bud/Fruit water content is in the apoplasm.
			Q_Petiole_Symp0 =0;
		}
		else
		{
			Q_Axil_Symp0 =0;
			Q_Axil_Apo0 =0;
			Q_Petiole_Symp0=0;
			Petiole_area =0;
		}
	}
	
	else if (FRACTAL==0 || FRACTAL==1 ) //then the total number of Axillary organ is N_Axil * Number_Branch
		{
			if (Type_Axil==3) //a flower, surface of a disc
			{
				Q_Axil_Symp0    = N_Axil *Number_Branch * WC_Axil * 3.1416 /4 * Diam_Axil * Diam_Axil;  // organ surface times surfacic water content  in Kg, WC in kg/m2
				Q_Petiole_Symp0 = N_Axil *Number_Branch * Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Symp_fraction * 1000; //total petiole Symp volume in Kg
				Q_Axil_Apo0     = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Apo_fraction * 1000;
				Petiole_area    = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * 3.1416;
			}
			else if (Type_Axil==2) //  a Fruit, volume of a sphere
			{
				Q_Axil_Symp0    = N_Axil *Number_Branch *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
				Q_Petiole_Symp0 = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Symp_fraction * 1000; //total Branch Symp volume in Kg
				Q_Axil_Apo0     = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_Apo_fraction * 1000;
				Petiole_area    = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * 3.1416;
			}
			
			else if (Type_Axil==1)  // a bud  volume of a sphere
			{
				Q_Axil_Symp0   = N_Axil *Number_Branch *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
				Q_Axil_Apo0    = N_Axil *Number_Branch *Q_Axil_Symp0/20;   // assume 5%  bud/Fruit water content is in the apoplasm.
				Q_Petiole_Symp0 =0;
			}
			else
			{
				Q_Axil_Symp0 =0;
				Q_Axil_Apo0 =0;
				Q_Petiole_Symp0=0;
				Petiole_area =0;
			}
		}
	
	Q_Root_Symp0   = Q_Root_Symp_FR;// * (Root_upper0 + Root_middle0 + Root_lower0); // always use fractal Roots
	Q_Root_Apo_t0  = Q_Root_Apo_FR *(Root_upper0 + Root_middle0 + Root_lower0);
	
	
	if (FRACTAL==1) //based on morphology
		{
			for (i=1;i<4;i++) Q_Branch_Symp0[i]= Branch_distri[i]* Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 * Branch_Symp_fraction * 1000; //total Branch Symp volume in Kg
			for (i=1;i<4;i++) Q_Branch_Apo0[i] = Branch_distri[i]* Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 * Branch_Apo_fraction * 1000; //total Branch apo volume in Kg
			Q_Trunk_Symp0 = Length_Trunk * Diam_Trunk * Diam_Trunk * 3.1416 / 4 * Trunk_Symp_fraction * 1000*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction; //in Kg of H2O in the sapwood
			Q_Trunk_Apo0  = Length_Trunk * Diam_Trunk * Diam_Trunk * 3.1416 / 4 * Trunk_Apo_fraction * 1000*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction;  //in Kg of H2O in the sapwood
		}
	else // Fractal tree
		{
			for (i=1;i<4;i++) Q_Branch_Symp0[i] = Branch_distri[i]* Q_Branch_Symp_FR;
			if (K_VAR==12) 
				{
				Q_Branch_Apo0[1]  = Branch_distri[1]* Q_Branch_Apo_FR*(100-K_VAR_P1)/100;
				Q_Branch_Apo0[2]  = Branch_distri[2]* Q_Branch_Apo_FR*(100-K_VAR_P2)/100;
				Q_Branch_Apo0[3]  = Branch_distri[3]* Q_Branch_Apo_FR*(100-K_VAR_P3)/100;
				}
			else for (i=1;i<4;i++) Q_Branch_Apo0[i]  = Branch_distri[i]* Q_Branch_Apo_FR;
				
			Q_Trunk_Symp0  = Q_Trunk_Sym_FR;
			Q_Trunk_Apo0   = Q_Trunk_Apo_FR;
		}
	Q_Leaf_Symp0[0]=Q_Branch_Symp0[0]=Q_Leaf_Apo0[0]=Q_Branch_Apo0[0]=0;
	for (i=1;i<4;i++) 
	{
		Q_Leaf_Symp0[0]+=Q_Leaf_Symp0[i];
		Q_Branch_Symp0[0]+=Q_Branch_Symp0[i];
		Q_Leaf_Apo0[0]+=Q_Leaf_Apo0[i];
		Q_Branch_Apo0[0]+=Q_Branch_Apo0[i];
	}
	
	Q_Plant_s0= (Q_Axil_Symp0 + Q_Petiole_Symp0 + Q_Leaf_Symp0[0]  + Q_Branch_Symp0[0] + Q_Trunk_Symp0 + Q_Root_Symp_FR*(Root_upper0 + Root_lower0 + Root_middle0 ));
	Q_Plant_a0= (Q_Leaf_Apo0[0] + Q_Axil_Apo0 + Q_Branch_Apo0[0] + Q_Trunk_Apo0 + Q_Root_Apo_FR*(Root_upper0 + Root_lower0 + Root_middle0 ));	
	Q_Plant0= Q_Plant_a0 + Q_Plant_s0;
	// printf(" \n%Le %Le  %Le %Le %Le %Le %Le %Le %Le\n", Q_Axil_Symp0, Q_Petiole_Symp0 , Q_Leaf_Symp0  ,Q_Branch_Symp0 , Q_Trunk_Symp0 , Q_Root_Symp_FR,Root_upper0, Root_lower0 , Root_middle0 );
	//convert from volume specific capacitance in kg MPa-1 L-1 to capacitance  in kg MPa-1 with water volume
	//if (LA_Var) C_Leaf_Apo   *=  (Succulence/1000 * Leaf_Apo_fraction * LA_max_init); else C_Leaf_Apo   *=  Q_Leaf_Apo0;
	
	for (i=1;i<4;i++) C_Leaf_Apo[i] = C_Leaf_Apo[0]  * Q_Leaf_Apo0[i];
	for (i=1;i<4;i++) C_Branch_Apo[i] = C_Branch_Apo[0]*  Q_Branch_Apo0[i];
	C_Trunk_Apo  *=  Q_Trunk_Apo0;
	C_Root_Apo   *=  Q_Root_Apo_t0;
	C_Axil_Apo   *=  Q_Axil_Apo0;
	
	//Diam_Root    = Diam_Root_FR;
	Length_Root = Length_Root_FR;
	Root_Area   = Root_Area_FR;
	
	if (FRACTAL==1)
	{
		for (i=1;i<4;i++) Branch_Area[i]   = Branch_distri[i]* Number_Branch *  Length_Branch * Diam_Branch * 3.1416;
		Trunk_Area    = Length_Trunk * Diam_Trunk * 3.1416;
	}
	else
	{
		for (i=1;i<4;i++) Branch_Area[i]   = Branch_distri[i]*Branch_Area_FR;
		Trunk_Area    = Trunk_Area_FR;
	}
	Branch_Area[0]=Branch_Area[1]+Branch_Area[2]+Branch_Area[3];
	
	if (Type_Axil==3) Axil_Area      = 3.1416*Diam_Axil*Diam_Axil/4;  // a flower, surface of a disc
	else              Axil_Area      = 3.1416*Diam_Axil*Diam_Axil;  // a bud or a Fruit, surface of a sphere
	if (FRACTAL==1)   Axil_Area*= (N_Axil*Number_Branch);
	else              Axil_Area*= N_Axil;
	
	for (i=1;i<4;i++) K_Leaf_Symp_0[i]  = K_Leaf_Symp_0[0]*Branch_distri[i]*LA_max_init;
	K_Root_Apo0    = K_Root_Apo_FR;
	K_Trunk_Symp0   *= Trunk_Area;
	K_Root_Symp1    = K_Root_Symp0 * Root_Area_fi;  // only the fine Roots absorb water
	K_Root_Symp2    = K_Root_Symp0 * Root_Area_FR;  // all the Root surface has a Symplasmic compartment
	//K_Leaf_Apo0    *= LA_max_Pheno;  // from Kl to k
	for (i=1;i<4;i++) K_Leaf_Apo0[i]    = K_Leaf_Apo0[0]*Branch_distri[i]*LA_max_init;  // from Kl to k
	
	if (FRACTAL==1)
	{
		for (i=1;i<4;i++) K_Branch_Symp0[i]  = K_Branch_Symp0[0]*Number_Branch * Branch_Area[i];
		for (i=1;i<4;i++) K_Branch_Apo0[i]   = K_Branch_Apo0[0]*Branch_distri[i]*Number_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 /Length_Branch; // from Ks to k
		K_Trunk_Apo0   *= Diam_Trunk * Diam_Trunk * 3.1416 / 4 /Length_Trunk*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction; // from Ks to k
		K_Axil_Apo0     *=(N_Axil*Number_Branch);
		K_Axil_Symp     *=(N_Axil*Number_Branch);
		K_Axil_Symp2    *=(N_Axil*Number_Branch);
	}
	else
	{
		for (i=1;i<4;i++) K_Branch_Symp0[i]  = K_Branch_Symp0[0]*Branch_Area[i];
		for (i=1;i<4;i++) K_Branch_Apo0[i]  =  Branch_distri[i]*K_Branch_Apo_FR; // from Ks to k
		K_Trunk_Apo0   =  K_Trunk_Apo_FR; // from Ks to k
		K_Axil_Apo0    *=N_Axil;
		K_Axil_Symp    *=N_Axil;
		K_Axil_Symp2   *=N_Axil;
	}
	
	C_Leaf_Apo[0]		*=   1000*1000/18; //convert Kg to mmol for capacitance & Q
	C_Branch_Apo[0]	*=   1000*1000/18;
	C_Trunk_Apo   	*=   1000*1000/18;
	C_Root_Apo    	*=   1000*1000/18;
	C_Axil_Apo    	*=   1000*1000/18;
	
	for (i=1;i<4;i++) 
	{
	Q_Leaf_Symp0[i]	*=   1000*1000/18;
	Q_Branch_Symp0[i]	*=   1000*1000/18;
	Q_Leaf_Apo0[i]   	*=   1000*1000/18;
	Q_Branch_Apo0[i] 	*=   1000*1000/18;
	C_Leaf_Apo[i]		*=   1000*1000/18; //convert Kg to mmol for capacitance & Q
	C_Branch_Apo[i]	*=   1000*1000/18;
	
	}
	Q_Trunk_Symp0 	*=   1000*1000/18;
	Q_Root_Symp0  	*=   1000*1000/18;
	Q_Axil_Symp0  	*=   1000*1000/18;
	Q_Petiole_Symp0	*=   1000*1000/18;	
	Q_Trunk_Apo0  	*=   1000*1000/18;
	Q_Root_Apo_t0	*=   1000*1000/18;
	Q_Axil_Apo0   	*=   1000*1000/18;
	
	Growth_Trunk=Q_Trunk_Symp0;

	Growth_Fruit=Q_Axil_Symp0;
		
}





void setparameters(void)
{
	Simul=para[0];DYNAMIC0=para[1]; dPLC_crit=para[2]; REW_crit=para[3]; debug=para[4]; TRANSIENT=para[5]; PRINT_SCREEN=para[6]; PRINT_GRAPH=para[7]; gs_cst=para[8]; K_VAR=para[9]; 
	CUT=para[10]; END_DEATH=para[11]; PLC_END=para[12]; days_simul=para[13]; dt_dyna=para[14]; dt_stat=para[15]; t_out=para[16]; CLIMAT=para[17]; TLEAF=para[18]; FLUID=para[19];
	SURFACE_TENSION=para[20]; T_OSMOTIC=para[21]; T_g_cuti=para[22]; GRAVITY=para[23]; CAPILLARITY=para[24]; Regul_gs=para[25]; Regul_gs_para1=para[26]; Regul_gs_para2=para[27]; 
	Regul_gs_para3=para[28]; Regul_gs_para4=para[29]; DOY_0=para[30]; Lat=para[31]; T_air_min=para[32]; T_air_max=para[33]; RH_air_min=para[34]; RH_air_max=para[35]; PAR_max=para[36]; Wind[0]=para[37]; 
	HW=para[38]; HW_day=para[39]; HW_duration=para[40]; HW_T=para[41]; Teta_s_1=para[42]; Teta_r_1=para[43]; alpha_1=para[44]; n_1=para[45]; K_sat_1=para[46]; L=para[47]; Rock_f1=para[48]; 
	PENMAN=para[49]; Penman_Coeff=para[50]; g_crown0=para[51]; INTERCEPTION=para[52]; Interception_min=para[53]; Threshold_rain=para[54]; Interception_factor=para[55]; IRRIGATE=para[56]; 
	RWC_Irr=para[57]; Daily_Irr=para[58]; IRR_DOY_S=para[59]; IRR_DOY_F=para[60]; REHYDRATE=para[61]; PLC_REHYD=para[62]; CONTINUOUS=para[63]; Extinction_Coeff=para[64]; Leaf_size=para[65]; 
	Leaf_angle[1]=para[66]; LA_Var=para[67]; LA_max_init=para[68]; LA_min=para[69]; LA_para1=para[70]; LA_para2=para[71]; LA_para3=para[72]; LA_para4=para[73]; Leaf_Fall=para[74]; 
	P50_Leaf_Fall=para[75]; Slope_Leaf_Fall=para[76]; Succulence=para[77]; Leaf_Apo_fraction=para[78]; LMA=para[79]; VcMax=para[80]; VjMax=para[81]; CO2_atm=para[82]; Rd25=para[83]; 
	Qye=para[84]; Kc25=para[85]; Ko25=para[86]; a_Res=para[87]; Number_Branch=para[88]; Length_Branch=para[89]; Diam_Branch=para[90]; Branch_Apo_fraction=para[91]; 
	Branch_Symp_fraction=para[92]; Density=para[93]; Length_Trunk=para[94]; Diam_Trunk=para[95]; Trunk_Apo_fraction=para[96]; Trunk_Symp_fraction=para[97]; Trunk_sapwood_fraction=para[98];
	Length_Root_fi=para[99]; Diam_Root=para[100]; Root_Apo_fraction=para[101]; Root_Symp_fraction=para[102]; FRACTAL=para[103]; Q_Branch_Symp_FR=para[104]; Q_Branch_Apo_FR=para[105]; 
	Q_Trunk_Sym_FR=para[106]; Q_Trunk_Apo_FR=para[107]; Q_Root_Symp_FR=para[108]; Q_Root_Apo_FR=para[109]; Branch_Area_FR=para[110]; Trunk_Area_FR=para[111]; Root_Area_FR=para[112]; 
	Root_Area_fi_0=para[113]; Length_Root_FR=para[114]; Diam_Root_FR=para[115]; K_Branch_Apo_FR=para[116]; K_Trunk_Apo_FR=para[117]; K_Root_Apo_FR=para[118]; gs_max[1]=para[119]; 
	gs_night=para[120]; Jarvis_PAR=para[121]; GS_MAX=para[122]; Tgs_optim=para[123]; Tgs_sens=para[124]; g_cuti_20=para[125]; TP=para[126]; Q10_1=para[127]; Q10_2=para[128]; 
	g_Branch=para[129]; g_Trunk=para[130]; g_Root[0]=para[131]; g_Soil0=para[132]; Extensibility_Trunk=para[133]; Yield_Trunk=para[134]; C_Leaf_Apo[0]=para[135]; C_Branch_Apo[0]=para[136]; 
	C_Trunk_Apo=para[137]; C_Root_Apo=para[138]; Soil_Depth=para[139]; Soil_Width=para[140]; Root_upper0=para[141]; Root_middle0=para[142]; Root_lower0=para[143]; gap=para[144]; 
	Epsilon_Leaf_Symp=para[145]; Epsilon_Branch_Symp=para[146]; Epsilon_Trunk_Symp=para[147]; Epsilon_Root_Symp=para[148]; Pi0_Leaf_Symp_0=para[149]; Pi0_Branch_Symp_0=para[150]; 
	Pi0_Trunk_Symp_0=para[151]; Pi0_Root_Symp=para[152]; K_Leaf_Apo0[0]=para[153]; K_Branch_Apo0[0]=para[154]; K_Trunk_Apo0=para[155]; K_Root_Apo0=para[156]; K_Leaf_Symp_0[0]=para[157]; 
	K_Branch_Symp0[0]=para[158]; K_Trunk_Symp0=para[159]; K_Root_Symp0=para[160]; P50_Leaf_Apo_0[1]=para[161]; P50_Branch_Apo_0[1]=para[162]; P50_Trunk_Apo_0=para[163]; P50_Root_Apo_0=para[164]; 
	Slope_Leaf_Apo[1]=para[165]; Slope_Branch_Apo[1]=para[166]; Slope_Trunk_Apo=para[167]; Slope_Root_Apo=para[168]; Type_Axil=para[169]; N_Axil=para[170]; K_Axil_Apo0=para[171]; 
	K_Axil_Symp=para[172]; K_Axil_Symp2=para[173]; Epsilon_Axil_Symp=para[174]; Pi0_Axil_Symp=para[175]; P50_Axil_Apo=para[176]; Slope_Axil_Apo=para[177]; g_Axil_min20=para[178]; 
	g_Axil_max=para[179]; g_Petiole=para[180]; C_Axil_Apo=para[181]; Diam_Axil=para[182]; Length_Petiole=para[183]; Diam_Petiole=para[184]; WC_Axil=para[185]; Extensibility_Fruit=para[186]; 
	Yield_Fruit=para[187]; REFILL=para[188]; P_REFILL=para[189]; leg_Leaf=para[190]; leg_Branch=para[191]; leg_Trunk=para[192]; leg_Root=para[193]; T_Soil_Crit=para[194]; 
	COMPET=para[195]; Rock_f2=para[196]; Rock_f3=para[197]; gs_CO2_sens=para[198]; HH1=para[199]; HH2=para[200]; K_VAR_P1=para[201]; K_VAR_P2=para[202]; K_VAR_P3=para[203]; gs_tc=para[204];
	RANDOMISE=para[205];LA_term_shoot=para[206];Tapering=para[207];Root_shoot_ratio=para[208];Root_ramif=para[209];Length_fine_Root=para[210];Length_term_shoot=para[211];
	ACCLIMATE=para[212];Acc_P1=para[213];Acc_P2=para[214];Acc_P3=para[215];Layer_1=para[216];Layer_2=para[217];Layer_3=para[218];Teta_s_2=para[219]; Teta_r_2=para[220]; alpha_2=para[221]; 
	n_2=para[222]; K_sat_2=para[223];Teta_s_3=para[224]; Teta_r_3=para[225]; alpha_3=para[226]; n_3=para[227]; K_sat_3=para[228];PhotoS_model=para[229]; 
	Branch_distri[1]=para[230]; Branch_distri[2]=para[231]; Branch_distri[3]=para[232]; P50_Leaf_Apo_0[2]=para[233]; P50_Leaf_Apo_0[3]=para[234];Slope_Leaf_Apo[2]=para[235];Slope_Leaf_Apo[3]=para[236];
	P50_Branch_Apo_0[2]=para[237];P50_Branch_Apo_0[3]=para[238];Slope_Branch_Apo[2]=para[239];Slope_Branch_Apo[3]=para[240];gs_max[2]=para[241]; gs_max[3]=para[242]; Leaf_angle[2]=para[243];Leaf_angle[3]=para[244];
	RWC_sens=para[245];crown_diam=para[246];PAR_att=para[247];bark=para[248];Rm_25=para[249];Rg_25=para[250];NSC=para[251];g_Axil_regul=para[252]; TP_Axil=para[253]; Q10_1_Axil=para[254];	Q10_2_Axil=para[255];
	Growth_para1=para[256];Growth_para2=para[257];Growth_control=para[258];PI0_Soil=para[259];Psoil_FC=para[260];Teta_ini_1=para[261];Teta_ini_2=para[262];Teta_ini_3=para[263];
}

void checkparameters(void)
{
	if (LA_day2<=LA_day1) LA_day2=LA_day1+1;
	if (LA_day3<=LA_day2) LA_day3=LA_day2+1;
	if (LA_day4<=LA_day3) LA_day4=LA_day3+1;	
	para[70]=LA_day1; para[71]=LA_day2; para[72]=LA_day3; para[73]=LA_day4;
	
}

void initialise(void)
{	
	FILE *transient,*transient_out,*File;
	int i,S,F;
		
		if (CLIMAT==7) //when daily climat file is passed as Simul
		{
			sprintf(filename_CLIM, "meteo/%.0Lf.txt",Simul); 
			CLIMAT=2;
		}
		
		
		else if (CLIMAT==8) //when hourly climat file is passed as Simul
		{
			sprintf(filename_CLIM, "%.0Lf.txt", Simul); 
			CLIMAT=1;
		}

		else if (CLIMAT==1 || CLIMAT==6 || CLIMAT==11)
		{
			sprintf(filename_CLIM, "climat_hour_in.txt");
		}  
		
		else if (CLIMAT==2 || CLIMAT==4)
		{
			sprintf(filename_CLIM, "climat_day_in.txt");

		}  

		if (debug==2) if (CLIMAT==2 || CLIMAT==4 || CLIMAT==7)
		{
		
		compute_climatic_stats();
			for (i=1;i<4;i++) 
				{
					PLC_Leaf_Apo[i]=0;
					PLC_Branch_Apo[i]=0;
				}
			PLC_Trunk_Apo=0;
			PLC_Root_Apo1=0;
			PLC_Root_Apo2=0;
			PLC_Root_Apo3=0;
		
			LA_day1=LA_para1;
			LA_day2=LA_para2;
			LA_day3=LA_para3;
			LA_day4=LA_para4;
			T_base1=LA_para1;
			S_GDD=LA_para2;
			//T_base2=LA_para3;
			LGE=LA_para3;
		}
		if (Regul_gs ==2 || Regul_gs ==15  )  
			{	
				
				Px_gs_0= Regul_gs_para2;

			}
		if (REFILL==2) 
		{
			SYMP_CAVIT=1;
			REFILL=0;
		}
		gs_0=0;
		
		Pgs_12=Regul_gs_para1;
		Pgs_88=Regul_gs_para2;
				
		if (DYNAMIC0==0)    //Steady
			{	
				printf("Entered dynamic 0");
				dt=dt_stat;
				DYNAMIC=0;
				if (PRINT_SCREEN) printf("quasi-STEADY\n");
				
			}
			
		else if (DYNAMIC0==1)    //Dynamic
			{	
				dt=dt_dyna;
				DYNAMIC=1;
				if (PRINT_SCREEN) 
				{	
					printf("DYNAMIC == 1\n");
				
				}
				
			}
		else if (DYNAMIC0>=2)    //Mix start with Steady
			{	
				dt=dt_stat;
				DYNAMIC=0;
				if (PRINT_SCREEN)  printf("MIX model\n");
			}
		
		

		// if -1 then use the Potential PAR to compute PAR_max
		if (PAR_max==-1) {
			PAR_max=Potential_PAR((12))*PAR_att;
		}  
		
		if (FRACTAL==1 || FRACTAL==2) Fractal();	
		g_Axil_min=g_Axil_min20;
		LA_max=LA_max_init;
		

		if (LA_Var) LA_max_Pheno=LA_min;
		else LA_max_Pheno=LA_max;
		Root_Area_fi=Root_Area_fi_0;
		
		Teta_fc_1= Teta_r_1+ (Teta_s_1-Teta_r_1)/(pow(1+pow(alpha_1*Psoil_FC*10,n_1),1-1/n_1));  	// soil humidity at field capacity = PIO_Soil cm pressure head 
		Teta_fc_2= Teta_r_2+ (Teta_s_2-Teta_r_2)/(pow(1+pow(alpha_2*Psoil_FC*10,n_2),1-1/n_2)); 
		Teta_fc_3= Teta_r_3+ (Teta_s_3-Teta_r_3)/(pow(1+pow(alpha_3*Psoil_FC*10,n_3),1-1/n_3)); 
		
		Teta_wp_1= Teta_r_1+ (Teta_s_1-Teta_r_1)/(pow(1+pow(alpha_1*15000,n_1),1-1/n_1));  	// soil humidity at field capacity = 15000cm pressure head 1.5MPa
		Teta_wp_2= Teta_r_2+ (Teta_s_2-Teta_r_2)/(pow(1+pow(alpha_2*15000,n_2),1-1/n_2)); 
		Teta_wp_3= Teta_r_3+ (Teta_s_3-Teta_r_3)/(pow(1+pow(alpha_3*15000,n_3),1-1/n_3)); 
	
		RWC_fc_1= (Teta_fc_1-Teta_r_1)/(Teta_s_1-Teta_r_1);                         		// soil RWC at field capacity
		RWC_fc_2= (Teta_fc_2-Teta_r_2)/(Teta_s_2-Teta_r_2); 
		RWC_fc_3= (Teta_fc_3-Teta_r_3)/(Teta_s_3-Teta_r_3);


		convert();
		
		PREM=1;
		init();
		init();  //do not know why it should be done twice???
		PREM=0;

		for (i=1;i<4;i++) TLeaf(i); // to init g_bl
		
		if (Type_Axil) TAxil();
		if (T_g_cuti) for (i=1;i<4;i++) Compute_g_cuti(i);
		if (!CUT) if (Regul_gs==1) Compute_Turgor_Ref();
		
		
		//if (debug==2) if (CLIMAT==2 || CLIMAT==4 || CLIMAT==7)  compute_climatic_stats();
		
		if (PREM1)
			{	

				//if (debug==2) if (CLIMAT==2 || CLIMAT==4 || CLIMAT==7)  compute_climatic_stats();
				if ((transient_out = fopen("transient_out.txt","r"))==NULL) //transient_init file not found; use default values
				{	
					
					Screen_out[0]=1;Screen_out[17]=1;Screen_out[30]=1;Screen_out[33]=1;Screen_out[78]=1;Screen_out[79]=1;Screen_out[104]=1;Screen_out[116]=1;Screen_out[119]=1;
					
					File_out[0]=1;File_out[17]=1;File_out[30]=1;File_out[33]=1;File_out[78]=1;File_out[79]=1;File_out[104]=1;File_out[116]=1;File_out[119]=1;
				}
				else
					while (!feof(transient_out))
					{
						
						fscanf(transient_out,"%d %d %d\n", &i, &S,&F);
						Screen_out[i]=S;
						File_out[i]=F;

						// This line creates and error in linux if put outside the {}
						fclose(transient_out);
					}
				
				if (TRANSIENT==11 || TRANSIENT==21) // use standard daily values
				{
					for (i=0;i<NPAR;i++) 
					{
						Screen_out[i]=0;
						File_out[i]=0;
					}
				Screen_out[0]=2;Screen_out[17]=1;Screen_out[30]=1;Screen_out[33]=1;Screen_out[40]=1;Screen_out[78]=1;Screen_out[79]=1;Screen_out[104]=1;Screen_out[186]=1;Screen_out[119]=1;
				File_out[0]=1;File_out[17]=1;File_out[30]=1;File_out[33]=1;File_out[40]=1;File_out[78]=1;File_out[79]=1;File_out[104]=1;File_out[186]=1;File_out[119]=1;	
				}
				
				if (TRANSIENT==3) // print only the test 
				{
					for (i=0;i<NPAR;i++) 
					{
						Screen_out[i]=0;
						File_out[i]=0;
					}
				Screen_out[0]=2;Screen_out[300]=1;Screen_out[301]=1;Screen_out[302]=1;Screen_out[303]=1;
				File_out[0]=1;File_out[300]=1;File_out[301]=1;File_out[302]=1;File_out[303]=1;
				}
				
				if (TRANSIENT==12 || TRANSIENT==22) // use standard yearly values
				{
					for (i=0;i<NPAR;i++) 
					{
						Screen_out[i]=0;
						File_out[i]=0;
					}
					Screen_out[0]=2;Screen_out[174]=1;Screen_out[147]=1;Screen_out[78]=1;Screen_out[79]=1;Screen_out[142]=1;Screen_out[120]=1;Screen_out[127]=1;Screen_out[129]=1;
					File_out[0]=1;	File_out[174]=1;	File_out[147]=1;	File_out[78]=1;File_out[79]=1;		File_out[142]=1;File_out[120]=1;	File_out[127]=1;File_out[129]=1;	
				}
				
				if ((File = fopen("annual_out.csv","r"))==NULL) //file does not exist then print headers in the file
				{
					if ((File = fopen("annual_out.csv","w"))==NULL) printf("\a\nCan't create file annual_out.csv !!");
					else
					{
						fprintf(File,"File;N;Simul;Year; ");
						for (i=0;i<NPAR;i++)  if (debug || File_out[i])  fprintf(File, "%s; ", Label[i]);
						fprintf(File,"\n");
						fclose(File);
					}
				}
				
				if (TRANSIENT)
				{
					if ((transient = fopen("transient_out.csv","r")) == NULL)  //if file doesn't exists then print headers
					{
						transient = fopen("transient_out.csv","w");
						fprintf(transient,"YEAR; ");
						if (debug)
							for (i=0;i<NPAR;i++)
							{
								fprintf(transient, "%s; ", Label[i]);
								File_out[i]=1;
							}
						else for (i=0;i<NPAR;i++)  if (File_out[i])  fprintf(transient, "%s; ", Label[i]);
						fprintf(transient,"\n");
					}
					else
					{
						transient = fopen("transient_out.csv","a");
						// fprintf(transient,"\n");
					}
					fclose(transient);
				}
				PREM1=0;
			}
		
		if (PRINT_GRAPH==3)
			{
				printf("Printing graphs only\n");
				
				exit(1);
			}
				
}

void Randomise(void)
{
	FILE *random_file;
	int  random1; //para number
	long double random2, random3, random4; //CV, min% max%
	long double proba,r, theta,u,v;

	if ((random_file = fopen("random_para.txt","r"))==NULL) printf("\a\nCan't find file random_para.txt!!");
	
	while (!feof(random_file)) 
	{
		fscanf(random_file,"%d %Le %Le %Le\n",&random1,&random2,&random3,&random4);
		if (RANDOMISE>0) //uniform
		{
			proba=(long double)(rand() % 10000)/10000; 
			para[random1]=para0[random1]*(1+random3/100+proba*(random4-random3)/100);  
		}
		if (RANDOMISE<0) //gaussian , transformation of Box-Muller
		{			
			u=(long double)(rand() % 10000)/10000; 
			v=(long double)(rand() % 10000)/10000; 
			theta= 2*3.1416*u ;
			if (v) r = pow(-2 *log(v), 0.5) ;
			else 
			{
				v =(long double)(rand() % 10000)/10000; 
				r = pow(-2 *log(v), 0.5) ;
			}
			proba= r*cos(theta) ;
			para[random1]=para0[random1]*(1+proba*random2/100);
		}
	} 
	fclose(random_file);	
}

void setup(void)  // load parameters for simulations and launch computation
{
	FILE *out,*transient,*para_file, *init_file;
	
	int parameter[NPAR];
	int i,j,N_para;
	
	if ((out = fopen(filename_OUT,"a"))==NULL) printf("\a\nCan't create file sureau.out!!");
	fclose(out);
	
	if ((para_file = fopen("sureau_para.txt","r"))==NULL) //para file not found; create a default one
	{
		printf("\a\nCan't find file sureau_para.txt!!");
		default_para_file(); 
		para_file = fopen("sureau_para.txt","r");
	}

	i=0;
	while (!feof(para_file)) 
	{
		fscanf(para_file,"%d\n",&parameter[i]);
		i++;
	}
	fclose(para_file);
	N_para=i-1;  //number of parameters
	
	if ((init_file = fopen(filename_IN,"r"))==NULL) //init file not found; create a default one
	{
		printf("\a\nCan't find file sureau_ini.txt!!");
		default_para(); 
		init_file = fopen(filename_IN,"r");
	}

		
	while (!feof(init_file))  //for all the simulations in sureau_ini.txt
		{
		
		for (j=0;j<=N_para;j++) 
		{
			if (!feof(init_file)) 
				{
				fscanf(init_file,"%Le",&para0[parameter[j]]);
				para[parameter[j]]=para0[parameter[j]];
				}			
			else exit(1);
		}
		setparameters();
	//	Reset();
		if (END_DEATH == 16) 
		{	
			RWC_END=PLC_END;
			PLC_END=99;
		}

		if (RANDOMISE) for(i=0;i<fabsl(RANDOMISE); i++) 
			{
				Randomise();
				
				setparameters();
				checkparameters();
				initialise();
				if (debug !=2) compute();
				N++;
				
			}
		else
			{
			initialise();

			if (debug !=2) compute(); //not climatic stats only
			N++;
			}
		} // end of ini file
		
	//if (TRANSIENT)
	//	{
	//	transient = fopen("transient_out.csv","a");
		//fprintf(transient,"\n");
	//	fclose(transient);
	//	}
	
	fclose (init_file);
}

int main(int argc, char * argv[])
{
	char buffer[40];
	int random1;
	long double random2;
	FILE *out, *random_file;
	printf("SurEau by H. Cochard UMR Piaf-INRAE version:%s\n",version) ;
	
	number=1;
	
	strcpy (filename_IN, "sureau_ini.txt");  // default values when no name is given as argument
	strcpy (filename_OUT,"sureau_out.csv");
	strcpy (filename_OUT2,"annual_out.csv");
	//strcpy (filename_CLIM, "climat_day_in.txt");
	strcpy (filename_TRANS,"transient_out.csv");
	
	if (argc==2) // passed with iXXX_sureau.ini OR cXXX_climat_day.ini as argument
	{
		strncpy(filenumber,argv[1],1);                                  	// first character of the argument
		if  (strcmp(filenumber,"i")  ==0) strcpy(filename_IN,argv[1]);  	// a ini file
		else strcpy(filename_CLIM,argv[1]);                             	// a clim file
		
		strncpy(filenumber,argv[1],4);                                  	//extract the file number from file name
		number=atoi(filenumber+1);
	
		strcpy(buffer,filenumber);											//print sureau_out
		strcat(buffer,"_sureau.out1");
		strcpy(filename_OUT,buffer);
		
		strcpy(buffer,filenumber);											//print annual_out
		strcat(buffer,"_annual.out2");
		strcpy(filename_OUT2,buffer);
		
		strcpy(buffer,filenumber);											//print transient_out
		strcat(buffer,"_transient.out3");
		strcpy(filename_TRANS,buffer);
	}
	
	if (argc==3) // passed with cXXX_climat_day.ini AND iXXX_sureau.ini as argument
	{
	   // strncpy(filenumber,argv[1],1);                  	          	// first character of the argument
		strcpy(filename_CLIM,argv[1]);                           			// clim file
		strcpy(filename_IN,argv[2]);  										// ini file
		
		strncpy(filenumber,argv[1],4);                               	// extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcat(buffer,filenumber);
		strcat(buffer,"_sureau.out1");
		strcpy(filename_OUT,buffer);
		
		strncpy(filenumber,argv[1],4);                                  // extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcpy(buffer,filenumber);
		strcat(buffer,"_annual.out2");
		strcpy(filename_OUT2,buffer);
		
		strncpy(filenumber,argv[1],4);                                  // extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcpy(buffer,filenumber);
		strcat(buffer,"_transient.out3");
		strcpy(filename_TRANS,buffer);
	}
	
	if ((out = fopen("sureau_out.csv","r"))==NULL) //file does not exist then print headers in the file
	{
		out = fopen("sureau_out.csv","w");
		fprintf(out, "N1;N2;Simul;YEAR;DOY;T_PLC_Leaf1;T_PLC_Leaf2;T_PLC_Leaf3;T_PLC_Axil;T_PLC_Branch1;T_PLC_Branch2;T_PLC_Branch3;T_PLC_Trunk;T_PLC_Root;T_PLC_Root1;T_RWC_Axil;T_REW_Soil; T_gs_regul;T_gs_50mmol;T_gs_close;T_budbreak;T_max_LAI;E_tot;A_net_tot;GPP;Radius;T_Leaf_max;PLC_Leaf;PLC_Leaf1;PLC_Leaf2;PLC_Leaf3;PLC_Branch;PLC_Branch1;PLC_Branch2;PLC_Branch3;REW_int1;REW_int2;RWC_int;RWC_min;P_soil_min;P_min_lf1;P_min_br1;RWC_br1;RWC_br2;RWC_br3;T_TLP_Leaf1;T_TLP_Axil;T_RWC_Br");;
	   
		if ((random_file = fopen("random_para.txt","r"))==NULL) 
		{
			printf("\a\nCan't find file random_para.txt!!");
			random_file = fopen("random_para.txt","w");
			fprintf(random_file,"79	10	-25	25");
			fclose(random_file);
			RANDOMISE=0;
		}
		else while (!feof(random_file)) 
			{
			fscanf(random_file,"%d %Le %Le %Le\n",&random1,&random2,&random2,&random2);
			fprintf(out,"para# %d;",random1);
			}
			//fprintf(out,"\n");
			fclose(random_file);
		//fclose(out);
		fprintf(out,"\n");
		fclose(out);
	}
	
	if (argc==2) srand((unsigned) number+time(NULL));
	
	else srand((unsigned) time(NULL));
	setup();
	
	
}
