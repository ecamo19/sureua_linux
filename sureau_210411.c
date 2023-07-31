// SurEau_C
// by Hervé Cochard, INRAE, UMR-PIAF, Clermont-Ferrand, France
// still uder developement...
// This model should not be distributed without prior notice to H Cochard
// thanks for reporting bugs...


#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

// Librairies
#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "time.h"
#include "stddef.h"
#include <string.h>

void print_screen(void);
void print_transient(void);
void Reset(void);
void Get_DATA(long double);
void Compute_ST(void);
void Randomise(void);
void setparameters(void);
void initialise(void);
void compute(void);
void setup(void);

char version[]=" 11-04-2021";

long double C3=1;
long double debug=0;  // if 1 then print ALL the data in the file out
long double dq_Soil1;
long double DYNAMIC0=1,DYNAMIC=1;
long double PRINT_GRAPH=1;
int WARNING=0;


long double para[250],para0[250];
char filenumber[5]="1";
int year=1,N_days=0,N_days2=0;
long indice=1;
long double indice_double=1;
long double FLUID=1;                // to account for the temperature dependance of K
long double SURFACE_TENSION=1;      // to account for the temperature dependance of P50
long double T_g_cuti=1;               // to account for the temperature dependance of g_min
long double T_OSMOTIC=1;            // to accout for the temperature dependance of osmotic potentials
long double TLEAF=1;                // if one then compute radiative leaf temperature
long double IRRIGATE=1, Irrigation, IRR_DOY_S, IRR_DOY_F;             // for irrigation
long double CONTINUOUS;             // if one the simulation is not reinitialised at the end of the year
long double FRACTAL;                // if one the tree is fractal
int Screen_out[250];
int File_out[250];
long double beta=0.5;
char *Label[]=    {  "  Days   ", "T_air      ", "T_leaf     ", "RH_air     ", "PAR        ", "VPD_Leaf   ", "VPD_Cut    ", "VPD_Axil   ", "VPD_Branch ", "VPD_Trunk  ",
	"VPD_Root1  ", "VPD_Root2  ", "VPD_Root3  ", "VPD_Soil   ", "PAR        ", "P_AIR      ", "E_clim_m2  ", "E_Leaf_m2  ", "E_cuti_m2  ", "E_Branch_m2", "E_Trunk_m2 ",
	"E_Axil_m2  ", "E_Root1_m2 ", "E_Root2_m2 ", "E_Root3_m2 ", "E_Soil_m2  ", "E_plant    ", "E_Plant_m2 ", "E_tot kg   ", "g_canopy   ", "g_s        ", "g_cuti     ",
	"P_evap_apo ", "P_leaf_sym ", "P_leaf_apo ", "turgor     ", "Turgor_leaf", "P_Axil_symp", "P_Axil_apo ", "P_branch_s ", "P_branch_a ", "P_trunk_s  ", "P_trunk_a  ",
	"Turgor     ", "P_root1_sy ", "P_endo1    ", "P_root1_ap ", "P_root2_sy ", "P_endo2    ", "P_root2_ap ", "P_root3_sy ", "P_endo3    ", "P_root3_a  ", "P_soil1    ", "P_soil2    ",
	"P_soil3    ", "RWC1       ", "RWC2       ", "RWC3       ", "K_Soil1    ", "K_Soil2    ", "K_Soil3    ", "K_Inter1   ", "K_Inter2   ", "K_Inter3   ", "K_leaf_s   ",
	"K_leaf_a   ", "K_leaf     ", "K_root1_s  ", "K_root1_a  ", "K_root2_s  ", "K_root2_a  ", "K_root3_s  ", "K_root3_a  ", "K_root     ", "K_trunk    ", "K_Plant    ",
	"K_tot      ", "PLC_leaf   ", "PLC_branch ", "PLC_Axil   ", "PLC_trunk  ", "PLC_root1  ", "PLC_root2  ", "PLC_root3  ", "Q_plant   ", "DQ_soil    ", "Q_evap_a   ",
	"Q_leaf_s   ", "Q_leaf_a   ", "Q_branch_s ", "Q_branch_a ", "Q_Axil_s   ", "Q_Axil_a   ", "Q_trunk_s  ", "Q_trunk_a  ", "Q_root_s   ", "Q_root_a   ", "Q_endo     ", "Q_soil1    ",
	"Q_soil2    ", "Q_soil3    ", "Growth_r   ", "Radius     ", "A_net      ", "A_net_tot_c", "WUE        ", "  Rm       ", "Reserve    ", "Leaf_Area  ", "Leaf_wc    ", "Branch_wc  ",
	"Br_Sy_wc   ", "Br_Apo_wc  ", "RWC_shoot  ", "RWC_Axil   ", "RWC_soil   ", "RWC_min    ", "RWC_int    ", "P_soil     ", "P_soil_min ", "P_soil_int ", "Irrigation ", "Sap_Flow_d ",
	"Sap_Flow_d2", "VPD_Air    ", "Fruit_Diam ", "ETP_Pen_t  ", "ETP_Pen    ", "Rain_tot   ", "Rain_soil  ", "VPD_air_t  ", "VPD_leaf_t ", "Rain_leaf_t", "ETP_leaf_t ",
	"T_air_an   ", "T_air_an_l ", "EvapoT     ", "EvapoT mm  ", "ETP_d, mm  ", "Sap_flow   ", "K_branch   ", "P_min_lf   ", "P_min_br   ", "K_tot2     ", "A_net_c    ",
	"Intercep   ", "GPP        ", "GPP_day    ", "ETR_d, mm  ", "E_day, mm  ", "T_Soil     ", "PAR_pot    ", "Cloud      ", "F_ls_le    ", "F_la_le    ", "F_ba_la    ", "F_ba_bs    ",
	"F_ta_ba    ", "F_ta_ts    ", "F_bua_bus  ", "F_ba_bua   ", "F_ra_ta    ", "F_re_ra    ", "F_rs_re    ", "F_s_re1    ", "F_s_12     ", "F_s_23     ", "A_net1     ", "A_net2     ",
	"g_Axil     ", "Root-Area  ", "Axil_WC    ", "RU_soil_t  ", "ETR mm     ", "Drainage   ", "A_gross    ", "A_net_tot  ", "A_gros_tot ", "Resp_tot   ", "Export     ", "Export_tot ",
	"T_Soil1    ", "T_Soil2    ", "T_Soil3    ", "RU_soil_wp ", "REW_t      ", "Q_pet_s    ", "E_Petiole  ", "P_Pet_s    ", "Pet_diam   ", "Q_Plant_a  ", "Q_Plant_s  ",
	"RWC_leaf   ", "RWC_leaf_s ", "RWC_leaf_a ", "RWC_Branch ", "RWC_Br_s   ", "RWC_Br_a   ", "RWC_Trunk  ", "RWC_Tr_s   ", "RWC_Tr_a   ", "RWC_Root   ", "RWC_Root_s ",
	"RWC_Root_a ", "RWC_Plant  ", "RWC_Plant_s", "RWC_Plant_a", "E_day mmol ", "P_min_lf_d ", "P_max_lf_d ", "gs_min_d   ", "gs_max_d   ", "SF_min_d   ", "SF_max_d   ", "gs_max_d2   ",
	"Fruit_gr   ", "Teta_Soil  ", "Teta_Soil1 ", "Teta_Soil2 ", "Teta_Soil3 ", "Ca  ppm    ", "K_Plant_20 ", "PLC_Plant  ", "P50_leaf   ", "P50_bran   ", "P50_trunk  ", "P50_root     ",
	"PI0_leaf   ", "Px_gs      ", "PI_leaf    ","","","","","","",""
	,"","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","",""};

char filename_IN[]="ixxxxxxx_sureau.ini";
char filename_OUT[]="ixxxxxxx_sureau.out1";
char filename_OUT2[]="ixxxxxxx_sureau.out2";
char filename_CLIM[]="cxxxxxxx_climat_day.ini";
char filename_TRANS[]="xxxxxxx_transient.out1";
long double DATA[250];
long TT=0;
int DEAD=0;
int END_CLIMAT=0;
int N=1;
long double CAPILLARITY=1;  // allow flow between soil layers by capillarity
int PREM=1, PREM1=1;
long double RANDOMISE=0;
long double TRANSIENT=0, gs_cst=0,PRINT_SCREEN=0;
long double YEAR1=2000, YEAR2=2001, dt,dt_dyna,dt_stat,T,T0,days_simul,T_max;     // dt time interval in seconds; T time of the day in seconds
long double CUT;                    // to simulate a cut branch (=1) or  trunk (=2)
long double CLIMAT;                 // if 1 then climatic data are from a data file; if 0 they are generated
//long j;
long double Lat=45;                     // latitude, between -90 and 90°
long double Day_length;              // length of day light, hours
long double DOY=180;                     // day of year 0-365
long double dPLC_crit=0.01,REW_crit=0.8;

//Atmosphere
long double cloud_cover, POTENTIAL_PAR, CO2_atm, T_air, T_air_min, T_air_max, Cum_T_air,Cum_T_air_l,T_air_an, T_air_an_l,RH_air, RH_air_min, RH_air_max, RH_air_min_0,P_air, VPD,VPD_Air_tot,VPD_Leaf_tot,PAR,PAR_max,Rain_1, Rain_2, Rain_tot,Rain_soil, Rain_leaf_tot, Proba_rain, Wind;
long double T_0, T_1, T_air_1, RH_air_1,PAR_1, T_2, T_air_2, RH_air_2,PAR_2,T_air_max_0, T_air_min_1, T_air_max_1, RH_air_min_1, RH_air_max_1,PAR_max_1,T_air_min_2, T_air_max_2, RH_air_min_2, RH_air_max_2,PAR_max_2,VPD_Air;
long double HW=1;                   //Heat Wave
long double HW_day, HW_duration, HW_T;
FILE *climat_in;
long double HH1=12.00;  // time of Tmax
long double HH2=0.0;    // time of RHmax
long double INTERCEPTION,    Interception_min,    Threshold_rain,    Interception_factor; // for rain interception model
long double ETP_Penman,ETP_Penman_tot,ETP_Penman_day,ETP_leaf_tot;
//SOIL
long double Rock_f1=0,Rock_f2=0,Rock_f3=0, Teta_Soil,Teta_Soil1,Teta_Soil2,Teta_Soil3, REW_t,  T_REW_Soil, gap, VPD_Soil,  T_Soil,T_Soil1,T_Soil2,T_Soil3, T_Soil_1, T_Soil_2,T_Soil_11, T_Soil_21, T_Soil_12, T_Soil_22,T_Soil_13, T_Soil_23,E_Soil, g_Soil,g_Soil0, dq_Soil,  Q_Soil1,Q_Soil2,Q_Soil3, Q_Soil01,Q_Soil02,Q_Soil03, K_Soil1,K_Soil2,K_Soil3,K_Soil_tot1,K_Soil_tot2,K_Soil_tot3, K_Interface1,K_Interface2,K_Interface3;
long double P_Soil1, P_Soil2, P_Soil3, Teta_s, Teta_fc, Teta_r,Teta_wp,alpha,m,n,K_sat,L,Volume_soil,Volume_soil1,Volume_soil2,Volume_soil3,Soil_Depth, Soil_Width, K_s,RWC1,RWC2,RWC3,RWC_Soil,RWC_min,RWC_int,RWC_Irr,RWC_fc, P_soil,P_soil_min,P_soil_int;
long double Surface_Soil, dq_Soil_Root_Symp,dq_Soil_12,dq_Soil_23,Fluidity_soil, Fluidity_soil1, Fluidity_soil2, Fluidity_soil3,Osmotic_TSoil,ST_Soil;
long double REHYDRATE=1, PLC_REHYD=99.99, Daily_Irr=1, Drainage ;
long double PENMAN=1, Penman_Coeff=0.75;
long double Teta_soil, teta_c=0.232;
long double PLC_LIMIT=0, T_LIMIT, T_Soil_Crit;
long double Q_soil1_init, Q_soil2_init, Q_soil3_init,COMPET;

//PLANT
FILE *File_fractal1,*File_fractal2;
long double Length_term_shoot,side[1001],N_fractal, Nb_segment[1001],Diam_segment[1001],Diam_sapwood[10001],Length_segment[1001],scale,angle,X[1001],Y[1001];
long double LA_term_shoot=100, Length_fine_root=50, Tapering=1.8,root_shoot_ratio=0.2,root_ramif=0.3;  //parameters for fractal
long double SF_min_d,SF_max_d,Q_Plant,Q_Plant_s,Q_Plant_a, Q_Plant0,Q_Plant_s0,Q_Plant_a0,K_Plant;
long double REFILL,P_REFILL; // to allow xylem refilling above P_REFILL
long double SYMP_CAVIT=0; // then water released by cavitation goes only into the adjacent symplasm; activated if REFILL=2
long double GPP, EvapoT_day,E_tot_day,GRAVITY, Pg_Leaf, Pg_branch, Pg_trunk,Extinction_Coeff, A_net_c, A_gross, A_gross_tot,A_net_tot_c,A_gross_tot_c,A_net_day_c,Resp_tot_c; // the gravimetric pressures
long double Regul_gs, Regul_gs_para1, Regul_gs_para2; // if 1 E regulated by leaf turgor, if 0 E regulated by P12
long double ST_Leaf=1, ST_Air=1; // Water surface tension at leaf and air temperatures relative to SF at 20 °C
long double Osmotic_TLeaf, Osmotic_TAir;
long double g_Canopy,gs_Jarvis,gs_max,gs_night, g_cuti, g_cuti_20, g_Branch, g_Trunk,g_Root,g_Root1,g_Root2,g_Root3, Jarvis_PAR,g_bl,g_crown,g_crown0;
long double Reserve,  Rm,Rg, Q_Wood;
long double K_tot,K_tot0; // Hydraulic conductance soil to leaf
long double Flow_save, Q_buffer;
long double Fluidity_leaf,Fluidity_air;
long double E_max=0, E_max_gs_close=0;
long double Px_gs; // the set xylem pressure for stomatal closure
long double turgor,K_to_Leaf_Apo,E_clim, E_t_B;
long double Slope_Leaf_Fall,P50_Leaf_Fall;
long double Tgs_optim=20, Tgs_sens=17; // variable to describe the effect of leaf temperature on gs_max;
long double Leaf_Fall=0;
long double K_VAR=0, K_VAR_P1,K_VAR_P2,K_VAR_P3; // variable K: 0 = Kleaf CST Kroot CST 1= Kleaf VAR Kroot CST 2: Kleaf CST Kroot VAR 3: Kleaf VAR Kroot VAR
long double END_DEATH=0;
long double PLC_END; // the PLC defining hydraulic failure and computation end
long double gs_tc=0; //  time constant of stomatal response in s
//Fractal tree
long double Q_Branch_Symp_FR, Q_Branch_Apo_FR, Q_Trunk_Sym_FR, Q_Trunk_Apo_FR, Q_Root_Symp_FR, Q_Root_Apo_FR, Branch_Area_FR, Trunk_Area_FR, Root_Area_FR, Root_Area_fi_0,Root_Area_fi, K_Branch_Apo_FR, K_Trunk_Apo_FR, K_Root_Apo_FR, Length_Root_FR, Diam_Root_FR ;

// variables for Plant
long double Alpha, Leaf_Area, LAI, Branch_Area, Trunk_Area, Root_Area, Root_Area1,Root_Area2,Root_Area3, E_Leaf, dq_cuti, dq_stomate, dq_Branch, dq_Trunk;
long double dq_Root1,dq_Root2,dq_Root3, E_cuti, E_Branch, E_Trunk, E_Root1,E_Root2,E_Root3, t_out, E_tot, EvapoT, EvapoT_tot, g;
long double P_min_leaf, P_min_stem; // minimal seasonal water potentials
long double K_Plant_20,K_Plant_20_0,PLC_Plant; // Whole plant conductance à 20°C. Use to compute PLC_Plant as 100*(1-K_Plant_20/K_Plant_0)

// variabes for legacy of PLC
long double leg_Leaf=0, leg_Branch=0.5, leg_Trunk=0.75, leg_Root=0;

// variables for growth
long double Growth_trunk, Growth_rate_trunk, Growth_fruit, Growth_rate_fruit, Extensibility_trunk, Yield_trunk, Extensibility_fruit, Yield_fruit;

// variables for Leaf
long double Leaf_rain=0, gs_min_d,gs_max_d,gs_max_d2;
long double Pgs_88,Pgs_12,P_min_lf_d ,P_max_lf_d;
long double gs_0, GS_MAX, E_Leaf_night, T_Leaf_max, TP, Q10_1, Q10_2, LMA, Succulence, Leaf_apo_fraction, PLF, g_s, Leaf_size, A_net, A_net1, A_net2, A_net_max, A_net_tot, A_net_day, Turgor_lp, VPD_Leaf, VPD_Cuti,T_Leaf,Px_Leaf_Apo;
long double PLCx, P12_Leaf_Apo,PLC_Leaf_Apo, Q_Leaf_Evap, Q_Leaf_Evap0, Q_Leaf_Symp, Q_Leaf_Apo, dQ_Leaf_Apo, Q_Leaf_Symp0, Q_Leaf_Apo0, Q_Leaf_Apo1, K_Leaf_Symp, K_Leaf_Symp_0, K_Leaf_Symp2;
long double K_Leaf_Apo, K_Leaf_Apo0, P_Leaf_Evap, P_Leaf_Symp, Turgor_Leaf_Symp,Turgor_Leaf_Symp_Ref, P_Leaf_Apo, C_Leaf_Evap, C_Leaf_Apo0,C_Leaf_Apo, Epsilon_Leaf_Symp, Pi0_Leaf_Symp;
long double P50_Leaf_Symp, Slope_Leaf_Symp, P50_Leaf_Apo, Slope_Leaf_Apo;
long double LA_Var, LA_max_Pheno, LA_min, LA_max, LA_max_init, LA_day1, LA_day2, LA_day3, LA_day4;
long double VcMax=110, VjMax=180, Kc25=40.4, Ko25=24800, Qye=0.12, Rd25=-0.37, Ca=400;
long double a_Res, Resp_tot;
long double Export, Export_tot=0;   //exportation of µmol of glucose /m2/s
long double leaf_angle=45;    //angle from horizontal    degrees
long double Gamma=0.8;       //an attenuation factor for E to reach a target potential
long double turgor_ref_factor;
long double gs_CO2_sens=0;

// variables for Axil: Axillary bud /Fruit/Flower
long double Type_Axil, T_RWC_Axil, T_PLC_Axil, RWC_Axil, E_Axil, E_Petiole, T_Axil, N_Axil, dq_Axil, K_Axil_Apo0,K_Axil_Apo, K_Axil_Symp, K_Axil_Symp2, Epsilon_Axil_Symp, Pi0_Axil_Symp;
long double Petiole_diam, P50_Axil_Apo, Slope_Axil_Apo, VPD_Axil, VPD_Petiole, g_Axil, g_Axil_min20, g_Axil_min,g_Axil_max, g_Petiole , C_Axil_Apo, Diam_Axil, Area_Axil, Petiole_area, WC_Axil;
long double PLC_Axil_Apo,Q_Axil_Symp, Q_Axil_Apo, Q_Axil_Symp0, Q_Petiole_Symp, Q_Petiole_Symp0, Q_Axil_Apo0,Q_Axil_Apo1, P_Axil_Symp, P_Petiole_Symp,P_Axil_Apo,Turgor_Axil_Symp, dq_fruit, Length_Petiole, Diam_Petiole;

// variables for Branch
long double RWC_shoot, Density, Branch_symp_fraction, Branch_apo_fraction,Length_Branch, Number_Branch, Diam_Branch, VPD_Branch,T_Branch,PLC_Branch_Apo,Q_Branch_Symp, Q_Branch_Apo, Q_Branch_Symp0;
long double Q_Branch_Apo0, Q_Branch_Apo1, K_Branch_Symp, K_Branch_Symp0,K_Branch_Apo, K_Branch_Apo0,P_Branch_Symp, P_Branch_Apo, C_Branch_Apo, C_Branch_Apo0,Epsilon_Branch_Symp, Pi0_Branch_Symp;
long double P50_Branch_Symp, Slope_Branch_Symp, P50_Branch_Apo, Slope_Branch_Apo;

// variables for Trunk
long double Radius_bark_trunk,Trunk_sapwood_fraction, Trunk_symp_fraction, Trunk_apo_fraction,Length_Trunk, Diam_Trunk,Turgor_Trunk_Symp,VPD_Trunk, T_Trunk,PLC_Trunk_Apo,Q_Trunk_Symp, Q_Trunk_Apo;
long double Q_Trunk_Symp0, Q_Trunk_Apo0,Q_Trunk_Apo1, K_Trunk_Symp, K_Trunk_Symp0,K_Trunk_Apo,K_Trunk_Apo0, P_Trunk_Symp, P_Trunk_Apo, C_Trunk_Apo,C_Trunk_Apo0, Epsilon_Trunk_Symp, Pi0_Trunk_Symp;
long double P50_Trunk_Symp, Slope_Trunk_Symp, P50_Trunk_Apo, Slope_Trunk_Apo;

// variables for Root
long double K_root_0, Root_symp_fraction, Root_apo_fraction,Length_Root, Length_Root_fi, Diam_Root,Q_Root_Endo0, Q_Root_Symp0, Q_Root_Apo_t0, Q_Root_Apo_t1,K_Root_Sympqq,K_Root_Symp1,K_Root_Symp2,K_Root_Symp0,K_Root_Apo,K_Root_Apo0, C_Root_Endo, C_Root_Apo,C_Root_Apo0;
long double Epsilon_Root_Symp, Pi0_Root_Symp, P50_Root_Symp, Slope_Root_Symp, P50_Root_Apo, Slope_Root_Apo, PRF;
long double VPD_Root1,PLC_Root_Apo1,Q_Root_Endo1,Q_Root_Endo01, Q_Root_Symp1, Q_Root_Apo1, Q_Root_Symp01, Q_Root_Apo01, Q_Root_Apo11, K_Root_Symp11,K_Root_Symp21, K_Root_Apo1,K_Root_Apo01;
long double P_Root_Endo1, P_Root_Symp1, P_Root_Apo1, C_Root_Endo1, C_Root_Apo1,C_Root_Apo01, Epsilon_Root_Symp1, Pi0_Root_Symp1, P50_Root_Symp1, Slope_Root_Symp1, P50_Root_Apo1, Slope_Root_Apo1;
long double VPD_Root2,PLC_Root_Apo2,Q_Root_Endo2,Q_Root_Endo02, Q_Root_Symp2, Q_Root_Apo2, Q_Root_Symp02, Q_Root_Apo02, Q_Root_Apo12, K_Root_Symp12,K_Root_Symp22, K_Root_Apo2,K_Root_Apo02;
long double P_Root_Endo2, P_Root_Symp2, P_Root_Apo2, C_Root_Endo2, C_Root_Apo2,C_Root_Apo02, Epsilon_Root_Symp2, Pi0_Root_Symp2, P50_Root_Symp2, Slope_Root_Symp2, P50_Root_Apo2, Slope_Root_Apo2;
long double VPD_Root3,PLC_Root_Apo3,Q_Root_Endo3,Q_Root_Endo03, Q_Root_Symp3, Q_Root_Apo3, Q_Root_Symp03, Q_Root_Apo03, Q_Root_Apo13, K_Root_Symp13,K_Root_Symp23, K_Root_Apo3,K_Root_Apo03;
long double P_Root_Apo, Q_Root_Apo_t, P_Root_Endo3, P_Root_Symp3, P_Root_Apo3, C_Root_Endo3, C_Root_Apo3,C_Root_Apo03, Epsilon_Root_Symp3, Pi0_Root_Symp3, P50_Root_Symp3, Slope_Root_Symp3, P50_Root_Apo3, Slope_Root_Apo3;
long double T_Root,Root_upper,Root_middle,Root_lower;

// variables for dq
long double dq_Leaf_Apo_Leaf_Symp,dq_Axil_Apo_Axil_Symp, dq_Leaf_Symp_Leaf_Evap,dq_Branch_Apo_Leaf_Apo,dq_Branch_Apo_Axil_Apo, dq_Branch_Apo_Branch_Symp, dq_Trunk_Apo_Branch_Apo;
long double dq_Trunk_Apo_Trunk_Symp, dq_Root_Apo_Trunk_Apo1, dq_Root_Symp_Root_Apo1,dq_Root_Apo_Trunk_Apo2, dq_Root_Symp_Root_Apo2,dq_Root_Apo_Trunk_Apo3, dq_Root_Symp_Root_Apo3;
long double dq_Leaf_Apo_Leaf_Evap,dq_Root_Endo_Root_Apo1, dq_Root_Symp_Root_Endo1, dq_Soil_Root_Endo1,dq_Root_Endo_Root_Apo2, dq_Root_Symp_Root_Endo2, dq_Soil_Root_Endo2;
long double dq_Root_Endo_Root_Apo3, dq_Root_Symp_Root_Endo3, dq_Soil_Root_Endo3, dq_Branch_Symp_Axil_Symp,dq_Axil_Apo_Petiole_Symp, dq_Petiole;
long double T_PLC_Leaf,T_PLC_Branch,T_PLC_Trunk,T_PLC_Root,T_PLC_Root1, T_gs_close,T_gs_regul;
long double dq_Root1_Root2,dq_Root1_Root3,dq_Root2_Root3, dq_Root_Apo_Trunk_Apo;

void init(void)  //intilialize a bunch of variables
{
	FILE *soil;
   // g_Axil_max= 67.39;
	long double K_Root, K_root1, K_root2, K_root3;
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
	LAI=Leaf_Area/Surface_Soil;   // LAI in m2/m2
	WARNING=0;
	T_air_an=0;
	T_air_an_l=0;
	Cum_T_air=0;
	Cum_T_air_l=0;
	N_days=0;
	N_days2=0;
	Rain_tot=0;
	Rain_leaf_tot=0;
	VPD_Air_tot=0;
	VPD_Leaf_tot=0;
	ETP_leaf_tot=0;
	A_net_day=0;
	E_tot_day=0;
	EvapoT_day=0;
	indice=1;
	indice_double=0;
	DEAD=0;
	END_CLIMAT=0;
	REW_t=1;
	T_PLC_Leaf=9999*3600*24;
	T_PLC_Branch=9999*3600*24;
	T_PLC_Trunk=9999*3600*24;
	T_PLC_Root=9999*3600*24;
	T_PLC_Root1=9999*3600*24;
	T_RWC_Axil=9999*3600*24;
	T_PLC_Axil=9999*3600*24;
	T_REW_Soil=9999*3600*24;
	T_gs_close=0;
	T_gs_regul=9999*3600*24;
	g_bl=1000;
	ST_Leaf=1.0;
	ST_Air=1.0;
	T_Leaf=T_air;
	g_cuti=g_cuti_20;
	turgor=1;
	E_max=0;
	E_max_gs_close=0;
	T_Leaf_max=0;
	E_Leaf_night=100;
	LA_max_Pheno=LA_max;
	Root_Area_fi=Root_Area_fi_0;
	Rain_soil=0;
	ETP_leaf_tot=0;
	ETP_Penman_tot=0;
	g_Soil=g_Soil0;
	Drainage=0;
	if (CUT==2)  // cut between branch and trunk
	{
		K_Trunk_Apo0=0;
		Turgor_Leaf_Symp_Ref=-Pi0_Leaf_Symp;
	}
	if (CUT==1) // cut between trunk and root
	{
		K_Branch_Apo0=0;
		Turgor_Leaf_Symp_Ref=-Pi0_Leaf_Symp;
	}
	Leaf_Area=LA_max_Pheno;
	Turgor_Leaf_Symp=-Pi0_Leaf_Symp;
	
	// valeur en attendant
	Q_Leaf_Evap0   = Q_Leaf_Apo/1000;
	C_Leaf_Evap    = C_Leaf_Apo;
	Q_Root_Endo01  = Q_Root_Apo_t0/1000;
	C_Root_Endo1   = C_Root_Apo;
	Q_Root_Endo02  = Q_Root_Apo_t0/1000;
	C_Root_Endo2   = C_Root_Apo;
	Q_Root_Endo03  = Q_Root_Apo_t0/1000;
	C_Root_Endo3   = C_Root_Apo;
	
   
	
	// initialise roots
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
	g_Root1=g_Root;
	g_Root2=g_Root;
	g_Root3=g_Root;
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
	K_Root_Symp11=K_Root_Symp1 * Root_upper;
	K_Root_Symp12=K_Root_Symp1 * Root_middle;
	K_Root_Symp13=K_Root_Symp1 * Root_lower;
	P50_Root_Apo1=P50_Root_Apo;
	P50_Root_Apo2=P50_Root_Apo;
	P50_Root_Apo3=P50_Root_Apo;
	Slope_Root_Apo1=Slope_Root_Apo;
	Slope_Root_Apo2=Slope_Root_Apo;
	Slope_Root_Apo3=Slope_Root_Apo;
	P_min_leaf = P_min_stem=0;
	Compute_ST();
	Px_Leaf_Apo= P50_Leaf_Apo*ST_Leaf+ 25/Slope_Leaf_Apo*log((100-PLCx)/PLCx); //compute the Pressure for x PLC
	A_net_tot=0;
	A_net_max=0;
	// Growth=0;
	// Reserve =50000;
	m=(1-1/n);
	Flow_save=0;
	
	if (CONTINUOUS==0 || PREM==1)
	{
		RWC1=RWC_fc;  // set RWC to field capacity
		RWC2=RWC_fc;
		RWC3=RWC_fc;
		RWC_Soil=RWC_fc;
		RWC_min=RWC_fc;
		RWC_int=RWC_fc;
		P_soil=-0.033;  // water potential at field capacity (Pf=2.5=-33kPa)
		P_soil_min=-0.033;
		P_soil_int=-0.033;
		P_Soil1=-0.033;
		P_Soil2=-0.033;
		P_Soil3=-0.033;
		Q_Soil01        =Teta_fc*Volume_soil1*1000*1000*1000/18;    // water in bulk soil in mmol at field capacity
		Q_Soil02        =Teta_fc*Volume_soil2*1000*1000*1000/18;
		Q_Soil03        =Teta_fc*Volume_soil3*1000*1000*1000/18;
		Q_Soil1        = Q_Soil01;
		Q_Soil2        = Q_Soil02;
		Q_Soil3        = Q_Soil03;
		Q_soil1_init   = Q_Soil01;
		Q_soil2_init   = Q_Soil02;
		Q_soil3_init   = Q_Soil03;
		if (COMPET)
		{
		soil = fopen("soil.trs","w");
		fprintf(soil,"%Lf %Lf %Lf",Q_Soil1,Q_Soil2,Q_Soil3);
		fclose(soil);
		}
		K_Soil1=1E8;
		K_Soil2=1E8;
		K_Soil3=1E8;
		K_Interface1= K_Soil1*10;
		K_Interface2= K_Soil2*10;
		K_Interface3= K_Soil3*10;
		P_Branch_Apo=0;
		P_Axil_Apo=0;
		P_Leaf_Apo=0;
		P_Trunk_Apo=0;
		P_Root_Apo1=0;
		P_Root_Apo2=0;
		P_Root_Apo3=0;
		Drainage=0;
	}
	
	T_Leaf=T_air_min; //leaf temperature at least equal to min air temperature to begin with
	E_tot=0;
	EvapoT_tot=0;
	dq_Branch_Apo_Branch_Symp=0;
	dq_Branch_Apo_Leaf_Apo=0;
	dq_Leaf_Apo_Leaf_Symp=0;
	dq_Root_Apo_Trunk_Apo1=0;
	dq_Root_Symp_Root_Apo1=0;
	dq_Root_Apo_Trunk_Apo2=0;
	dq_Root_Symp_Root_Apo2=0;
	dq_Root_Apo_Trunk_Apo3=0;
	dq_Root_Symp_Root_Apo3=0;
	dq_Trunk_Apo_Branch_Apo=0;
	dq_Trunk_Apo_Trunk_Symp=0;
	dq_Branch=E_Branch*Branch_Area*dt;
	dq_Trunk=E_Trunk*Trunk_Area*dt;
	dq_Root1=E_Root1*Root_Area1*dt;
	dq_Root2=E_Root2*Root_Area2*dt;
	dq_Root3=E_Root3*Root_Area3*dt;
	K_Leaf_Symp = K_Leaf_Symp_0;
	K_Leaf_Symp2=1*K_Leaf_Symp;                   // assume the resitance of the symplasmic sap pathway is the same as the resitance to the symplasmis reservoir
	K_Root_Symp21=K_Root_Symp2 * Root_upper;
	K_Root_Symp22=K_Root_Symp2 * Root_middle;
	K_Root_Symp23=K_Root_Symp2 * Root_lower;
	
	K_Leaf_Apo=K_Leaf_Apo0;
	K_Axil_Apo=K_Axil_Apo0;
	K_Branch_Apo=K_Branch_Apo0;
	K_Trunk_Apo=K_Trunk_Apo0;
	K_Root_Apo1=K_Root_Apo01;
	K_Root_Apo2=K_Root_Apo02;
	K_Root_Apo3=K_Root_Apo03;
	
	if (K_Root_Symp11 && K_Root_Apo1) K_root1=1/(1/K_Root_Symp11 + 1/K_Root_Apo1); else K_root1=0;
	if (K_Root_Symp12 && K_Root_Apo2) K_root2=1/(1/K_Root_Symp12 + 1/K_Root_Apo2); else K_root2=0;
	if (K_Root_Symp13 && K_Root_Apo3) K_root3=1/(1/K_Root_Symp13 + 1/K_Root_Apo3); else K_root3=0;
	K_Root=K_root1+K_root2+K_root3;
	K_root_0=K_Root;
	if (K_Leaf_Apo && K_Leaf_Symp && K_Branch_Apo && K_Trunk_Apo && K_Root)            
			K_Plant_20_0=  1/(1/K_Leaf_Apo +1/K_Leaf_Symp + 1/K_Branch_Apo+ 1/K_Trunk_Apo + 1/K_Root);  
	else  K_Plant_20_0=0;
	
	P12_Leaf_Apo=P50_Leaf_Apo+50/Slope_Leaf_Apo;
	if (CLIMAT!=2)
	{
		PLC_Leaf_Apo=0;
		PLC_Branch_Apo=0;
		PLC_Trunk_Apo=0;
		PLC_Root_Apo1=0;
		PLC_Root_Apo2=0;
		PLC_Root_Apo3=0;
	}
	
	Surface_Soil=Soil_Width*Soil_Width; //surface of the soil in m2 assuming a cube
	//Volume_soil  = Surface_Soil*Soil_Depth*(1-1/3*(Rock_f1+Rock_f2+Rock_f3)); //rock in the soil contain no available water for roots
	Volume_soil1 = Surface_Soil*Soil_Depth*(1-Rock_f1)/3;
	Volume_soil2 = Surface_Soil*Soil_Depth*(1-Rock_f2)/3;
	Volume_soil3 = Surface_Soil*Soil_Depth*(1-Rock_f3)/3;
	Volume_soil  = Volume_soil1 + Volume_soil2 +Volume_soil3;
	
	Irrigation     =0;
	C_Leaf_Apo0    =C_Leaf_Apo;
	C_Branch_Apo0  =C_Branch_Apo;
	C_Trunk_Apo0   =C_Trunk_Apo;
	C_Root_Apo01   =C_Root_Apo1;
	C_Root_Apo02   =C_Root_Apo2;
	C_Root_Apo03   =C_Root_Apo3;
	Turgor_lp= Pi0_Leaf_Symp*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp+Epsilon_Leaf_Symp);
	
	
	if (GRAVITY==1)
	{
		Pg_Leaf  = -9.81*(Length_Branch+Length_Trunk)/1000;
		Pg_trunk = -9.81*Length_Trunk/1000;                     // gravimetric pressure drop at the top leaves, trunk and branches in MPa
		Pg_branch= -9.81*(Length_Branch+Length_Trunk)/1000;
	}
	else
	{
		Pg_Leaf=0;
		Pg_trunk =0;
		Pg_branch=0;
	}
	P_Leaf_Evap=Pg_branch;        // the water potential at the site of evaporation in the leaf
	P_Leaf_Symp=Pg_branch;         // leaf symplasmic water potential
	P_Leaf_Apo=Pg_branch;
	P_Axil_Apo=Pg_branch;
	P_Axil_Symp=Pg_branch;
	P_Branch_Apo=Pg_branch;
	P_Branch_Symp=Pg_branch;
	P_Trunk_Apo=Pg_trunk;
	P_Trunk_Symp=Pg_trunk;
	P_Root_Apo1=0;
	P_Root_Symp1=0;
	P_Root_Apo2=0;
	P_Root_Symp2=0;
	P_Root_Apo3=0;
	P_Root_Symp3=0;
	
	if(GRAVITY==0 || GRAVITY==3)
	{
		Q_Leaf_Evap    = Q_Leaf_Evap0;
		Q_Leaf_Symp    = Q_Leaf_Symp0;
		Q_Leaf_Apo     = Q_Leaf_Apo0;
		Q_Leaf_Apo1    = Q_Leaf_Apo0;
		Q_Axil_Symp    = Q_Axil_Symp0;
		Q_Axil_Apo     = Q_Axil_Apo0;
		Q_Axil_Apo1    = Q_Axil_Apo0;
		Q_Branch_Symp  = Q_Branch_Symp0;
		Q_Branch_Apo   = Q_Branch_Apo0;
		Q_Branch_Apo1  = Q_Branch_Apo0;
		Q_Trunk_Symp   = Q_Trunk_Symp0;
		Q_Trunk_Apo    = Q_Trunk_Apo0;
		Q_Trunk_Apo1   = Q_Trunk_Apo0;
	}
	
	else // compute steady-state Q at Pg with C (apo) or PV curves (Symp); assume potential is above turgor loss point
	{
		long double discri, tlp;
		Q_Leaf_Evap    = Q_Leaf_Evap0   + C_Leaf_Evap   * P_Leaf_Evap ;
		Q_Leaf_Apo     = Q_Leaf_Apo0    + C_Leaf_Apo    * P_Leaf_Apo ;
		Q_Leaf_Apo1    = Q_Leaf_Apo ;
		Q_Axil_Apo     = Q_Axil_Apo0   + C_Axil_Apo   * P_Axil_Apo  ;
		Q_Axil_Apo1    = Q_Axil_Apo;
		Q_Branch_Apo   = Q_Branch_Apo0  + C_Branch_Apo  * P_Branch_Apo ;
		Q_Branch_Apo1  = Q_Branch_Apo;
		Q_Trunk_Apo    = Q_Trunk_Apo0   + C_Trunk_Apo   * P_Trunk_Apo;
		Q_Trunk_Apo1   = Q_Trunk_Apo;
		
		// LEAF
		tlp=Pi0_Leaf_Symp*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp+Epsilon_Leaf_Symp);
		if (P_Leaf_Symp>tlp)
		{
			discri=pow(Pi0_Leaf_Symp*Osmotic_TLeaf + Epsilon_Leaf_Symp + P_Leaf_Symp,2)-4*Pi0_Leaf_Symp*Osmotic_TLeaf*Epsilon_Leaf_Symp;
			Q_Leaf_Symp = ((Pi0_Leaf_Symp*Osmotic_TLeaf + Epsilon_Leaf_Symp + P_Leaf_Symp) + pow(discri,0.5))/(2*Epsilon_Leaf_Symp/Q_Leaf_Symp0);
		}
		else   Q_Leaf_Symp = Q_Leaf_Symp0*Pi0_Leaf_Symp*Osmotic_TLeaf/P_Leaf_Symp;
		
		// BUD Flower
		if (Type_Axil)
		{
			tlp=Pi0_Axil_Symp*Epsilon_Axil_Symp/(Pi0_Axil_Symp+Epsilon_Axil_Symp);
			if (P_Axil_Symp>tlp)
			{
				discri=pow(Pi0_Axil_Symp*Osmotic_TLeaf + Epsilon_Axil_Symp + P_Axil_Symp,2)-4*Pi0_Axil_Symp*Osmotic_TLeaf*Epsilon_Axil_Symp;
				Q_Axil_Symp = ((Pi0_Axil_Symp*Osmotic_TLeaf + Epsilon_Axil_Symp + P_Axil_Symp) + pow(discri,0.5))/(2*Epsilon_Axil_Symp/Q_Axil_Symp0);
			}
			else Q_Axil_Symp = Q_Axil_Symp0*Pi0_Axil_Symp*Osmotic_TLeaf/P_Axil_Symp;
			if (Type_Axil>=2)  // petiole
			{
				tlp=Pi0_Axil_Symp*Epsilon_Axil_Symp/(Pi0_Axil_Symp+Epsilon_Axil_Symp);
				if (P_Petiole_Symp>tlp)
				{
					discri=pow(Pi0_Axil_Symp*Osmotic_TLeaf + Epsilon_Axil_Symp + P_Petiole_Symp,2)-4*Pi0_Axil_Symp*Osmotic_TLeaf*Epsilon_Axil_Symp;
					Q_Petiole_Symp = ((Pi0_Axil_Symp*Osmotic_TLeaf + Epsilon_Axil_Symp + P_Petiole_Symp) + pow(discri,0.5))/(2*Epsilon_Axil_Symp/Q_Petiole_Symp0);
				}
				else Q_Petiole_Symp = Q_Petiole_Symp0*Pi0_Axil_Symp*Osmotic_TLeaf/P_Petiole_Symp;
			}
		}
		// BRANCH
		tlp=Pi0_Branch_Symp*Epsilon_Branch_Symp/(Pi0_Branch_Symp+Epsilon_Branch_Symp);
		if (P_Branch_Symp>tlp)
		{
			discri=pow(Pi0_Branch_Symp*Osmotic_TAir + Epsilon_Branch_Symp + P_Branch_Symp,2)-4*Pi0_Branch_Symp*Osmotic_TAir*Epsilon_Branch_Symp;
			Q_Branch_Symp = ((Pi0_Branch_Symp*Osmotic_TAir + Epsilon_Branch_Symp + P_Branch_Symp) + pow(discri,0.5))/(2*Epsilon_Branch_Symp/Q_Branch_Symp0);
		}
		else Q_Branch_Symp = Q_Branch_Symp0*Pi0_Branch_Symp*Osmotic_TAir/P_Branch_Symp;
		
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
	
}

void default_para(void) // create a init file with default values when init file is missing; 
{
	FILE *new_ini;
	
	printf("SurEau_ini.txt not found, using default values\n");
	new_ini = fopen("sureau_ini.txt","w");
	fprintf(new_ini,"1	0.001	0.9	0	11	1	0	0	0	0	0	99.9	0	0.01	60.000	1	0	1	1	1	1	1	1	0	1	0.80	-1.50	-1.50	-2.00	0	45.000	15.000	30.000	30.000	80.000	1500.000	1	0.000	3.000	3.000	15	0.28	0.1	0.0005	2	5.00	0.5	0	0	0.7	75	1	20	2	1	0	0.9	10	0	192	0	90	0	0.5	50	45	0	10.5	0	123	140	285	325	0	-0.5	200	100	0.250	45.849	110.000	180	2.000	-0.370	0.120	41.140	27350.000	-1.815	6.000	2.500E+00	3.000E-03	0.400	0.200	700.000	10	0.086	0.400	0.200	0.331	1.114E+03	1.000E-03	0.400	0.200	1	2.4256E+00	4.8512E+00	6.3808E+00	1.2762E+01	1.9759E+00	3.9517E+00	5.7839E+00	2.6912E+00	1.8043E+01	3.5000E+00	3.8358E+03	1.4713E-03	8.2538E+01	3.9387E+02	3.5680E+03	170	20.000	0.006	1	25.000	17.000	4.00	37.500	1.200	4.800	2.000	1.000	1.000	3.00E+01	0.0000E+00	9.0000E-01	5.0000E-05	5.0000E-06	2.0000E-05	5.0000E-04	1.1195	1.3229	2.5000	0.4000	0.1000	0.0000	10.0000	10.0000	10.0000	10.0000	-2.1000	-2.1000	-2.1000	-2.1000	3.0000	5.00E+03	5.00E+04	2.50E+04	2.0000	1.0000	1.0000	2.0000	-3.4	-3.4	-3.4	-3.4	60	60	60	60	0	0.0000E+00	2.0000E+00	5.0000E-03	2.0000E-02	1.000E+01	-2.100	-3.400	60.000	1.000	1.0000E+00	1.0000E+00	5.000E-02	0.001	0.010	0.002	0.500	0.000	0.500	0.000	-1.000	0.000	0.000	0.000	0.000	0.00E+00	0.000	0	0	0	12	0");
	fclose(new_ini);
}

void default_para_file(void)
{
	FILE *new_para;
	
	printf("SurEau_para.txt not found, using default values\n");
	new_para = fopen("sureau_para.txt","w");
	fprintf(new_para,"1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	189	190	191	192	193	194	195	196	197	198	199	200");
	fclose(new_para);
}

void Graph_ps(void) // Print graphs in postcript files
{
	FILE *File[250], *transient, *GRAPH;
	char FileName[]="                        ";
	long i,j,k;
	long double X0=20,Y0=150,DATAX0=100,DATAY0[250],XMAX=216,YMAX=279,Xscale,Yscale;
	long double yearqq;
	char LABEL[1000],c;
	long double DATAmin[250],DATAmax[250];
	long double color=0,R=0,G=0,B=0.4;
	long double Ncurves=1;
	long double NbDATA=0;
	
	YMAX=279;
	Xscale=0.8*XMAX;
	Yscale=Xscale/2; // a rectangular graph
	
	//Look for min max values in transient_out.csv
	transient = fopen(filename_TRANS,"r+");
	fscanf(transient, "%[^\n]", LABEL); //first text line
	fscanf(transient,"%Le %c",&yearqq,&c); // year
	for (i=0;i<250;i++)  if (File_out[i])
	{
		fscanf(transient, "%Le%c", &DATA[i],&c);
		DATAmin[i]=DATA[i];
		DATAmax[i]=DATA[i];
	}
	NbDATA+=1;
	while   (!feof(transient))
	{
		fscanf(transient,"%Le %c",&yearqq,&c);
		for (i=0;i<250;i++)  if (File_out[i])
		{
			fscanf(transient, "%Le%c", &DATA[i],&c);
			if (DATA[i]<DATAmin[i]) DATAmin[i]=DATA[i];
			if (DATA[i]>DATAmax[i]) DATAmax[i]=DATA[i];
		}
		if (DATA[0]< DATAX0) Ncurves++;  // a new simulation
		DATAX0 = DATA[0];
		NbDATA+=1;
	}
	fclose(transient);
	for (i=0;i<250;i++)
		if (File_out[i])
			if (DATAmin[i]==DATAmax[i])
			{
				DATAmin[i]*=0.9;
				DATAmax[i]*=1.1;
			}
	if (PRINT_GRAPH==2) //continuous
	{
		DATAmin[0]=0;
		DATAmax[0]=NbDATA;
	}
	
	//create a different .ps file for each variable in transient_out.csv; assume YEAR and Days are first two variables
	//first, create all .ps files with headers and axis
	for (i=0;i<250;i++)  if (File_out[i])
	{
		strcpy (FileName,"Graph_");
		strcat (FileName,Label[i]);
		strcat (FileName,".ps") ;
		File[i] = fopen(FileName, "w+" );
		
		fputs("%!PS-Adobe-2.0\n%%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File[i]);
		fputs("%%DocumentFonts: Courier\n%%EndProlog\n%%Page: 1 1\n",File[i]);
		fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File[i]);
		fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File[i]);
		
		//draw axis
		fputs("\n%Tracé du graphe \ngsave \n  .3 M setlinewidth\n  0 0 0 setrgbcolor\n newpath \n",File[i]);  //black
		fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",File[i]);
		fputs("\n/Times-Italic findfont\n  5 M scalefont\n  setfont\n\n",File[i]);
		fprintf(File[i],"\n%Lf M %Lf M moveto \n (%s) show\n ",X0-7+Xscale/2,Y0+10+Yscale,Label[i]);// Graph Title: variable
		fprintf(File[i],"\n%Lf M %Lf M moveto \n (%s) show\n ",X0-7+Xscale/2,Y0-10,Label[0]);       // x axis Title: Days
		fputs("\n/Times-Italic findfont\n  2 M scalefont\n  setfont\n\n",File[i]);
		fprintf(File[i],"\n150 M 5 M moveto \n (SurEau.c  ver.%s (C) Herve Cochard, INRAE) show\n ",version);
		fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",File[i]);
		if (DATAmin[i]<0 && DATAmax[i]>0)                                                           //draw Y=0 line
		{
			fputs("\n  .1 M setlinewidth\n  0 0 0 setrgbcolor\n",File[i]);
			fprintf(File[i],"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0+(0-DATAmin[i])/(DATAmax[i]-DATAmin[i])*Yscale,X0+Xscale,Y0+(0-DATAmin[i])/(DATAmax[i]-DATAmin[i])*Yscale);
			fprintf(File[i],"\n%Lf M %Lf M moveto \n (0) show\n",X0-3, Y0+(0-DATAmin[i])/(DATAmax[i]-DATAmin[i])*Yscale);
		}
		fputs("\n  .2 M setlinewidth\n  0 0 0 setrgbcolor\n",File[i]);
		fprintf(File[i],"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0+Xscale,Y0);  //x axis
		fprintf(File[i],"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0,Y0+Yscale);  //y axis
		for (j=0;j<=10;j++)
		{
			fprintf(File[i],"%Lf M %Lf M %Lf M %Lf M T\n",X0+j*Xscale/10,Y0,X0+j*Xscale/10,Y0-1); //x axis ticks
			fprintf(File[i],"%Lf M %Lf M moveto \n (%.1Lf) show\n ",X0-4+j*Xscale/10,Y0-5,DATAmin[0]+j* (DATAmax[0]-DATAmin[0])/10);
			fprintf(File[i],"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0+j*Yscale/10,X0-1,Y0+j*Yscale/10); //Y axis ticks
			fprintf(File[i],"%Lf M %Lf M moveto \n (%.2Le) show\n ",X0-18,Y0-1+j*Yscale/10,DATAmin[i]+j* (DATAmax[i]-DATAmin[i])/10);
		}
		fprintf(File[i],"  stroke \n.2 M setlinewidth\n%Lf %Lf %Lf setrgbcolor \n newpath\n",R,G,B);
		fprintf(File[i],"%Lf M %Lf M moveto \n (curve 1) show\n ",X0+Xscale,Yscale+70+color*40);
	}
	
	transient = fopen(filename_TRANS,"r+");
	//first line with text labels
	fscanf(transient, "%[^\n]", LABEL);
	NbDATA=0;
	//first line of data
	fscanf(transient,"%Le %c",&yearqq,&c); //Year
	fscanf(transient, "%Le %c", &DATA[0],&c); //Days
	if (PRINT_GRAPH==2) DATAX0=0;
	else                DATAX0=DATA[0];
	for (i=1;i<250;i++)  if (File_out[i]) //the remaining variables
	{
		fscanf(transient, "%Le %c", &DATA[i],&c);
		DATAY0[i]=DATA[i];
	}
	k=0;
	//all the other lines
	while (!feof(transient))
	{
		
		NbDATA+=1;
		if (PRINT_GRAPH==2) DATAX0=NbDATA-1;
		else                DATAX0=DATA[0];
		fscanf(transient,"%Le %c",&yearqq,&c); //Year
		fscanf(transient, "%Le %c", &DATA[0],&c); //Days
		if (PRINT_GRAPH==2) DATA[0]=NbDATA;
		if (DATA[0]< DATAX0)  // a new simulation
		{
			k++;
			if (Ncurves-1)  color+=(1/(Ncurves-1));
			else            color=0;
			if (color>1) color=1;
			if (color<0) color=0;
			if (color<=0.12){R=0;G=0;B=0.4+color*5.0;}
			else if (color<=0.32){R=0;G=(color-0.12)*5;B=1.0;}
			else if (color<=0.52){R=0;G=1.0;B=1.0-(color-0.32)*5;}
			else if (color<=0.72){R=(color-0.52)*5;G=1.0;B=0;}
			else if (color<=0.92){R=1.0;G=1.0-(color-0.72)*5;B=0.0;}
			else if (color<=1){R=1.0-(color-0.92)*5;G=0.0;B=0.0;}
			
			for (i=1;i<250;i++) if (File_out[i])
			{
				fscanf(transient, "%Le %c", &DATA[i],&c);
				fprintf(File[i],"  stroke \n.2 M setlinewidth\n%Lf %Lf %Lf setrgbcolor \n newpath\n",R,G,B);
				fprintf(File[i],"%Lf M %Lf M moveto \n (curve %ld) show\n ",X0+Xscale,Yscale+70+3*k,k+1);
			}
			
		}
		else for (i=1;i<250;i++) if (File_out[i])
		{
			DATAY0[i]=DATA[i];
			fscanf(transient, "%Le %c", &DATA[i],&c);
			fprintf(File[i],"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i]-DATAmin[i])/(DATAmax[i]-DATAmin[i])*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i]-DATAmin[i])/(DATAmax[i]-DATAmin[i])*Yscale);
		}
	}
	for (i=0;i<250;i++)  if (File_out[i])
	{
		fputs("\nstroke\nshowpage\n%%Trailer\n",File[i]);
		fclose(File[i]);
	}
	
	// all PLC data in one graph
	GRAPH= fopen("GRAPH_PLC.ps", "w+");
	fputs("%!PS-Adobe-2.0\n%%SurEau model Creator: Herve Cochard INRAE-PIAF\n",GRAPH);
	fputs("%%DocumentFonts: Courier\n%%EndProlog\n%%Page: 1 1\n",GRAPH);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",GRAPH);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",GRAPH);
	fputs("\n%Tracé du graphe \ngsave \n  .3 M setlinewidth\n  0 0 0 setrgbcolor\n newpath \n",GRAPH);  //black
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n/Times-Italic findfont\n  5 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (PLC) show\n ",X0-7+Xscale/2,Y0+10+Yscale);// Graph Title: variable
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (%s) show\n ",X0-5+Xscale/2,Y0-10,Label[0]);       // x axis Title: Days
	fputs("\n/Times-Italic findfont\n  2 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n150 M 5 M moveto \n (SurEau.c  ver.%s (C) Herve Cochard, INRAE) show\n ",version);
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n  .2 M setlinewidth\n  0 0 0 setrgbcolor\n",GRAPH);
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0+Xscale,Y0);  //x axis
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0,Y0+Yscale);  //y axis
	for (j=0;j<=10;j++)
	{
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0+j*Xscale/10,Y0,X0+j*Xscale/10,Y0-1); //x axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.1Lf) show\n ",X0-4+j*Xscale/10,Y0-5,DATAmin[0]+j* (DATAmax[0]-DATAmin[0])/10);
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0+j*Yscale/10,X0-1,Y0+j*Yscale/10); //Y axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.2ld) show\n ",X0-7,Y0-1+j*Yscale/10,j*10);
	}
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0 0 0.4 setrgbcolor \n newpath\n");
	if (File_out[78]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Leaf) show\n ",X0+Xscale-5,Yscale+22+3*6);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 0.114286 1.000000 setrgbcolor \n newpath\n");
	if (File_out[79])fprintf(GRAPH,"%Lf M %Lf M moveto \n (Branch) show\n ",X0+Xscale-5,Yscale+22+3*5);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 0.828571 1.000000 setrgbcolor \n newpath\n");
	if (File_out[80]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Axil) show\n ",X0+Xscale-5,Yscale+22+3*4);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 0.4000000 0.457143 setrgbcolor \n newpath\n");
	if (File_out[81]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Trunk) show\n ",X0+Xscale-5,Yscale+22+3*3);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.257143 1.000000 0.000000 setrgbcolor \n newpath\n");
	if (File_out[82]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root1) show\n ",X0+Xscale-5,Yscale+22+3*2);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.971429 0.800000 0.000000 setrgbcolor \n newpath\n");
	if (File_out[83]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root2) show\n ",X0+Xscale-5,Yscale+22+3*1);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6 0 0 setrgbcolor \n newpath\n");
	if (File_out[84]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root3) show\n ",X0+Xscale-5,Yscale+22+3*0);
	
	
	transient = fopen(filename_TRANS,"r+");
	//first line with text labels
	fscanf(transient, "%[^\n]", LABEL);
	NbDATA=0;
	//first line of data
	fscanf(transient,"%Le %c",&yearqq,&c); //Year
	fscanf(transient, "%Le %c", &DATA[0],&c); //Days
	DATAX0=DATA[0];
	for (i=1;i<250;i++)  if (File_out[i]) //the remaining variables
	{
		fscanf(transient, "%Le %c", &DATA[i],&c);
		DATAY0[i]=DATA[i];
	}
	k=0;
	//all the other lines
	while (!feof(transient))
	{
		NbDATA+=1;
		DATAX0=DATA[0];
		fscanf(transient,"%Le %c",&yearqq,&c); //Year
		fscanf(transient, "%Le %c", &DATA[0],&c); //Days
		for (i=1;i<250;i++) if (File_out[i])
		{
			DATAY0[i]=DATA[i];
			fscanf(transient, "%Le %c", &DATA[i],&c);
			if (i==78)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0 0 0.4 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==79)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 0.114286 1.000000 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==80)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 0.828571 1.000000 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==81)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.000000 1.000000 0.457143 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==82)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.257143 1.000000 0.000000 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==83)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.971429 1.000000 0.000000 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
			if (i==84)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6 0 0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(100)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(100)*Yscale);
			}
		}
	}
	
		fputs("\nstroke\nshowpage\n%%Trailer\n",GRAPH);
		fclose(GRAPH);
	
	// all RWC data in one graph
	GRAPH= fopen("GRAPH_RWC.ps", "w+");
	fputs("%!PS-Adobe-2.0\n%%SurEau model Creator: Herve Cochard INRAE-PIAF\n",GRAPH);
	fputs("%%DocumentFonts: Courier\n%%EndProlog\n%%Page: 1 1\n",GRAPH);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",GRAPH);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",GRAPH);
	fputs("\n%Tracé du graphe \ngsave \n  .3 M setlinewidth\n  0 0 0 setrgbcolor\n newpath \n",GRAPH);  //black
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n/Times-Italic findfont\n  5 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (RWC) show\n ",X0-7+Xscale/2,Y0+10+Yscale);// Graph Title: variable
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (%s) show\n ",X0-5+Xscale/2,Y0-10,Label[0]);       // x axis Title: Days
	fputs("\n/Times-Italic findfont\n  2 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n150 M 5 M moveto \n (SurEau.c  ver.%s (C) Herve Cochard, INRAE) show\n ",version);
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n  .2 M setlinewidth\n  0 0 0 setrgbcolor\n",GRAPH);
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0+Xscale,Y0);  //x axis
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0,Y0+Yscale);  //y axis
	for (j=0;j<=10;j++)
	{
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0+j*Xscale/10,Y0,X0+j*Xscale/10,Y0-1); //x axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.1Lf) show\n ",X0-4+j*Xscale/10,Y0-5,DATAmin[0]+j* (DATAmax[0]-DATAmin[0])/10);
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0+j*Yscale/10,X0-1,Y0+j*Yscale/10); //Y axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.2ld) show\n ",X0-7,Y0-1+j*Yscale/10,j*10);
	}
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0 0 0.4 setrgbcolor \n newpath\n");
	if (File_out[193]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Leaf) show\n ",X0+Xscale-5,Yscale+22+3*16);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    0.7 setrgbcolor \n newpath\n");
	if (File_out[194])fprintf(GRAPH,"%Lf M %Lf M moveto \n (Leaf_s) show\n ",X0+Xscale-5,Yscale+22+3*15);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    1 setrgbcolor \n newpath\n");
	if (File_out[195]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Leaf_a) show\n ",X0+Xscale-5,Yscale+22+3*14);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.35    1 setrgbcolor \n newpath\n");
	if (File_out[196]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Branch) show\n ",X0+Xscale-5,Yscale+22+3*13);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.65    1 setrgbcolor \n newpath\n");
	if (File_out[197]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Branch_s) show\n ",X0+Xscale-5,Yscale+22+3*12);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    1 setrgbcolor \n newpath\n");
	if (File_out[198]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Branch_a) show\n ",X0+Xscale-5,Yscale+22+3*11);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.7 setrgbcolor \n newpath\n");
	if (File_out[199]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Trunk) show\n ",X0+Xscale-5,Yscale+22+3*10);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.35 setrgbcolor \n newpath\n");
	if (File_out[200]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Trunk_s) show\n ",X0+Xscale-5,Yscale+22+3*9);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.05 setrgbcolor \n newpath\n");
	if (File_out[201]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Trunk_a) show\n ",X0+Xscale-5,Yscale+22+3*8);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.4    1    0 setrgbcolor \n newpath\n");
	if (File_out[202]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root) show\n ",X0+Xscale-5,Yscale+22+3*7);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.75    1    0 setrgbcolor \n newpath\n");
	if (File_out[203]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root_s) show\n ",X0+Xscale-5,Yscale+22+3*6);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.8    0 setrgbcolor \n newpath\n");
	if (File_out[204]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root_a) show\n ",X0+Xscale-5,Yscale+22+3*5);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.45    0 setrgbcolor \n newpath\n");
	if (File_out[205]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Plant) show\n ",X0+Xscale-5,Yscale+22+3*4);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.15    0 setrgbcolor \n newpath\n");
	if (File_out[206]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Plant_s) show\n ",X0+Xscale-5,Yscale+22+3*3);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.95    0    0 setrgbcolor \n newpath\n");
	if (File_out[207]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Plant_a) show\n ",X0+Xscale-5,Yscale+22+3*2);
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6    0    0 setrgbcolor \n newpath\n");
	if (File_out[186]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Soil) show\n ",X0+Xscale-5,Yscale+22+3*1);
	
	transient = fopen(filename_TRANS,"r+");
	//first line with text labels
	fscanf(transient, "%[^\n]", LABEL);
	NbDATA=0;
	//first line of data
	fscanf(transient,"%Le %c",&yearqq,&c); //Year
	fscanf(transient, "%Le %c", &DATA[0],&c); //Days
	DATAX0=DATA[0];
	for (i=1;i<250;i++)  if (File_out[i]) //the remaining variables
	{
		fscanf(transient, "%Le %c", &DATA[i],&c);
		DATAY0[i]=DATA[i];
	}
	k=0;
	//all the other lines
	while (!feof(transient))
	{
		NbDATA+=1;
		DATAX0=DATA[0];
		fscanf(transient,"%Le %c",&yearqq,&c); //Year
		fscanf(transient, "%Le %c", &DATA[0],&c); //Days
		for (i=1;i<250;i++) if (File_out[i])
		{
			DATAY0[i]=DATA[i];
			fscanf(transient, "%Le %c", &DATA[i],&c);
			if (i==193)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0 0 0.4 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==194)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    0.7 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==195)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    1 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==196)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.35    1 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==197)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.65    1 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==198)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.35 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==199)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.7 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==200)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.05 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==201)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.4    1    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==202)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.75    1    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==203)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.8    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==204)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.45    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==205)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.15    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==206)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.95    0    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==207)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6    0    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			if (i==186)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.4    0    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
		}
	}
	
	
	// all RWC_s data in one graph
	GRAPH= fopen("GRAPH_RWC_symp.ps", "w+");
	fputs("%!PS-Adobe-2.0\n%%SurEau model Creator: Herve Cochard INRAE-PIAF\n",GRAPH);
	fputs("%%DocumentFonts: Courier\n%%EndProlog\n%%Page: 1 1\n",GRAPH);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",GRAPH);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",GRAPH);
	fputs("\n%Tracé du graphe \ngsave \n  .3 M setlinewidth\n  0 0 0 setrgbcolor\n newpath \n",GRAPH);  //black
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n/Times-Italic findfont\n  5 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (RWC_Symp) show\n ",X0-7+Xscale/2,Y0+10+Yscale);// Graph Title: variable
	fprintf(GRAPH,"\n%Lf M %Lf M moveto \n (%s) show\n ",X0-5+Xscale/2,Y0-10,Label[0]);       // x axis Title: Days
	fputs("\n/Times-Italic findfont\n  2 M scalefont\n  setfont\n\n",GRAPH);
	fprintf(GRAPH,"\n150 M 5 M moveto \n (SurEau.c  ver.%s (C) Herve Cochard, INRAE) show\n ",version);
	fputs("\n/Courier findfont\n  3 M scalefont\n  setfont\n\n",GRAPH);
	fputs("\n  .2 M setlinewidth\n  0 0 0 setrgbcolor\n",GRAPH);
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0+Xscale,Y0);  //x axis
	fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0,X0,Y0+Yscale);  //y axis
	for (j=0;j<=10;j++)
	{
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0+j*Xscale/10,Y0,X0+j*Xscale/10,Y0-1); //x axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.1Lf) show\n ",X0-4+j*Xscale/10,Y0-5,DATAmin[0]+j* (DATAmax[0]-DATAmin[0])/10);
		fprintf(GRAPH,"%Lf M %Lf M %Lf M %Lf M T\n",X0,Y0+j*Yscale/10,X0-1,Y0+j*Yscale/10); //Y axis ticks
		fprintf(GRAPH,"%Lf M %Lf M moveto \n (%.2ld) show\n ",X0-7,Y0-1+j*Yscale/10,j*10);
	}
   
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    0.7 setrgbcolor \n newpath\n");
	if (File_out[194])fprintf(GRAPH,"%Lf M %Lf M moveto \n (Leaf_s) show\n ",X0+Xscale-5,Yscale+22+3*15);
	
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.65    1 setrgbcolor \n newpath\n");
	if (File_out[197]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Branch_s) show\n ",X0+Xscale-5,Yscale+22+3*12);
	
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.35 setrgbcolor \n newpath\n");
	if (File_out[200]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Trunk_s) show\n ",X0+Xscale-5,Yscale+22+3*9);
	
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.75    1    0 setrgbcolor \n newpath\n");
	if (File_out[203]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Root_s) show\n ",X0+Xscale-5,Yscale+22+3*6);
	
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.15    0 setrgbcolor \n newpath\n");
	if (File_out[206]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Plant_s) show\n ",X0+Xscale-5,Yscale+22+3*3);
	
	fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6    0    0 setrgbcolor \n newpath\n");
	if (File_out[186]) fprintf(GRAPH,"%Lf M %Lf M moveto \n (Soil) show\n ",X0+Xscale-5,Yscale+22+3*1);
	
	transient = fopen(filename_TRANS,"r+");
	//first line with text labels
	fscanf(transient, "%[^\n]", LABEL);
	NbDATA=0;
	//first line of data
	fscanf(transient,"%Le %c",&yearqq,&c); //Year
	fscanf(transient, "%Le %c", &DATA[0],&c); //Days
	DATAX0=DATA[0];
	for (i=1;i<250;i++)  if (File_out[i]) //the remaining variables
	{
		fscanf(transient, "%Le %c", &DATA[i],&c);
		DATAY0[i]=DATA[i];
	}
	k=0;
	//all the other lines
	while (!feof(transient))
	{
		NbDATA+=1;
		DATAX0=DATA[0];
		fscanf(transient,"%Le %c",&yearqq,&c); //Year
		fscanf(transient, "%Le %c", &DATA[0],&c); //Days
		for (i=1;i<250;i++) if (File_out[i])
		{
			DATAY0[i]=DATA[i];
			fscanf(transient, "%Le %c", &DATA[i],&c);
		   
			if (i==194)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0    0.7 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
		   
		   
			if (i==197)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    0.65    1 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
	   
			if (i==200)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0    1    0.35 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			
			if (i==203)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.75    1    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
			
			if (i==206)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 1    0.15    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
		   
			if (i==186)
			{
				fprintf(GRAPH,"  stroke \n.2 M setlinewidth\n 0.6    0    0 setrgbcolor \n newpath\n");
				fprintf(GRAPH,"%Le M %Le M %Le M %Le M T\n",X0+(DATAX0-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale,Y0+(DATAY0[i])/(1)*Yscale,X0+(DATA[0]-DATAmin[0])/(DATAmax[0]-DATAmin[0])*Xscale, Y0+(DATA[i])/(1)*Yscale);
			}
		}
	}
	
	fputs("\nstroke\nshowpage\n%%Trailer\n",GRAPH);
	fclose(GRAPH);
}

void Leaf_C_budget(long double dt_long)
{
	long double export0=0.0;             //0.74477;   //exportation of 0.833 µmol of glucose /m2/s
	long double Reserve0=50000;          // maximum reserve in µmol Assuming 10% mg/g de NSC in leave and a LMA of 90g/m2 then 9g NSC/m2 or Reserve0=50 mmol equ glucose with MW=180
	
	// Reserve in µmol of glucose
	
	Export=-Turgor_Leaf_Symp/Pi0_Leaf_Symp*export0;
	Export_tot+=Export*dt_long*dt;
	//if (g_s) Reserve+=(A_gross*dt_long*dt/6-export*dt_long*dt);   //in µmol of glucose; 6 CO2 for 1 glucose
	//else
	Reserve+=(Rd25*exp(18.72-46390/(8.314*(T_Leaf+273.15)))*dt_long*dt/6);  // if the carbon loss is due only to the photosynthesis respiration if CO2 from Glucose is recycled
	if (Reserve>Reserve0) Reserve=Reserve0;
}

long double Respiration(void)  //under development;
{
	long double  Tref=25,  RTref, Resp;
	// long double temp, Q10=2, Rm_25=100;
	
	// Q_Wood= (Q_Branch_Symp + Q_Branch_Apo + Q_Trunk_Symp + Q_Trunk_Apo + Q_Root_Symp1 + Q_Root_Apo1 + Q_Root_Endo1+ Q_Root_Symp2 + Q_Root_Apo2 + Q_Root_Endo2+ Q_Root_Symp3 + Q_Root_Apo3 + Q_Root_Endo3)*18/1000/1000/1000;
	// temp=T_air;
	// Rm= Rm_25*pow(Q10,(temp-Tb)/10); //Basic Q10 respiration rate in µmol CO2/s/m3
	// Rg=0;
	
	// Only leaves Rm is included so far; from Heskel et al PNAS 2016
	RTref=exp(a_Res+0.1012*Tref-0.0005*Tref*Tref);
	Resp=RTref*exp(0.1012*(T_Leaf-Tref)-0.0005*(T_Leaf*T_Leaf-Tref*Tref));
	return Resp;
}

long double Net_Photosynthesis_C3 (void)  // From Tyree spreadsheet Farquhar et al model
{
	// Cp CO2 Compensation point [umol.s-1.m-2]
	long double temp, VcMaxT, VjMaxT , Rd, Kc,Ko, KmT, Jv, gg,Cp=3.7;
	long double a,b,c,delta;
	long double PAR_apparent=PAR*cos(leaf_angle*3.1416/180);
	long double A_LEAF;
	
	temp=T_Leaf;
	gg=g_s/1.6; // only the stomatal conductance involved here
	
	if (gg<0) gg=0;
	   //  Temperature effects on parameters
	//  Maximum velocity of carboxilation (Rubisco-limited) [umol.s-1.m-2]
	VcMaxT = VcMax * (1 + 0.0505 * (temp - 25) - 0.000248 * pow((temp - 25), 2) - 8.09E-05 * pow((temp - 25),3));
	if (VcMaxT<0) VcMaxT=0;
	
	//  Maximum velocity of RuBP regeneration [umol.s-1.m-2]
	VjMaxT = VjMax * (1 + 2.9741E-02 * (temp - 25) - 1.2413E-03 * pow((temp - 25),2) - 3.3085E-05 * pow((temp - 25),3));
	if (VjMaxT<0) VjMaxT=0;
	
	//  Apparent michaelis constant [Pa] in Rubisco-limited situation
	Kc = Kc25 * exp(79300 * 0.000404 * (temp - 25) / (temp + 273.15));
	Ko = Ko25 * exp(36400 * 0.000404 * (temp - 25) / (temp + 273.15));
	KmT = Kc * (1 + 20900 / Ko);
	Rd= Rd25*exp(18.72-46390/(8.314*(temp+273.15)));
	
	//  Velocity of RuBP regeneration [umol.s-1.m-2]
	if (VjMaxT) Jv = PAR_apparent * Qye / pow((1 + pow((Qye * PAR_apparent / VjMaxT),2)), 0.5);
	else        Jv=0;
	
	// Quadratic solution for Rubisco-limited rate of photosynthesis
	a= -1000;
	b= gg*(Ca+10*KmT) + 1000*(VcMaxT+ Rd);
	c= gg*VcMaxT*(10*Cp-Ca) - Rd*gg*(Ca+10*KmT);
	delta= b*b - 4*a*c;
	if (delta>0) A_net1= (-b+pow(delta,0.5))/(2*a); 
	 else A_net1=0;
	
	//  Quadratic solution for Light-limited rate of photosynthesis
	a= -4000;
	b= gg*(4*Ca+80*Cp) + 1000*Jv + 4000*Rd;
	c= gg*Jv*(10*Cp-Ca) - Rd*gg*(4*Ca+80*Cp);
	delta= b*b - 4*a*c;
	if (delta>=0) A_net2 = (-b+pow(delta,0.5))/(2*a);
	else          A_net2=0;
	if (delta>0) A_net2 = (-b+pow(delta,0.5))/(2*a); else A_net2=0;
	
	//  Net assimilation rate [umol.s-1.m-2]
	if (A_net1<A_net2)  A_LEAF=A_net1; else A_LEAF=A_net2;
	return A_LEAF;
}

long double Arrhenius(long double T2, long double Ea, long double Hd, long double DS)
{
	long double woo;
	long double Rgas = 8.314;  // forgotten in inital R code... appears to be the ideal gas law coefficient
	woo = exp(Ea * (T2 - 298) / (298 * Rgas * T2)) * (1 + exp((298 * DS - Hd) / (298 * Rgas))) / (1 + exp((T2 * DS - Hd) / (T2 * Rgas)));
	return woo;
}

long double Net_Photosynthesis_C4 (void)  //not tested yet
{
	long double Ci=500, VPMAX25=120, JMAX25=400, Vcmax=60, Vpr=80, Alpha2=0.0, gbs=3e-3;
	long double Vpmax=1, Jmax=1, low_gammastar=1.93e-4;
	long double Kc, Kp, Ko, K, Rd, Rm2;
	long double a_c, b_c, c_c, A_enzyme;
	long double Qp2, J, a_j, b_j, c_j, A_light, shape2;
	long double Vp, A_net3, Ac, Ad, Aj;
	long double A_LEAF, O2=210, x=0.4, THETA=0.7, Q10 = 2.3, RD0=1, RTEMP=25, TBELOW=0, DAYRESP=1,  FRM=0.5; //Q10F=2
	long double TK = T_Leaf + 273.15;  // init TK... change later
	long double res1, res2, res3;
	long double foo;
	/*  description and implementation of the A-Ci curve for C4 plants, based on von Caemmerer et al. (2000)
 R code written by Rhys Whitley (plantecophys R package)... converted to C by Sean Gleason */

// Temperature effects on Vcmax, Vpmax and Jmax (Massad et al. 2007)
// This function returns value between 0 and 1.
	// Michaelis-Menten coefficients for CO2 (Kc, mu mol mol-1) and O (Ko, mmol mol-1) and combined (K)
	Kc = 650 * pow(Q10, ((T_Leaf-25)/10));
	Kp = 80 * pow(Q10, ((T_Leaf-25)/10));
	Ko = 450 * pow(Q10, ((T_Leaf-25)/10));
	K = Kc*(1+O2/Ko);
	
	// T effects according to Massad et al. (2007)
	res1 = Arrhenius(TK, 67294, 144568, 472);
	res2 = Arrhenius(TK, 70373, 117910, 376);
	res3 = Arrhenius(TK, 77900, 191929, 627);
	
	Vcmax = Vcmax * res1;
	Vpmax = VPMAX25 * res2;
	Jmax = JMAX25 * res3;
	
	// Day leaf respiration, umol m-2 s-1
	if (T_Leaf > TBELOW)	Rd = RD0 * pow(Q10, ((T_Leaf-RTEMP)/10)) * DAYRESP;
	else                	Rd = 0.0;
	Rm2 =  FRM*Rd;
	
	// PEP carboxylation rate
	//  Vp = pmin(Ci*Vpmax/(Ci+Kp),Vpr) // origional R code
	foo = Ci * Vpmax / (Ci + Kp);
	if(foo < Vpr) Vp = foo;
	else (Vp = Vpr);
	
	
	// Quadratic solution for enzyme limited C4 assimilation
	a_c = 1 - (Alpha2 * Kc) / (0.047 * Ko);
	b_c = -((Vp - Rm + gbs * Ci) + (Vcmax - Rd) + gbs * K + Alpha2 * low_gammastar / 0.047 * (low_gammastar * Vcmax+Rd * Kc / Ko));
	c_c = (Vcmax - Rd) * (Vp - Rm2 + gbs * Ci) - (Vcmax * gbs * low_gammastar * O2 + Rd * gbs * K);
	A_enzyme = (-b_c - sqrt(pow(b_c, 2) - 4 * a_c * c_c)) / (2 * a_c);
	
	// # Non-rectangular hyperbola describing light effect on electron transport rate (J)
	Qp2 = PAR * (1 - 0.15) / 2;
	J = (1 / (2 * THETA)) * (Qp2 + Jmax - sqrt(pow((Qp2 + Jmax), 2) - 4 * THETA * Qp2 * Jmax));
	
	// Quadratic solution for light-limited C4 assimilation
	a_j = 1 - 7 * low_gammastar * Alpha2 / (3*0.047);
	b_j = -( (x * J / 2 - Rm2 + gbs * Ci) + ((1 - x) * J / 3 - Rd) + gbs * (7 * low_gammastar * O2 / 3) +Alpha2 * low_gammastar / 0.047 * ((1 - x) * J / 3 + Rd));
	c_j = ( (x * J / 2 - Rm2 + gbs * Ci) * ((1 - x) * J / 3 - Rd) - gbs * low_gammastar * O2 * ((1 - x) * J / 3 - 7 * Rd / 3));
	A_light = (-b_j - sqrt(pow(b_j, 2) - 4 * a_j * c_j)) / (2 * a_j);
	
	//  Actual assimilation rate
	// A_net3 <- pmin(A.enzyme,A.light) // origional R code
	if(A_enzyme < A_light)  A_net3 = A_enzyme;
	else                    A_net3 = A_light;
	
	Ac = A_enzyme;
	Aj = A_light;
	
	//  Hyperbolic minimum (Buckley), to avoid discontinuity at transition from Ac to Aj
	// can turn this off... just makes the curve look pretty
	shape2 = 0.999;
	Ad = (Ac + Aj - sqrt(pow((Ac+Aj), 2) - 4 * shape2 * Ac * Aj)) / (2 * shape2) - Rd;
	Ac = Ac - Rd;
	Aj = Aj - Rd;
	A_LEAF=Ad;
	return A_LEAF;
}

long double Declination (void)
{
	long double c1 = 0.398749068925246; // =Sin(23.5*pi/180), 23.5 = Earth declination
	long double c2  = 2 * 3.1416 / 365;
	long double c3  = 80; // date of spring
	long double x;
	
	x = c1 * sin((DOY - c3) * c2);
	return atan(x / pow((1 - x*x),0.5));;
}

long double Potential_PAR(long double timeOfDay)
{
	long double decl,pn,pz,hRad, se, sn, sz, alt ,azi , pfd ,dpfd;
   // long double timeOfDay=(T/3600-DOY*24);
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

long double Flux (long double TTT,long double gb,long double hh, long double SWRa,long double GsIn,long double RH,long double Tair,long double em_leaf,long double em_air)
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
	//... Leaf emissivity: em_leaf
	//... Air emissivity: em_air
	//... Output: Heat balance (W/m2)
	
	return SWRa+BB(Tair,em_air)-BB(TTT,em_leaf)-Transpi(TTT,gb,RH,GsIn)-Convec(TTT,hh,Tair);
}

long double brentq( long double xa, long double xb, long double xtol, long double rtol, long double gb, long double hh,long double SWRa,long double GsIn,long double RH,long double Tair,long double em_leaf,long double em_air, int maxiter)
{
	// brentq method taken from scipy.optimize library
	
	double xpre = xa, xcur = xb;
	double xblk = 0.0, fpre, fcur, fblk = 0.0, spre = 0.0, scur = 0.0, sbis, tol;
	double stry, dpre, dblk;
	int i;
	fpre = Flux(xpre,gb,hh,SWRa,GsIn,RH,Tair,em_leaf,em_air);
	fcur = Flux(xcur,gb,hh,SWRa,GsIn,RH,Tair,em_leaf,em_air);
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
		fcur = Flux(xcur,gb,hh,SWRa,GsIn,RH,Tair,em_leaf,em_air);
	}
	return xcur;
}

void Tleaf(void)          //Energy budget from Ecofiz_Tleaf_K2_v3 www.landflux.org/r
{
	long double SWR;   					// short-wave radiation    W m-2
	long double WS;    					// windspeed    m s-1
	long double Tair;   					// air temperature    oC
	long double RH;    					// relative humidity    %
	long double aSWR=0.5;  				// absorptance to SWR     %
	long double em_leaf=0.97;    			// emissivity    none
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
	long double rbl ;  					// leaf boundary-layer resistance    s m-1
	long double Delta_T  ; 				// leaf-to-air temperature difference    oC
	long double TLeaf, Tleaf_NonLinear;   // leaf temperature    oC
	int maxiter = 50;
	long double em_air;
	
	
	if (POTENTIAL_PAR) cloud_cover=PAR/POTENTIAL_PAR; else cloud_cover=0;
	if (cloud_cover>1)	cloud_cover=1;
	if (CLIMAT==5)		cloud_cover=1;
	d=Leaf_size;
	Tair=T_air;
	RH= RH_air;
	
	if (g_s+g_cuti) rst=1/(g_s+g_cuti)*1000*40;   // conversion fromm mmol/s/m2 to s/m 
	else     rst=9999.99;
	SWR=    PAR*0.5495;    // from µmol/m²/s to Watts/m²
	WS=     Wind;
	esat=   a*exp(b*Tair/(Tair+z)); //kPa
	ea=     esat*(RH/100);
	s=      esat*b*z/(pow((Tair+z),2));
	em_air =((1-0.84*cloud_cover)*1.31*pow(10*ea/(Tair+273.15),0.14285714)+0.84*cloud_cover);
	VPDx=   esat-ea;
	SWRabs= aSWR*cos(leaf_angle*3.1416/180)*SWR;
	LWRin=  em_air*SB*pow(Tair+273.15,4);  // for clear and cloudy sky
	LWRouti=em_leaf*SB*pow(Tair+273.15,4);
	Rni=    SWRabs+LWRin-LWRouti;
	rr=     p*Cp/(4*em_leaf*SB*pow(Tair+273.15,3));
	
	if (Leaf_size > 3 ) rbl=1/(1.5*gflat*(pow(WS,jflat)/pow(d/1000,(1-jflat))));     //  a flat leaf if > 3mm
	else                rbl=1/(1.5*gcyl*(pow(WS,jcyl)/pow(d/1000,(1-jcyl))));                      // a needle, formula for a cylinder
	g_bl=1/rbl*1000*40;     //leaf boundary layer conductance in mmol/s/m2
	rblr=1/(1/rbl+1/rr);
	ym=y*(rst/rblr);
	// compute Tleaf with linear approximation
	if (TLEAF==1)
	{
		Delta_T= (ym*Rni*rblr/(p*Cp)-VPDx)/(s+ym);
		TLeaf=   Tair+Delta_T;
		T_Leaf=  TLeaf;
	}
	// compute non linear Tleaf
	else if (TLEAF==2)
	{
		Tleaf_NonLinear =  brentq(Tair-100.0, Tair+100.0, 1e-16, 1e-16, 1/rbl,p*Cp/rbl,SWRabs,1/rst,RH,Tair,em_leaf,em_air,maxiter);
		T_Leaf= Tleaf_NonLinear;
	}
	else T_Leaf=T_air;
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
	else if ((IRRIGATE==3 || IRRIGATE==7)&& T>IRR_DOY_S*3600*24 && T<IRR_DOY_F*3600*24) //automatic daily irrigation of Daily_irr mm only between DOY_start and DOY_end
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}
	else if (IRRIGATE==5 && T>IRR_DOY_S*3600*24 && T<IRR_DOY_F*3600*24 && REW_t<RWC_Irr) //automatic daily irrigation of Daily_irr mm only between DOY_start and DOY_end when below RWC_irr
	{
		Q_Soil1+=   Daily_Irr *Surface_Soil*1000*1000/18; //add rain from the past interval to soil
		Irrigation+=Daily_Irr;
		fill_soil_layers();
	}

}

void Climat (long double dt_long)
{
	long double e_sat, e,  e_sat_air, e_air, e_air_sol1, e_air_sol2, e_air_sol3, slope, tangente,TTTT,day;
	if      (CO2_atm==0)  Ca=  400;
	else if (CO2_atm==1)                   // RCP 2.6
	{
		if (YEAR1<2000) 	Ca=  255.1911 + 282.7456/pow(1+exp(-(YEAR1-2058.0442)/11.3900),0.1781); 
		else 				Ca= -0.00000001516648304*pow(YEAR1,5) + 0.0001586265578*pow(YEAR1,4) - 0.663190082*pow(YEAR1,3) + 1385.434635*pow(YEAR1,2) - 1446186.605*YEAR1 + 603458228.5;
	}
	else if (CO2_atm==2)  Ca=  255.1911 + 282.7456/pow(1+exp(-(YEAR1-2058.0442)/11.3900),0.1781);                 // RCP 4.5
	else if (CO2_atm==3)  Ca=  282.8391 + 902.2988/pow(1+exp(-(YEAR1-2104.6814)/17.3348),0.3910);                 // RCP 8.5
	else                  Ca=  CO2_atm;
	

	if (CLIMAT==2 || CLIMAT==3 || CLIMAT==4) // except for interpolated values
	{
		tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOY-80)/365*2*3.1416)));
		TTTT=(T/3600-DOY*24); // time of day
		if (tangente<-1)     Day_length=0;
		else if (tangente>1) Day_length=24;
		else                 Day_length=    24 -  acos(tangente)*7.6394194;
		if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR=0;
		else PAR = PAR_max*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
		if (PAR<0) PAR=0;
	}
	
	if (CLIMAT==0)
	{
		Day_length=12;
		TTTT=T/3600-24*(float)(int)(T/3600/24);
		if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR=0;
		else PAR = PAR_max*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
		if (PAR<0) PAR=0;
		// climatic variable are computed from cos functions with constant values
	
		T_Soil=20;
		T_air  = (T_air_max+T_air_min)/2+(T_air_max-T_air_min)/2*cos(3.14159265359/12*(T/3600-HH1));    //daily cos variation between T_min et T_max et T_max at HH1
		//if (T_air<0.5) T_air= 0.5;
		RH_air = (RH_air_max+RH_air_min)/2+(RH_air_max-RH_air_min)/2*cos(3.14159265359/12*(T/3600-HH2)); //daily cos variation between HR_min et T_max et HR_max at 0h00
		if (HW) // a Heat Wave
		{
			day= T/24/3600;
			if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
			
			 else T_air+=HW_T;
			}
		}
		
	 //   if (IRRIGATE==1 ||IRRIGATE==4 || IRRIGATE==5) Irrigate();  // for cases 2 & 3 irrigate only once at midnight
		
	}
	
	if (CLIMAT==5)
	{
		Day_length=12;
		T_Soil=20;
		T_air  = T_air_max;
		RH_air = RH_air_min;
		PAR=PAR_max;
		//if (IRRIGATE==1 ||IRRIGATE==4 || IRRIGATE==5) Irrigate(); // for cases 2 & 3 irrigate only once at midnight
		if (HW) // a Heat Wave
		{
			day= T/24/3600;
			if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
			 else T_air+=HW_T;
			}
		}
	}
	
	
	if (CLIMAT==1)  // climatic variable interpolated from 2 hourly values
	{
		tangente=tan(Lat*3.1416/180)*tan(asin(sin(23.44*3.1416/180)*sin((DOY-80)/365*2*3.1416)));
		if (tangente<-1)Day_length=0;
		if (tangente>1) Day_length=24;
		else Day_length=    24 -  acos(tangente)*7.6394194;
		T_Soil=20;
		T_air  = T_air_1 + (T-T_1) * (T_air_2-T_air_1)/(T_2-T_1);
		RH_air = RH_air_1 + (T-T_1) * (RH_air_2-RH_air_1)/(T_2-T_1);
		PAR    = PAR_1 + (T-T_1) /(T_2-T_1) * (PAR_2-PAR_1);
	}
	
	if (CLIMAT==2 || CLIMAT==3 || CLIMAT==4)  // climatic variable are computed from cos functions with daily values from a file
	{
		//if (IRRIGATE) Irrigate();  //refill soil water content by irrigation
		// T_air
		if (sin(3.14159265359/12*(T/3600-HH1))<0)   // from Tmim_n to Tmax_n starting at 24-HH1 to HH1
			T_air  = (T_air_max+T_air_min)/2+(T_air_max-T_air_min)/2*cos(3.14159265359/12*(T/3600-HH1));    //daily cos variation between T_min et T_max et T_max at 14h00
		else   //from T_max_n to T_min_n+1 starting at HH1 to 24-HH1
		{
			if (sin(3.14159265359/12*(T/3600-12))>-0.0001)  
					T_air  = (T_air_max+T_air_min_2)/2+(T_air_max-T_air_min_2)/2*cos(3.14159265359/12*(T/3600-HH1));   // before midnigh from T_max_n to T_min_n+1
			else    T_air  = (T_air_max_0+T_air_min)/2+(T_air_max_0-T_air_min)/2*cos(3.14159265359/12*(T/3600-HH1));  // after midnigh from T_max_n-1 to T_min_n
		}
		 
		// T_soil
		if (CLIMAT==4) T_Soil= T_Soil_1 + (T_Soil_2-T_Soil_1)*(T/3600-DOY*24)/24; else T_Soil=20;       // if not measured assumed T_Soil is constant and == 15°C
		
		// RH_air
		if (sin(3.14159265359/12*(T/3600-HH2))>=0)  // from HRmax_n to HRmin_n
			RH_air = (RH_air_max+RH_air_min)/2+(RH_air_max-RH_air_min)/2*cos(3.14159265359/12*(T/3600-HH2)); //daily cos variation between HR_min et T_max et HR_max at 4h00
		else //   from HRmin to HRmax
			if (sin(3.14159265359/12*(T/3600-12))>-0.001)  // before midnigh from T_max_n to T_min_n+1
					RH_air = (RH_air_max_2+RH_air_min)/2+(RH_air_max_2-RH_air_min)/2*cos(3.14159265359/12*(T/3600-HH2)); //daily cos variation between HR_min et T_max et HR_max at 4h00
			else    RH_air = (RH_air_max+RH_air_min_0)/2+(RH_air_max-RH_air_min_0)/2*cos(3.14159265359/12*(T/3600-HH2));
	}
	
	if (IRRIGATE==1 ||IRRIGATE==4 || IRRIGATE==5) Irrigate();  // for cases 2,3,6,7 irrigate only once at midnight or 19h00
	
	if (RH_air<0.1) RH_air=0.1;
	if (RH_air>100) RH_air=99.999;
	if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
				 e_air=e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
	P_air =     0.4609418*(273.15+T_air)*log(e_air/e_sat_air);              // air water potential, MPa
	VPD_Air=    (e_sat_air - e_air)/1000;                                  // vpd in kPa
	if (VPD_Air<0) VPD_Air=0;
	VPD_Air_tot+=VPD_Air*dt_long*dt/3600/24;
	slope=4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
	ETP_Penman= 0.5625*(0.408*slope*PAR*0.5495*3.6e-3+0.066*37*Wind*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind));  //ETP in mm/h  0.5625 to fit observed ETP
	ETP_Penman/=3600;  // ETP in mm/s
	ETP_Penman=ETP_Penman/18*1000*1000; // ETP in mmol/s  to be multiplied by soil area to obtain the volume
	POTENTIAL_PAR=Potential_PAR((T/3600-DOY*24));
	// LEAF
	if (!TLEAF) T_Leaf=T_air;             //leaf temperature is then assumed to equal the air temperature
	if (T_Leaf>0) e_sat=611.21*exp((18.678-T_Leaf/234.5)*T_Leaf/(257.14+T_Leaf));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else          e_sat=611.15*exp((23.036-T_Leaf/333.7)*T_Leaf/(279.82+T_Leaf));
	
	e=        e_sat*exp(P_Leaf_Evap*2.16947115/(T_Leaf+273.15));
	if (Leaf_Area)  VPD_Leaf= (e-e_air)/1000; //vpd between leaf and air in kPa
	else            VPD_Leaf=0;
	if (VPD_Leaf<0) VPD_Leaf=0;  // to avoid back flow
	VPD_Leaf_tot+=  VPD_Leaf*dt_long*dt/24/3600;
	
	e=        e_sat*exp(P_Leaf_Symp*2.16947115/(T_Leaf+273.15));
	VPD_Cuti= (e-e_air)/1000; //vpd between leaf and air in kPa
	if (VPD_Cuti<0) VPD_Cuti=0;  // to avoid back flow
	
	// BUD
	if (Type_Axil)
	{
		T_Axil=T_air;  // could also be T_leaf, or T_bud!
		if (T_Axil>0) e_sat=611.21*exp((18.678-T_Axil/234.5)*T_Axil/(257.14+T_Axil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else          e_sat=611.15*exp((23.036-T_Axil/333.7)*T_Axil/(279.82+T_Axil));
		e=            e_sat*exp(P_Axil_Symp*2.16947115/(T_Axil+273.15));
		VPD_Axil=(e-e_air)/1000; //vpd between branch and air in mmol/m3
		if (VPD_Axil<0) VPD_Axil=0;
		if (Type_Axil>=2) //petiole
		{
		T_Axil=T_air;  // could also be T_leaf, or T_bud!
		if (T_Axil>0) e_sat=611.21*exp((18.678-T_Axil/234.5)*T_Axil/(257.14+T_Axil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else          e_sat=611.15*exp((23.036-T_Axil/333.7)*T_Axil/(279.82+T_Axil));
		e=            e_sat*exp(P_Petiole_Symp*2.16947115/(T_Axil+273.15));
		VPD_Petiole=(e-e_air)/1000; //vpd between branch and air in mmol/m3
		if (VPD_Petiole<0) VPD_Petiole=0;
		}

	}
	// BRANCH
	T_Branch=T_air;
	if (T_Branch>0) e_sat=611.21*exp((18.678-T_Branch/234.5)*T_Branch/(257.14+T_Branch));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Branch/333.7)*T_Branch/(279.82+T_Branch));
	e=              e_sat*exp(P_Branch_Symp*2.16947115/(T_Branch+273.15));
	VPD_Branch=(e-e_air)/1000; //vpd between branch and air in mmol/m3
	if (VPD_Branch<0) VPD_Branch=0;
	
	// TRUNK
	T_Trunk=T_air;
	if (T_Trunk>0)  e_sat=611.21*exp((18.678-T_Trunk/234.5)*T_Trunk/(257.14+T_Trunk));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Trunk/333.7)*T_Trunk/(279.82+T_Trunk));
	e=              e_sat*exp(P_Trunk_Symp*2.16947115/(T_Trunk+273.15));
	VPD_Trunk=(e-e_air)/1000; //vpd between Trunk and air in mmol/m3
	if (VPD_Trunk<0) VPD_Trunk=0;
	
	// ROOT
	
	if (CLIMAT !=4) T_Root=T_Soil;
	else T_Root=T_Soil1;
	if (T_Root>0)   e_sat=611.21*exp((18.678-T_Root/234.5)*T_Root/(257.14+T_Root));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Root/333.7)*T_Root/(279.82+T_Root));
	
	e_air_sol1=     e_sat*exp(P_Soil1 *2.16947115/(T_Soil+273.15));               // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp1*2.16947115/(T_Root+273.15));         // root vapour pressure at root water potential
	VPD_Root1=     (e-e_air_sol1)/1000; //vpd between Root   and air in mmol/m3
	if (VPD_Root1<0) VPD_Root1=0;
	
	if (CLIMAT ==4) T_Root=T_Soil2;
	e_air_sol2=     e_sat*exp(P_Soil2 *2.16947115/(T_Soil+273.15));           // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp2*2.16947115/(T_Root+273.15));
	VPD_Root2=     (e-e_air_sol2)/1000; //vpd between Root and air in mmol/m3
	if (VPD_Root2<0) VPD_Root2=0;
	
	if (CLIMAT ==4)T_Root=T_Soil3;
	e_air_sol3=     e_sat*exp(P_Soil3 *2.16947115/(T_Soil+273.15));           // soil vapour pressure at soil water potential
	e=              e_sat*exp(P_Root_Symp3*2.16947115/(T_Root+273.15));
	VPD_Root3=     (e-e_air_sol3)/1000; //vpd between Root and air in mmol/m3
	if (VPD_Root3<0) VPD_Root3=0;
	
	// SOIL
	if (CLIMAT !=4)T_Root=T_Soil;
	else T_Root=T_Soil1;
	if (T_Soil>0)   e_sat=611.21*exp((18.678-T_Soil/234.5)*T_Soil/(257.14+T_Soil));    // saturation vapour water pressure at Tair in Pa from Buck's equation
	else            e_sat=611.15*exp((23.036-T_Soil/333.7)*T_Soil/(279.82+T_Soil));
	
	e=              e_sat*exp(P_Soil1*2.16947115/(T_Soil+273.15));
	VPD_Soil=      (e-e_air)/1000; //vpd between soil and air in mmol/m3
	if (VPD_Soil<0) VPD_Soil=0;
	if (CLIMAT!=4)
	{
		T_Soil1=T_Soil;
		T_Soil2=T_Soil;
		T_Soil3=T_Soil;
	}
}

void Interception (void) //rain interception from foliage
{
	long double In; //interception in mm
	
	if (Leaf_Area)
	{
		if (Leaf_rain<0) Leaf_rain=0;
		if (Rain_1>Threshold_rain)  In=Rain_1*(Interception_min+(100- Interception_min)/(exp(Interception_factor*(Rain_1 - Threshold_rain))))/100;  // then compute the interception in mm
		else                        In=Rain_1;
		Leaf_rain+=In;                                                                                                  // fill leaf reservoir with IN, max is Interception_min
		if (Leaf_rain>Threshold_rain)                                                                                   // reservoir is saturated, the excess is not intercepted
		{
			In -= (Leaf_rain-Threshold_rain);
			Leaf_rain=Threshold_rain ;
		}
		Rain_1 = Rain_1-In;
	}
	Rain_soil+=Rain_1;
}

void LA_acclimatation(void)  //Leaf area acclimatation from year to year; work only in continuuous mode!
{
	if (LA_Var==2) LA_max=LA_max_init*(100-PLC_Branch_Apo)/100;
	if (LA_Var==3) LA_max=LA_max_init*(100-PLC_Leaf_Apo)/100;
}

void purge(void)
{
//	printf("skip end of year \n");
	while (T_2>T_1 && !END_CLIMAT) 
	{
		++N_days;
		if (Leaf_Area) ++N_days2;
			
		if (!feof(climat_in))
			{
			if (CLIMAT==2)   fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
			else             fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
			}
			else END_CLIMAT=1;
			T_2=T_2*3600*24;
			Cum_T_air+= (T_air_min_2+T_air_max_2)/2;
			if (Leaf_Area)   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
			Rain_tot+=Rain_2;
			if (Leaf_Area) Rain_leaf_tot+=Rain_2;
	} 
	
	DEAD=0;
	if (PRINT_SCREEN) print_screen();
	if (TRANSIENT) print_transient();
	DOY=T_2/3600/24;
	if (LA_Var>1) LA_acclimatation();
	T0=T_2;
	T=T0;
	T_1=T_0;
	indice=0;
	indice_double=0;
	if (!END_CLIMAT) Reset(); 
}

void next_climat (void)  // load new set of climatic data
{
	long double buffer, day , e_sat_air, e_air;
	if (CLIMAT==1)  //data on a hour basis
	{
		if (T>T_2)      // reach the end of the time interval; load new values
		{
			if (floor(T_2/3600/24)>floor(T_1/3600/24))  //a new day
			{
				ETP_Penman_day=0; // reset at midnight
				A_net_day=0;
				E_tot_day=0;
				EvapoT_day=0;
				P_min_lf_d=0;
				P_max_lf_d=-1000;
				gs_min_d=10000;
				gs_max_d=0;
				gs_max_d2=0;
				SF_min_d=10000;
				SF_max_d=0;
			}
			Rain_tot+=Rain_1;
			if (INTERCEPTION)Interception();
			else Rain_soil+=Rain_1;
			if (Leaf_Area) Rain_leaf_tot+=Rain_1;
			Q_Soil1+=Rain_1*Surface_Soil*1000*1000/18; //add rain from the past interval to soil
			fill_soil_layers();
	
			T_1=T_2;
			T_air_1=T_air_2;
			RH_air_1=RH_air_2;
			PAR_1=PAR_2;
			T_Soil_1=T_Soil_2;
			/*if (CLIMAT==4)
			{
				T_Soil_11=T_Soil_21;
				T_Soil_12=T_Soil_22;
				T_Soil_13=T_Soil_23;
				if (T_Soil_11<T_Soil_Crit || T_Soil_12<T_Soil_Crit || T_Soil_13<T_Soil_Crit) T_LIMIT=1;
				else if (T_Soil_11>T_Soil_Crit+0.1 && T_Soil_12>T_Soil_Crit+0.1 && T_Soil_13>T_Soil_Crit+0.1) T_LIMIT=0;
				else T_LIMIT=1;
			}*/
		//    if (IRRIGATE) Irrigate();  //refill soil water content at midnight
			Rain_1=Rain_2;
			if (!feof(climat_in))
			{
				YEAR1=YEAR2;
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR2, &T_2, &PAR_2,&T_air_2, &RH_air_2,&Rain_2, &Wind);
				T_2=T_2*3600*24;  //time in sec
			   // if (T_air_2 < 0.5) T_air_2=0.5; // otherwise leaf energy budget doesn't work
				if (HW) // a Heat Wave
				{
					day= T/24/3600;
					if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
						else T_air_2+=HW_T; 
					}
				}
					
				
		
				if (T_2<T_1) // then it is a new run or a new year. Reset simulation to zero
				{
					DOY=T_2/3600/24;
					if (LA_Var>1) LA_acclimatation();
					 T0=T_2;
					 T=T0;
					 T_1=T_0;
					 indice=0;
					 indice_double=0;
					 Reset();
				}
			}
			else END_CLIMAT=1;
		}
	}
	
	if (CLIMAT==2 || CLIMAT==4) // data on a day basis
	{
		if (T>=T_2) //end of day or DEAD 
		{
			++N_days;
			if (Leaf_Area) ++N_days2;
			Rain_tot+=Rain_1;
			if (Leaf_Area) Rain_leaf_tot+=Rain_1;
			if (INTERCEPTION) Interception();
			else Rain_soil+=Rain_1;
			ETP_Penman_day=0; // reset at midnight
			A_net_day=0;
			E_tot_day=0;
			EvapoT_day=0;
			P_min_lf_d=0;
			P_max_lf_d=-1000;
			gs_min_d=10000;
			gs_max_d=0;
			gs_max_d2=0;
			SF_min_d=10000;
			SF_max_d=0;
			Q_Soil1+=Rain_1*Surface_Soil*1000*1000/18; //add rain from the past interval to soil
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
			PAR_max=PAR_max_2;
			Rain_1=Rain_2;
			T_Soil=T_Soil_2;
			if (CLIMAT==4)
			{
				T_Soil1=T_Soil_21;
				T_Soil2=T_Soil_22;
				T_Soil3=T_Soil_23;
				if (T_Soil1<T_Soil_Crit || T_Soil2<T_Soil_Crit || T_Soil3<T_Soil_Crit) T_LIMIT=1;
				else if (T_Soil1>T_Soil_Crit+0.1 && T_Soil2>T_Soil_Crit+0.1 && T_Soil3>T_Soil_Crit+0.1) T_LIMIT=0;
				else T_LIMIT=1;
			}
			T_1=T_2;
		   
			if (!feof(climat_in))
			{
				YEAR1=YEAR2;
				DOY=T_2/3600/24;
				if (CLIMAT==2) fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
				else           fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_21,&T_Soil_22,&T_Soil_23, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
			
				T_2=T_2*3600*24;
				if (HW) // a Heat Wave
				{
					day= T/24/3600;
					if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
						else 
						{
						T_air_max_2+=HW_T;
						T_air_min_2+=HW_T;
						}
					}
				}
			  
				Cum_T_air += (T_air_min_2+T_air_max_2)/2;
				if (Leaf_Area)   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
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
				if (Wind<0.1) Wind=0.1;
				
				if (T_2<T_1) // then it is a new run or a new year. Reset simulation to zero
				{
					DOY=T_2/3600/24;
					if (LA_Var>1) LA_acclimatation();
					 T0=T_2;
					 T=T0;
					 T_1=T_0;
					 indice=0;
					 indice_double=0;
					 Reset();
				}
				
				if (T_2<T_1) // then it is a new run or a new year. Reset simulation to zero
				{
					DOY=T_2/3600/24;
					if (LA_Var>1) LA_acclimatation();
					 T0=T_2;
					 T=T0;
					 T_1=T_0;
					 indice=0;
					 indice_double=0;
					 Reset();
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
			if (Leaf_Area) ++N_days2;
			Cum_T_air += (T_air_min_1+T_air_max_1)/2;
			if (Leaf_Area)   Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
			
			proba=(long double)(rand() % 1000)/1000;                 // probability of rain this day
			if (proba<Proba_rain) rain=Rain_1/Proba_rain;       // rain this this
			else rain=0;
			Rain_tot+=rain;
			if (INTERCEPTION) Interception();
			else Rain_soil+=rain;
			Q_Soil1+=rain*Surface_Soil*1000*1000/18; //add rain from the past interval to soil
			fill_soil_layers();
			
	   //     if (IRRIGATE) Irrigate();  //refill soil water content by irrigation
		}
		
		if (T>T_2)                // a new month
		{
			T_1=T_2;
			if (!feof(climat_in))
			{
				fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Proba_rain, &Wind);
				
				if (HW) // a Heat Wave
				{
					day= T/24/3600;
					if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
				if (Leaf_Area)   Cum_T_air_l += (T_air_min_2+T_air_max_2)/2;
				T_air_min=T_air_min_2;  //set values for the first day
				T_air_max=T_air_max_2;
				T_air_max_0=T_air_max_2;
				T_air_min_1=T_air_min_2;
				T_air_max_1=T_air_max_2;
				
				RH_air_min=RH_air_min_2;
				RH_air_max=RH_air_max_2;
				PAR_max=PAR_max_2;
				RH_air_min_1=RH_air_min_2;
				RH_air_max_1=RH_air_max_2;
				PAR_max_1=PAR_max_2;
				Rain_1=Rain_2;
			}
			else END_CLIMAT=1;
			T=T_1;
		}
	}
	if (Wind<=0) Wind=0.1;
	if (CLIMAT!=4)
	{
		T_Soil1=T_Soil;
		T_Soil2=T_Soil;
		T_Soil3=T_Soil;
	}
}

void load_climat (void)  // load the first set of climatic data from .txt files
{
	long double day, e_sat_air, e_air;
	if (CLIMAT==1) // data on a hour basis
	{
	if ((climat_in = fopen("climat_hour_in.txt","r+"))==NULL) CLIMAT=0;
	else
		{
		fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR1,&T_1, &PAR_1,&T_air_1, &RH_air_1,&Rain_1, &Wind);
		fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le\n", &YEAR2,&T_2, &PAR_2,&T_air_2, &RH_air_2,&Rain_2, &Wind);
		if (HW) // a Heat Wave
		{
			day= T/24/3600;
			if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
		
		T_1=T_1*24*3600;        // time of the day in secondes
		T_2=T_2*24*3600;
		T0=T_1;
		T=T_1;
		T_air=T_air_1;
		RH_air=RH_air_1;
		PAR=PAR_1;
	  }
	}
	
	if (CLIMAT==2 || CLIMAT==4) // data on a day basis
	{
		if ((climat_in = fopen(filename_CLIM,"r+"))==NULL) CLIMAT=0;
		else
		{
		++N_days;
		if (Leaf_Area) ++N_days2;
	   
		if (CLIMAT==2)
		{
			fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &T_air_min_1, &T_air_max_1, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Wind);
			fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
	   }
		else   // CLIMAT==4 inludes soil temperature in the climat file
		{
			fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR1, &T_1, &T_air_min_1, &T_air_max_1, &T_Soil_11,&T_Soil_12,&T_Soil_13, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Wind);
			fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR2, &T_2, &T_air_min_2, &T_air_max_2, &T_Soil_21, &T_Soil_22, &T_Soil_23, &RH_air_min_2, &RH_air_max_2, &PAR_max_2, &Rain_2, &Wind);
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
			if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
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
		if (Leaf_Area)  Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
		T_air_min=      T_air_min_1;  //set values for the first day
		T_air_max=      T_air_max_1;
		T_air_max_0=    T_air_max_1;
		RH_air_min=     RH_air_min_1;
		RH_air_max=     RH_air_max_1;
		PAR_max=        PAR_max_1;
		DOY=T_1;
		T_1=T_1*24*3600;            // time of the day in secondes
		T_2=T_2*24*3600;
		T=T_1;
		}
	}
	
	if (CLIMAT==3) // monthly average of daily values (including for precipitation)
	{
		if ((climat_in = fopen("climat_month_in.txt","r+"))==NULL) CLIMAT=0;
		else 
		{
		++N_days;
		if (Leaf_Area) ++N_days2;
	   
		fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&T_1, &T_air_min_1, &T_air_max_1, &RH_air_min_1, &RH_air_max_1, &PAR_max_1, &Rain_1, &Proba_rain, &Wind);
		if (HW) // a Heat Wave
		{
			day= T/24/3600;
			if (day>HW_day && day<(HW_day+HW_duration)) // a Heat wave starting at day HW_day, lasting HW_duration days with a temperature increase of HW_T °C
			{
				if (HW==1) RH_air_min_1*= exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1))/exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+T_air_max_1+HW_T));  //new RH_air assuming constant e_air
				 if (HW==3) RH_air_min_1= 100-(611.21*exp((18.678-T_air_max_1/234.5)*T_air_max_1/(257.14+T_air_max_1)))/(611.21*exp((18.678-(T_air_max_1+HW_T)/234.5)*(T_air_max_1+HW_T)/(257.14+(T_air_max_1+HW_T))))*(100-RH_air_min_1);    
				   
				 T_air_max_1+=HW_T;
				T_air_min_1+=HW_T;
			}
		}
		Cum_T_air+= (T_air_min_1+T_air_max_1)/2;
		if (Leaf_Area)   Cum_T_air_l += (T_air_min_1+T_air_max_1)/2;
		T_air_min=T_air_min_1;  //set values for the first day
		T_air_max=T_air_max_1;
		T_air_max_0=T_air_max;
		RH_air_min=RH_air_min_1;
		RH_air_max=RH_air_max_1;
		PAR_max=PAR_max_1;
		
		T_air_min_2=T_air_min_1;
		T_air_max_2=T_air_max_1;
		RH_air_min_2=RH_air_min_1;
		RH_air_max_2=RH_air_max_1;
		PAR_max_2=PAR_max_1;
		Rain_2=Rain_1;
		T_2=(T_1)*24*3600*365/12;
		T_1=0 ;          // time of the day in secondes
		T=T_1;
		}
	}
	if (Wind<=0) Wind=0.1;
	
}

void E_day(long double dt_court, long double dt_long) //compute plant transpiration based on ETP, gs and VPD on mmol/m2/s
{
	long double E_Soil1, E_Soil2,E_stomata;
	//FILE *out;
	long double K_root, K_root1, K_root2, K_root3, K_to_branch, K_to_leaf, P_Soil;
	long double gs_max2, gs_night2, PL_gs, slope_gs, P50_gs;
	long double SLmax;
	
	gs_0=g_s;
	g_crown=g_crown0*pow(Wind,0.6);             // wind effect on canopy conductance
	
	if (GS_MAX==1)      gs_max2 = gs_max/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));                // temperature effect on gs_max
	else if (GS_MAX==2) gs_max2=  gs_max*300/Ca;                                                // CO2 effect on gs_max; assume 300ppm in 1950
	else if (GS_MAX==3) gs_max2 = gs_max*300/Ca/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));         // CO2+T effects on gs_max
	else if (GS_MAX==4) gs_max2=  gs_max*pow(300/Ca,0.5);  
	else if (GS_MAX==5) gs_max2 = gs_max*pow(300/Ca,0.5)/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2)); 
	else if (GS_MAX==6) gs_max2 = gs_max*(1 + gs_CO2_sens/100*(Ca-300)/100);
	else if (GS_MAX==7) gs_max2 =(gs_max*(1 + gs_CO2_sens/100*(Ca-300)/100))/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));
	else                gs_max2=  gs_max;  
	if (gs_max2<0) gs_max2=0;                                                     // gs_max does not vary with temperature and CO2
	
	if (GS_MAX==1)      gs_night2 = gs_night/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));                // temperature effect on gs_min
	else if (GS_MAX==2) gs_night2=  gs_night*300/Ca;                                                // CO2 effect on gs_night; assume 300ppm in 1950
	else if (GS_MAX==3) gs_night2 = gs_night*300/Ca/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));         // CO2+T effects on gs_night
	else if (GS_MAX==4) gs_night2=  gs_night*pow(300/Ca,0.5);                                                // CO2 effect on gs_night; assume 300ppm in 1950
	else if (GS_MAX==5) gs_night2 = gs_night*pow(300/Ca,0.5)/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));         // CO2+T effects on gs_night
	else if (GS_MAX==6) gs_night2 = gs_night* (1 + gs_CO2_sens/100*(Ca-300)/100);
	else if (GS_MAX==7) gs_night2 = (gs_night* (1 + gs_CO2_sens/100*(Ca-300)/100))/(1+pow((T_Leaf-Tgs_optim)/Tgs_sens,2));
	else                gs_night2=  gs_night;                                                       // gs_night does not vary with temperature and CO2
	if (gs_night2<0) gs_night2=0;
	
	if (!gs_cst) gs_Jarvis=   (gs_night2 + (gs_max2-gs_night2)*(1-exp(-Jarvis_PAR*PAR)));           // limitation of gs by light following Jarvis
	else         gs_Jarvis=   (gs_night2 + (gs_max2-gs_night2)*(1-exp(-Jarvis_PAR*PAR_max)));       // gs does not respond to PAR
	
	if (Leaf_Area && gs_Jarvis)  E_clim =    1/ (1/g_bl + 1/(gs_Jarvis+g_cuti)  + 1/g_crown) * VPD_Leaf/101.3;                //leaf  evaporation rate due to climatic demand in mmol s-1 m-2
	else E_clim=0;
	
	if (PENMAN)
	{
		if (Penman_Coeff==0) Penman_Coeff= -0.006*LAI*LAI +0.134*LAI + 0.036;           // Empirical formula from Granier
		if ((E_clim*Leaf_Area) >(ETP_Penman*Penman_Coeff*Surface_Soil)) E_clim=ETP_Penman*Penman_Coeff*Surface_Soil/Leaf_Area;        //Evaporation limited to ETP_Penman*coeff following Granier. soil evaporation is not included!
	}
	
	if (Leaf_Area && g_cuti) 	E_cuti =    1/( 1/g_bl + 1/g_cuti + 1/g_crown)  * VPD_Cuti/101.3;                //evaporation through leaf cuticle due to climatic demand
	else                     	E_cuti = 0;
	if (Type_Axil)
	{
	//g_Axil=(g_Axil_max-g_Axil_min)*(P_Axil_Symp + 2) +g_Axil_min;
	g_Axil=g_Axil_max*(0.76*(P_Axil_Symp-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
	if (g_Axil<g_Axil_min) g_Axil=g_Axil_min;
	if (g_Axil>g_Axil_max) g_Axil=g_Axil_max;
	}
	else g_Axil=0;
	E_Branch=   g_Branch*VPD_Branch/101.3;        //evaporation through branch cuticle due to climatic demand
	E_Trunk=    g_Trunk*VPD_Trunk/101.3;          //evaporation through trunk cuticle due to climatic demand
	E_Axil =    g_Axil*VPD_Axil/101.3;            //for all the axillary organs
	if (Type_Axil >= 2) E_Petiole=  g_Petiole*VPD_Petiole/101.3;  else E_Petiole=0;
	E_Root1=    g_Root1*VPD_Root1/101.3;          //evaporation through root cuticle due to climatic demand
	E_Root2=    g_Root2*VPD_Root2/101.3;
	E_Root3=    g_Root3*VPD_Root3/101.3;
	
	g_Soil=g_Soil0*RWC1;
	E_Soil1=    g_Soil*VPD_Soil/101.3;                      //VPD effect
	if (g_Soil0) E_Soil2=    g_Soil/g_Soil0*ETP_Penman*exp(-0.5*LAI);    // limitation by ETP depending on radiation reaching the soil
	else         E_Soil2=0;
	E_Soil=     min(E_Soil1,E_Soil2);
	if (T_Soil<0 || T_Soil1<0) E_Soil=0;                       //no evaporation from frozen soil

	dq_Branch=   E_Branch*Branch_Area*dt_court;           // flows during dt in mmol
	dq_Axil=     E_Axil*Area_Axil*dt_court;               // for all the axillary organs
	if (Type_Axil >= 2) dq_Petiole=  E_Petiole*Petiole_area*dt_court;   else dq_Petiole=0;     // for all petioles
	dq_cuti=     Leaf_Area*E_cuti*dt_court;              // Plant evaporation through cuticle only during dt in mmol on the time interval dt; use one LA !
	if (Type_Axil == 2) dq_fruit=    Growth_rate_fruit*dt_court;  else dq_fruit=0;           // water flow due to fruit growth
	dq_Trunk=    E_Trunk*Trunk_Area*dt_court;
	dq_Root1=    E_Root1*Root_Area1*dt_court;
	dq_Root2=    E_Root2*Root_Area2*dt_court;
	dq_Root3=    E_Root3*Root_Area3*dt_court;
	dq_Soil=     E_Soil*Soil_Width*Soil_Width*dt_court;
 
	if (K_Soil1 && K_Interface1 &&K_Root_Symp11 && K_Root_Apo1) K_root1=1/(1/K_Soil1 + 1/K_Interface1 + 1/K_Root_Symp11 + 1/K_Root_Apo1); else K_root1=0;
	if (K_Soil2 && K_Interface2 &&K_Root_Symp12 && K_Root_Apo2) K_root2=1/(1/K_Soil2 + 1/K_Interface2 + 1/K_Root_Symp12 + 1/K_Root_Apo2); else K_root2=0;
	if (K_Soil3 && K_Interface3 &&K_Root_Symp13 && K_Root_Apo3) K_root3=1/(1/K_Soil3 + 1/K_Interface3 + 1/K_Root_Symp13 + 1/K_Root_Apo3); else K_root3=0;
	K_root=K_root1+K_root2+K_root3;
	if (Leaf_Area && K_root && K_Trunk_Apo && K_Branch_Apo && K_Leaf_Apo && K_Leaf_Symp)  K_tot=1/(1/K_root+1/K_Trunk_Apo+1/K_Branch_Apo+1/K_Leaf_Apo +1/K_Leaf_Symp)/Leaf_Area;
	else K_tot=0;

	switch((int)Regul_gs) // if gs is computed then E is derived and vice versa
	{
		
		case 0: // gs is not regulated by water stress but still respond to PAR
			E_Leaf=E_clim;
			break;
			
		case 1:            // gs is decreased proportionnaly to the decrease in turgor pressure; ref is midday turgor of control
			if (Turgor_Leaf_Symp<Turgor_Leaf_Symp_Ref) turgor=Turgor_Leaf_Symp/Turgor_Leaf_Symp_Ref;
			else  turgor=1;
			g_s=gs_max2*turgor;
			break;
			
		case 2:  // gs regulation by turgor_max=-PI0*turgor_ref_factor
			if (Turgor_Leaf_Symp<(-Pi0_Leaf_Symp*turgor_ref_factor))  turgor=-Turgor_Leaf_Symp/(Pi0_Leaf_Symp*turgor_ref_factor);
			else  turgor=1;
			g_s=gs_max2*turgor;
			break;
			
		case 3:     // E regulated at P_branch_Apo, MPa
			K_to_branch=1/(1/K_root+1/K_Trunk_Apo+1/K_Branch_Apo);
			if (K_root) P_Soil= (P_Soil1*K_root1 + P_Soil2*K_root2 + P_Soil3*K_root3)/K_root;
			else        P_Soil=Px_gs+Pg_Leaf;
			if (Leaf_Area) E_Leaf=Gamma*(K_to_branch*(P_Soil-Px_gs+Pg_Leaf))/Leaf_Area;
			else           E_Leaf=0;
			 E_stomata=E_Leaf-Gamma*E_cuti;
			if (VPD_Leaf) g_Canopy=E_stomata/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown and g_bl
			if (g_Canopy) g_s= 1/(1/g_Canopy - 1/g_bl - 1/g_crown); else g_s=0;
			break;
			
		case 4:     // E regulated at P_leaf_Apo MPa
			K_to_Leaf_Apo=1/(1/K_root+1/K_Trunk_Apo+1/K_Branch_Apo+1/K_Leaf_Apo);
			if (K_root) P_Soil= (P_Soil1*K_root1 + P_Soil2*K_root2 + P_Soil3*K_root3)/K_root;
			else        P_Soil=Px_gs+Pg_Leaf;
			if (Leaf_Area) E_Leaf=Gamma*(K_to_Leaf_Apo*(P_Soil-Px_gs+Pg_Leaf))/Leaf_Area;
			else           E_Leaf=0;
			E_stomata= E_Leaf - Gamma*E_cuti;
			if (VPD_Leaf) g_Canopy=E_stomata/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown and g_bl
			if (g_Canopy) g_s= 1/(1/g_Canopy - 1/g_bl - 1/g_crown); else g_s=0;
			break;
			
		case 5:     // E regulated at P_leaf_symp, MPa
			K_to_leaf=1/(1/K_root+1/K_Trunk_Apo+1/K_Branch_Apo+1/K_Leaf_Apo +1/K_Leaf_Symp);
			if (K_root) P_Soil= (P_Soil1*K_root1 + P_Soil2*K_root2 + P_Soil3*K_root3)/K_root;
			else P_Soil=Px_gs+Pg_Leaf;
			if (Leaf_Area) E_Leaf=Gamma*(K_to_leaf*(P_Soil-Px_gs+Pg_Leaf))/Leaf_Area;
			else           E_Leaf=0;
			E_stomata= E_Leaf - Gamma*E_cuti;
			if (VPD_Leaf) g_Canopy=E_stomata/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown and g_bl
			if (g_Canopy) g_s= 1/(1/g_Canopy - 1/g_bl - 1/g_crown); else g_s=0;
			break;
			
		case 6:   // E regulated at PLCx_leaf, %
			K_to_Leaf_Apo=1/(1/K_root+1/K_Trunk_Apo+1/K_Branch_Apo+1/K_Leaf_Apo);
			if (K_root) P_Soil= (P_Soil1*K_root1 + P_Soil2*K_root2 + P_Soil3*K_root3)/(K_root);
			else        P_Soil=Px_Leaf_Apo+Pg_Leaf;
			if (Leaf_Area) E_Leaf=Gamma*(K_to_Leaf_Apo*(P_Soil-Px_Leaf_Apo+Pg_Leaf))/Leaf_Area;
			else           E_Leaf=0;
			 E_stomata= E_Leaf - Gamma*E_cuti;
			if (E_stomata<0) E_stomata=0;
			if (VPD_Leaf) g_Canopy=E_stomata/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown and g_bl
			if (g_Canopy) g_s= 1/(1/g_Canopy - 1/g_bl - 1/g_crown); else g_s=0;
			break;
			
		case 8: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Soil
			g_s=gs_max2*(0.76*(P_soil-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
			
		case 9: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Leaf_Symp
			g_s=gs_max2*(0.76*(P_Leaf_Symp-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
		
		case 10: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Leaf_Symp
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Symp-P50_gs)));
			g_s = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 11: //   a power function of  g_s = Regul_gs_para1*pow(-P_Leaf_Apo,Regul_gs_para2)
			if (P_Leaf_Apo<0) g_s = Regul_gs_para1*pow(-P_Leaf_Apo,Regul_gs_para2); else g_s=gs_max2;
			break;
			
		case 12: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Leaf_Apo
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Apo-P50_gs)));
			g_s = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 13: // gs regulated by Pgs_12 and Pgs_88 with a sigmoidal fit to P_Branch_Apo
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Branch_Apo-P50_gs)));
			g_s = gs_max2 *(100-PL_gs)/100;
			break;
			
		case 14: // gs regulated by Pgs_12 and Pgs_88 with a linear fit to P_Branch_Apo
			g_s=gs_max2*(0.76*(P_Branch_Apo-Pgs_88)/(Pgs_12-Pgs_88)+0.12);
			break;
			
		case 15: //like 2 but with a sigmoid fit between P_leaf_symp at 0.12 *max_turgor*para1 and 0.88*max_turgor*para1
			Pgs_12=-Pi0_Leaf_Symp*turgor_ref_factor*Osmotic_TLeaf*0.88+Epsilon_Leaf_Symp*Pi0_Leaf_Symp*Osmotic_TLeaf/(Epsilon_Leaf_Symp+Pi0_Leaf_Symp*Osmotic_TLeaf-Pi0_Leaf_Symp*turgor_ref_factor*Osmotic_TLeaf*0.88);
			Pgs_88=-Pi0_Leaf_Symp*turgor_ref_factor*Osmotic_TLeaf*0.12+Epsilon_Leaf_Symp*Pi0_Leaf_Symp*Osmotic_TLeaf/(Epsilon_Leaf_Symp+Pi0_Leaf_Symp*Osmotic_TLeaf-Pi0_Leaf_Symp*turgor_ref_factor*Osmotic_TLeaf*0.12);
			P50_gs = (Pgs_12 + Pgs_88)/2;
			slope_gs = 100/(Pgs_12-Pgs_88);
			PL_gs = 100/(1+exp(slope_gs/25*(P_Leaf_Symp-P50_gs)));
			g_s = gs_max2 *(100-PL_gs)/100;
			break;

		default: E_Leaf=E_clim;
			break;
	}
	
	if (g_s<0) g_s=0;
	if (!gs_tc) if (g_s > gs_Jarvis) g_s = gs_Jarvis;
	if (g_s > gs_max2) g_s = gs_max2;
	if (g_s) g_Canopy=1/ (1/g_bl + 1/g_s + 1/g_crown); else g_Canopy=0;
	E_stomata=g_Canopy*VPD_Leaf/101.3;
	E_Leaf = E_stomata + E_cuti;
	if (E_Leaf >= E_clim) 
	{
		E_Leaf = E_clim;      // E is less than climatic demand
		g_s=gs_Jarvis;
	}
	else if (T<T_gs_regul) T_gs_regul=(T-T0);    // define the timing of onset of stomatal regulation
   

   if (gs_tc) // if the time constant of stomatal response is not null 
	{
	SLmax= (g_s-gs_0)/(gs_tc*2.718);
	g_s=gs_0+SLmax*dt_long*dt;
	if (g_s) g_Canopy=1/ (1/g_bl + 1/g_s + 1/g_crown); else g_Canopy=0;
	E_stomata=g_Canopy*VPD_Leaf/101.3;
	E_Leaf = E_stomata+ E_cuti;
	} 
	 
	dq_stomate=  Leaf_Area*(E_Leaf-E_cuti)*dt_court;     // Plant evaporation through stomata only during dt in mmol
	if (VPD_Leaf) g_Canopy=E_stomata/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown, and g_bl not g_cuti so that g_s can be computed
	if (g_Canopy) g_s= 1/(1/g_Canopy - 1/g_bl - 1/g_crown); else g_s=0;
	if (VPD_Leaf) g_Canopy=E_Leaf/VPD_Leaf*101.3;  else g_Canopy = 0; //total canopy vapor conductance including g_s, g_crown, g_cuti and g_bl
  
	if (g_s<0) g_s=0; // this is to eliminate very negative values when 100% PLC in the leaves
	if (!Leaf_Area) g_s=0;
	if ((g_s + g_cuti)>gs_max_d) gs_max_d = g_s + g_cuti;
	if ((g_s + g_cuti)<gs_min_d) gs_min_d = g_s + g_cuti;
  
 // STATS  
	if (g_s == gs_night2) if (E_Leaf < E_Leaf_night) E_Leaf_night = E_Leaf;
	if (E_Leaf <= (E_cuti*1.1))  //when sigmoid fit E never reaches E_cuti so 1.1 is necessary
	{
		E_Leaf = E_cuti;                                  	 // E cannot be less than E_cuti
		if (E_Leaf>E_max_gs_close) E_max_gs_close=E_Leaf;    // Max E_leaf when stomata are closed
	}
	else T_gs_close=(T-T0);         							// define the timing of total of stomatal regulation
	if (E_Leaf > E_max) E_max=E_Leaf;  // grand Max E_leaf
 
 if (!Leaf_Area)
	{
	gs_max_d  = 0;
	gs_min_d  = 0;
	gs_max_d2 = 0;
	}
	
 if (Leaf_rain)                                 // evaporate first the water on the leaves until dry
	{
	// set all E to 0
	E_Branch=0; E_Trunk=0; E_Axil=0; E_cuti=0; E_Leaf=0; E_clim=0; E_Petiole=0;
	// set all dq to 0
	dq_Branch=0; dq_Axil=0; dq_Petiole=0;dq_cuti=0; dq_stomate=0; dq_fruit=0; dq_Trunk=0;
	if (DYNAMIC==0) Leaf_rain-= ETP_Penman*18/1000/1000*dt;  //water loss by ETP in mm during dt_stat;
	else Leaf_rain-= ETP_Penman*18/1000/1000*dt_long*dt;
	if (Leaf_rain<0) Leaf_rain=0;
	}

}

void soil(long double time)   //pedotransfer functions
{
	RWC1=(Q_Soil1/1000/1000/1000*18/Volume_soil1- Teta_r)/(Teta_s - Teta_r);         // Relative volumetric water content
	if(RWC1>RWC_fc) RWC1=RWC_fc;
	if (RWC1) P_Soil1=-(pow(((pow((1/RWC1),(1/m)))-1),(1/n)))/alpha/10000;  else P_Soil1=-9999;
	K_s=1000*K_sat*2*3.1416*Length_Root_fi*Root_upper/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_upper/Volume_soil1,0.5)/(Diam_Root/2))*Fluidity_soil; //based on total root system
	K_Soil1=K_s*pow(RWC1,L)*pow(1-pow(1-pow(RWC1,1/m),m),2);
	if (K_Soil1) K_Interface1=K_Soil1*10*pow(Q_Root_Symp1/Q_Root_Symp01,gap);    // root-soil interface Conductance linked to root diameter
	else         K_Interface1=0;
	K_Soil_tot1=1000*K_sat/Soil_Depth*3*pow(RWC1,L)*pow(1-pow(1-pow(RWC1,1/m),m),2)*Fluidity_soil;
	
	RWC2=(Q_Soil2/1000/1000/1000*18/Volume_soil2- Teta_r)/(Teta_s - Teta_r);
	if(RWC2>RWC_fc) RWC2=RWC_fc;
	if (RWC2) P_Soil2=-(pow(((pow((1/RWC2),(1/m)))-1),(1/n)))/alpha/10000; else P_Soil2=-9999;
	K_s=1000*K_sat*2*3.1416*Length_Root_fi*Root_middle/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_middle/Volume_soil2,0.5)/(Diam_Root/2))*Fluidity_soil;
	K_Soil2=K_s*pow(RWC2,L)*pow(1-pow(1-pow(RWC2,1/m),m),2);
	if (K_Soil2)  K_Interface2=K_Soil2*10*pow(Q_Root_Symp2/Q_Root_Symp02,gap);
	else          K_Interface2 =0;
	K_Soil_tot2=1000*K_sat/Soil_Depth*3*pow(RWC2,L)*pow(1-pow(1-pow(RWC2,1/m),m),2)*Fluidity_soil;
	
	RWC3=(Q_Soil3/1000/1000/1000*18/Volume_soil3- Teta_r)/(Teta_s - Teta_r);
	if(RWC3>RWC_fc) RWC3=RWC_fc;
	if (RWC3) P_Soil3=-(pow(((pow((1/RWC3),(1/m)))-1),(1/n)))/alpha/10000; else P_Soil3=-9999;
	K_s=1000*K_sat*2*3.1416*Length_Root_fi*Root_lower/Surface_Soil/log(1/pow(3.1416*Length_Root_fi*Root_lower/Volume_soil3,0.5)/(Diam_Root/2))*Fluidity_soil;
	K_Soil3=K_s*pow(RWC3,L)*pow(1-pow(1-pow(RWC3,1/m),m),2);
	if (K_Soil3)  K_Interface3=K_Soil3*10*pow(Q_Root_Symp3/Q_Root_Symp03,gap);
	else          K_Interface3 =0;
	K_Soil_tot3=1000*K_sat/Soil_Depth*3*pow(RWC3,L)*pow(1-pow(1-pow(RWC3,1/m),m),2)*Fluidity_soil;
	
	Teta_Soil=(Q_Soil1+Q_Soil2+Q_Soil3)/1000/1000/1000*18/Volume_soil;
	RWC_Soil=(Teta_Soil- Teta_r)/(Teta_s - Teta_r); //mean soil RWC

	REW_t=   (Teta_Soil- Teta_r)/(Teta_fc-Teta_r);  //REW based on Teta_r
	if (RWC_Soil<RWC_min) RWC_min=RWC_Soil;  // min soil RWC
	if (RWC_Soil>RWC_fc) RWC_Soil=RWC_fc;
	RWC_int+=(1-RWC_Soil)*time/3600/24; // cumulated RWC deficit
	
	if (RWC_Soil>0) P_soil=-(pow(((pow((1/RWC_Soil),(1/m)))-1),(1/n)))/alpha/10000; else P_soil=-9999;//mean soil water potential
	P_soil_int +=P_soil*time/3600/24;
	if (P_soil<P_soil_min) P_soil_min=P_soil; // min soil WP
}



void Compute_Growth (void)  //under development. Lockhart model
{
	long double T1,T2,T3,sigma;
	
	// effect of temperature on wall extensibility: max at T2, zero at T1 and T3
	T1=5;
	T2=40;
	T3=50;
	if (T_air<T1) sigma=0;                          // a test to include a temperature dependance of growth rate
	else if (T_air<T2) sigma=(T_air-T1)/(T2-T1);
	else if (T_air<T3) sigma=(T_air-T3)/(T2-T3);
	else sigma=0;
	
	//sigma=1;
	
	if (Turgor_Trunk_Symp > Yield_trunk)  Growth_rate_trunk = Q_Trunk_Symp0*sigma*Extensibility_trunk * (Turgor_Trunk_Symp - Yield_trunk);  // volume growth rate of the branch symplasm in mmol/s
	else                                  Growth_rate_trunk = 0;
	if (Turgor_Axil_Symp > Yield_fruit)   Growth_rate_fruit = Q_Axil_Symp0*sigma*Extensibility_fruit * (Turgor_Axil_Symp - Yield_fruit);  // volume growth rate of the fruit symplasm in mmol/s
	else                                  Growth_rate_fruit = 0;
}

void  Compute_g_cuti(void)  //temperature dependance of g_cuti
{
	long double g_cuti_tp, g_Axil_tp;
	
	g_cuti_tp= g_cuti_20*pow(Q10_1,(TP-20)/10);
	g_Axil_tp=g_Axil_min20*pow(Q10_1,(TP-20)/10);
	if (T_Leaf<TP)
	{
		g_cuti= g_cuti_20*pow(Q10_1,(T_Leaf-20)/10);
		g_Axil_min= g_Axil_min20*pow(Q10_1,(T_Leaf-20)/10);
	}
	else
	{
		g_cuti= g_cuti_tp*pow(Q10_2,(T_Leaf-TP)/10);
		g_Axil_min= g_Axil_tp*pow(Q10_2,(T_Leaf-TP)/10);
	}
	
}

void Compute_C (void)   //The apoplasmic capacitance if function of the apoplaspic water content;
{
	long double a,b;
	a=0.0 ;  // fraction to prevent null C at 100 PLC
	b=1;
	
	C_Leaf_Apo    = b*C_Leaf_Apo0   *(a+(1-a)*(100-PLC_Leaf_Apo)   /100);
	C_Branch_Apo  = b*C_Branch_Apo0 *(a+(1-a)*(100-PLC_Branch_Apo) /100);
	C_Trunk_Apo   = b*C_Trunk_Apo0  *(a+(1-a)*(100-PLC_Trunk_Apo)  /100);
	C_Root_Apo1   = b*C_Root_Apo01  *(a+(1-a)*(100-PLC_Root_Apo1)  /100);
	C_Root_Apo2   = b*C_Root_Apo02  *(a+(1-a)*(100-PLC_Root_Apo2)  /100);
	C_Root_Apo3   = b*C_Root_Apo03  *(a+(1-a)*(100-PLC_Root_Apo3)  /100);
}

void Compute_Q_steady(void)
{
	long double Rs, Qold;
	//Leaf
	Qold=Q_Leaf_Symp;
	Rs=(-(P_Leaf_Symp+Pi0_Leaf_Symp-Epsilon_Leaf_Symp)-pow((pow((P_Leaf_Symp+Pi0_Leaf_Symp-Epsilon_Leaf_Symp),2)+ 4*Epsilon_Leaf_Symp*P_Leaf_Symp),0.5))/(2*Epsilon_Leaf_Symp);
	if (Rs<(1-Pi0_Leaf_Symp/P_Leaf_Symp)) Rs=1-Pi0_Leaf_Symp/P_Leaf_Symp;
	Turgor_Leaf_Symp=-Pi0_Leaf_Symp*Osmotic_TLeaf - Epsilon_Leaf_Symp*Rs;
	if (Turgor_Leaf_Symp<0) Turgor_Leaf_Symp=0;
	Q_Leaf_Symp=Q_Leaf_Symp0*(1-Rs);
	Q_Soil3+= Qold-Q_Leaf_Symp;
	
	//Branch
	Qold=Q_Branch_Symp;
	Rs=(-(P_Branch_Symp+Pi0_Branch_Symp-Epsilon_Branch_Symp)-pow((pow((P_Branch_Symp+Pi0_Branch_Symp-Epsilon_Branch_Symp),2)+ 4*Epsilon_Branch_Symp*P_Branch_Symp),0.5))/(2*Epsilon_Branch_Symp);
	if (Rs<(1-Pi0_Branch_Symp/P_Branch_Symp)) Rs=1-Pi0_Branch_Symp/P_Branch_Symp;
	Q_Branch_Symp=Q_Branch_Symp0*(1-Rs);
	Q_Soil3+= Qold-Q_Branch_Symp;
	
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
		if (Type_Axil>=2) // a flower or a fruit
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
	long double F_Leaf, F_Branch, F_Trunk, F_Axil, F_Petiole, F_Cuti, F_Root1,F_Root2,F_Root3;
	long double LIMIT=1e-5;
	
	if (K_Soil1 && K_Interface1 && K_Root_Apo1 && K_Root_Symp11) K_Root1=1/(1/K_Soil1+1/K_Interface1+1/K_Root_Symp11+1/K_Root_Apo1);  else K_Root1=LIMIT;
	if (K_Soil2 && K_Interface2 && K_Root_Apo2 && K_Root_Symp12) K_Root2=1/(1/K_Soil2+1/K_Interface2+1/K_Root_Symp12+1/K_Root_Apo2);  else K_Root2=LIMIT;
	if (K_Soil3 && K_Interface3 && K_Root_Apo3 && K_Root_Symp13) K_Root3=1/(1/K_Soil3+1/K_Interface3+1/K_Root_Symp13+1/K_Root_Apo3);  else K_Root3=LIMIT;
	
	if (K_Soil1 && K_Interface1 && K_Root_Symp11) K_Root1b=1/(1/K_Soil1+1/K_Interface1+1/K_Root_Symp11); else K_Root1b=LIMIT;
	if (K_Soil2 && K_Interface2 && K_Root_Symp12) K_Root2b=1/(1/K_Soil2+1/K_Interface2+1/K_Root_Symp12); else K_Root2b=LIMIT;
	if (K_Soil3 && K_Interface3 && K_Root_Symp13) K_Root3b=1/(1/K_Soil3+1/K_Interface3+1/K_Root_Symp13); else K_Root3b=LIMIT;
	
	F_Leaf=  (E_Leaf-E_cuti)*Leaf_Area; //flow only through stomata
	F_Cuti=   E_cuti*Leaf_Area;
	F_Branch= E_Branch*Branch_Area; // for all branches
	F_Trunk=  E_Trunk*Trunk_Area;

	if (Type_Axil) F_Axil=   E_Axil*Area_Axil;  else F_Axil=0;   // for all flowers
	if (Type_Axil>=2) F_Petiole=E_Petiole*Petiole_area; else F_Petiole=0; //for one petioles
	
	FE_Root1= E_Root1*Root_Area1;
	FE_Root2= E_Root2*Root_Area2;
	FE_Root3= E_Root3*Root_Area3;
	
	P_Root_Apo1= -((F_Leaf+F_Branch+F_Trunk+F_Axil+F_Petiole+F_Cuti) - K_Root1*P_Soil1  - K_Root2*P_Soil2 - K_Root3*P_Soil3)/(K_Root1+K_Root2+K_Root3);
	P_Root_Apo2=P_Root_Apo1;
	P_Root_Apo3=P_Root_Apo1;
	F_Root1= (P_Soil1-P_Root_Apo1)*K_Root1;
	F_Root2= (P_Soil2-P_Root_Apo2)*K_Root2;
	F_Root3= (P_Soil3-P_Root_Apo3)*K_Root3;
	
	Q_Soil1-=F_Root1*dt_court;
	Q_Soil2-=F_Root2*dt_court;
	Q_Soil3-=F_Root3*dt_court;
	
	P_Root_Endo1= P_Soil1 - F_Root1/K_Root1b;
	P_Root_Endo2= P_Soil2 - F_Root2/K_Root2b;
	P_Root_Endo3= P_Soil3 - F_Root3/K_Root3b;
	
	if (K_Root_Symp21)  P_Root_Symp1= P_Root_Endo1 -FE_Root1/(K_Root_Symp21); else P_Root_Symp1= P_Root_Endo1;
	if (K_Root_Symp22)  P_Root_Symp2= P_Root_Endo2 -FE_Root2/(K_Root_Symp22); else P_Root_Symp2= P_Root_Endo2;
	if (K_Root_Symp23)  P_Root_Symp3= P_Root_Endo3 -FE_Root3/(K_Root_Symp23); else P_Root_Symp3= P_Root_Endo3;
	
	if (K_Trunk_Apo) P_Trunk_Apo= Pg_trunk+ P_Root_Apo1-(F_Leaf+F_Branch+F_Trunk+F_Axil+F_Petiole+F_Cuti)/K_Trunk_Apo; else P_Trunk_Apo=P_air;
					 P_Trunk_Symp=P_Trunk_Apo-F_Trunk/K_Trunk_Symp;

	if (K_Branch_Apo) P_Branch_Apo=Pg_branch- Pg_trunk + P_Trunk_Apo-(F_Leaf+F_Cuti+F_Branch+F_Axil+F_Petiole)/K_Branch_Apo; else P_Branch_Apo=P_air;
					  P_Branch_Symp=P_Branch_Apo-F_Branch/K_Branch_Symp;
	
	if (Type_Axil)
	{
		if (K_Axil_Apo) P_Axil_Apo=P_Branch_Apo-(F_Axil+F_Petiole)/(K_Axil_Apo); else P_Axil_Apo=P_air; // no symplasmic connexion here !
		P_Axil_Symp=P_Axil_Apo-F_Axil/(K_Axil_Symp);
		P_Petiole_Symp=P_Axil_Apo-F_Petiole/(K_Axil_Symp2);
	}
	
	if (K_Leaf_Apo) P_Leaf_Apo=Pg_Leaf-Pg_branch + P_Branch_Apo-(F_Leaf+F_Cuti)/(K_Leaf_Apo); else P_Leaf_Apo=P_air;
	if (K_Leaf_Symp) P_Leaf_Evap=P_Leaf_Apo-(F_Leaf+F_Cuti)/(K_Leaf_Symp); else P_Leaf_Evap=P_Leaf_Apo;
	if (K_Leaf_Symp2) P_Leaf_Symp=P_Leaf_Evap-F_Cuti/(K_Leaf_Symp2); else P_Leaf_Symp=P_Leaf_Apo;
		
}



void Fluidity (void)   // compute water fluidity variation with temperature
{
	if (FLUID)
	{
		Fluidity_leaf =   1.01212E-04*T_Leaf * T_Leaf  + 2.04152E-02*T_Leaf + 5.51781E-01;
		Fluidity_air  =   1.01212E-04*T_air  * T_air   + 2.04152E-02*T_air  + 5.51781E-01;
		Fluidity_soil  =  1.01212E-04*T_Soil * T_Soil  + 2.04152E-02*T_Soil + 5.51781E-01;
		Fluidity_soil1  = 1.01212E-04*T_Soil1* T_Soil1 + 2.04152E-02*T_Soil1+ 5.51781E-01;
		Fluidity_soil2  = 1.01212E-04*T_Soil2* T_Soil2 + 2.04152E-02*T_Soil2+ 5.51781E-01;
		Fluidity_soil3  = 1.01212E-04*T_Soil3* T_Soil3 + 2.04152E-02*T_Soil3+ 5.51781E-01;
	}
	else
	{
		Fluidity_leaf = 1;
		Fluidity_air  = 1;
		Fluidity_soil = 1;
	}
}

void Compute_K(void) // changes of K due to cavitation
{
	long double K,Cortex_Gap,R;
	K_Leaf_Apo    = K_Leaf_Apo0       *(100-PLC_Leaf_Apo)/100   * (100-PLF)/100 * Fluidity_leaf;
	K_Axil_Apo    = K_Axil_Apo0       *(100-PLC_Axil_Apo)/100   * Fluidity_air;
	K_Branch_Apo  = K_Branch_Apo0     *(100-PLC_Branch_Apo)/100 * Fluidity_air;
	K_Trunk_Apo   = K_Trunk_Apo0      *(100-PLC_Trunk_Apo)/100  * Fluidity_air;
	K_Root_Apo1   = K_Root_Apo01      *(100-PLC_Root_Apo1)/100  * Fluidity_soil ;  // T_soil assumed constant
	K_Root_Apo2   = K_Root_Apo02      *(100-PLC_Root_Apo2)/100  * Fluidity_soil ;
	K_Root_Apo3   = K_Root_Apo03      *(100-PLC_Root_Apo3)/100  * Fluidity_soil ;
	K_Branch_Symp = K_Branch_Symp0 * Fluidity_air;
	K_Trunk_Symp  = K_Trunk_Symp0  * Fluidity_air;
	
	// Leaf
	if (K_VAR==1)  // K_leaf is variable
	{
		//K= 6.83 + 81.4* exp (7.56 * P_Leaf_Symp);
		K= K_VAR_P1 + K_VAR_P2* exp (K_VAR_P3 * P_Leaf_Symp) ;
		if (Leaf_Area) K_Leaf_Symp= 1/(1/K -1/(K_Leaf_Apo0/Leaf_Area)) * Fluidity_leaf; else K_Leaf_Symp=99999;
		if (K_Leaf_Symp > 106) K_Leaf_Symp=106;
		K_Leaf_Symp*=Leaf_Area;
	}
	else
	{
		K_Leaf_Symp  = K_Leaf_Symp_0/LA_max*Leaf_Area*(100-PLF)/100 * Fluidity_leaf;
		K_Leaf_Symp2 = K_Leaf_Symp;
	}
	
	// Root
	if (K_VAR>1)  //K is variable
	{
		if (K_VAR==2) // Arabidopsis Scoffoni
		{
			//K=  4.25+495*exp(5.85*P_Leaf_Symp);
			 K=K_VAR_P1+K_VAR_P2*exp(K_VAR_P3*P_Leaf_Symp);
			K*=LA_max;
			K_Root_Symp1= 1/(1/K -1/(K_Root_Apo0)) ;
			if (K_Root_Symp1>1.6) K_Root_Symp1=1.6;
		}
		
	   else if (K_VAR==3) // Effet température sol pour Doux-Glace
		{
			if (T_Soil1 >20) K_Root_Symp11 = K_Root_Symp0 * Fluidity_soil1 * Root_Area_fi * Root_upper;                  //  follows fluidity
			else if (T_Soil1<T_Soil_Crit) K_Root_Symp11 = 0;                                                             //  zero when below T_Soil_Crit
			else K_Root_Symp11 = K_Root_Symp0 * Root_Area_fi * Root_upper * T_Soil1/(20-T_Soil_Crit);                    //  between T_Soil_Crit and 20 linear with T_soil
			
			if (T_Soil2 >20) K_Root_Symp12 = K_Root_Symp0 * Fluidity_soil2 * Root_Area_fi * Root_middle;
			else if (T_Soil2<T_Soil_Crit) K_Root_Symp12 = 0;
			else K_Root_Symp12 = K_Root_Symp0 * Root_Area_fi * Root_middle * T_Soil2/(20-T_Soil_Crit);
			
			if (T_Soil3 >20) K_Root_Symp13 = K_Root_Symp0 * Fluidity_soil3 * Root_Area_fi * Root_lower;
			else if (T_Soil3<T_Soil_Crit) K_Root_Symp13 = 0;
			else K_Root_Symp13 = K_Root_Symp0 * Root_Area_fi * Root_lower * T_Soil3/(20-T_Soil_Crit);
			
		}
		
	  else  if (K_VAR==4) // air gap in the cortex following a sigmoidal curve K_VAR_P1=P50 K_VAR_P2=slope with P_root_symp; non reversible loss of K
		{
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp1-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi*Root_upper*(100-Cortex_Gap)/100;
			if (K<K_Root_Symp11)  K_Root_Symp11=K;
			//K_Root_Symp21= 1 * K_Root_Symp11;
			
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp2-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi* Root_middle *(100-Cortex_Gap)/100;
			if (K<K_Root_Symp12)  K_Root_Symp12=K;
			//K_Root_Symp22= 1 * K_Root_Symp12;
			
			Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Root_Symp3-K_VAR_P1)));
			K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi*Root_lower*(100-Cortex_Gap)/100;
			if (K<K_Root_Symp13)  K_Root_Symp13=K;
			//K_Root_Symp23= 1 * K_Root_Symp13;
		}
		
	  else  if (K_VAR==5) // For Maddy and Tim
	  {
		  R=pow(K_Root_Symp0/K_VAR_P1,-1/K_VAR_P2);
		  K_Root_Symp11=K_VAR_P1*pow(R-P_Soil1,K_VAR_P2)*Fluidity_soil * Root_Area_fi* Root_upper;
		  K_Root_Symp12=K_VAR_P1*pow(R-P_Soil2,K_VAR_P2)*Fluidity_soil * Root_Area_fi* Root_middle;
		  K_Root_Symp13=K_VAR_P1*pow(R-P_Soil3,K_VAR_P2)*Fluidity_soil * Root_Area_fi* Root_lower;
	  }
	  
	  else  if (K_VAR==6) // air gap in the cortex following a sigmoidal curve with P_soil; non reversible loss of K
	  {
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil1-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi*Root_upper*(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp11)  K_Root_Symp11=K;
		 // K_Root_Symp21= 1 * K_Root_Symp11;
		  
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil2-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi* Root_middle *(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp12)  K_Root_Symp12=K;
		//  K_Root_Symp22= 1 * K_Root_Symp12;
		  
		  Cortex_Gap=100/(1+exp(K_VAR_P2/25*(P_Soil3-K_VAR_P1)));
		  K = K_Root_Symp0 * Fluidity_soil * Root_Area_fi*Root_lower*(100-Cortex_Gap)/100;
		  if (K<K_Root_Symp13)  K_Root_Symp13=K;
	   //   K_Root_Symp23= 1 * K_Root_Symp13;
	  }
	
	else if (K_VAR==7) //K_branch_symp variable with PLC
	{
		K_Branch_Symp=K_Branch_Symp0*(100-PLC_Branch_Apo)/100 * Fluidity_air;
		K_Trunk_Symp=K_Trunk_Symp0*(100-PLC_Trunk_Apo)/100 * Fluidity_air;
	}
	
	else if (K_VAR==8) //For Carola and Tim with K_VAR_P1=P50  K_VAR_P2=slope K_VAR_P3=Kmax fraction to have Kmin
	{
		K_Root_Symp11=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp1-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil * Root_Area_fi* Root_upper;
		K_Root_Symp12=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp2-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil * Root_Area_fi* Root_middle;
		K_Root_Symp13=((K_Root_Symp0*(1-K_VAR_P3))*(100-100/(1+exp(K_VAR_P2/25*(P_Root_Symp3-K_VAR_P1))))/100+K_Root_Symp0*K_VAR_P3)*Fluidity_soil * Root_Area_fi* Root_lower;
	}
	
	else // K_root is constant
	{
		K_Root_Symp1 =  K_Root_Symp0 * Fluidity_soil * Root_Area_fi; //follows fluidity
		K_Root_Symp11=  K_Root_Symp1 * Root_upper;
		K_Root_Symp12=  K_Root_Symp1 * Root_middle;
		K_Root_Symp13=  K_Root_Symp1 * Root_lower;
		K_Root_Symp21=  K_Root_Symp2 * Root_upper;
		K_Root_Symp22=  K_Root_Symp2 * Root_middle;
		K_Root_Symp23=  K_Root_Symp2 * Root_lower;
	}
   }
}

void Compute_Leaf_Fall(void)   //force leaves to fall according to leaf water potential
{
	long double LA, PLF_i;
	
	PLF_i=100/(1+exp(Slope_Leaf_Fall/25*(P_Leaf_Symp-P50_Leaf_Fall)));
	LA = LA_max_Pheno *(100-PLF_i)/100;
	if (LA==0) LA=0.001;
	if (LA<Leaf_Area)
	{
		if (Leaf_Area)  Q_Leaf_Apo1*=(LA/Leaf_Area);
		if (Leaf_Area)  Q_Leaf_Apo*=(LA/Leaf_Area);
		if (Leaf_Area)  Q_Leaf_Symp0*=(LA/Leaf_Area);
		if (Leaf_Area)  Q_Leaf_Symp*=(LA/Leaf_Area);
		if (Leaf_Area)  Q_Leaf_Evap*=(LA/Leaf_Area);
		if (Leaf_Area)  Q_Leaf_Evap0*=(LA/Leaf_Area);
		Leaf_Area=LA;
		PLF=PLF_i;
	}
}

void Compute_Root_Fall(void)   //force root to die according to root water potential
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

void Phenology (void)  //computes seasonal changes in leaf area
{
	//if (LA_min==0) LA_min=0.0001;
	if (T/24/3600<LA_day1) LA_max_Pheno=LA_min;
	else if (T/24/3600>LA_day1 && T/24/3600<LA_day2) LA_max_Pheno=LA_min + (LA_max - LA_min)*(T/24/3600-LA_day1)/(LA_day2-LA_day1);
	else if (T/24/3600>LA_day2 && T/24/3600<LA_day3) LA_max_Pheno=LA_max;
	else if (T/24/3600>LA_day3 && T/24/3600<LA_day4) LA_max_Pheno=LA_max- (LA_max - LA_min)*(T/24/3600-LA_day3)/(LA_day4-LA_day3);
	else if (T/24/3600>LA_day4) LA_max_Pheno=LA_min;
	
	Leaf_Area=LA_max_Pheno;
	if (Leaf_Area) K_Leaf_Symp=K_Leaf_Symp_0/LA_max*Leaf_Area; else K_Leaf_Symp=10;
	if (Leaf_Area) K_Leaf_Apo=K_Leaf_Apo0/LA_max*Leaf_Area;    else K_Leaf_Apo=10;
	if (Leaf_Area) K_Leaf_Symp2=K_Leaf_Symp;                   else K_Leaf_Symp2=10;
	LAI=Leaf_Area/Surface_Soil;
}

void Compute_T_Osmotic(void) // temperature dependance of osmotic potentials
{
	if (T_OSMOTIC)
	{
		Osmotic_TLeaf =  (T_Leaf + 273.16)/293.16;
		Osmotic_TAir  =  (T_air  + 273.16)/293.16;
		Osmotic_TSoil=   (T_Soil + 273.16)/293.16;
	}
	else
	{
		Osmotic_TLeaf = 1.0;
		Osmotic_TAir  = 1.0;
		Osmotic_TSoil = 1.0;
	}
}

void Compute_ST(void)  // temperature dependance of water surface tension
{
	if (SURFACE_TENSION)
	{
		ST_Leaf =  (75.6986 - (2.64569E-4) * T_Leaf * T_Leaf - 1.42361E-1 * T_Leaf )/72.7455;
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

void Compute_Cavitation(void)
{
	long double PLC=0.0,dQ,KP;

   	
	Compute_ST();
	//LEAF
	PLC=100/(1+exp(Slope_Leaf_Apo/25*(P_Leaf_Apo-P50_Leaf_Apo*ST_Leaf)));
	if (Regul_gs ==6) Px_Leaf_Apo= P50_Leaf_Apo*ST_Leaf + 25/Slope_Leaf_Apo*log((100-PLCx)/PLCx); 
			
	if ((PLC-PLC_Leaf_Apo)>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Leaf_Apo>P_REFILL) || (PLC>PLC_Leaf_Apo))
	{
		dQ=(PLC-PLC_Leaf_Apo)/100*Q_Leaf_Apo0;
		if (SYMP_CAVIT) 
		{
			if (DYNAMIC && dQ>0) Q_Leaf_Symp+=dQ; 
		}
		else if (DYNAMIC && dQ>0)     // water released is distributed to the connected symplasmic reservoirs except for Refilling (dQ<0)
		{
		 KP=K_Leaf_Symp*P_Leaf_Symp + K_Leaf_Apo*P_Branch_Apo;
		 if (KP)
			{
			Q_Leaf_Symp+=     dQ*(K_Leaf_Symp*P_Leaf_Symp)/KP;  //water released is distributed to the connected reservoirs according to the F=K*P ratio
			Q_Branch_Apo+=    dQ*(K_Leaf_Apo*P_Branch_Apo)/KP;
			}
		}
		else Q_Soil3+=  dQ;                         // water released is distributed the deeper soil layer
		Q_Leaf_Apo1-=   dQ;
		Q_Leaf_Apo-=    dQ;
		PLC_Leaf_Apo=   PLC;
	}
	
	
	//Axil
	if (Type_Axil)
	{
	PLC=100/(1+exp(Slope_Axil_Apo/25*(P_Axil_Apo-P50_Axil_Apo*ST_Air)));
	if ((PLC-PLC_Axil_Apo)>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Axil_Apo>P_REFILL) || (PLC>PLC_Axil_Apo))
	{
		dQ=(PLC-PLC_Axil_Apo)/100*Q_Axil_Apo0;
		if (DYNAMIC && dQ>0)
		{
			KP=K_Axil_Symp*P_Axil_Symp  + N_Axil *K_Axil_Apo*P_Branch_Apo;                  // the water released by cavitation is distributed
			if (KP)
			{
				Q_Branch_Apo+=      N_Axil * dQ*(K_Axil_Apo*P_Branch_Apo)/KP;                                            // according to the relative potential flows = K*pressure gradient
				if (Type_Axil>=2) Q_Petiole_Symp+=dQ*(K_Axil_Symp*P_Axil_Symp)/KP;                                     // into the connected apo and symplasmic compartements
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
	PLC=100/(1+exp(Slope_Branch_Apo/25*(P_Branch_Apo-P50_Branch_Apo*ST_Air)));
	if ((PLC-PLC_Branch_Apo)>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Branch_Apo>P_REFILL) || (PLC>PLC_Branch_Apo))
	{
		dQ=(PLC-PLC_Branch_Apo)/100*Q_Branch_Apo0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0)  Q_Branch_Symp+=dQ;  
		else if (DYNAMIC && dQ>0)  // into the connected apo and symplasmic compartements
		{                                
			if (Type_Axil) KP=K_Branch_Symp*P_Branch_Symp + K_Leaf_Apo*P_Leaf_Apo+ K_Axil_Apo*P_Axil_Apo+ K_Branch_Apo*P_Trunk_Apo;   // the water released by cavitation is distributed
			else KP=K_Branch_Symp*P_Branch_Symp + K_Leaf_Apo*P_Leaf_Apo+K_Branch_Apo*P_Trunk_Apo;   // the water released by cavitation is distributed
			if (KP)
			{
				Q_Branch_Symp+=   dQ*K_Branch_Symp*P_Branch_Symp/KP ;                                     // into the connected apo and symplasmic compartements
				Q_Leaf_Apo+=      dQ*K_Leaf_Apo*P_Leaf_Apo/KP;                                            // according to the relative potential flows = K*pressure gradient
				Q_Trunk_Apo+=     dQ*K_Branch_Apo*P_Trunk_Apo/KP;
				if (Type_Axil) 	Q_Axil_Apo+=      dQ*K_Axil_Apo*P_Axil_Apo/KP;
			}
		}
		else Q_Soil3+=dQ;
		Q_Branch_Apo-=dQ;
		Q_Branch_Apo1-=dQ;
		PLC_Branch_Apo=PLC;
	}
	
	//TRUNK
	PLC=100/(1+exp(Slope_Trunk_Apo/25*(P_Trunk_Apo-P50_Trunk_Apo*ST_Air)));
	if ((PLC-PLC_Trunk_Apo)>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Trunk_Apo>P_REFILL) || (PLC>PLC_Trunk_Apo))
	{
		dQ=(PLC-PLC_Trunk_Apo)/100*Q_Trunk_Apo0;
		if (SYMP_CAVIT && DYNAMIC && dQ>0) Q_Trunk_Symp+=dQ;
		else if (DYNAMIC && dQ>0) 
		{
			KP=K_Trunk_Symp*P_Trunk_Symp + K_Branch_Apo*P_Branch_Apo + K_Trunk_Apo*P_Root_Apo;
			if (KP)
			{
			Q_Trunk_Symp+=  	dQ*K_Trunk_Symp*P_Trunk_Symp/KP;
			Q_Branch_Apo+= 	dQ*K_Branch_Apo*P_Branch_Apo/KP;
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
	if ((PLC-PLC_Root_Apo1)>dPLC_crit) PLC_LIMIT=1;
	if ((REFILL && P_Root_Apo1>P_REFILL) || (PLC>PLC_Root_Apo1))
	{
		dQ=(PLC-PLC_Root_Apo1)/100*Q_Root_Apo01;
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
	if ((REFILL && P_Root_Apo2>P_REFILL) || (PLC>PLC_Root_Apo2))
	{
		dQ=(PLC-PLC_Root_Apo2)/100*Q_Root_Apo02;
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
	if ((REFILL && P_Root_Apo3>P_REFILL) || (PLC>PLC_Root_Apo3))
	{
		dQ=(PLC-PLC_Root_Apo3)/100*Q_Root_Apo03;
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
	
  //  Get_DATA(dt);
	if (PLC_Leaf_Apo>PLC_END     && T_PLC_Leaf>(T-T0))      T_PLC_Leaf= T-T0;                                            // time to leaf hydraulic failure
	if (PLC_Branch_Apo>PLC_END   && T_PLC_Branch>(T-T0))    T_PLC_Branch= T-T0;                                        // time to branch hydraulic failure
	if (PLC_Trunk_Apo>PLC_END    && T_PLC_Trunk>(T-T0))     T_PLC_Trunk= T-T0;                                          // time to trunk hydraulic failure
	if ((PLC_Root_Apo1>PLC_END   && PLC_Root_Apo2>PLC_END && PLC_Root_Apo3>PLC_END) && T_PLC_Root>(T-T0)) T_PLC_Root= T-T0;   // time to ALL roots hydraulic failure
	if (PLC_Root_Apo1>PLC_END  && T_PLC_Root>(T-T0)) T_PLC_Root1= T-T0;
	if (PLC_Axil_Apo>PLC_END     && T_PLC_Axil>(T-T0))      T_PLC_Axil= T-T0;
	  
}

void init_Cavitation(void)
{
	long double PLC=0.0;
	Compute_ST();
	PLC = 100/(1+exp(Slope_Leaf_Apo  /25*(P_Leaf_Apo  -P50_Leaf_Apo  * ST_Leaf)));  PLC_Leaf_Apo  = PLC;
	PLC = 100/(1+exp(Slope_Axil_Apo  /25*(P_Axil_Apo  -P50_Axil_Apo  * ST_Air )));  PLC_Axil_Apo  = PLC;
	PLC = 100/(1+exp(Slope_Branch_Apo/25*(P_Branch_Apo-P50_Branch_Apo* ST_Air )));  PLC_Branch_Apo= PLC;
	PLC = 100/(1+exp(Slope_Trunk_Apo /25*(P_Trunk_Apo -P50_Trunk_Apo * ST_Air )));  PLC_Trunk_Apo = PLC;
	PLC = 100/(1+exp(Slope_Root_Apo1 /25*(P_Root_Apo1 -P50_Root_Apo1 * ST_Soil)));  PLC_Root_Apo1 = PLC;
	PLC = 100/(1+exp(Slope_Root_Apo2 /25*(P_Root_Apo2 -P50_Root_Apo2 * ST_Soil)));  PLC_Root_Apo2 = PLC;
	PLC = 100/(1+exp(Slope_Root_Apo3 /25*(P_Root_Apo3 -P50_Root_Apo3 * ST_Soil)));  PLC_Root_Apo3 = PLC;
}

void Compute_Turgor_Ref(void)
{
	// compute the midday leaf turgor for a well watered plant to use this value to compute the regulation of E by leaf turgor;
	// assume no drougth no cavitation and an ETP=4mm and 25°C
	long double ETP_ref0 = 4;  // 5mm per day
	long double ETP_ref;     // ETP in mmol/s/m2 of Leaf area
	long double K,K_Root;           // whole plant hydraulic LS conductance in mmol/s/m2/MPA
	long double LWP;         //Leaf Water Potential in MPa
	long double tlp, discri, Rs,Tp,Q=0;
	long double T1,T2,T3;
	
	T1=T_air;
	T2=T_Leaf;
	T3=T_Soil;
	T_air=25;
	T_Leaf=25;
	T_Soil=20;
	LAI=Leaf_Area/Surface_Soil;   // LAI in m2/m2
	if (Penman_Coeff==0) Penman_Coeff= -0.006*LAI*LAI +0.134*LAI + 0.036;
	if (Leaf_Area) ETP_ref= (ETP_ref0*1.6106-0.5616)*Soil_Width*Soil_Width*Penman_Coeff/Leaf_Area;  else ETP_ref=0;//ETP at midday in mmol/s/m2
	Fluidity();
	Compute_K();
	Compute_T_Osmotic();
	K_Root=(1/(1/K_Root_Symp11+1/K_Root_Apo1) + 1/(1/K_Root_Symp12+1/K_Root_Apo2) + 1/(1/K_Root_Symp13+1/K_Root_Apo3));
	if (Leaf_Area) K_tot= (1/( 1/K_Leaf_Apo +1/K_Leaf_Symp + 1/K_Branch_Apo+ 1/K_Trunk_Apo + 1/K_Root)/Leaf_Area);
	else K=0;
	LWP=-ETP_ref/K_tot;
	if (GRAVITY==1) LWP-=(Length_Trunk+Length_Branch)*0.01;
	
	tlp=Pi0_Leaf_Symp*Osmotic_TLeaf*Epsilon_Leaf_Symp/(Pi0_Leaf_Symp*Osmotic_TLeaf+Epsilon_Leaf_Symp);
	if (LWP>tlp)
	{
		discri=pow(Pi0_Leaf_Symp*Osmotic_TLeaf + Epsilon_Leaf_Symp + LWP,2)-4*Pi0_Leaf_Symp*Osmotic_TLeaf*Epsilon_Leaf_Symp;
		Q = ((Pi0_Leaf_Symp*Osmotic_TLeaf + Epsilon_Leaf_Symp + LWP) + pow(discri,0.5))/(2*Epsilon_Leaf_Symp/Q_Leaf_Symp0);
	}
	else   Turgor_Leaf_Symp_Ref=0;
	
	Rs=(Q_Leaf_Symp0-Q)/Q_Leaf_Symp0;
	Tp=-Pi0_Leaf_Symp*Osmotic_TLeaf - Epsilon_Leaf_Symp*Rs;
	if (Tp<0)Tp=0;
	Turgor_Leaf_Symp_Ref =Tp;
	if (PRINT_SCREEN) printf("LWP=%Lf TLP=%Lf Turgor_ref=%Lf", LWP, tlp, Tp);
	T_air=T1;           //restore previous values
	T_Leaf=T2;
	T_Soil=T3;
}

void Get_DATA(long double dt_court)
{
	long double Radius_bark_trunk0, E_Plant,  K_Root, K_tot2, Diam_Fruit, Specific_Leaf_WC,Specific_Br_WC, Specific_Br_Symp_WC, Specifc_Br_Apo_WC, Volume_Branch,Sap_Flow_d,Sap_Flow_d2,Axil_WC;
	long double K_Rhizo,K_root1, K_root2, K_root3, P_Soil;
	long double RWC_leaf, RWC_leaf_s,RWC_leaf_a,RWC_Branch,RWC_Branch_s,RWC_Branch_a,RWC_Trunk,RWC_Trunk_s,RWC_Trunk_a,RWC_Root,RWC_Root_s,RWC_Root_a,RWC_Plant,RWC_Plant_s,RWC_Plant_a;
	
	if (FRACTAL==1)       Volume_Branch= Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4;                                                 // branch volume in m3
	else                Volume_Branch=  Q_Branch_Apo_FR/Branch_apo_fraction;                                                                                           // if Fractal then sapwood volume in m3
	if (Leaf_Area)      Specific_Leaf_WC= (Q_Leaf_Symp   + Q_Leaf_Apo0   *(100-PLC_Leaf_Apo)/100) /1000*18 / (LMA * Leaf_Area); else  Specific_Leaf_WC=0;    // g H20 per g leaf dry mass
	Specific_Br_WC=     (Q_Branch_Symp + Q_Branch_Apo0 *(100-PLC_Branch_Apo)/100) /1000/1000*18 / (Density * Volume_Branch);                            // kg H20 per kg total branch dry mass
	Specific_Br_Symp_WC= Q_Branch_Symp/1000/1000*18 / (Density * Volume_Branch* Branch_symp_fraction);                                                  // kg H20 per kg symplasm branch dry mass
	Specifc_Br_Apo_WC=  (Q_Branch_Apo0 *(100-PLC_Branch_Apo)/100) /1000/1000*18 / (Density * Volume_Branch* Branch_apo_fraction);                       // kg H20 per kg apoplasm branch dry mass
   
	Radius_bark_trunk0=  1000*(sqrt((Q_Trunk_Symp0*18/1000/1000/1000/Length_Trunk + 3.1416 *Diam_Trunk * Diam_Trunk/ 4)/3.1416) -Diam_Trunk/2);  //equivalent bark thickness in mm
	Radius_bark_trunk =  1000*(-Radius_bark_trunk0+ 1000*(sqrt((Growth_trunk*18/1000/1000/1000/Length_Trunk + 3.1416 *Diam_Trunk * Diam_Trunk/ 4)/3.1416) -Diam_Trunk/2)); //equivalent variation in bark thickness in µm
	//printf(" %.2Le %.2Le %.2Le %.2Le",Radius_bark_trunk0 , Radius_bark_trunk,Q_Trunk_Symp0,Length_Trunk );
	Petiole_diam =      1000*1000*pow(4/3.1416*Q_Petiole_Symp*18/1000/1000/Length_Petiole/Branch_symp_fraction/1000,0.5);                                                                    //in µm
	Diam_Fruit   =      pow(Q_Axil_Symp/1000/1000/1000*18*6/N_Axil/(WC_Axil*3.1416),0.333333333333333);
	Growth_fruit =      pow(Q_Axil_Symp0/1000/1000/1000*18*6/N_Axil/(WC_Axil*3.1416),0.333333333333333);
	
	if      (CUT==1)    E_Plant= (E_Leaf*Leaf_Area + E_Branch*Branch_Area + E_Axil*Area_Axil + E_Petiole*Petiole_area);                                     //cut base of branch
	else if (CUT==2)    E_Plant= (E_Leaf*Leaf_Area + E_Branch*Branch_Area + E_Axil*Area_Axil + E_Petiole*Petiole_area + E_Trunk*Trunk_Area     );           //cut base of trunk
	else                E_Plant= (E_Leaf*Leaf_Area + E_Branch*Branch_Area + E_Axil*Area_Axil + E_Petiole*Petiole_area + E_Trunk*Trunk_Area + E_Root1*Root_Area1 + E_Root2*Root_Area2 + E_Root3*Root_Area3 ); //total plant water loss in mmol/s
	EvapoT =            E_Plant + E_Soil*Soil_Width*Soil_Width;                                                                                             //total plant + soil water loss in mmol/s
	if      (CUT==1)
	{
		Q_Plant_s=(Q_Axil_Symp +Q_Petiole_Symp + Q_Leaf_Symp  + Q_Branch_Symp)*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo0*(100-PLC_Leaf_Apo) + Q_Axil_Apo0*(100-PLC_Axil_Apo) + Q_Branch_Apo0*(100-PLC_Branch_Apo))*18/1000/1000/100 ;
	}
	else if (CUT==2)
	{
		Q_Plant_s=(Q_Axil_Symp +Q_Petiole_Symp + Q_Leaf_Symp  + Q_Branch_Symp + Q_Trunk_Symp )*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo0*(100-PLC_Leaf_Apo)/100 + Q_Axil_Apo0*(100-PLC_Axil_Apo)/100 + Q_Branch_Apo0*(100-PLC_Branch_Apo)/100 + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100)*18/1000/1000 ;
	}
	else
	{
		Q_Plant_s=(Q_Axil_Symp +Q_Petiole_Symp + Q_Leaf_Symp  + Q_Branch_Symp + Q_Trunk_Symp +Q_Root_Symp1 +Q_Root_Symp2 +Q_Root_Symp3 )*18/1000/1000;
		Q_Plant_a= (Q_Leaf_Apo0*(100-PLC_Leaf_Apo) + Q_Axil_Apo0*(100-PLC_Axil_Apo) + Q_Branch_Apo0*(100-PLC_Branch_Apo) + Q_Trunk_Apo0*(100-PLC_Trunk_Apo) + Q_Root_Apo01*(100-PLC_Root_Apo1) + Q_Root_Apo02*(100-PLC_Root_Apo2) + Q_Root_Apo03*(100-PLC_Root_Apo3))*18/1000/1000/100 ;
	}
	Q_Plant= Q_Plant_a+Q_Plant_s;
	
	if (K_Soil1 && K_Interface1 && K_Root_Symp11 && K_Root_Apo1) K_root1=1/(1/K_Soil1 + 1/K_Interface1 + 1/K_Root_Symp11 + 1/K_Root_Apo1); else K_root1=0;
	if (K_Soil2 && K_Interface2 && K_Root_Symp12 && K_Root_Apo2) K_root2=1/(1/K_Soil2 + 1/K_Interface2 + 1/K_Root_Symp12 + 1/K_Root_Apo2); else K_root2=0;
	if (K_Soil3 && K_Interface3 && K_Root_Symp13 && K_Root_Apo3) K_root3=1/(1/K_Soil3 + 1/K_Interface3 + 1/K_Root_Symp13 + 1/K_Root_Apo3); else K_root3=0;
	K_Rhizo=K_root1+K_root2+K_root3;
	if (K_Root_Symp11 && K_Root_Apo1) K_root1=1/(1/K_Root_Symp11 + 1/K_Root_Apo1); else K_root1=0;
	if (K_Root_Symp12 && K_Root_Apo2) K_root2=1/(1/K_Root_Symp12 + 1/K_Root_Apo2); else K_root2=0;
	if (K_Root_Symp13 && K_Root_Apo3) K_root3=1/(1/K_Root_Symp13 + 1/K_Root_Apo3); else K_root3=0;
	K_Root=K_root1+K_root2+K_root3;
	if (K_Leaf_Apo && K_Leaf_Symp && K_Branch_Apo && K_Trunk_Apo && K_Root)            
	{
		K_Plant_20=  1/(1/K_Leaf_Apo*Fluidity_leaf +1/K_Leaf_Symp*Fluidity_leaf + 1/K_Branch_Apo*Fluidity_air+ 1/K_Trunk_Apo*Fluidity_air + 1/K_Root*Fluidity_soil);  
		K_Plant=  1/(1/K_Leaf_Apo +1/K_Leaf_Symp + 1/K_Branch_Apo+ 1/K_Trunk_Apo + 1/K_Root);             
	}
	else 
		{
			K_Plant=0;
			K_Plant_20=0;
		}
	if (K_Plant_20_0) PLC_Plant=100*(1-K_Plant_20/K_Plant_20_0); else PLC_Plant=0;
	if (Leaf_Area && K_Rhizo && K_Plant)  K_tot  =  1/(1/K_Rhizo+1/K_Plant)/Leaf_Area; else K_tot=0;
	if (K_Rhizo)   P_Soil= (P_Soil1*K_root1 + P_Soil2*K_root2 + P_Soil3*K_root3)/K_Rhizo; else P_Soil=-10;
	if (Leaf_Area) K_tot2= E_Plant/(P_Soil - P_Leaf_Evap + Pg_Leaf )/Leaf_Area; else K_tot2=0;
	
	RWC_leaf= (Q_Leaf_Apo0*(100-PLC_Leaf_Apo)/100+Q_Leaf_Symp)/(Q_Leaf_Apo0 + Q_Leaf_Symp0);
	RWC_leaf_s= Q_Leaf_Symp/Q_Leaf_Symp0;
	RWC_leaf_a= Q_Leaf_Apo0*(100-PLC_Leaf_Apo)/Q_Leaf_Apo0/100;
	RWC_Branch= (Q_Branch_Apo0*(100-PLC_Branch_Apo)/100 + Q_Branch_Symp)/(Q_Branch_Apo0 + Q_Branch_Symp0);
	RWC_Branch_s= Q_Branch_Symp/ Q_Branch_Symp0;
	RWC_Branch_a= Q_Branch_Apo0*(100-PLC_Branch_Apo)/Q_Branch_Apo0/100;
	RWC_Trunk= (Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Trunk_Symp)/(Q_Trunk_Apo0 + Q_Trunk_Symp0);
	RWC_Trunk_s= Q_Trunk_Symp/ Q_Trunk_Symp0;
	RWC_Trunk_a= Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/Q_Trunk_Apo0/100;
	RWC_Root=(Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 )/(Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03+Q_Root_Symp01+Q_Root_Symp02+Q_Root_Symp03);
	RWC_Root_s= (Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3)/(Q_Root_Symp01+Q_Root_Symp02+Q_Root_Symp03);
	RWC_Root_a= (Q_Root_Apo01*(100-PLC_Root_Apo1) + Q_Root_Apo02*(100-PLC_Root_Apo2) + Q_Root_Apo03*(100-PLC_Root_Apo3))/(Q_Root_Apo01+Q_Root_Apo02+Q_Root_Apo03)/100;
	RWC_Plant= (Q_Leaf_Apo0*(100-PLC_Leaf_Apo)/100+Q_Leaf_Symp +Q_Branch_Apo0*(100-PLC_Branch_Apo)/100 + Q_Branch_Symp + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Trunk_Symp + Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp + Q_Axil_Apo*(100-PLC_Axil_Apo)/100)/(Q_Leaf_Apo0+Q_Leaf_Symp0 +Q_Branch_Apo0 + Q_Branch_Symp0 + Q_Trunk_Apo0 + Q_Trunk_Symp0 + Q_Root_Apo01 + Q_Root_Apo02 + Q_Root_Apo03 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp + Q_Axil_Apo);
	RWC_Plant_s=(Q_Leaf_Symp  + Q_Branch_Symp + + Q_Trunk_Symp + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp)/(Q_Leaf_Symp0 + Q_Branch_Symp0 +Q_Trunk_Symp0 + Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp3 + Q_Axil_Symp + Q_Petiole_Symp);
	RWC_Plant_a=(Q_Leaf_Apo0*(100-PLC_Leaf_Apo)/100+Q_Branch_Apo0*(100-PLC_Branch_Apo)/100 + Q_Trunk_Apo0*(100-PLC_Trunk_Apo)/100 + Q_Root_Apo01*(100-PLC_Root_Apo1)/100 + Q_Root_Apo02*(100-PLC_Root_Apo2)/100 + Q_Root_Apo03*(100-PLC_Root_Apo3)/100 + Q_Axil_Apo*(100-PLC_Axil_Apo)/100)/(Q_Leaf_Apo0+Q_Branch_Apo0 +  Q_Trunk_Apo0 +  Q_Root_Apo01 + Q_Root_Apo02 + Q_Root_Apo03 + Q_Axil_Apo);
	RWC_shoot=  (Q_Leaf_Symp + Q_Leaf_Apo0   *(100-PLC_Leaf_Apo)/100 + Q_Branch_Symp + Q_Branch_Apo0 *(100-PLC_Branch_Apo)/100 ) / (Q_Leaf_Symp0 + Q_Leaf_Apo0   + Q_Branch_Symp0 + Q_Branch_Apo0);
	Teta_Soil=(Q_Soil1+Q_Soil2+Q_Soil3)/1000/1000/1000*18/Volume_soil;
	Teta_Soil1=(Q_Soil1)/1000/1000/1000*18/Volume_soil1;
	Teta_Soil2=(Q_Soil2)/1000/1000/1000*18/Volume_soil2;
	Teta_Soil3=(Q_Soil3)/1000/1000/1000*18/Volume_soil3;
	GPP=A_net_tot/1000/1000/Surface_Soil*12 ;
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
	if (DYNAMIC)    Sap_Flow_d=(dq_Trunk_Apo_Branch_Apo+dq_Trunk_Apo_Trunk_Symp)/1000*18/dt_court/(Diam_Trunk * Diam_Trunk * 3.1416 / 4 *(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction); //sap flow density in the trunk, in g/s/m2 of sapwood
	else            Sap_Flow_d=(E_Plant)/1000*18/(Diam_Trunk * Diam_Trunk * 3.1416 / 4 *(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction);
					Sap_Flow_d2=Sap_Flow_d/100000*3600; //sap flow density in the trunk, in dm3/dm2/h of sapwood
	if (Sap_Flow_d2<SF_min_d) SF_min_d=Sap_Flow_d2;
	if (Sap_Flow_d2>SF_max_d) SF_max_d=Sap_Flow_d2;
	if (CLIMAT==0 || CLIMAT==1 || CLIMAT==3 || CLIMAT==5) DATA[0]= T/3600/24; else DATA[0]= T/3600/24-1;
	DATA[1]=  T_air;
	DATA[2]=  T_Leaf;
	DATA[3]=  RH_air;
	DATA[4]=  PAR;
	DATA[5]=  VPD_Leaf;
	DATA[6]=  VPD_Cuti;
	DATA[7]=  VPD_Axil;
	DATA[8]=  VPD_Branch;
	DATA[9]=  VPD_Trunk;
	DATA[10]= VPD_Root1;
	DATA[11]= VPD_Root2;
	DATA[12]= VPD_Root3;
	DATA[13]= VPD_Soil;
	DATA[14]= PAR;
	DATA[15]= P_air;
	DATA[16]= E_clim;
	DATA[17]= E_Leaf;
	DATA[18]= E_cuti;
	DATA[19]= E_Branch;
	DATA[20]= E_Trunk;
	DATA[21]= E_Axil;
	DATA[22]= E_Root1;
	DATA[23]= E_Root2;
	DATA[24]= E_Root3;
	DATA[25]= E_Soil;
	DATA[26]= E_Plant;
	if (Leaf_Area) DATA[27]= E_Plant/Leaf_Area; else DATA[27]=0;
	DATA[28]= E_tot*18/1000/1000;
	DATA[29]= g_Canopy;
	DATA[30]= g_s;
	DATA[31]= g_cuti;
	DATA[32]= P_Leaf_Evap;
	DATA[33]= P_Leaf_Symp;
	DATA[34]= P_Leaf_Apo;
	DATA[35]= turgor;
	DATA[36]= Turgor_Leaf_Symp;
	DATA[37]= P_Axil_Symp;
	DATA[38]= P_Axil_Apo;
	DATA[39]= P_Branch_Symp;
	DATA[40]= P_Branch_Apo;
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
	if (Leaf_Area) DATA[65]= K_Leaf_Symp/Leaf_Area;
	if (Leaf_Area) DATA[66]= K_Leaf_Apo/Leaf_Area;
	if (Leaf_Area) DATA[67]= 1/(1/K_Leaf_Apo+1/K_Leaf_Symp)/Leaf_Area;
	DATA[68]= K_Root_Symp11;
	DATA[69]= K_Root_Apo1;
	DATA[70]= K_Root_Symp12;
	DATA[71]= K_Root_Apo2;
	DATA[72]= K_Root_Symp13;
	DATA[73]= K_Root_Apo3;
	if (Leaf_Area) DATA[74]= K_Root/Leaf_Area;
	if (Leaf_Area) DATA[75]= K_Trunk_Apo/Leaf_Area;
	DATA[76]= K_Plant;
	DATA[77]= K_tot;
	DATA[78]= PLC_Leaf_Apo;
	DATA[79]= PLC_Branch_Apo;
	DATA[80]= PLC_Axil_Apo;
	DATA[81]= PLC_Trunk_Apo;
	DATA[82]= PLC_Root_Apo1;
	DATA[83]= PLC_Root_Apo2;
	DATA[84]= PLC_Root_Apo3;
	DATA[85]= Q_Plant;
	DATA[86]= (Q_Soil01-Q_Soil1+Q_Soil02-Q_Soil2+Q_Soil03-Q_Soil3)*18/1000/1000;
	DATA[87]= Q_Leaf_Evap*18/1000/1000;
	DATA[88]= Q_Leaf_Symp*18/1000/1000;
	DATA[89]= Q_Leaf_Apo*18/1000/1000;
	DATA[90]= Q_Branch_Symp*18/1000/1000;
	DATA[91]= Q_Branch_Apo*18/1000/1000;
	DATA[92]= Q_Axil_Symp*18/1000/1000;
	DATA[93]= Q_Axil_Apo*18/1000/1000;
	DATA[94]= Q_Trunk_Symp*18/1000/1000;
	DATA[95]= Q_Trunk_Apo*18/1000/1000;
	DATA[96]= (Q_Root_Symp1+Q_Root_Symp2+Q_Root_Symp2)*18/1000/1000;
	//DATA[97]= (Q_Root_Apo1+Q_Root_Apo2+Q_Root_Apo3)*18/1000/1000;
	DATA[97]= Q_Root_Apo_t*18/1000/1000;
	DATA[98]= (Q_Root_Endo1+Q_Root_Endo2+Q_Root_Endo3)*18/1000/1000;
	DATA[99]=  Q_Soil1*18/1000/1000;
	DATA[100]= Q_Soil2*18/1000/1000;
	DATA[101]= Q_Soil3*18/1000/1000;
	DATA[102]= Growth_rate_trunk;
	DATA[103]= Radius_bark_trunk;
	DATA[104]= A_net;
	DATA[105]= A_net_tot_c/1000/1000;  //Annual plant canopy A_net in mol of CO2
	DATA[106]= A_net_tot/E_tot/1000;
	DATA[107]= Rm;
	DATA[108]= Reserve/1000;
	DATA[109]= Leaf_Area;
	DATA[110]= Specific_Leaf_WC;
	DATA[111]= Specific_Br_WC;
	DATA[112]= Specific_Br_Symp_WC;
	DATA[113]= Specifc_Br_Apo_WC;
	DATA[114]= RWC_shoot;
	DATA[115]= RWC_Axil;
	DATA[116]= RWC_Soil;
	DATA[117]= RWC_min;
	DATA[118]= RWC_int;
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
	DATA[133]= Rain_leaf_tot;
	DATA[134]= ETP_leaf_tot;
	if (N_days) DATA[135]= Cum_T_air/(long double)N_days; else DATA[135]=0;
	if (N_days2) DATA[136]= Cum_T_air_l/(long double)N_days2; else DATA[136]=0;
	DATA[137]= EvapoT;
	DATA[138]= EvapoT_tot*18/1000/1000/Surface_Soil;  // cumul plant + soil evaporation in mm
	DATA[139]= ETP_Penman_day;           // cumul ETP per day in mm
	DATA[140]= (dq_Trunk_Apo_Branch_Apo+dq_Trunk_Apo_Trunk_Symp)/dt_court;  //sap flow in the trunk in mmol/s; includes trunk bark transpiration
	if (Leaf_Area) DATA[141]=K_Branch_Apo/Leaf_Area;
	DATA[142]= P_min_leaf;
	DATA[143]= P_min_stem;
	DATA[144]= K_tot2;
	DATA[145]= A_net_c;                              // Canopy A_net per ground area µmol
	DATA[146]= Rain_tot-Rain_soil ;					// Canopy interception, mm
	DATA[147]= A_net_tot_c/1000/1000/Surface_Soil*12 ; //GPP in g of C per soil surface per year
	DATA[148]= A_net_day_c/1000/1000/Surface_Soil*12 ; //GPP in g of C per soil surface per day
	DATA[149]= EvapoT_day*18/1000/1000/Surface_Soil;
	DATA[150]= E_tot_day*18/1000/1000/Surface_Soil;
	DATA[151]= T_Soil;
	DATA[152]= POTENTIAL_PAR;
	DATA[153]= cloud_cover;
	DATA[154]= -dq_Leaf_Symp_Leaf_Evap/dt_court;
	DATA[155]= dq_Leaf_Apo_Leaf_Evap/dt_court;
	DATA[156]= dq_Branch_Apo_Leaf_Apo/dt_court;
	DATA[157]= dq_Branch_Apo_Branch_Symp/dt_court;
	DATA[158]= dq_Trunk_Apo_Branch_Apo/dt_court;
	DATA[159]= dq_Trunk_Apo_Trunk_Symp/dt_court;
	DATA[160]= dq_Axil_Apo_Axil_Symp/dt_court;
	DATA[161]= dq_Branch_Apo_Axil_Apo/dt_court;
	DATA[162]= (dq_Root_Apo_Trunk_Apo1+dq_Root_Apo_Trunk_Apo2+dq_Root_Apo_Trunk_Apo3)/dt_court;
	DATA[163]= (dq_Root_Endo_Root_Apo1+dq_Root_Endo_Root_Apo2+dq_Root_Endo_Root_Apo3)/dt_court;
	DATA[164]= (dq_Root_Symp_Root_Endo1+dq_Root_Symp_Root_Endo2+dq_Root_Symp_Root_Endo3)/dt_court;
	DATA[165]= (dq_Soil_Root_Endo1+dq_Soil_Root_Endo2+dq_Soil_Root_Endo3)/dt_court;
	DATA[166]= dq_Soil_12/dt_court;
	DATA[167]= dq_Soil_23/dt_court;
	DATA[168]= A_net1;
	DATA[169]= A_net2;
	DATA[170]= g_Axil;
	DATA[171]= Root_Area_fi;
	DATA[172]= Axil_WC;
	DATA[173]= ((Q_Soil1+Q_Soil2+Q_Soil3)*18/1000/1000/1000/Volume_soil- Teta_r)*Soil_Depth*1000; //RU tot
	DATA[174]= E_tot*18/1000/1000/Surface_Soil;
	DATA[175]= Drainage*18/1000/1000/Surface_Soil;
	DATA[176]= A_net-Rm;                 //A_gross in µmol
	DATA[177]= A_net_tot/1000;           //Annual A_net in mmol
	DATA[178]= A_gross_tot/1000;         //Annual A_gross in mmol
	DATA[179]= Resp_tot/1000;            //Annual Leaf Respiration in mmol
	DATA[180]= Export;
	DATA[181]= Export_tot/1000;
	DATA[182]= T_Soil1;
	DATA[183]= T_Soil2;
	DATA[184]= T_Soil3;
	DATA[185]= ((Q_Soil1+Q_Soil2+Q_Soil3)*18/1000/1000/1000/Volume_soil- Teta_wp)*Soil_Depth*1000;
	DATA[186]= REW_t; //RWC based on Teta_r
	DATA[187]= Q_Petiole_Symp*18/1000/1000;   // in Kg
	DATA[188]= E_Petiole;
	DATA[189]= P_Petiole_Symp;
	DATA[190]= Petiole_diam;
	DATA[191]= Q_Plant_a;
	DATA[192]= Q_Plant_s;
	DATA[193]= RWC_leaf;
	DATA[194]= RWC_leaf_s;
	DATA[195]= RWC_leaf_a;
	DATA[196]= RWC_Branch;
	DATA[197]= RWC_Branch_s;
	DATA[198]= RWC_Branch_a;
	DATA[199]= RWC_Trunk;
	DATA[200]= RWC_Trunk_s;
	DATA[201]= RWC_Trunk_a;
	DATA[202]= RWC_Root;
	DATA[203]= RWC_Root_s;
	DATA[204]= RWC_Root_a;
	DATA[205]= RWC_Plant;
	DATA[206]= RWC_Plant_s;
	DATA[207]= RWC_Plant_a;
	if (Leaf_Area) DATA[208]= E_tot_day/3600/24/Leaf_Area; else DATA[208]=0;  // in mmol/s/m2
	DATA[209]= P_min_lf_d;
	DATA[210]= P_max_lf_d;
	DATA[211]= gs_min_d;
	DATA[212]= gs_max_d;
	DATA[213]= SF_min_d;
	DATA[214]= SF_max_d;
	DATA[215]= gs_max_d2;   //gs_max at Pmin
	DATA[216]= Growth_fruit*1000; //fruit diameter in mm without elastic variation
	DATA[217]= Teta_Soil; 
	DATA[218]= Teta_Soil1; 
	DATA[219]= Teta_Soil2; 
	DATA[220]= Teta_Soil3; 
	DATA[221]= Ca; 
	DATA[222]= K_Plant_20;
	DATA[223]= PLC_Plant;
}

void compute_climatic_stats(void) //computes some stats on the annual climatic data when given on a daily basis
{
	long double YEAR, DOI,DOI0,Tmin,Tmax,Tsoil,RHmin,RHmax,PAR_day,Rain,Wind_day;
	long double T_annual_mean=0,T_summer_mean=0,T_annual_min=0,T_summer_min=0,T_annual_max=0,T_summer_max=0,HR_annual_mean=0,HR_summer_mean=0,HR_annual_min=0,HR_summer_min=0,Wind_annual_mean=0,Wind_summer_mean=0;
	long double PAR_annual_max=0,PAR_summer_max=0,Rain_annual=0,Rain_summer=0;
	long double Days=0,days_summer=0,tangente,TTTT;
	long double e_sat_air,e_air,slope,vpd,ETP,VPD_annual=0,VPD_summer=0,ETP_annual=0,ETP_summer=0;
	FILE *out;
	int i;
	printf("Computing climatic data\n");
	out = fopen("Stat_Clim.csv","w");
	fprintf(out,"YEAR;T_annual_mean;T_annual_min;T_annual_max;HR_annual_mean;HR_annual_min;Wind_annual_mean;PAR_annual_max;Rain_annual;VPD_annual;ETP_annual;;T_summer_mean;T_summer_min;T_summer_max;HR_summer_mean;HR_summer_min;Wind_summer_mean;PAR_summer_max;Rain_summer;VPD_summer;ETP_summer\n");
	climat_in = fopen(filename_CLIM,"r+");
	if (CLIMAT==4) fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI0, &Tmin, &Tmax, &Tsoil, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
	else           fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",    &YEAR, &DOI0, &Tmin, &Tmax, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
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
		if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR=0;
		else PAR = PAR_day*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
		
		if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
		else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
		
		e_air=e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
		vpd+= (e_sat_air - e_air)/1000;
		slope=4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
		ETP+= 0.5625*(0.408*slope*PAR*0.5495*3.6e-3+0.066*37*Wind_day*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind_day));  //ETP in mm/h  0.5625 to fit observed ETP
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
	
	if (DOI0>171 && DOI0<265)  //summer period
	{
		days_summer++;
		T_summer_mean=      (Tmin+Tmax)/2;
		T_summer_min=       Tmin;
		T_summer_max=       Tmax;
		HR_summer_mean=     (RHmin+RHmax)/2;
		HR_summer_min=      RHmin;
		Wind_summer_mean=   Wind_day;
		PAR_summer_max=     PAR_day;
		Rain_summer=        Rain;
		VPD_summer=         vpd;
		ETP_summer=         ETP;
	}
	
	while (!feof(climat_in))
	{
		if (CLIMAT==4) fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI, &Tmin, &Tmax, &Tsoil, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
		else fscanf(climat_in,"%Le %Le %Le %Le %Le %Le %Le %Le %Le\n",&YEAR, &DOI, &Tmin, &Tmax, &RHmin, &RHmax, &PAR_day, &Rain, &Wind_day);
		if (DOI<DOI0) // a new year
		{
			fprintf(out,"%.0Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;; ", YEAR-1,T_annual_mean,T_annual_min,T_annual_max,HR_annual_mean,HR_annual_min,Wind_annual_mean,PAR_annual_max,Rain_annual,VPD_annual,ETP_annual);
			fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;\n", T_summer_mean,T_summer_min,T_summer_max,HR_summer_mean,HR_summer_min,Wind_summer_mean,PAR_summer_max,Rain_summer,VPD_summer,ETP_summer);
			T_annual_mean=T_annual_min=T_annual_max=HR_annual_mean=HR_annual_min=Wind_annual_mean=PAR_annual_max=Rain_annual=VPD_annual=ETP_annual=0;
			T_summer_mean=T_summer_min=T_summer_max=HR_summer_mean=HR_summer_min=Wind_summer_mean=PAR_summer_max=Rain_summer=VPD_summer=ETP_summer=0;
			Days=0;
			days_summer=0;
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
			if ((TTTT < (12 - Day_length/2)) || (TTTT > (12 + Day_length/2))) PAR=0;
			else PAR = PAR_day*cos(3.14159265359*(TTTT-12)/Day_length);   //daily cos variation of incident PAR between 6am to 6pm, max at 12h00
			
			if (T_air>0) e_sat_air=611.21*exp((18.678-T_air/234.5)*T_air/(257.14+T_air));    // saturation vapour water pressure at Tair in Pa from Buck's equation
			else         e_sat_air=611.15*exp((23.036-T_air/333.7)*T_air/(279.82+T_air));
			
			e_air=   e_sat_air*RH_air/100;                                         // vapour water pressure at Tair and RHair
			vpd+=   (e_sat_air - e_air)/1000;
			slope=  4098*0.6018*exp(17.27*T_air/(T_air+237.3))/(pow(T_air+237.3,2));
			ETP+=   0.5625*(0.408*slope*PAR*0.5495*3.6e-3+0.066*37*Wind_day*(e_sat_air-e_air)/1000/(T_air+273))/(slope+0.066*(1+0.34*Wind_day));  //ETP in mm/h  0.5625 to fit observed ETP
			//  printf("%Lf T %Lf RH %Lf PAR %Lf vpd %Lf ETP %Lf \n", TTTT, T_air,RH_air, PAR,vpd,ETP);
			
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
		
		if (DOI>171 && DOI<265)  //summer period
		{
			days_summer++;
			T_summer_mean=  (T_summer_mean*(days_summer-1)+(Tmin+Tmax)/2)/days_summer;
			T_summer_min=   (T_summer_min*(days_summer-1)+Tmin)/days_summer;
			T_summer_max=   (T_summer_max*(days_summer-1)+Tmax)/days_summer;
			HR_summer_mean= (HR_summer_mean*(days_summer-1)+(RHmin+RHmax)/2)/days_summer;
			HR_summer_min=  (HR_summer_min*(days_summer-1) + RHmin)/days_summer;
			Wind_summer_mean=(Wind_summer_mean*(days_summer-1)+Wind_day)/days_summer;
			PAR_summer_max= (PAR_summer_max*(days_summer-1)+PAR_day)/days_summer;
			//VPD_summer=(VPD_summer*(days_summer-1)+vpd)/days_summer;
			//ETP_summer=(ETP_summer*(days_summer-1)+ETP)/days_summer;
			ETP_summer+=ETP;
			VPD_summer+=vpd;
			Rain_summer+=Rain;
		}
	}
	fprintf(out,"%.0Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;; ", YEAR,T_annual_mean,T_annual_min,T_annual_max,HR_annual_mean,HR_annual_min,Wind_annual_mean,PAR_annual_max,Rain_annual,VPD_annual,ETP_annual);
	fprintf(out,"%Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf; %Lf;\n", T_summer_mean,T_summer_min,T_summer_max,HR_summer_mean,HR_summer_min,Wind_summer_mean,PAR_summer_max,Rain_summer,VPD_summer,ETP_summer);
	fclose(climat_in);
	fclose(out);
}

void print_screen(void)
{
	int i;
	
	Get_DATA(dt);
	printf("  %5.0Lf",YEAR2);
	if (PRINT_SCREEN==2) printf("%7.2Lf %7.2Lf", DATA[0],PLC_Leaf_Apo);
	else
		for (i=0;i<250;i++)
			if (Screen_out[i])
			{
				if (i==0) printf("%7.2Lf ", DATA[0]);
				else printf("%9.2Le  ", DATA[i]);
			}
	// printf("%d %Le", DEAD,dt);
	if (WARNING)
	{
		printf("  WARNING %Lf", dt);
		WARNING=0;
	}
	  
	printf("\n"); 
	
   // printf("dq:%Le %Le\n",K_Soil1, K_Soil_tot1 );
}

void print_transient(void)
{
	FILE *transient;
	int i;
	
	transient = fopen(filename_TRANS,"a");
	Get_DATA(dt);
	fprintf(transient,"%.0Lf; ",YEAR1);
	for (i=0;i<250;i++)
		if (debug || File_out[i] )
		{
			if (i==0) fprintf(transient, "%.4Lf; ", DATA[i]);
			else fprintf(transient, "%.8Le; ", DATA[i]);
		}
	fprintf(transient,"\n");
	fclose(transient);
}

void Soil_compet (void)  //when plants compet for water in the same soil volume
{
	long double  Q_soil1_saved, Q_soil2_saved, Q_soil3_saved;
	long double  dQ_soil1, dQ_soil2, dQ_soil3;
	FILE *soil;
	
	dQ_soil1=Q_soil1_init-Q_Soil1;
	dQ_soil2=Q_soil2_init-Q_Soil2;
	dQ_soil3=Q_soil3_init-Q_Soil3;
	
	soil = fopen("soil.trs","r+");  
	fscanf(soil,"%Lf %Lf %Lf",&Q_soil1_saved,&Q_soil2_saved,&Q_soil3_saved);
	rewind (soil);
	Q_soil1_init= Q_soil1_saved - dQ_soil1;
	Q_soil2_init= Q_soil2_saved - dQ_soil2;
	Q_soil3_init= Q_soil3_saved - dQ_soil3;
	fprintf(soil,"%Lf %Lf %Lf",Q_soil1_init,Q_soil2_init,Q_soil3_init);
	fclose(soil);

	Q_Soil1=Q_soil1_init;
	Q_Soil2=Q_soil2_init;
	Q_Soil3=Q_soil3_init;
}

void Gravity_test(void)
{
	if(T>(DOY+6.25)*3600*24)
		{
		Pg_Leaf  = -9.81*(Length_Branch+Length_Trunk*4)/1000;
		Pg_trunk = -9.81*Length_Trunk*4/1000;                     // gravimetric pressure drop at the top leaves, trunk and branches in MPa
		Pg_branch= -9.81*(Length_Branch+Length_Trunk*4)/1000;
		P_Trunk_Apo=Pg_trunk;
		}
	else if(T>(DOY+4.25)*3600*24)
		{
		Pg_Leaf  = -9.81*(Length_Branch+Length_Trunk*3)/1000;
		Pg_trunk = -9.81*Length_Trunk*3/1000;                     // gravimetric pressure drop at the top leaves, trunk and branches in MPa
		Pg_branch= -9.81*(Length_Branch+Length_Trunk*3)/1000;
		P_Trunk_Apo=Pg_trunk;
		}
	else if(T>(DOY+2.25)*3600*24)
		{
		Pg_Leaf  = -9.81*(Length_Branch+Length_Trunk*2)/1000;
		Pg_trunk = -9.81*Length_Trunk*2/1000;                     // gravimetric pressure drop at the top leaves, trunk and branches in MPa
		Pg_branch= -9.81*(Length_Branch+Length_Trunk*2)/1000;
		P_Trunk_Apo=Pg_trunk;
		}
	else if(T>(DOY+0.25)*3600*24)
		{
		Pg_Leaf  = -9.81*(Length_Branch+Length_Trunk)/1000;
		Pg_trunk = -9.81*Length_Trunk/1000;                     // gravimetric pressure drop at the top leaves, trunk and branches in MPa
		Pg_branch= -9.81*(Length_Branch+Length_Trunk)/1000;
		P_Trunk_Apo=Pg_trunk;
		}
}

//***************************//
// SHORT TIME STEP FUNCTIONS //
//***************************//
void Compute_dq(long double dt_court)  // in mmol H20 exchanged between compartments during dt
{
	
								dq_Leaf_Symp_Leaf_Evap       	=   dt_court * K_Leaf_Symp2   *   (P_Leaf_Symp     -  P_Leaf_Evap);
	if (PLC_Leaf_Apo<100)   	dq_Leaf_Apo_Leaf_Evap       	=   dt_court * K_Leaf_Symp    *   (P_Leaf_Apo      -  P_Leaf_Evap); else dq_Leaf_Apo_Leaf_Evap=0;
	if (PLC_Branch_Apo<100) 	dq_Branch_Apo_Leaf_Apo      	=   dt_court * K_Leaf_Apo     *   (P_Branch_Apo    -  Pg_branch - P_Leaf_Apo + Pg_Leaf); else dq_Branch_Apo_Leaf_Apo=0;
	if (PLC_Branch_Apo<100) 	dq_Branch_Apo_Branch_Symp    	=   dt_court * K_Branch_Symp  *   (P_Branch_Apo    -  P_Branch_Symp); else dq_Branch_Apo_Branch_Symp =0;
	if (PLC_Trunk_Apo<100)  	dq_Trunk_Apo_Branch_Apo       	=   dt_court * K_Branch_Apo   *   (P_Trunk_Apo     -  Pg_trunk - P_Branch_Apo + Pg_branch); else dq_Trunk_Apo_Branch_Apo=0;
	if (PLC_Trunk_Apo<100)  	dq_Trunk_Apo_Trunk_Symp    	=   dt_court * K_Trunk_Symp   *   (P_Trunk_Apo     -  P_Trunk_Symp); else dq_Trunk_Apo_Trunk_Symp=0;
								dq_Root_Apo_Trunk_Apo      	=   dt_court*K_Trunk_Apo      *   (P_Root_Apo      -  P_Trunk_Apo + Pg_trunk);

	if (Type_Axil) // flow for one axillary organ
	{
		if (PLC_Axil_Apo<100)   dq_Axil_Apo_Axil_Symp       =   dt_court*K_Axil_Symp   *   (P_Axil_Apo      -  P_Axil_Symp);      else dq_Axil_Apo_Axil_Symp=0;
								 dq_Branch_Symp_Axil_Symp    =   dt_court*K_Axil_Symp2  *   (P_Branch_Symp   -  P_Axil_Symp);
		if (PLC_Branch_Apo<100) dq_Branch_Apo_Axil_Apo      =   dt_court*K_Axil_Apo    *   (P_Branch_Apo    -  P_Axil_Apo);       else dq_Branch_Apo_Axil_Apo =0;
		if (PLC_Axil_Apo<100)   dq_Axil_Apo_Petiole_Symp    =   dt_court*K_Axil_Symp2  *   (P_Axil_Apo      -  P_Petiole_Symp);   else dq_Axil_Apo_Petiole_Symp=0;
	}

	if (Root_upper)
	{
		if (PLC_Root_Apo1<100) dq_Root_Apo_Trunk_Apo1      =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_trunk); else dq_Root_Apo_Trunk_Apo1 =0;
		if (PLC_Root_Apo1<100) dq_Root_Endo_Root_Apo1      =   dt_court*K_Root_Apo1    *   (P_Root_Endo1     -  P_Root_Apo); else  dq_Root_Endo_Root_Apo1=0;
		dq_Root_Symp_Root_Endo1     =   dt_court*K_Root_Symp21  *   (P_Root_Symp1     -  P_Root_Endo1);
		dq_Soil_Root_Endo1          =   dt_court*(1/(1/K_Soil1 +1/K_Root_Symp11 +1/K_Interface1))       *  (P_Soil1 -  P_Root_Endo1);
	}
	
	if (Root_middle)
	{
		if (PLC_Root_Apo2<100) dq_Root_Apo_Trunk_Apo2      =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_trunk); else dq_Root_Apo_Trunk_Apo2 =0;
		if (PLC_Root_Apo2<100) dq_Root_Endo_Root_Apo2      =   dt_court*K_Root_Apo2    *   (P_Root_Endo2     -  P_Root_Apo); else dq_Root_Endo_Root_Apo2=0;
		dq_Root_Symp_Root_Endo2     =   dt_court*K_Root_Symp22  *   (P_Root_Symp2     -  P_Root_Endo2);
		dq_Soil_Root_Endo2          =   dt_court*(1/(1/K_Soil2 +1/K_Root_Symp12 +1/K_Interface2))       *  (P_Soil2 -  P_Root_Endo2);
	}
	
	if (Root_lower)
	{
	   if (PLC_Root_Apo3<100) dq_Root_Apo_Trunk_Apo3      =   dt_court*K_Trunk_Apo    *   (P_Root_Apo      -  P_Trunk_Apo + Pg_trunk); else dq_Root_Apo_Trunk_Apo3 =0;
		if (PLC_Root_Apo3<100) dq_Root_Endo_Root_Apo3      =   dt_court*K_Root_Apo3    *   (P_Root_Endo3     -  P_Root_Apo); else dq_Root_Endo_Root_Apo3=0;
		dq_Root_Symp_Root_Endo3     =   dt_court*K_Root_Symp23  *   (P_Root_Symp3     -  P_Root_Endo3);
		dq_Soil_Root_Endo3          =   dt_court*(1/(1/K_Soil3 +1/K_Root_Symp13 +1/K_Interface3))       *  (P_Soil3 -  P_Root_Endo3);
	}
	
	if (CAPILLARITY)
	{
		dq_Soil_12                  =   dt_court*K_Soil_tot1*K_Soil_tot2/(K_Soil_tot1+K_Soil_tot2) *(P_Soil1-P_Soil2);
		dq_Soil_23                  =   dt_court*K_Soil_tot2*K_Soil_tot3/(K_Soil_tot2+K_Soil_tot3) *(P_Soil2-P_Soil3);
	}
	else dq_Soil_12 = dq_Soil_23 =0.0;
}

void Compute_Q(void)  // changes in water content of the different compartments
{
	Q_Leaf_Evap+=           dq_Leaf_Apo_Leaf_Evap       +   dq_Leaf_Symp_Leaf_Evap          - dq_stomate;
	Q_Leaf_Symp+=        -  dq_Leaf_Symp_Leaf_Evap      -   dq_cuti  ;
	Q_Leaf_Apo+=         -  dq_Leaf_Apo_Leaf_Evap       +   dq_Branch_Apo_Leaf_Apo;
	
	if (Type_Axil)				// if there is an axillary organ
	{
		if (Type_Axil==3) 		// a flower
		{
		Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp       -   dq_Axil;
		Q_Petiole_Symp+=       dq_Axil_Apo_Petiole_Symp    -   dq_Petiole;
		Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo      -   dq_Axil_Apo_Axil_Symp  -  dq_Axil_Apo_Petiole_Symp;
		Q_Branch_Symp+=        dq_Branch_Apo_Branch_Symp   -   dq_Branch;
		Q_Branch_Apo+=       - dq_Branch_Apo_Leaf_Apo      -   dq_Branch_Apo_Branch_Symp       + dq_Trunk_Apo_Branch_Apo -  dq_Branch_Apo_Axil_Apo;
		}
		else if (Type_Axil==2) // a fruit
		{
		Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp      -   dq_Axil   + dq_fruit;
		Q_Petiole_Symp+=       dq_Axil_Apo_Petiole_Symp   -   dq_Petiole;
		Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo     -   dq_Axil_Apo_Axil_Symp         - dq_fruit;
		Q_Branch_Symp+=        dq_Branch_Apo_Branch_Symp  -   dq_Branch;
		Q_Branch_Apo+=       - dq_Branch_Apo_Leaf_Apo     -   dq_Branch_Apo_Branch_Symp     + dq_Trunk_Apo_Branch_Apo -  dq_Branch_Apo_Axil_Apo;
		}
		else if (Type_Axil==1) // a bud
		{
		Q_Axil_Symp+=          dq_Axil_Apo_Axil_Symp      +   dq_Branch_Symp_Axil_Symp       - dq_Axil  ;
		Q_Axil_Apo+=           dq_Branch_Apo_Axil_Apo     -   dq_Axil_Apo_Axil_Symp        ;
		Q_Branch_Symp+=        dq_Branch_Apo_Branch_Symp  -   dq_Branch_Symp_Axil_Symp       - dq_Branch;
		Q_Branch_Apo+=       - dq_Branch_Apo_Leaf_Apo     -   dq_Branch_Apo_Branch_Symp       + dq_Trunk_Apo_Branch_Apo -   dq_Branch_Apo_Axil_Apo;
		}
	}
	
	else //no organ
	{
		Q_Branch_Apo+=       -  dq_Branch_Apo_Leaf_Apo      -   dq_Branch_Apo_Branch_Symp       + dq_Trunk_Apo_Branch_Apo;
		Q_Branch_Symp+=         dq_Branch_Apo_Branch_Symp   -   dq_Branch;
	}
	
	Q_Trunk_Symp+=          dq_Trunk_Apo_Trunk_Symp     -   dq_Trunk;
	Q_Trunk_Apo+=        -  dq_Trunk_Apo_Branch_Apo     -   dq_Trunk_Apo_Trunk_Symp      + dq_Root_Apo_Trunk_Apo;
	Q_Root_Apo1+=        -  dq_Root_Apo_Trunk_Apo1      +   dq_Root_Endo_Root_Apo1 		;
	Q_Root_Apo2+=        -  dq_Root_Apo_Trunk_Apo2      +   dq_Root_Endo_Root_Apo2 		;
	Q_Root_Apo3+=        -  dq_Root_Apo_Trunk_Apo3      +   dq_Root_Endo_Root_Apo3 		;
	Q_Root_Apo_t+=		  -  dq_Root_Apo_Trunk_Apo	  	   +   dq_Root_Endo_Root_Apo1 			+  dq_Root_Endo_Root_Apo2 +   dq_Root_Endo_Root_Apo3; 
	Q_Root_Endo1+=       -  dq_Root_Endo_Root_Apo1      +   dq_Root_Symp_Root_Endo1         + dq_Soil_Root_Endo1;
	Q_Root_Endo2+=       -  dq_Root_Endo_Root_Apo2      +   dq_Root_Symp_Root_Endo2         + dq_Soil_Root_Endo2;
	Q_Root_Endo3+=       -  dq_Root_Endo_Root_Apo3      +   dq_Root_Symp_Root_Endo3         + dq_Soil_Root_Endo3;
	Q_Root_Symp1+=       -  dq_Root_Symp_Root_Endo1     -   dq_Root1;
	Q_Root_Symp2+=       -  dq_Root_Symp_Root_Endo2     -   dq_Root2;
	Q_Root_Symp3+=       -  dq_Root_Symp_Root_Endo3     -   dq_Root3;
	Q_Soil1+=            -  dq_Soil_Root_Endo1          -   dq_Soil                         - dq_Soil_12;
	Q_Soil2+=            -  dq_Soil_Root_Endo2          +   dq_Soil_12                      - dq_Soil_23;
	Q_Soil3+=            -  dq_Soil_Root_Endo3          +   dq_Soil_23;
}

void Compute_P_dynamic(void)   // compute water potentials
{
	long double Tp, Pi, Rs;
	
	//leaf
	P_Leaf_Evap=(Q_Leaf_Evap-Q_Leaf_Evap0)/C_Leaf_Evap + Pg_Leaf;
	Rs=(Q_Leaf_Symp0-Q_Leaf_Symp)/Q_Leaf_Symp0;
	Tp=-Pi0_Leaf_Symp*Osmotic_TLeaf - Epsilon_Leaf_Symp*Rs;
	if (Tp<0)Tp=0;
	if (1-Rs) Pi=Pi0_Leaf_Symp*Osmotic_TLeaf/(1-Rs); else Pi=-1000;
	P_Leaf_Symp=Tp+Pi;
	P_Leaf_Apo=(Q_Leaf_Apo-Q_Leaf_Apo1)/C_Leaf_Apo+Pg_Leaf;
	Turgor_Leaf_Symp=Tp;
	
	//branch
	Rs=(Q_Branch_Symp0-Q_Branch_Symp)/Q_Branch_Symp0;
	Tp=-Pi0_Branch_Symp*Osmotic_TAir - Epsilon_Branch_Symp*Rs;
	if (Tp<0)Tp=0;
	if (1-Rs) Pi=Pi0_Branch_Symp*Osmotic_TAir/(1-Rs); else Pi=-1000;
	P_Branch_Symp=Tp+Pi;
	P_Branch_Apo=(Q_Branch_Apo-Q_Branch_Apo1)/C_Branch_Apo+Pg_branch;
	
	//Axil
	if (Type_Axil)
	{
		Rs=(Q_Axil_Symp0-Q_Axil_Symp)/Q_Axil_Symp0;
		Tp=-Pi0_Axil_Symp*Osmotic_TAir - Epsilon_Axil_Symp*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Axil_Symp* Osmotic_TAir /(1-Rs); else Pi=-1000;
		P_Axil_Symp=Tp+Pi;
		P_Axil_Apo=(Q_Axil_Apo-Q_Axil_Apo1)/C_Axil_Apo+Pg_branch;
		Turgor_Axil_Symp=Tp;
		
		//Axil Petiole
		if (Type_Axil>=2) // for flower or fruit with a petiole only
		{
			Rs=(Q_Petiole_Symp0-Q_Petiole_Symp)/Q_Petiole_Symp0;
			Tp=-Pi0_Axil_Symp*Osmotic_TAir - Epsilon_Axil_Symp*Rs;
			if (Tp<0)Tp=0;
			if (1-Rs) Pi=Pi0_Axil_Symp* Osmotic_TAir /(1-Rs); else Pi=-1000;
			P_Petiole_Symp=Tp+Pi;
		}
	}
	//trunk
	Rs=(Q_Trunk_Symp0-Q_Trunk_Symp)/Q_Trunk_Symp0;
	Tp=-Pi0_Trunk_Symp*Osmotic_TAir - Epsilon_Trunk_Symp*Rs;
	if (Tp<0)Tp=0;
	if (1-Rs) Pi=Pi0_Trunk_Symp*Osmotic_TAir/(1-Rs); else Pi=-1000;
	P_Trunk_Symp=Tp+Pi;
	P_Trunk_Apo=(Q_Trunk_Apo-Q_Trunk_Apo1)/C_Trunk_Apo+Pg_trunk;
	Turgor_Trunk_Symp=Tp;
	
	//roots
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
	}
	if (Root_middle)
	{
		Rs=(Q_Root_Symp02-Q_Root_Symp2)/Q_Root_Symp02;
		Tp=-Pi0_Root_Symp2*Osmotic_TSoil - Epsilon_Root_Symp2*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Root_Symp2*Osmotic_TSoil/(1-Rs); else Pi=-1000;
		P_Root_Symp2=Tp+Pi;
		P_Root_Apo2=(Q_Root_Apo2-Q_Root_Apo12)/C_Root_Apo2;
		P_Root_Endo2=(Q_Root_Endo2-Q_Root_Endo02)/C_Root_Endo2;
	}
	if (Root_lower)
	{
		Rs=(Q_Root_Symp03-Q_Root_Symp3)/Q_Root_Symp03;
		Tp=-Pi0_Root_Symp3*Osmotic_TSoil - Epsilon_Root_Symp3*Rs;
		if (Tp<0)Tp=0;
		if (1-Rs) Pi=Pi0_Root_Symp3*Osmotic_TSoil/(1-Rs); else Pi=-1000;
		P_Root_Symp3=Tp+Pi;
		P_Root_Apo3=(Q_Root_Apo3-Q_Root_Apo13)/C_Root_Apo3;
		P_Root_Endo3=(Q_Root_Endo3-Q_Root_Endo03)/C_Root_Endo3;
	}
	 P_Root_Apo=(Q_Root_Apo_t-Q_Root_Apo_t1)/C_Root_Apo;
	 P_Root_Apo1=P_Root_Apo2=P_Root_Apo3=P_Root_Apo;	

		 
	 
}

void print_annual(void)
{
	FILE *out;
	int i;
	
	if ((out = fopen(filename_OUT2,"a"))==NULL) printf("\nCan't create file annual.out!!");
	else
	{
		Get_DATA(dt);
		fprintf(out, "%s; %d;  %.0Lf; ", filenumber, N, YEAR1);
		for (i=0;i<250;i++)  if (File_out[i])  fprintf(out, "%.5Le; ", DATA[i]);
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
		fprintf(out,"%s; %d; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le; %Le;",filenumber,  N, YEAR1, T_PLC_Leaf/3600/24,T_PLC_Axil/3600/24,T_PLC_Branch/3600/24,T_PLC_Trunk/3600/24,T_PLC_Root/3600/24,T_PLC_Root1/3600/24,T_RWC_Axil/3600/24,T_REW_Soil/3600/24, T_gs_regul/3600/24, T_gs_close/3600/24, A_net_tot/1000/1000, A_net_max,  E_max, T_Leaf_max, PLC_Leaf_Apo, PLC_Branch_Apo, PLC_Plant, GPP,Q_Plant0,Q_Plant_s0,Q_Plant_a0, K_Plant_20_0,K_root_0,Q_Branch_Symp0*18/1000/1000,Q_Branch_Apo0*18/1000/1000,Q_Trunk_Symp0*18/1000/1000,Q_Trunk_Apo0*18/1000/1000,Q_Root_Symp0*18/1000/1000*(Root_upper + Root_middle + Root_lower) ,Q_Root_Apo_t0*18/1000/1000, Branch_Area,Trunk_Area,Root_Area);
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



//***************************//
//        MAIN FUNCTION      //
//***************************//

void compute (void) //the main function where the job is done
{
	long i=0;
	FILE *transient;
	long double Q1,Q2, dt_long;
	printf("\nSimulation %d Year %d", N, year++);
	Q1=Q_Trunk_Symp0;
	Q2=Q_Axil_Symp0;
	if (CLIMAT!=0 && CLIMAT!=5) load_climat();
	if (CLIMAT !=1) T0=DOY*24*3600;
	T=T0;
	
	if (DYNAMIC==1) dt_long=60/dt;            // 1 min 
	else            dt_long=600/dt;           // 10 min
	
	Climat(dt_long);
	if (LA_Var) Phenology();
	if (Leaf_Fall==1) Compute_Leaf_Fall();
	if (Leaf_Fall==2) Compute_Root_Fall();
	if (TLEAF) Tleaf();
	Fluidity();
	Compute_ST();
	Compute_T_Osmotic();
	if (T_g_cuti) Compute_g_cuti();

	if (Leaf_Area)
	{
		if (C3) A_net=Net_Photosynthesis_C3();
		else    A_net=Net_Photosynthesis_C4();
	}
	else        A_net=0;
	Rm=Respiration();
  //  E_day(dt,dt_long);
	Compute_Cavitation();
	Respiration();
	soil(dt_long);
	init_Cavitation();
	Compute_K();
	
	
	if (PRINT_SCREEN)
	{
		printf("   \n   Year ");
		if (PRINT_SCREEN==2)  printf("   Days  PLC_Leaf");
		else for (i=0;i<250;i++)  if (Screen_out[i])  printf("%s", Label[i]);
		printf("\n");
	}
	
	while(!DEAD && !END_CLIMAT)
	{
		indice++;
		indice_double+=1;
		T=T0+indice_double*dt;
		

		if (DYNAMIC==1)     // Dynamic model
		{
			Compute_dq(dt);
			Compute_Q();
			Compute_P_dynamic();
			//if (GRAVITY==3)   Gravity_test();// simulates a suden change in trunk xylem potential
			if (DYNAMIC0==2)         // mix model
			{
				if (!PLC_LIMIT && !T_LIMIT)
				{
					DYNAMIC=0;
					dt=dt_stat;
					dt_long=600/dt;
					 indice=indice*dt_dyna/dt_stat+1;
					indice_double=indice_double*dt_dyna/dt_stat+1;
					//if (PRINT_SCREEN) printf("STATIC\n");
				}
				if (REW_t>=REW_crit || !Leaf_Area) PLC_LIMIT=0;

			}
		}
		
		else if (DYNAMIC==0)  // quasi-STEADY
		{
			Climat(dt_long);
			Compute_P_steady(dt);
			Compute_Q_steady();
			Compute_K();
			E_day(dt,dt_long);
			soil(dt);
			Compute_Cavitation();
			if (DYNAMIC0==2)         // mix model
			{
				if (PLC_LIMIT || T_LIMIT)
				{
					DYNAMIC=1;
					dt=dt_dyna;
					dt_long=60/dt;
					 indice=indice*dt_stat/dt_dyna+1;
					indice_double=indice_double*dt_stat/dt_dyna+1;
				  //  if (PRINT_SCREEN) printf("DYNAMIC\n");
					E_day(dt, dt_long);
				}
			}			
			
			if (Leaf_Area)
			{
				if (C3) A_net=Net_Photosynthesis_C3();
				else    A_net=Net_Photosynthesis_C4();
			}
			else A_net=0;
			Rm=Respiration();
		}
		
		//  compute these data only once per mn for Dynamic and once per 10min for Steady
		if(!(indice%(unsigned long)(dt_long)))
		{
			if (DYNAMIC==1)     //long loop for Dynamic model
			{
				Compute_K();
				Climat(dt_long);   // compute climatic variables
				E_day(dt, dt_long);
				Compute_Cavitation();
				soil(dt_long);
				//Top Leaf Photosynthesis rates
				if (Leaf_Area)
				{
					if (C3) A_net=Net_Photosynthesis_C3();
					else    A_net=Net_Photosynthesis_C4();
				}
				else A_net=0;
				Rm=Respiration();
			}
			
			Fluidity();
			Compute_T_Osmotic();
			if (TLEAF) Tleaf();
			if (LA_Var) Phenology();
			if (Leaf_Fall==1) Compute_Leaf_Fall();
			if (Leaf_Fall==2) Compute_Root_Fall();
			if (T_g_cuti)  Compute_g_cuti();
			if (Leaf_Area) Compute_Growth();        // compute trunk and fruit growth only when leafy
			else
			{
				Growth_rate_fruit=0;
				Growth_rate_trunk=0;
			}
			
			if (T_Leaf > T_Leaf_max) T_Leaf_max=T_Leaf;
			if (P_Leaf_Symp < P_min_leaf) P_min_leaf=P_Leaf_Symp;
			if (P_Leaf_Symp < P_min_lf_d)
			{
				P_min_lf_d=P_Leaf_Symp;
				gs_max_d2= g_cuti+g_s;
			}
			if (P_Leaf_Symp > P_max_lf_d) P_max_lf_d=P_Leaf_Symp;
			if (P_Branch_Apo < P_min_stem) P_min_stem=P_Branch_Apo;
			
			if (CUT==1)     E_tot+=dt_long*dt*(E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Axil*Area_Axil + E_Petiole*Petiole_area);  // in mmol
			else if (CUT==2)E_tot+=dt_long*dt*(E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Trunk*Trunk_Area+E_Axil*Area_Axil + E_Petiole*Petiole_area);  // in mmol
			else
			{
				E_tot+=         dt_long*dt*(E_Axil*Area_Axil + E_Petiole*Petiole_area+E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3);  // in mmol
				E_tot_day+=     dt_long*dt*(E_Axil*Area_Axil + E_Petiole*Petiole_area+E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3);  // in mmol
				EvapoT_tot+=    dt_long*dt*(E_Axil*Area_Axil + E_Petiole*Petiole_area+E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3+E_Soil*Soil_Width*Soil_Width);
				EvapoT_day+=    dt_long*dt*(E_Axil*Area_Axil + E_Petiole*Petiole_area+E_Leaf*Leaf_Area+E_Branch*Branch_Area+E_Trunk*Trunk_Area+E_Root1*Root_Area1+E_Root2*Root_Area2+E_Root3*Root_Area3+E_Soil*Soil_Width*Soil_Width);
			}
			
			//Top Leaf Photosynthesis rates
			if (Leaf_Area)
			{
				if (C3) A_net=Net_Photosynthesis_C3();
				else    A_net=Net_Photosynthesis_C4();
			}
			else A_net=0;
			Rm=Respiration();
			if (A_net>A_net_max) A_net_max=A_net;
			
			A_net_tot+=     A_net*dt_long*dt;        //total leaf annual A_net in µmol
			A_net_day+=     A_net*dt_long*dt;        //total leaf daily  A_net in µmol
			Resp_tot+=      Rm*dt_long*dt;            //total lad annual respiration in µmol
			A_gross=        A_net - Rm;                // gross A
			A_gross_tot+=   A_gross*dt_long*dt;    //total leaf annual A_gross in µmol
			
			//Canopy Photosynthesis data
			if (Leaf_Area)  A_net_c=A_net/Extinction_Coeff*(1-exp(-Extinction_Coeff*LAI)); //Canopy photosynthesis per ground area
			else            A_net_c=0;
		
			A_net_tot_c +=A_net_c*Surface_Soil*dt_long*dt;  //total plant canopy annual A_net in µmol
			A_net_day_c +=A_net_c*Surface_Soil*dt_long*dt;  //total plant canopy daily  A_net in µmol
			Resp_tot_c  +=Rm*Leaf_Area*dt_long*dt;          //total plant canopy annual  Respiration in µmol
			
			A_gross_tot_c=A_net_tot_c-Resp_tot_c;
			
			Leaf_C_budget(dt_long);
			//Reserve += (A_net*Leaf_Area - (Rm + Rg)*Q_Wood)*dt_long*dt;  // in µmol
			ETP_Penman_tot+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP totale in mm or Kg
			ETP_Penman_day+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP per day in mm or Kg
			
			if (Leaf_Area)  ETP_leaf_tot+=ETP_Penman*18/1000/1000*dt_long*dt;  //ETP totale in mm or Kg
			if (Type_Axil==2)
			{
				Q_Axil_Symp0+= Growth_rate_fruit*dt_long*dt;                            //plastic growth for a fruit
			  //  Growth_fruit+= (Q_Axil_Symp-Q2);                                        //elastic growth
			}
			Growth_trunk+= (Q_Trunk_Symp-Q1+ Growth_rate_trunk*dt_long*dt);                                                         // cumulative growth in mmol
			Q1=Q_Trunk_Symp;
			Q2=Q_Axil_Symp;
			
			
			//refill soil at XXX%
			if (REHYDRATE)
				if (PLC_Branch_Apo>=PLC_REHYD)
				{
					if (REHYDRATE==1)
					{
						Q_Soil1        =Teta_fc*Volume_soil1*1000*1000*1000/18;
						Q_Soil2        =Teta_fc*Volume_soil2*1000*1000*1000/18;
						Q_Soil3        =Teta_fc*Volume_soil3*1000*1000*1000/18;
						if (!REFILL) PLC_Branch_Apo-=0.01; 							// to stop irrigation after rehydration and allow a new dehydration
					}
					if (REHYDRATE==2)
					{
						Q_Soil1+=Daily_Irr*dt_long*dt/(24*60*60)*Surface_Soil*1000*1000/18; //rehydrate top soil layer with daily irrigation in mm
						Irrigation+=Daily_Irr*dt_long*dt/(24*60*60);
						fill_soil_layers();
					}
				}
			
			// test for END
			if      (END_DEATH== 1 && PLC_Leaf_Apo>=PLC_END) DEAD=1;                                                     // dead when leaf xylem is fully embolised
			else if (END_DEATH== 2 && PLC_Branch_Apo>=PLC_END) DEAD=1;                                              // dead when branch xylem is fully embolised
			else if (END_DEATH== 3 && PLC_Trunk_Apo>=PLC_END) DEAD=1;                                               // dead when trunk xylem is fully embolised
			else if (END_DEATH== 4 && (PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END)) DEAD=1;        // dead when ALL root xylem is fully embolised
			else if (END_DEATH== 5 && (PLC_Leaf_Apo>=PLC_END  ||  PLC_Branch_Apo>=PLC_END ||  PLC_Axil_Apo>=PLC_END || PLC_Trunk_Apo>=PLC_END || PLC_Root_Apo1>=PLC_END || PLC_Root_Apo2>=PLC_END || PLC_Root_Apo3>=PLC_END)) DEAD=1; // dead when ONE organ is fully embolised
			else if (END_DEATH== 6 && (PLC_Leaf_Apo>=PLC_END  &&  PLC_Branch_Apo>=PLC_END &&  PLC_Axil_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
			else if (END_DEATH== 7 && (REW_t<=0.2)) DEAD=1;  // Stop when RWC<0.1
			else if (END_DEATH== 8 && PLC_Axil_Apo>=PLC_END) DEAD=1;
			else if (END_DEATH== 9 && (PLC_Axil_Apo>=PLC_END && PLC_Leaf_Apo>=PLC_END)) DEAD=1;
			else if (END_DEATH==10 && (RWC_Axil<=0.02)) DEAD=1;  // Stop when 2% RWC
			else if (END_DEATH==11 && Type_Axil  && (RWC_Axil<=0.02) && (PLC_Leaf_Apo>=PLC_END  &&  PLC_Branch_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END  &&  PLC_Axil_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END) ) DEAD=1;  // Stop when 2% RWC         
			else if (END_DEATH==11 && !Type_Axil && (PLC_Leaf_Apo>=PLC_END  &&  PLC_Branch_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END && PLC_Root_Apo2>=PLC_END && PLC_Root_Apo3>=PLC_END) ) DEAD=1;  // Stop when 2% RWC
			 else if (END_DEATH==12 && Type_Axil  && (PLC_Leaf_Apo>=PLC_END  &&  PLC_Branch_Apo>=PLC_END &&  PLC_Axil_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
			 else if (END_DEATH==12 && !Type_Axil &&(PLC_Leaf_Apo>=PLC_END  &&  PLC_Branch_Apo>=PLC_END && PLC_Trunk_Apo>=PLC_END && PLC_Root_Apo1>=PLC_END)) DEAD=1; // dead when ALL organs are fully embolised
		
			if (CLIMAT!=0 && CLIMAT!=5)  next_climat();
			if (END_DEATH && DEAD && CLIMAT==2)  purge(); // if PLC=100% skip the end of the year if END_DEATH not zero
 
			if (COMPET) Soil_compet();
			if (REW_t<0.2 && T_REW_Soil>(T-T0))                  	   T_REW_Soil= T-T0;
			if (RWC_Axil<0.02 && T_RWC_Axil>(T-T0))                 T_RWC_Axil= T-T0;
		//	}  //of long loop
		
		if (t_out==24 && (TRANSIENT==1 || TRANSIENT==12 || TRANSIENT==22) ) // if T_out is 24h then print only midday max values at the time of T_max HH1
		{
			if (!(indice%(unsigned long)(3600*t_out/dt)-(long)(3600*HH1/dt)))
			{
				if (PRINT_SCREEN) print_screen();
				if (TRANSIENT) print_transient();
			}
		}
		else if(!(indice%(unsigned long)(3600*t_out/dt)))
		{
			if (PRINT_SCREEN) print_screen();
			if (TRANSIENT) print_transient();
		}
		if (days_simul && T>(DOY+days_simul)*3600*24) DEAD=1; //if days=0 never stops, otherwise compute for x days
		
		if (CLIMAT==0 || CLIMAT==5)  if (!(indice%(unsigned long)(3600*24/dt)))  // a new day
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
		}
		if (!(indice%(unsigned long)(3600*24/dt))) 	if (IRRIGATE == 2 || IRRIGATE == 3) Irrigate(); // irrigate at 24h00
		if (!(indice%(unsigned long)(3600*19/dt)))	if (IRRIGATE == 6 || IRRIGATE == 7) Irrigate(); // irrigate at 19h00
		if (!(indice%(unsigned long)(3600*24/dt)))    indice=0;  // a new day
	   }//of long loop
	if (isnan(E_Leaf)|| isnan(g_s) || isnan(P_Leaf_Apo) || isnan(P_Trunk_Apo) || isnan(P_Branch_Apo) || isnan(P_Root_Apo1) ||isnan(P_Root_Apo2) || isnan(P_Root_Apo3)) exit(1);
	if(P_Leaf_Apo>0.1 || P_Trunk_Apo>0.1 || P_Branch_Apo>0.1 || P_Root_Apo1>0.1 || P_Root_Apo2>0.1 || P_Root_Apo3>0.1)
		{
		 printf(" Overflow");
		 DEAD=1;
		 if (RANDOMISE)setup(); 
		 else exit(1);
		}
	} //of while not dead loop
	
	Get_DATA(dt);
	print_final();
	print_annual();
	 
	if (TRANSIENT==1)
	{ 	
		print_transient();
		transient = fopen(filename_TRANS,"a");
		fprintf(transient,"\n");
		fclose(transient);
	}
	if (PRINT_SCREEN)
	{
		print_screen();
	//	printf("Days to leaf HF=\t %.2Lf \nDays to Branch HF=\t %.2Lf \nDays to Trunk HF=\t %.2Lf \nDays to Root HF=\t %.2Lf \nT_gs_regul= %.2Lf T_gs_close= %.2Lf \nA_net_tot=%Lf E_day_max=%Lf E_day_min=%Lf E_night=%Lf T_max= %Lf\n",T_PLC_Leaf/3600/24,T_PLC_Branch/3600/24,T_PLC_Trunk/3600/24,T_PLC_Root/3600/24,T_gs_regul/3600/24, T_gs_close/3600/24,A_net_tot/1000/1000, E_max, E_max_gs_close, E_Leaf_night, T_Leaf_max);
	}

	if (END_CLIMAT) fclose (climat_in);
	PLC_Leaf_Apo=  0;
	PLC_Branch_Apo=0;
	PLC_Axil_Apo=  0;
	PLC_Trunk_Apo= 0;
	PLC_Root_Apo1= 0;
	PLC_Root_Apo2= 0;
	PLC_Root_Apo3= 0;
	
	year=1;
}

void convert(void)  //compute different parameters and convert Kg to mmol for capacitance Q
{
	Q_Leaf_Apo0   = Succulence/1000 * Leaf_apo_fraction * LA_max_Pheno;  // in Kg
	Q_Leaf_Symp0  = Succulence/1000 * (1- Leaf_apo_fraction) * LA_max_Pheno; // in Kg
	
	// compute the total Q for all organs
	if (FRACTAL==2) //then the total number of axillary organs is N_Axil
	{
		if (Type_Axil==3) //a flower, surface of a disc
		{
			Q_Axil_Symp0    = N_Axil * WC_Axil * 3.1416 /4 * Diam_Axil * Diam_Axil;  // organ surface times surfacic water content  in Kg, WC in kg/m2
			Q_Petiole_Symp0 = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_symp_fraction * 1000; //total petiole symp volume in Kg
			Q_Axil_Apo0     = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_apo_fraction * 1000;
			Petiole_area    = N_Axil *Length_Petiole * Diam_Petiole * 3.1416;
		}
		else if (Type_Axil==2) //  a fruit, volume of a sphere
		{
			Q_Axil_Symp0    = N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
			Q_Petiole_Symp0 = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_symp_fraction * 1000; //total branch symp volume in Kg
			Q_Axil_Apo0     = N_Axil *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_apo_fraction * 1000;
			Petiole_area    = N_Axil *Length_Petiole * Diam_Petiole * 3.1416;
		}
		
		else if (Type_Axil==1)  // a bud  volume of a sphere
		{
			Q_Axil_Symp0   = N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
			Q_Axil_Apo0    = N_Axil *Q_Axil_Symp0/20;   // assume 5%  bud/fruit water content is in the apoplasm.
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
	
	if (FRACTAL==0 || FRACTAL==1) //then the total number of axillary organ is N_Axil * Number_Branch
		{
			if (Type_Axil==3) //a flower, surface of a disc
			{
				Q_Axil_Symp0    = N_Axil *Number_Branch * WC_Axil * 3.1416 /4 * Diam_Axil * Diam_Axil;  // organ surface times surfacic water content  in Kg, WC in kg/m2
				Q_Petiole_Symp0 = N_Axil *Number_Branch * Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_symp_fraction * 1000; //total petiole symp volume in Kg
				Q_Axil_Apo0     = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_apo_fraction * 1000;
				Petiole_area    = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * 3.1416;
			}
			else if (Type_Axil==2) //  a fruit, volume of a sphere
			{
				Q_Axil_Symp0    = N_Axil *Number_Branch *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
				Q_Petiole_Symp0 = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_symp_fraction * 1000; //total branch symp volume in Kg
				Q_Axil_Apo0     = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * Diam_Petiole * 3.1416 / 4 * Branch_apo_fraction * 1000;
				Petiole_area    = N_Axil *Number_Branch *Length_Petiole * Diam_Petiole * 3.1416;
			}
			
			else if (Type_Axil==1)  // a bud  volume of a sphere
			{
				Q_Axil_Symp0   = N_Axil *Number_Branch *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil;  // organ volume times bud volumetric water content  in Kg; WC in m3/m3
				Q_Axil_Apo0    = N_Axil *Number_Branch *Q_Axil_Symp0/20;   // assume 5%  bud/fruit water content is in the apoplasm.
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
	
	Q_Root_Symp0   = Q_Root_Symp_FR; // always use fractal roots
	Q_Root_Apo_t0  = Q_Root_Apo_FR;
	
	if (FRACTAL==1) //based on morphology
	{
		Q_Branch_Symp0= Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 * Branch_symp_fraction * 1000; //total branch symp volume in Kg
		Q_Branch_Apo0 = Number_Branch *Length_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 * Branch_apo_fraction * 1000; //total branch apo volume in Kg
		Q_Trunk_Symp0 = Length_Trunk * Diam_Trunk * Diam_Trunk * 3.1416 / 4 * Trunk_symp_fraction * 1000*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction; //in Kg of H2O in the sapwood
		Q_Trunk_Apo0  = Length_Trunk * Diam_Trunk * Diam_Trunk * 3.1416 / 4 * Trunk_apo_fraction * 1000*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction;  //in Kg of H2O in the sapwood
	}
	else // Fractal tree
	{
		Q_Branch_Symp0 = Q_Branch_Symp_FR;
		Q_Branch_Apo0  = Q_Branch_Apo_FR;
		Q_Trunk_Symp0  = Q_Trunk_Sym_FR;
		Q_Trunk_Apo0   = Q_Trunk_Apo_FR;
	}
	Q_Plant_s0= (Q_Axil_Symp0 + Q_Petiole_Symp0 + Q_Leaf_Symp0  + Q_Branch_Symp0 + Q_Trunk_Symp0 + Q_Root_Symp_FR*(Root_upper + Root_lower + Root_middle ));
	Q_Plant_a0= (Q_Leaf_Apo0 + Q_Axil_Apo0 + Q_Branch_Apo0 + Q_Trunk_Apo0 + Q_Root_Apo_FR*(Root_upper + Root_lower + Root_middle ));	
	Q_Plant0= Q_Plant_a0 + Q_Plant_s0;
	
	//convert from volume specific capacitance in kg MPa-1 L-1 to capacitance  in kg MPa-1 with water volume
	C_Leaf_Apo   *=  Q_Leaf_Apo0;
	C_Branch_Apo *=  Q_Branch_Apo0;
	C_Trunk_Apo  *=  Q_Trunk_Apo0;
	C_Root_Apo   *=  Q_Root_Apo_t0;
	C_Axil_Apo   *=  Q_Axil_Apo0;
	
	Diam_Root    = Diam_Root_FR;
	Length_Root = Length_Root_FR;
	Root_Area   = Root_Area_FR;
	
	if (FRACTAL==1)
	{
		Branch_Area   = Number_Branch *  Length_Branch * Diam_Branch * 3.1416;
		Trunk_Area    = Length_Trunk * Diam_Trunk * 3.1416;
	}
	else
	{
		Branch_Area   = Branch_Area_FR;
		Trunk_Area    = Trunk_Area_FR;
	}
	
	if (Type_Axil==3) Area_Axil      = 3.1416*Diam_Axil*Diam_Axil/4;  // a flower, surface of a disc
	else              Area_Axil      = 3.1416*Diam_Axil*Diam_Axil;  // a bud or a fruit, surface of a sphere
	if (FRACTAL==2)      Area_Axil*= N_Axil;
	else              Area_Axil*= (N_Axil*Number_Branch);
	
	K_Leaf_Symp_0  *= LA_max_Pheno;
	K_Root_Apo0    = K_Root_Apo_FR;
	K_Trunk_Symp0   *= Trunk_Area;
	K_Root_Symp1    = K_Root_Symp0 * Root_Area_fi;  // only the fine roots absorb water
	K_Root_Symp2    = K_Root_Symp0 * Root_Area_FR;  // all the root surface has a symplasmic compartment
	K_Leaf_Apo0    *= LA_max_Pheno;  // from Kl to k
	
	if (FRACTAL==1)
	{
		K_Branch_Symp0  *= Number_Branch * Branch_Area;
		K_Branch_Apo0  *= Number_Branch * Diam_Branch * Diam_Branch * 3.1416 / 4 /Length_Branch; // from Ks to k
		K_Trunk_Apo0   *= Diam_Trunk * Diam_Trunk * 3.1416 / 4 /Length_Trunk*(2-Trunk_sapwood_fraction)*Trunk_sapwood_fraction; // from Ks to k
		K_Axil_Apo0     *=(N_Axil*Number_Branch);
		K_Axil_Symp     *=(N_Axil*Number_Branch);
		K_Axil_Symp2    *=(N_Axil*Number_Branch);
	}
	else
	{
		K_Branch_Symp0  *= Branch_Area;
		K_Branch_Apo0  =  K_Branch_Apo_FR; // from Ks to k
		K_Trunk_Apo0   =  K_Trunk_Apo_FR; // from Ks to k
		K_Axil_Apo0    *=N_Axil;
		K_Axil_Symp    *=N_Axil;
		K_Axil_Symp2   *=N_Axil;
	}
	
	C_Leaf_Apo    *=   1000*1000/18; //convert Kg to mmol for capacitance & Q
	C_Branch_Apo  *=   1000*1000/18;
	C_Trunk_Apo   *=   1000*1000/18;
	C_Root_Apo    *=   1000*1000/18;
	C_Axil_Apo    *=   1000*1000/18;
	
	Q_Leaf_Symp0  *=   1000*1000/18;
	Q_Branch_Symp0*=   1000*1000/18;
	Q_Trunk_Symp0 *=   1000*1000/18;
	Q_Root_Symp0  *=   1000*1000/18;
	Q_Axil_Symp0  *=   1000*1000/18;
	Q_Petiole_Symp0*=   1000*1000/18;
	
	Q_Leaf_Apo0   *=    1000*1000/18;
	Q_Branch_Apo0 *=    1000*1000/18;
	Q_Trunk_Apo0  *=    1000*1000/18;
	Q_Root_Apo_t0   *=    1000*1000/18;
	Q_Axil_Apo0   *=    1000*1000/18;
	
	Growth_trunk=Q_Trunk_Symp0;

	Growth_fruit=Q_Axil_Symp0;
}

void Reset (void) // When it is a new run or a new year. Reset simulation to zero
{
	FILE *transient;

	if (DYNAMIC0==0) dt=dt_stat;
	if (DYNAMIC0==2) dt=dt_stat;
	if (DYNAMIC0==1) dt=dt_dyna;  
	if (TRANSIENT)
	{
		transient = fopen(filename_TRANS,"a");
		fprintf(transient, "\n");
		fclose(transient);
	}
	
	print_annual();
	indice=1;
	T_1=T_2;
	T_air_min=T_air_min_2;  //set values for next day loaded before
	T_air_max=T_air_max_2;
	RH_air_min=RH_air_min_2;
	RH_air_max=RH_air_max_2;
	PAR_max=PAR_max_2;
	Rain_1=Rain_2;
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
		

	PLC_Leaf_Apo*=leg_Leaf;
	PLC_Axil_Apo*=leg_Leaf;  //same as for leaf
	PLC_Branch_Apo*=leg_Branch;
	PLC_Trunk_Apo*=leg_Trunk;
	PLC_Root_Apo1*=leg_Root;
	PLC_Root_Apo2*=leg_Root;
	PLC_Root_Apo3*=leg_Root;
	ETP_Penman_tot=0;
	
	Q_Leaf_Apo=Q_Leaf_Apo0;
	Q_Axil_Apo=Q_Axil_Apo0;
	Q_Axil_Apo1=Q_Axil_Apo0;
	Q_Branch_Apo=Q_Branch_Apo0;
	Q_Trunk_Apo=Q_Trunk_Apo0;
	Q_Root_Apo1=Q_Root_Apo01;
	Q_Root_Apo2=Q_Root_Apo02;
	Q_Root_Apo3=Q_Root_Apo03;
	Q_Leaf_Apo1=Q_Leaf_Apo0;
	Q_Branch_Apo1=Q_Branch_Apo0;
	Q_Trunk_Apo1=Q_Trunk_Apo0;
	Q_Root_Apo11=Q_Root_Apo01;
	Q_Root_Apo12=Q_Root_Apo02;
	Q_Root_Apo13=Q_Root_Apo03;
	Q_Root_Apo_t=Q_Root_Apo_t0;
	
	E_tot=0;EvapoT_tot=0;A_net_tot=0;A_net_tot_c=0;ETP_Penman_tot=0;Irrigation=0;RWC_int=0;P_soil_int=0;P_soil_min=0;P_min_stem=P_min_leaf=0;
	RWC_min=1;Rain_tot=0;Rain_soil=0;VPD_Air_tot=0;VPD_Leaf_tot=0;Rain_leaf_tot=0;ETP_leaf_tot=0;Cum_T_air_l=0;g_Soil=g_Soil0;A_net_max=0;
	Drainage=0;A_gross_tot=0;A_gross_tot_c=0;Resp_tot=0;Resp_tot_c=0;A_net_day=0;E_tot_day=0;EvapoT_day=0;
	Q_Axil_Symp0=N_Axil *WC_Axil * 3.1416 /6* 1000 * Diam_Axil * Diam_Axil * Diam_Axil*1000*1000/18;
	Q_Axil_Symp=Q_Axil_Symp0;
	init();
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

void Fractal_loop_root (int ordre,int i)
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
			if ((rand()%1000)<((N_fractal-1)*1000)) Fractal_loop_root(ordre+1,j);
		}
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",Diam_Root*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_fine_root/100*scale*cos(side[ordre-1]*angle+3.1416/2),y-Length_fine_root/100*scale*sin(side[ordre-1]*angle+3.1416/2));	
	
	//fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_fine_root/100*scale,y-Length_fine_root/100*scale);	
	
}

void Fractal(void)  //compute a fractal tree for the whole plant when FRACTAL=2 or for the root only when FRACTAL=0;
{
	long double SEGMENTS=100; //discretisation of the shoot and root systems 100 is enough
	long double Diam_term_sapwood, Diam_trunk_sapwood;
	long double A_term_root, Length_term_root,Diam_taproot, Length_root, Length_taproot, Diam_taproot_sapwood,N_root, Mean_diam_root;
	long double Nb_term_shoot,Nb_term_root,Fineroot_area, Taproot_area, Volume_fineroot, Volume_taproot, Volume_fineroot_sapwood, Volume_taproot_sapwood, R_fineroot, R_taproot;
	long double F,Delta1,Delta2;
	long double Discriminant, Node_base_crown;
	long double Volume_branch ,Volume_branch_sapwood, Volume_trunk ,Volume_trunk_sapwood, R_branch, R_trunk; 
	long double l_branch,leaf; //side1,side2,side3,side4,side5,side6,side7,side8;
	int i,j;
	long double X0=108,Y0=100,XMAX=216,YMAX=279,Y_crown,x,y;
	
	
	
	for (i=0;i<=1000;i++) side[i]=1;
	
// for shoots
if(FRACTAL==2) //shoot is fractal
{
	Diam_term_sapwood=Diam_Branch;
	Diam_trunk_sapwood=pow(Diam_Trunk*Trunk_sapwood_fraction*(2*Diam_Trunk-Diam_Trunk*Trunk_sapwood_fraction),0.5);
	if (!Length_term_shoot) Length_term_shoot=(Length_Branch+Length_Trunk)/SEGMENTS; else Length_term_shoot=Length_term_shoot/100;
	Nb_term_shoot= LA_max_init / LA_term_shoot*10000;
	Nb_segment[0]=Nb_term_shoot;
	Diam_segment[0]=Diam_Branch;
	Length_segment[0]=Length_term_shoot;
	Branch_Area_FR=Length_term_shoot*3.1416*Diam_segment[0]*Nb_segment[0];
	Trunk_Area_FR=0;
	Volume_branch=Nb_segment[0]*Diam_segment[0]*Diam_segment[0]*3.1416/4*Length_segment[0];
	Volume_branch_sapwood=Nb_segment[0]*Diam_term_sapwood*Diam_term_sapwood*3.1416/4*Length_segment[0];
	Volume_trunk=0;
	Volume_trunk_sapwood=0;

	Q_Branch_Symp_FR=Volume_branch_sapwood*Branch_symp_fraction*1000; //in Kg
	Q_Branch_Apo_FR=Volume_branch_sapwood*Branch_apo_fraction*1000;
	Q_Trunk_Sym_FR=Volume_trunk_sapwood*Trunk_symp_fraction*1000;
	Q_Trunk_Apo_FR=Volume_trunk_sapwood*Trunk_apo_fraction*1000;
	R_branch=Length_segment[0]/(Diam_term_sapwood*Diam_term_sapwood*3.1416/4)*1/K_Branch_Apo0/Nb_term_shoot*pow((Diam_term_sapwood/Diam_term_sapwood),Tapering);
	R_trunk=0;
	
	F=2*(Length_Branch+Length_Trunk-Length_term_shoot*(SEGMENTS+1))/((SEGMENTS+1)*SEGMENTS);
	Discriminant=pow((F+2*Length_term_shoot),2)+8*F*Length_Branch;
	Node_base_crown=(-(F+2*Length_term_shoot)+pow(Discriminant,0.5))/(2*F);
	N_fractal=pow(Nb_term_shoot,1/Node_base_crown);
	Delta1=Node_base_crown*log(N_fractal)/log(Diam_Trunk/Diam_Branch);
	Delta2=Node_base_crown*log(N_fractal)/log(Diam_trunk_sapwood/Diam_term_sapwood);
	
	for(i=1;i<(int)SEGMENTS;i++)
	{		
		Nb_segment[i]=max(Nb_segment[i-1]/N_fractal,1);
		Length_segment[i]=Length_term_shoot+i*F;
		if ((long double)i<Node_base_crown) // for the branches
		{
				Diam_segment[i]=(Diam_Branch)*pow(N_fractal,(i/Delta1));
				Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(i/Delta2));
				Volume_branch+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
				Volume_branch_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
				Branch_Area_FR+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
				R_branch+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Branch_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);	
		}
		else // for the trunk
		{
			Diam_segment[i]=(Diam_Branch)*pow(N_fractal,(Node_base_crown/Delta1));
			Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(Node_base_crown/Delta2));
			Trunk_Area_FR+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
			Volume_trunk+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
			Volume_trunk_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
			R_trunk+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Branch_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
		}		
	}
	
	Q_Branch_Symp_FR=Volume_branch_sapwood*Branch_symp_fraction*1000;  //in Kg
	Q_Branch_Apo_FR=Volume_branch_sapwood*Branch_apo_fraction*1000;
	Q_Trunk_Sym_FR=Volume_trunk_sapwood*Trunk_symp_fraction*1000;
	Q_Trunk_Apo_FR=Volume_trunk_sapwood*Trunk_apo_fraction*1000;
	K_Branch_Apo_FR=1/R_branch;
	K_Trunk_Apo_FR=1/R_trunk;	
}

// Fractal map1
    if (PRINT_GRAPH && FRACTAL==2) //fractal shoot
	{
	srand((unsigned)time(NULL));
	YMAX=200;
	angle=3.1416/8;
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	leaf=pow(LA_term_shoot,0.5)/100;
	File_fractal1 = fopen("fractal_map1.ps", "w+" );   
	fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal1);
	fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal1);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal1);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal1);

	//soil limit
	fprintf(File_fractal1,"stroke 0.5 M setlinewidth 1 0 0 setrgbcolor  newpath ");
	fprintf(File_fractal1,"0 M %Le M %Le M %Le M T\n",Y[0],XMAX,Y[0]);

	//trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal1);
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //trunk
	{
		fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_segment[i]*scale);
		Y[0]+=Length_segment[i]*scale;
	}
	Y_crown=Y[0];

	//branch wood
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
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //trunk
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
	
	l_branch=0;
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(Node_base_crown);i>=0;i--) //branches
	{
		for (j=1;j<=Nb_segment[i];j++)
		{		
		x=X[0]+l_branch*scale*cos(j*3.1416/Nb_segment[i]);
		y=Y_crown+l_branch*scale*sin(j*3.1416/Nb_segment[i]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[i]*scale*cos(j*3.1416/Nb_segment[i]),y+Length_segment[i]*scale*sin(j*3.1416/Nb_segment[i]));
		}
		l_branch+=Length_segment[i];
	}
	
	for (j=0;j<Nb_segment[0];j++)//leaves
	{	
		x=X[0]+(l_branch)*scale*cos(j*3.1416/Nb_segment[0]);
		y=Y_crown+(l_branch)*scale*sin(j*3.1416/Nb_segment[0]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",leaf*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+leaf*scale*cos(j*3.1416/Nb_segment[0]),y+leaf*scale*sin(j*3.1416/Nb_segment[0]));	
	}
			
}

if (PRINT_GRAPH && FRACTAL==1) //shoot is morphometric
{
	YMAX=200;	
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	leaf=pow(LA_max_init/Number_Branch,0.5);
	File_fractal1 = fopen("fractal_map1.ps", "w+" );   
	fputs("%!PS-Adobe-2.0\n%SurEau model Creator: Herve Cochard INRAE-PIAF\n",File_fractal1);
	fputs("%DocumentFonts: Courier\n%%EndProlog\n%Page: 1 1\n",File_fractal1);
	fputs("/M {\n  0.3527 div\n  } def\n\n 1 setlinejoin  %coins arrondis\n",File_fractal1);
	fputs("/T { newpath\n  moveto\n  lineto\n  stroke\n } def \n\n",File_fractal1);
	
	//soil limit

	fprintf(File_fractal1,"stroke 0.5 M setlinewidth 1 0 0 setrgbcolor  newpath ");
	fprintf(File_fractal1,"0 M %Le M %Le M %Le M T\n",Y[0],XMAX,Y[0]);
	//trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal1);
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_Trunk*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
	Y[0]+=Length_Trunk*scale;
	

	//heartwood
	fputs("\n%%Heartwood\n",File_fractal1);
	X[0]=X0;Y[0]=Y0;	
	fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.67 0 0 setrgbcolor newpath ",Diam_Trunk*Trunk_sapwood_fraction*scale);
	fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
			
	//branch wood
	Y[0]+=Length_Trunk*scale;
	fputs("%Branches\n",File_fractal1);
	for(i=1;i<=Number_Branch;i++)  //ordre 1
		{
			angle=3.1416/(Number_Branch+1);
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_Branch*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0]+Length_Branch*scale*cos(3.1416-i*angle),Y[0]+Length_Branch*scale*sin(3.1416-i*angle));
			//leaf
			fprintf(File_fractal1,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",leaf*scale);
			fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0]+Length_Branch*scale*cos(3.1416-i*angle),Y[0]+Length_Branch*scale*sin(3.1416-i*angle),X[0]+(Length_Branch+leaf)*scale*cos(3.1416-i*angle),Y[0]+(Length_Branch+leaf)*scale*sin(3.1416-i*angle));	

		}
	fputs("\nstroke\nshowpage\n%%Trailer\n",File_fractal1);
	fclose(File_fractal1);
	}	
		
// for one root
	A_term_root = LA_max_init/3; //assume the area of fine roots is equal to the leaf area cm2
	Length_root=(Length_Branch+Length_Trunk)*root_shoot_ratio;
	Length_taproot= Length_root*(1-root_ramif);
	
	Diam_term_sapwood=Diam_Root ;
	Diam_taproot=Diam_Trunk/pow(3,0.5); //assume the cross section of the 3 taproots=dbh
	Diam_taproot_sapwood=pow(Diam_taproot*Trunk_sapwood_fraction*(2*Diam_taproot-Diam_taproot*Trunk_sapwood_fraction),0.5);  //assume same sapwood fraction in the taproot that in the trunk
	//Length_term_root=Length_root/SEGMENTS; 
	Length_term_root=Length_fine_root/100;
	Nb_term_root= A_term_root / (Length_fine_root/100*Diam_Root*3.1416);
	Nb_segment[0]=Nb_term_root;
	N_root=Nb_segment[0];
	
	Diam_segment[0]=Diam_Root;
	Mean_diam_root=Diam_segment[0]*Nb_segment[0];
	Length_segment[0]=Length_term_root;
	Fineroot_area=Length_term_root*3.1416*Diam_segment[0]*Nb_segment[0];
	Taproot_area=0;
	Volume_fineroot=Nb_segment[0]*Diam_segment[0]*Diam_segment[0]*3.1416/4*Length_segment[0];
	Volume_fineroot_sapwood=Nb_segment[0]*Diam_term_sapwood*Diam_term_sapwood*3.1416/4*Length_segment[0];
	Volume_taproot=0;
	Volume_taproot_sapwood=0;
	Q_Root_Symp_FR=Volume_fineroot_sapwood*Root_symp_fraction*1000+Volume_taproot_sapwood*Trunk_symp_fraction*1000;  //in Kg
	Q_Root_Apo_FR=Volume_fineroot_sapwood*Root_apo_fraction*1000 + Volume_taproot_sapwood*Trunk_apo_fraction*1000;
	R_fineroot=Length_segment[0]/(Diam_term_sapwood*Diam_term_sapwood*3.1416/4)*1/K_Root_Apo0/Nb_term_root*pow((Diam_term_sapwood/Diam_term_sapwood),Tapering);
	R_taproot=0;
	
	F=2*(Length_root-Length_term_root*(SEGMENTS+1))/((SEGMENTS+1)*SEGMENTS);
	Discriminant=pow((F+2*Length_term_root),2)+8*F*(Length_root-Length_taproot);
	Node_base_crown=(-(F+2*Length_term_root)+pow(Discriminant,0.5))/(2*F);
	N_fractal=pow(Nb_term_root,1/Node_base_crown);
	Delta1=Node_base_crown*log(N_fractal)/log(Diam_taproot/Diam_Root);
	Delta2=Node_base_crown*log(N_fractal)/log(Diam_taproot_sapwood/Diam_term_sapwood);
	Length_Root_FR=Nb_segment[0]*Length_segment[0];
	
	for(i=1;i<(int)SEGMENTS;i++)
	{
		Nb_segment[i]=max(Nb_segment[i-1]/N_fractal,1);
		N_root+=Nb_segment[i];
		Length_segment[i]=Length_term_root+i*F;
		if ((long double)i<Node_base_crown) // for the fineroot
		{
				Diam_segment[i]=(Diam_Root)*pow(N_fractal,(i/Delta1));
				Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(i/Delta2));
				Volume_fineroot+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
				Volume_fineroot_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
				Fineroot_area+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
				R_fineroot+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Root_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
				Length_Root_FR+=Nb_segment[i]*Length_segment[i];
				Mean_diam_root+=Diam_segment[i]*Nb_segment[i];
		}
		else // for the taproot
		{
			Diam_segment[i]=(Diam_Root)*pow(N_fractal,(Node_base_crown/Delta1));
			Diam_sapwood[i]=(Diam_term_sapwood)*pow(N_fractal,(Node_base_crown/Delta2));
			Taproot_area+=Length_segment[i]*3.1416*Diam_segment[i]*Nb_segment[i];
			Volume_taproot+=Nb_segment[i]*Diam_segment[i]*Diam_segment[i]*3.1416/4*Length_segment[i];
			Volume_taproot_sapwood+=Nb_segment[i]*Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4*Length_segment[i];
			R_taproot+=Length_segment[i]/(Diam_sapwood[i]*Diam_sapwood[i]*3.1416/4)*1/K_Root_Apo0/Nb_segment[i]*pow((Diam_term_sapwood/Diam_sapwood[i]),Tapering);
			Length_Root_FR+=Nb_segment[i]*Length_segment[i];
			Mean_diam_root+=Diam_segment[i]*Nb_segment[i];
		}	
	}
	Q_Root_Symp_FR=Volume_fineroot_sapwood*Root_symp_fraction*1000+Volume_taproot_sapwood*Trunk_symp_fraction*1000;  //in Kg
	Q_Root_Apo_FR=Volume_fineroot_sapwood*Root_apo_fraction*1000 + Volume_taproot_sapwood*Trunk_apo_fraction*1000;
	K_Root_Apo_FR=1/(R_fineroot+R_taproot);
	Root_Area_FR=Fineroot_area+Taproot_area;
	Diam_Root_FR=Mean_diam_root/N_root;
	Root_Area_fi_0=A_term_root;
	Length_Root_fi=Root_Area_fi_0/(3.1416*Diam_Root);

if (PRINT_GRAPH) // Fractal map2
{
	if (FRACTAL==2) //print fractal root in File_fractal1
	{
	angle=3.1416/8;
	X[0]=X0;Y[0]=Y0;
	fputs("%TapRoot\n",File_fractal1);
	fprintf(File_fractal1,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //tap root
	{
		fprintf(File_fractal1,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal1,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]-Length_segment[i]*scale);
		Y[0]-=Length_segment[i]*scale;
	}
	Y_crown=Y[0];
	fputs("%Fine roots\n",File_fractal1);
	l_branch=Length_segment[0];
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
				Fractal_loop_root(1,i);
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
//root	
	fputs("%TapRoot\n",File_fractal2);
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
	for(i=(int)(SEGMENTS-1);i>(int)Node_base_crown;i--) //tap root
	{
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]-Length_segment[i]*scale);
		Y[0]-=Length_segment[i]*scale;
	}
	Y_crown=Y[0];
	fputs("%Fine roots\n",File_fractal2);
	l_branch=Length_segment[0];
	fprintf(File_fractal2,"stroke 0.5 0 0 setrgbcolor  newpath ");
		for(i=(int)(Node_base_crown);i>=0;i--) //fine roots
		{
			for (j=1;j<=Nb_segment[i];j++)
			{		
			x=X[0]+l_branch*scale*cos(j*3.1416/Nb_segment[i]);
			y=Y_crown-l_branch*scale*sin(j*3.1416/Nb_segment[i]);
			fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_segment[i]*scale);
			fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_segment[i]*scale*cos(j*3.1416/Nb_segment[i]),y-Length_segment[i]*scale*sin(j*3.1416/Nb_segment[i]));
			}
			l_branch+=Length_segment[i];
		}
	fprintf(File_fractal2,"stroke  0 0.5 0 setrgbcolor newpath ");	
	for (j=0;j<Nb_term_root;j++)//terminal roots
	{	
		x=X[0]+(l_branch)*scale*cos(j*3.1416/Nb_segment[0]);
		y=Y_crown-(l_branch)*scale*sin(j*3.1416/Nb_segment[0]);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth newpath ",Diam_Root*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",x,y,x+Length_term_root*scale*cos(j*3.1416/Nb_segment[0]),y-Length_term_root*scale*sin(j*3.1416/Nb_segment[0]));	
	}
	

if (FRACTAL==1)
	{
	YMAX=200;	
	scale=0.5*YMAX/(Length_Branch+Length_Trunk);
	X[0]=X0;Y[0]=Y0;
	leaf=pow(LA_max_init/Number_Branch,0.5);

	//trunk 
	fputs("%!PS-Adobe-2.0\n%%Trunk\n",File_fractal2);
	fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor  newpath ",Diam_Trunk*scale);
	fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
	Y[0]+=Length_Trunk*scale;
	

	//heartwood
	fputs("\n%%Heartwood\n",File_fractal1);
	X[0]=X0;Y[0]=Y0;	
	fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.67 0 0 setrgbcolor newpath ",Diam_Trunk*Trunk_sapwood_fraction*scale);
	fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0],Y[0]+Length_Trunk*scale);
			
	//branch wood
	Y[0]+=Length_Trunk*scale;
	fputs("%Branches\n",File_fractal2);
	for(i=1;i<=Number_Branch;i++)  //ordre 1
	{
		angle=3.1416/(Number_Branch+1);
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0.5 0 0 setrgbcolor newpath ",Diam_Branch*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0],Y[0],X[0]+Length_Branch*scale*cos(3.1416-i*angle),Y[0]+Length_Branch*scale*sin(3.1416-i*angle));
		//leaf
		fprintf(File_fractal2,"stroke %Lf M setlinewidth 0 0.5 0 setrgbcolor newpath ",leaf*scale);
		fprintf(File_fractal2,"%Le M %Le M %Le M %Le M T\n",X[0]+Length_Branch*scale*cos(3.1416-i*angle),Y[0]+Length_Branch*scale*sin(3.1416-i*angle),X[0]+(Length_Branch+leaf)*scale*cos(3.1416-i*angle),Y[0]+(Length_Branch+leaf)*scale*sin(3.1416-i*angle));	
	}
}	

	fputs("\nstroke\nshowpage\n%%Trailer\n",File_fractal2);
	fclose(File_fractal2);
	}
}

void setparameters(void)
{
	DYNAMIC0=para[1]; dPLC_crit=para[2]; REW_crit=para[3]; debug=para[4]; TRANSIENT=para[5]; PRINT_SCREEN=para[6]; PRINT_GRAPH=para[7]; gs_cst=para[8]; K_VAR=para[9]; 
	CUT=para[10]; END_DEATH=para[11]; PLC_END=para[12]; days_simul=para[13]; dt_dyna=para[14]; dt_stat=para[15]; t_out=para[16]; CLIMAT=para[17]; TLEAF=para[18]; FLUID=para[19];
	SURFACE_TENSION=para[20]; T_OSMOTIC=para[21]; T_g_cuti=para[22]; GRAVITY=para[23]; CAPILLARITY=para[24]; Regul_gs=para[25]; Regul_gs_para1=para[26]; Regul_gs_para2=para[27]; 
	Pgs_12=para[28]; Pgs_88=para[29]; DOY=para[30]; Lat=para[31]; T_air_min=para[32]; T_air_max=para[33]; RH_air_min=para[34]; RH_air_max=para[35]; PAR_max=para[36]; Wind=para[37]; 
	HW=para[38]; HW_day=para[39]; HW_duration=para[40]; HW_T=para[41]; Teta_s=para[42]; Teta_r=para[43]; alpha=para[44]; n=para[45]; K_sat=para[46]; L=para[47]; Rock_f1=para[48]; 
	PENMAN=para[49]; Penman_Coeff=para[50]; g_crown0=para[51]; INTERCEPTION=para[52]; Interception_min=para[53]; Threshold_rain=para[54]; Interception_factor=para[55]; IRRIGATE=para[56]; 
	RWC_Irr=para[57]; Daily_Irr=para[58]; IRR_DOY_S=para[59]; IRR_DOY_F=para[60]; REHYDRATE=para[61]; PLC_REHYD=para[62]; CONTINUOUS=para[63]; Extinction_Coeff=para[64]; Leaf_size=para[65]; 
	leaf_angle=para[66]; LA_Var=para[67]; LA_max_init=para[68]; LA_min=para[69]; LA_day1=para[70]; LA_day2=para[71]; LA_day3=para[72]; LA_day4=para[73]; Leaf_Fall=para[74]; 
	P50_Leaf_Fall=para[75]; Slope_Leaf_Fall=para[76]; Succulence=para[77]; Leaf_apo_fraction=para[78]; LMA=para[79]; VcMax=para[80]; VjMax=para[81]; CO2_atm=para[82]; Rd25=para[83]; 
	Qye=para[84]; Kc25=para[85]; Ko25=para[86]; a_Res=para[87]; Number_Branch=para[88]; Length_Branch=para[89]; Diam_Branch=para[90]; Branch_apo_fraction=para[91]; 
	Branch_symp_fraction=para[92]; Density=para[93]; Length_Trunk=para[94]; Diam_Trunk=para[95]; Trunk_apo_fraction=para[96]; Trunk_symp_fraction=para[97]; Trunk_sapwood_fraction=para[98];
	Length_Root_fi=para[99]; Diam_Root=para[100]; Root_apo_fraction=para[101]; Root_symp_fraction=para[102]; FRACTAL=para[103]; Q_Branch_Symp_FR=para[104]; Q_Branch_Apo_FR=para[105]; 
	Q_Trunk_Sym_FR=para[106]; Q_Trunk_Apo_FR=para[107]; Q_Root_Symp_FR=para[108]; Q_Root_Apo_FR=para[109]; Branch_Area_FR=para[110]; Trunk_Area_FR=para[111]; Root_Area_FR=para[112]; 
	Root_Area_fi_0=para[113]; Length_Root_FR=para[114]; Diam_Root_FR=para[115]; K_Branch_Apo_FR=para[116]; K_Trunk_Apo_FR=para[117]; K_Root_Apo_FR=para[118]; gs_max=para[119]; 
	gs_night=para[120]; Jarvis_PAR=para[121]; GS_MAX=para[122]; Tgs_optim=para[123]; Tgs_sens=para[124]; g_cuti_20=para[125]; TP=para[126]; Q10_1=para[127]; Q10_2=para[128]; 
	g_Branch=para[129]; g_Trunk=para[130]; g_Root=para[131]; g_Soil0=para[132]; Extensibility_trunk=para[133]; Yield_trunk=para[134]; C_Leaf_Apo=para[135]; C_Branch_Apo=para[136]; 
	C_Trunk_Apo=para[137]; C_Root_Apo=para[138]; Soil_Depth=para[139]; Soil_Width=para[140]; Root_upper=para[141]; Root_middle=para[142]; Root_lower=para[143]; gap=para[144]; 
	Epsilon_Leaf_Symp=para[145]; Epsilon_Branch_Symp=para[146]; Epsilon_Trunk_Symp=para[147]; Epsilon_Root_Symp=para[148]; Pi0_Leaf_Symp=para[149]; Pi0_Branch_Symp=para[150]; 
	Pi0_Trunk_Symp=para[151]; Pi0_Root_Symp=para[152]; K_Leaf_Apo0=para[153]; K_Branch_Apo0=para[154]; K_Trunk_Apo0=para[155]; K_Root_Apo0=para[156]; K_Leaf_Symp_0=para[157]; 
	K_Branch_Symp0=para[158]; K_Trunk_Symp0=para[159]; K_Root_Symp0=para[160]; P50_Leaf_Apo=para[161]; P50_Branch_Apo=para[162]; P50_Trunk_Apo=para[163]; P50_Root_Apo=para[164]; 
	Slope_Leaf_Apo=para[165]; Slope_Branch_Apo=para[166]; Slope_Trunk_Apo=para[167]; Slope_Root_Apo=para[168]; Type_Axil=para[169]; N_Axil=para[170]; K_Axil_Apo0=para[171]; 
	K_Axil_Symp=para[172]; K_Axil_Symp2=para[173]; Epsilon_Axil_Symp=para[174]; Pi0_Axil_Symp=para[175]; P50_Axil_Apo=para[176]; Slope_Axil_Apo=para[177]; g_Axil_min20=para[178]; 
	g_Axil_max=para[179]; g_Petiole=para[180]; C_Axil_Apo=para[181]; Diam_Axil=para[182]; Length_Petiole=para[183]; Diam_Petiole=para[184]; WC_Axil=para[185]; Extensibility_fruit=para[186]; 
	Yield_fruit=para[187]; REFILL=para[188]; P_REFILL=para[189]; leg_Leaf=para[190]; leg_Branch=para[191]; leg_Trunk=para[192]; leg_Root=para[193]; T_Soil_Crit=para[194]; 
	COMPET=para[195]; Rock_f2=para[196]; Rock_f3=para[197]; gs_CO2_sens=para[198]; HH1=para[199]; HH2=para[200]; K_VAR_P1=para[201]; K_VAR_P2=para[202]; K_VAR_P3=para[203]; gs_tc=para[204];
	RANDOMISE=para[205];LA_term_shoot=para[206];Tapering=para[207];root_shoot_ratio=para[208];root_ramif=para[209];Length_fine_root=para[210];Length_term_shoot=para[211];
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
	
		if (REFILL==2) 
		{
			SYMP_CAVIT=1;
			REFILL=0;
		}

		if (Regul_gs ==6) 
			{
				Gamma= Regul_gs_para1;			
				PLCx = Regul_gs_para2;
			}
		if (Regul_gs ==2 || Regul_gs ==3 || Regul_gs ==4 || Regul_gs ==5  )  
			{
				Gamma= Regul_gs_para1;
				Px_gs= Regul_gs_para2;
			}
		if (Regul_gs ==2 || Regul_gs ==15 ) turgor_ref_factor=Regul_gs_para1;		
			
		if (DYNAMIC0==0)    //Steady
			{
				dt=dt_stat;
				DYNAMIC=0;
			//	if (PRINT_SCREEN) printf("quasi-STEADY\n");
			}
			
		else if (DYNAMIC0==1)    //Dynamic
			{
				dt=dt_dyna;
				DYNAMIC=1;
			//	if (PRINT_SCREEN) printf("DYNAMIC\n");
			}
		else if (DYNAMIC0==2)    //Mix start with Steady
			{
				dt=dt_dyna;
				DYNAMIC=1;
		//		if (PRINT_SCREEN)  printf("MIX model\n");
			}
		if (PAR_max==-1) PAR_max=Potential_PAR((12)); // if -1 then use the Potential PAR to compute PAR_max
		if (FRACTAL) Fractal();	
		g_Axil_min=g_Axil_min20;
		LA_max=LA_max_init;
		LA_max_Pheno=LA_max;
		Root_Area_fi=Root_Area_fi_0;
		Teta_fc= Teta_r+ (Teta_s-Teta_r)/(pow(1+pow(alpha*330,n),1-1/n));  	// soil humidity at field capacity = 330cm pressure head 33kPa
		Teta_wp= Teta_r+ (Teta_s-Teta_r)/(pow(1+pow(alpha*15000,n),1-1/n));  // soil humidity at field capacity = 15000cm pressure head 1.5MPa
		RWC_fc= (Teta_fc-Teta_r)/(Teta_fc-Teta_r);                         	// soil RWC at field capacity
		convert();
		PREM=1;
		init();
		init();  //do not know why it should be done twice???
		PREM=0;
		Tleaf(); // to init g_bl
		if (T_g_cuti) Compute_g_cuti();
		if (!CUT) if (Regul_gs==1) Compute_Turgor_Ref();
			
		if (PREM1)
			{
				if (debug==2) if (CLIMAT==2 || CLIMAT==4)  compute_climatic_stats();
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
					}
				fclose(transient_out);
				
				if (TRANSIENT==11 || TRANSIENT==21) // use standard daily values
				{
					for (i=0;i<250;i++) 
					{
						Screen_out[i]=0;
						File_out[i]=0;
					}
				Screen_out[0]=1;Screen_out[17]=1;Screen_out[30]=1;Screen_out[33]=1;Screen_out[78]=1;Screen_out[79]=1;Screen_out[104]=1;Screen_out[186]=1;Screen_out[119]=1;
				File_out[0]=1;File_out[17]=1;File_out[30]=1;File_out[33]=1;File_out[78]=1;File_out[79]=1;File_out[104]=1;File_out[186]=1;File_out[119]=1;	
				}
				
				if (TRANSIENT==12 || TRANSIENT==22) // use standard yearly values
				{
					for (i=0;i<250;i++) 
					{
						Screen_out[i]=0;
						File_out[i]=0;
					}
					Screen_out[0]=1;Screen_out[174]=1;Screen_out[147]=1;Screen_out[78]=1;Screen_out[79]=1;Screen_out[142]=1;Screen_out[120]=1;Screen_out[127]=1;Screen_out[129]=1;
					File_out[0]=1;	File_out[174]=1;	File_out[147]=1;	File_out[78]=1;File_out[79]=1;		File_out[142]=1;File_out[120]=1;	File_out[127]=1;File_out[129]=1;	
				}
				
				if ((File = fopen("annual_out.csv","r"))==NULL) //file does not exist then print headers in the file
				{
					if ((File = fopen("annual_out.csv","w"))==NULL) printf("\a\nCan't create file annual_out.csv !!");
					else
					{
						fprintf(File,"File;Simul#;Year; ");
						for (i=0;i<250;i++)  if (File_out[i])  fprintf(File, "%s; ", Label[i]);
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
							for (i=0;i<250;i++)
							{
								fprintf(transient, "%s; ", Label[i]);
								File_out[i]=1;
							}
						else for (i=0;i<250;i++)  if (File_out[i])  fprintf(transient, "%s; ", Label[i]);
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
				Graph_ps();
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
	int parameter[250];
	int i,j,N_para;
	
	
	if ((out = fopen(filename_OUT,"a"))==NULL) printf("\a\nCan't create file sureau.out!!");
	fclose(out);
	if ((para_file = fopen("sureau_para.txt","r"))==NULL) //para file not found; create a default one
	{
		printf("\a\nCan't find file sureau_para.txt!!");
		default_para_file(); 
		para_file = fopen("sureau_para.txt","r");
	}
	i=1;
	while (!feof(para_file)) 
	{
		fscanf(para_file,"%d\n",&parameter[i]);
		i++;
	}
	fclose(para_file);
	N_para=i-1;
	
	if ((init_file = fopen(filename_IN,"r"))==NULL) //init file not found; create a default one
	{
		printf("\a\nCan't find file sureau_ini.txt!!");
		default_para(); 
		init_file = fopen(filename_IN,"r");
	}

	while (!feof(init_file))  //for all the simulations in sureau_ini.txt
		{
		for (j=1;j<=N_para;j++) 
		{
			if (!feof(init_file)) 
				{
				fscanf(init_file,"%Le",&para0[parameter[j]]);
				para[parameter[j]]=para0[parameter[j]];
				}			
			else break;
		}
		setparameters();
			
		if (RANDOMISE) for(i=0;i<fabsl(RANDOMISE); i++) 
			{
				Randomise();
				setparameters();
				checkparameters();
				initialise();
				compute();
				N++;
			}
			else
			{
			initialise();
			compute();
		//	printf("\n");
			N++;
			}
		} // end of ini file
		
	if (TRANSIENT)
		{
		transient = fopen("transient_out.csv","a");
		fprintf(transient,"\n");
		fclose(transient);
		}
	
	fclose (init_file);
}

int main(int argc, char * argv[])
{
	char buffer[40];
	unsigned number=1;
	int random1;
	long double random2;
	FILE *out, *random_file;
	printf("SurEau by H. Cochard UMR Piaf-INRAE version:%s \n",version) ;
	//printf("Simulation 1 Year 1\n");
	
	strcpy (filename_IN, "sureau_ini.txt");  // default values when no name is given as argument
	strcpy (filename_OUT,"sureau_out.csv");
	strcpy (filename_OUT2,"annual_out.csv");
	strcpy (filename_CLIM, "climat_day_in.txt");
	strcpy (filename_TRANS,"transient_out.csv");
	
	if (argc==2) // passed with iXXX_sureau.ini OR cXXX_climat_day.ini as argument
	{
		strncpy(filenumber,argv[1],1);                                  // first character of the argument
		if  (strcmp(filenumber,"i")  ==0) strcpy(filename_IN,argv[1]);  //a ini file
		else strcpy(filename_CLIM,argv[1]);                             //a clim file
		
		strncpy(filenumber,argv[1],4);                                  //extract the file number from file name
		number=atoi(filenumber+1);
	
		strcpy(buffer,filenumber);
		strcat(buffer,"_sureau.out1");
		strcpy(filename_OUT,buffer);
		
		strcpy(buffer,filenumber);
		strcat(buffer,"_annual.out2");
		strcpy(filename_OUT2,buffer);
		
		strcpy(buffer,filenumber);
		strcat(buffer,"_transient.out3");
		strcpy(filename_TRANS,buffer);
	}
	if (argc==3) // passed with cXXX_climat_day.ini AND iXXX_sureau.ini as argument
	{
	   // strncpy(filenumber,argv[1],1);                                  // first character of the argument
		strcpy(filename_CLIM,argv[1]);                           	//clim file
		strcpy(filename_IN,argv[2]);  									//ini file
		
		strncpy(filenumber,argv[1],4);                                  //extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcat(buffer,filenumber);
		strcat(buffer,"_sureau.out1");
		strcpy(filename_OUT,buffer);
		
		strncpy(filenumber,argv[1],4);                                  //extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcpy(buffer,filenumber);
		strcat(buffer,"_annual.out2");
		strcpy(filename_OUT2,buffer);
		
		strncpy(filenumber,argv[1],4);                                  //extract the file number from file name
		strcpy(buffer,filenumber);
		strncpy(filenumber,argv[2],4);
		strcpy(buffer,filenumber);
		strcat(buffer,"_transient.out3");
		strcpy(filename_TRANS,buffer);
	}
	
	if ((out = fopen("sureau_out.csv","r"))==NULL) //file does not exist then print headers in the file
	{
		out = fopen("sureau_out.csv","w");
		fprintf(out, "N1;N2;YEAR;T_PLC_Leaf;T_PLC_Axil;T_PLC_Branch;T_PLC_Trunk;T_PLC_Root;T_PLC_Root1;T_RWC_Axil;T_REW_Soil; T_gs_regul;T_gs_close;A_net_tot;A_net_max;E_max;T_Leaf_max;PLC_leaf;PLC_branch;PLC_Plant;GPP;Q_Plant;Q_Plant_s;Q_Plant_a; K_Plant;K_root;Q_Branch_Symp;Q_Branch_Apo;Q_Trunk_Sym;Q_Trunk_Apo;Q_Root_Symp;Q_Root_Apo; Branch_Area;Trunk_Area;Root_Area;");
	   
		if ((random_file = fopen("random_para.txt","r"))==NULL) printf("\a\nCan't find file random_para.txt!!");
			else while (!feof(random_file)) 
			{
			fscanf(random_file,"%d %Le %Le %Le\n",&random1,&random2,&random2,&random2);
			fprintf(out,"para# %d;",random1);
			}
			fprintf(out,"\n");
			fclose(random_file);
		
		
		fclose(out);
		fprintf(out,"\n");
		fclose(out);
	}
	if (argc==2) srand((unsigned) number+time(NULL));
	else srand((unsigned) time(NULL));
	setup();
	if (PRINT_GRAPH) Graph_ps();
}
